
close all;
clear all;
%%
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS15','JS14','ER1','KL8'}; % animal prefix
animaltestday = 1;% animal experimental day
eps = 1:2:17;% epochs
%%
%---- set parameters ----%
minstd = 3; % min std above mean for ripple detection
minrip = 1; % min number of tetrode detected ripples
% minenergy = 0; % min energy threshold
velfilter = 4; %velocity <= 4cm for ripple detection
matcheegtime = 0; % match EEG time?
mismatch = 0;
day = 1;

analyzeData = 1;
indAmp = [];
coordAmp = [];
indFreq = [];
coordFreq = [];
indLength = [];
coordLength = [];
%%
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    % set animal directory
    animaldir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    if analyzeData == 1
        load(sprintf('%s/%sctxripple_amp_freq_%02d.mat', animaldir, animalprefix, day));
        for ep = 1:length(eps)
            e = eps(ep);
            tmpdata = rippleampfreq{day}{e};
            if ~isempty(tmpdata)
                indAmp = [indAmp; tmpdata(find(tmpdata(:,3)==1),4)];
                coordAmp = [coordAmp; tmpdata(find(tmpdata(:,3)==2),4)];
                indFreq = [indFreq; tmpdata(find(tmpdata(:,3)==1),5)];
                coordFreq = [coordFreq; tmpdata(find(tmpdata(:,3)==2),5)];
                indLength = [indLength; tmpdata(find(tmpdata(:,3)==1),2)-tmpdata(find(tmpdata(:,3)==1),1)];
                coordLength = [coordLength; tmpdata(find(tmpdata(:,3)==2),2)-tmpdata(find(tmpdata(:,3)==2),1)];
            end
        end
    else
        eegdir = [animaldir,'EEG/'];% EEG directory
        %%
        tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo'); % get tetrode info
        %%
        % day loop
        for day = animaltestday
            load(sprintf('%s/%sctxrippletime_noncoordSWS%02d.mat', animaldir, animalprefix, day));
            nc_rip = ctxripple; clear ripple
            load(sprintf('%s/%sctxrippletime_coordSWS%02d.mat', animaldir, animalprefix, day));
            load(sprintf('%s/%sspikes%02d.mat', animaldir, animalprefix, day));
            c_rip = ctxripple; clear ripple


            % loop for each epoch per day
            for id = 1:length(eps)
                d = day;%day
                e = eps(id);%epoch
                nc_riptimes = [nc_rip{day}{e}.starttime nc_rip{day}{e}.endtime];
                nc_riptimes(:,3) = 1;

                c_riptimes = [c_rip{day}{e}.starttime c_rip{day}{e}.endtime];
                c_riptimes(:,3) = 2;

                %combine riptimes
                riptimes = sortrows([nc_riptimes; c_riptimes],1);
                disp(['Animal: ',animalprefix,' Epoch:',num2str(e)])% display current animal, day and epoch

                % calculate using riptet only
                tetfilter = 'isequal($descrip, ''ctxriptet'')';% use riptet only
                tetlist =  evaluatefilter(tetinfo{d}{e}, tetfilter);
                tetlist = unique(tetlist(:,1))';
                if isequal(animalprefix,'JS14')
                    tetlist(find(tetlist == 17)) = [];
                end

                % get the mean and std for each frequency channel for z-scoring
                % note that only run it for the first time
                baselinespecgram_forref(animalprefix, d, e, tetlist, 'fpass',[0 400])% high frequency range, 0-400 Hz, for all tets

                ripples = loaddatastruct(animaldir, animalprefix, 'ctxripples', d);% load ripple info

                r = ripples{d}{e}{tetlist(1)};
                % time range, 10ms bin
                times = r.timerange(1):0.001:r.timerange(end);
                %reset
                nrip = zeros(size(times));
                nstd=[];
                ripplestd = zeros(size(times));

                % tetrode loop
                if ~isempty(riptimes)
                    for t = 1:length(tetlist)
                        tmprip = ripples{d}{e}{tetlist(t)};
                        % get the indeces for the ripples with energy above minenergy
                        % and maxthresh above minstd
%                         rvalid = find((tmprip.energy >= minenergy) & (tmprip.maxthresh >= minstd));
                        rvalid = find(tmprip.maxthresh >= minstd);
                        rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
                        tmpripplestd = [tmprip.maxthresh(rvalid) tmprip.maxthresh(rvalid)];
                        % create another parallel vector with bordering times for zeros
                        nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
                        rtimes = reshape(rtimes', length(rtimes(:)), 1);
                        rtimes(:,2) = 1;
                        tmpriplestd = [rtimes(:,1) tmpripplestd(:)];
                        nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
                            length(nrtimes(:)), 1) ; r.timerange(2)];
                        nrtimes(:,2) = 0;
                        % create a new list with all of the times in it
                        tlist = sortrows([rtimes ; nrtimes]);
                        [junk, ind] = unique(tlist(:,1));
                        tlist = tlist(ind,:);

                        stdlist = sortrows([tmpriplestd ; nrtimes]);
                        stdlist =stdlist(ind,:);
                        nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
                        nstd(t,:) = interp1(stdlist(:,1), stdlist(:,2), times, 'nearest');  % carry forward amplitude of ripple
                    end

                    %find the start and end borders of each ripple
                    inripple = (nrip >= minrip);
                    startrippleind = find(diff(inripple) == 1)+1;
                    endrippleind = find(diff(inripple) == -1)+1;
                    ripplestdout = [];

                    if (endrippleind(1) < startrippleind(1))
                        endrippleind = endrippleind(2:end);
                    end
                    if (endrippleind(end) < startrippleind(end))
                        startrippleind = startrippleind(1:end-1);
                    end
                    startripple = times(startrippleind);
                    endripple = times(endrippleind);
                    %----- measure amplitude of each ripple-----%
                    % Get amplitude of "global" ripple: maximum across tetrodes
                    [max_nstd,tetid] = max(nstd,[],1);
                    ampripple = max_nstd(startrippleind);
                    riptet = tetid(startrippleind);

                    out = [startripple(:) endripple(:) ampripple(:)]; % amplitude of ripple
                    riptimes(:,4) = 0;
                    for r = 1:length(riptimes(:,1))
                        ripmidtmp = riptimes(r,1) + ((riptimes(r,2) - riptimes(r,1))/2);
                        idx = find((ripmidtmp > out(:,1)) & (ripmidtmp < out(:,2)));
                        if ~isempty(idx)
                            amptmp = out(idx,3);
                            riptimes(r,4) = amptmp;
                            riptetnew(r) = riptet(idx);
                        else
                            amptmp = NaN;
                            riptimes(r,4) = amptmp;
                            riptetnew(r) = NaN;
                            mismatch = mismatch+1;
                        end
                    end
                    riptimes(:,5) = 0;

                    %----- measure frequncy of each ripple-----%
                    for r = 1:length(riptimes(:,1))
                        riptime = riptimes(r,1:2);
                        if ~isnan(riptetnew(r))
                            riptimes(r,5) = ripple_frequency_fun(animalprefix, d, e, tetlist(riptetnew(r)), riptime);
                        else
                            riptimes(r,5) = NaN;
                        end
                    end
                else
                    riptimes = [];
                end
                rippleampfreq{d}{e} = riptimes;% save result
                clear out; clear freqripple; clear riptimes
            end
            save(sprintf('%s/%sctxripple_amp_freq_%02d.mat', animaldir, animalprefix, d), 'rippleampfreq');% save files
            clear rippleampfreq;
        end
    end
end

if analyzeData == 1
    
    meanAmpInd = nanmean(indAmp)
    semAmpInd = nanstd(indAmp)./sqrt(length(find(~isnan(indAmp))))
    medianAmpInd = nanmedian(indAmp)

    meanAmpCoord = nanmean(coordAmp)
    semAmpCoord = nanstd(coordAmp)./sqrt(length(find(~isnan(coordAmp))))
    medianAmpCoord = nanmedian(coordAmp)

    meanFreqInd = nanmean(indFreq)
    semFreqInd = nanstd(indFreq)./sqrt(length(find(~isnan(indFreq))))
    medianFreqInd = nanmedian(indFreq)

    meanFreqCoord = nanmean(coordFreq)
    semFreqcoord = nanstd(coordFreq)./sqrt(length(find(~isnan(coordFreq))))
    medianFreqCoord = nanmedian(coordFreq)

    meanLengthInd = nanmean(indLength)
    semLengthInd = nanstd(indLength)./sqrt(length(find(~isnan(indLength))))
    medianLengthInd = nanmedian(indLength)

    meanLengthCoord = nanmean(coordLength)
    semLengthcoord = nanstd(coordLength)./sqrt(length(find(~isnan(coordLength))))
    medianLengthCoord = nanmedian(coordLength)

    [p h] = ranksum(indAmp,coordAmp)
    datacombinedRipAmp = [indAmp; coordAmp];
    g1 = repmat({'Ind'},length(indAmp),1);
    g2 = repmat({'Coord'},length(coordAmp),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipAmp,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    % ylim([-0.02 0.2])
    title(['CA1 - Ripple amplitude-p = ' num2str(p)])
    ylim([2 12])
    set(gcf, 'renderer', 'painters')

    [p2 h2] = ranksum(indFreq,coordFreq)
    datacombinedRipFreq = [indFreq; coordFreq];
    g1 = repmat({'Ind'},length(indFreq),1);
    g2 = repmat({'Coord'},length(coordFreq),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipFreq,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    % ylim([-0.02 0.2])
    title(['CA1 - Ripple frequency-p = ' num2str(p2)])
    set(gcf, 'renderer', 'painters')

    [p3 h3] = ranksum(indLength,coordLength)
    datacombinedRipLength = [indLength; coordLength];
    g1 = repmat({'Ind'},length(indLength),1);
    g2 = repmat({'Coord'},length(coordLength),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipLength,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    % ylim([-0.02 0.2])
    title(['CA1 - Ripple lengths-p = ' num2str(p3)])
    set(gcf, 'renderer', 'painters')
end

