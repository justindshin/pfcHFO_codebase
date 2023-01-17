function jds_replay_decoding_CA1(animalprefix,day,ep,cellcountthresh)
% addpath('C:\Users\wbtang\Desktop\HMM\Src_Matlab\Functions');
% addpath(genpath('C:\Users\wbtang\Desktop\HMM\Src_Matlab'))
% savedir = ('D:\Single_Day_Experiment\ER1_ProcessedData\');
% ripdir = ('D:\Single_Day_Experiment\ER1_NEW_direct\');

savedata = 1;

if strcmp(animalprefix,'ZT2')
    savedir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
elseif strcmp(animalprefix,'JS34')
    savedir = ('/Volumes/JUSTIN/SingleDay/JS34_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/JS34_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/JS34_direct/');
elseif strcmp(animalprefix,'BG1')
    savedir = ('/Volumes/JUSTIN/SingleDay/BG1_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/BG1_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/BG1_direct/');
elseif strcmp(animalprefix,'JS17')
    savedir = ('/Volumes/JUSTIN/SingleDay/JS17_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/JS17_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/JS17_direct/');
elseif strcmp(animalprefix,'JS21')
    savedir = ('/Volumes/JUSTIN/SingleDay/JS21_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/JS21_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/JS21_direct/');
elseif strcmp(animalprefix,'JS14')
    savedir = ('/Volumes/JUSTIN/SingleDay/JS14_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/JS14_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/JS14_direct/');
elseif strcmp(animalprefix,'JS15')
    savedir = ('/Volumes/JUSTIN/SingleDay/JS15_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/JS15_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/JS15_direct/');
elseif strcmp(animalprefix,'KL8')
    savedir = ('/Volumes/JUSTIN/SingleDay/KL8_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/KL8_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/KL8_direct/');
elseif strcmp(animalprefix,'ER1')
    savedir = ('/Volumes/JUSTIN/SingleDay/ER1_direct/');
    ripdir = ('/Volumes/JUSTIN/SingleDay/ER1_direct/');
    dir = ('/Volumes/JUSTIN/SingleDay/ER1_direct/');
end

% if animalprefix == 'KL8'
%     if ep == [3:2:11, 15 17];
%         eprun = (ep-1);% using current epoch for run epoch
%     elseif ep == [2 4 6 8 10 12 14 16];
%         eprun = [2 4 6 8 10 12 14 16];% using the last run epoch for sleep epoch
%     end
% else
%     if ep == [4 6 8 10 12 14 16];
%         eprun = [4 6 8 10 12 14 16];
%     else
%         eprun = [4 6 8 10 12 14 16];
%     end
% end

if mod(ep,2) == 0
    eprun = ep;
elseif ep == 1
    eprun = 2;
elseif mod(ep,2) ~= 0
    eprun = (ep-1);
end


switch animalprefix
    case 'ZT2'
        animal_order = 1;
    case 'JS34'
        animal_order = 2;
    case 'BG1'
        animal_order = 3;
    case 'JS17'
        animal_order = 4;
    case 'JS21'
        animal_order = 5;
    case 'JS14'
        animal_order = 6;
    case 'JS15'
        animal_order = 7;
    case 'KL8'
        animal_order = 8;
    case 'ER1'
        animal_order = 9;
end

%%
%-----match neurons across epochs-----%
[ctxidx, hpidx] = matchidx_across2ep_singleday(dir, animalprefix, day, [eprun ep], []); %(tet, cell)
ctxnum = length(ctxidx(:,1));
hpnum = length(hpidx(:,1));
tBinSz = 15; %default temporal bin in ms [typically around 15ms]
wellcutoff = 0; %cm; make 0 for sleep
%%
%-----create the ratemaps [nPosBin x nHPCells]-----%
rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
cellidxm = [];
load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
for i = 1:hpnum
    cind = hpidx(i,:);
    if (length(linfields{day}{eprun})>= cind(1))
        if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
            linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
        else
            linfield1 =[];
        end
    else
        linfield1=[];
    end

    if ~isempty(linfield1)
        linfield_hp = [];
        lintrack_hp = [];
        pos_hp = [];
        for track = 1:4 %concatenating all linearized firing rates
            temp1 = linfield1{track};
            pos1 = temp1(:,1);
            lintrack1 = ones(size(pos1))*track; %get length of lin traj and convert to vector of traj number
            occnormrate1 = temp1(:,5);
            linfield_hp = [linfield_hp;occnormrate1];
            pos_hp = [pos_hp;pos1];
            lintrack_hp = [lintrack_hp;lintrack1];
        end
        if (max(linfield_hp) >= 3) % peak firing rate max larger than 3 Hz
            rm = [rm;linfield_hp']; %horizontal vector of firing rates (each row, diff cell)
            pm = [pm;pos_hp']; %horizontal vector of pos
            tm = [tm;lintrack_hp'];
            cellidxm = [cellidxm; cind];
        end
    end
end
rm = rm'; %[nPosBin x nHPCells]
pm = pm';
tm = tm';
for i = 1:4 %Cutt off reward well positions
    pm_traj = pm(find(tm == i));
    maxpos = max(max(pm_traj));
    rm(find(tm == i & pm <= wellcutoff)) = 0;
    rm(find(tm == i & pm >= maxpos-wellcutoff)) = 0;
end
rm = rm+ (eps.^8); %Add a small number so there are no zeros
expecSpk =rm.*tBinSz./1000; %[nPos x nCells] Expected number of spikes per bin
hpnum = length(rm(1,:)); % update (only cells with >3hz FR)

%%
%-----create the event matrix during SWRs-----%
spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
load(sprintf('%s/%srippletime_noncoordSWS%02d.mat', dir, animalprefix, day));
nc_rip = ripple; clear ripple
load(sprintf('%s/%srippletime_coordSWS%02d.mat', dir, animalprefix, day));
load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
c_rip = ripple; clear ripple

nc_riptimes = [nc_rip{day}{ep}.starttime nc_rip{day}{ep}.endtime];
nc_riptimes(:,3) = 1;

c_riptimes = [c_rip{day}{ep}.starttime c_rip{day}{ep}.endtime];
c_riptimes(:,3) = 2;

%combine riptimes
riptimes = sortrows([nc_riptimes; c_riptimes],1);

rip_starttime = 1000*riptimes(:,1);  % in ms
rip_endtime = 1000*riptimes(:,2);  % in ms
%
% % Find ripples separated by at least a second -- sj_getpopulationevents2
%
% iri = diff(rip_starttime);
% keepidx = [1;find(iri>=1000)+1];
%
riplength = rip_endtime - rip_starttime;
keepidx2 = find(riplength >= 50);% use the ripple last for more than 50 ms
% keepidx = intersect(keepidx,keepidx2);
% %
riptimes = riptimes(keepidx2,:);

ripnum = size(riptimes,1);

%%
if ripnum > 1
    celldata = [];
    spikecounts = [];
    for cellcount = 1:hpnum %get spikes for each cell
        index = [day,ep,cellidxm(cellcount,:)] ;
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
            spikebins = spikebins(validspikes);
            tmpcelldata = [spiketimes spikebins];
        end
        if ~isempty(spiketimes)
            tmpcelldata(:,3) = cellcount; %keep count of how many cells active during rip event
        else
            tmpcelldata = [0 0 cellcount];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount]; %concatenating num spikes per cell, per event
    end
    cellcounts = sum((spikecounts > 0));
    eventindex = find(cellcounts >= cellcountthresh); %Find event indices that have more than cellthresh # of cells
    if ~isempty(eventindex)
        for event = 1:length(eventindex)
            event
            rippletype = riptimes(eventindex(event),3);
            cellsi = celldata(find(celldata(:,2)==eventindex(event)),3);
            [cellsi,ia] = unique(cellsi,'first');
            [~,sortorder] = sort(ia);
            event_cellSeq = cellsi(sortorder); %sort in order of spikes
            tmpind = find(celldata(:,2) == eventindex(event));
            spiketimes = celldata(tmpind,1); %spike times during ripple event
            cellindex = celldata(tmpind,3);
            %-----create the event matrix during SWRs (spkT{cells}.spiketimes) -----%
            for cell = event_cellSeq'
                validspikeidx = find(cellindex == cell);
                spkT{cell} = spiketimes(validspikeidx).*1000;
            end
            %--------calculate the posterior probability on an event-by-event basis--------%
            startevent = riptimes(eventindex(event),1).*1000;
            endevent = riptimes(eventindex(event),2).*1000;
            timebins = startevent:tBinSz:endevent; % timebins are the binedges
            nTBin = length(timebins)-1;
            nCell = hpnum;
            spkPerBin = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
            for nn  = 1:hpnum
                cellInd = nn; %current cell
                if length(spkT) >= cellInd
                    if ~isempty(spkT{cellInd})
                        %                     spkPerBin(1,:,cellInd) = histcounts(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
                        temp = histc(spkT{cellInd}, timebins); %[1 x nTBin x nCell] Number of spikes per timebing per cell
                        spkPerBin(1,:,cellInd) = temp(1:end-1);
                    end
                end
            end
            nSpkPerTBin = squeeze(sum(spkPerBin,3)); %[nTBin x 1] number of spikes in tBin
            nPBin = size(rm,1); %N positional bin
            expecSpk  = reshape(expecSpk, [nPBin,1, nCell]); %[nPos x 1 x nCell]
            expon = exp(-expecSpk); %Exponent of equation. Need explanation for this
            factSpkPerBin = factorial(spkPerBin); %Factorial to divide by

            wrking = bsxfun(@power, expecSpk, spkPerBin); %[nPos x nTbin x nCell]
            wrking = bsxfun(@rdivide, wrking, factSpkPerBin); %[nPos x nTbin x nCell]
            wrking = bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
            post = prod(wrking,3); %Non normalised prob [nPos x Tbin]
            post(:,nSpkPerTBin==0)  =0; % FO changed this 27/5/15 so the posterior matrix can be smoothed.
            post(isnan(post)) = 0;
            trajinfo = mean(tm,2);% trajactory number
            jumpDistTmp = nan(1,4);
            for traj = 1:4 % four trajactory
                trajidx = find(trajinfo == traj);
                pMat = post(trajidx,:);% create a posterior matrix for each traj (split post by traj)
                distvector = pm(trajidx,1)';

                szPM1 = size(pMat,1);
                szPM2 = size(pMat,2);
                for i = 1:szPM2
                    if (sum(pMat(:,i))>0)
                        pMat(:,i) = pMat(:,i)./sum(pMat(:,i)); % normalized across positions to 1 for each time bin
                    end;
                end
                %get jump distances
                jDistBins = [];
                for j = 1:szPM2-1
                    tmp = pMat(:,j);
                    tmp2 = pMat(:,j+1);
                    if (sum(tmp) ~= 0) && (sum(tmp2) ~= 0)
                        [M1,I1] = max(tmp);
                        [M2,I2] = max(tmp2);
                        jDist = abs(I1-I2)*2;
                        jDistBins = [jDistBins; jDist];
                    end
                end
                jumpDistTmp(traj) = mean(jDistBins);
                nonzerobins = find(nSpkPerTBin > 0);
                rvalues = [];
                slopes = [];
                entropy = [];
                interc = [];
                totalsamples = 10000;
                for rloop = 1:500
                    tBinPicks = distsample(totalsamples,nSpkPerTBin);
                    regressdata = [];
                    for i = 1:length(nonzerobins)
                        if (nSpkPerTBin(nonzerobins(i)) > 0)
                            tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                            if ~isempty(find(pMat(:,nonzerobins(i)) ~= 0))
                                distpicks = distvector(distsample(tmpnumsamples,pMat(:,nonzerobins(i))))';
                                entropy_loop(i) = -nansum((hist(distpicks,0:5:200)./length(distpicks)).*log(hist(distpicks,0:5:200)./length(distpicks)));
                                distpicks(:,2) = nonzerobins(i);
                                regressdata = [regressdata; distpicks];
                            end
                        end
                    end
                    regressdata(:,3) = 1;
                    [b,bint,r,rint,stats] = regress(regressdata(:,1),[regressdata(:,3),regressdata(:,2)]);
                    rvalues = [rvalues; stats(1)];
                    slopes = [slopes; b(2)];
                    interc = [interc; b(1)];

                    entropy = [entropy; mean(entropy_loop)];
                end
                Res(event,traj) =  mean(rvalues);
                YInterc(event,traj) = mean(interc);
                Spd(event,traj) = mean(slopes);
                Entropy(event,traj) = mean(entropy);
                pMat_cell{event}{traj}.pMat = pMat;
                pMat_cell{event}{traj}.timevec = 1:szPM2;
                pMat_cell{event}{traj}.posvec = distvector;
                pMat_cell{event}{traj}.timebinsz = tBinSz;
%                 % debug stuff...
%                 imagesc(1:szPM2,distvector,pMat);
%                 colormap(jet)
%                 hold on
%                 plot([1, szPM2], [0, Spd(event,traj) * (szPM2-1)] + (YInterc(event,traj) +1), 'w');
%                 hold off
                %-------Shuffling to get the pvalue for each traj------%
                permbins = nonzerobins;
                srvalues = [];
                for iteration = 1:1500
                    permbins = permbins(randperm(length(permbins)));
                    tmpspkPerBin = zeros(size(spkPerBin));
                    tmpspkPerBin(:,permbins,:) = spkPerBin(:,nonzerobins,:);
                    tmpfactSpkPerBin = factorial(tmpspkPerBin); %Factorial to divide by
                    tmpnSpkPerTBin = squeeze(sum(tmpspkPerBin,3)); %[nTBin x 1] number of spikes in tBin

                    tmpexpecSpk = expecSpk(trajidx,:,:);
                    tmpexpon = expon(trajidx,:,:);
                    wrking = bsxfun(@power, tmpexpecSpk, tmpspkPerBin); %[nPos x nTbin x nCell]
                    wrking = bsxfun(@rdivide, wrking, tmpfactSpkPerBin); %[nPos x nTbin x nCell]
                    wrking = bsxfun(@times,wrking, tmpexpon); %[nPos x nTbin x nCell]
                    tmppMat = prod(wrking,3); %Non normalised prob [nPos x Tbin]
                    tmppMat(:,tmpnSpkPerTBin==0)  =0; % FO changed this 27/5/15 so the posterior matrix can be smoothed.
                    tmppMat(isnan(tmppMat)) = 0;
                    for i = 1:szPM2; if (sum(tmppMat(:,i))>0)
                            tmppMat(:,i) = tmppMat(:,i)./sum(tmppMat(:,i));end;end % normalized across positions to 1 for each time bin
                    clear wrking tmpfactSpkPerBin tmpexpon
                    tBinPicks = distsample(totalsamples,nSpkPerTBin);
                    regressdata = [];
                    for i = 1:length(permbins)
                        if (nSpkPerTBin(nonzerobins(i)) > 0)
                            tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                            if ~isempty(find(pMat(:,nonzerobins(i)) ~= 0))
                                distpicks = distvector(distsample(tmpnumsamples,tmppMat(:,permbins(i))))';
                                distpicks(:,2) = permbins(i);
                                regressdata = [regressdata; distpicks];
                            end
                        end
                    end
                    regressdata(:,3) = 1;
                    warning('off','all')
                    [b,bint,r,rint,stats] = regress(regressdata(:,1),[regressdata(:,3),regressdata(:,2)]);
                    srvalues = [srvalues; stats(1)];
                end
                pvalue(event,traj) =sum(Res(event,traj) < srvalues)/length(srvalues);
                shuffle_rvalues{event}{traj} = srvalues;
            end
            [minP,tidx] = min(pvalue(event,:));
            if minP < 0.05
                decode_traj(event) = tidx;
            else
                decode_traj(event) = 0;% no significant traj
            end
            activecell{event} = cellsi;
            activecellidx{event} = cellidxm(cellsi,:);
            jumpDistance{event} = jumpDistTmp;
            if rippletype == 1
                eventripple(event) = 1;
            elseif rippletype == 2
                eventripple(event) = 2;
            end
            clear wrking factSpkPerBin expon
        end

        %%
        replaytraj.pMat = pMat_cell;
        replaytraj.eventtimes = riptimes(eventindex,[1 2]);
        replaytraj.rsquare = Res;
        replaytraj.slopes = Spd;
        replaytraj.YInterc = YInterc;
        replaytraj.Entropy = Entropy;
        replaytraj.eventidx = eventindex;
        replaytraj.besttraj = decode_traj;
        replaytraj.shuffle_rsquare = shuffle_rvalues;
        replaytraj.jumpdistances = jumpDistance;
        replaytraj.ripType = eventripple;
        replaytraj.pvalue = pvalue;
        replaytraj.activecell = activecell;
        replaytraj.activecellidx = activecellidx;
        replaytraj.sigeventprc = length(find(decode_traj~=0))./length(decode_traj);
        replaytraj.sigeventnum = length(find(decode_traj~=0));
        replaytraj.candeventnum = length(decode_traj);
        replaytraj.allripples = riptimes;

        replaytrajactory{day}{ep} = replaytraj;
        %%
        if savedata
            save(sprintf('%s%sreplaydecode%02d_%02d.mat', savedir,animalprefix,day,ep), 'replaytraj');
        end
    end
end