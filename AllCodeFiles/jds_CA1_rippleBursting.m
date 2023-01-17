%Get firing rate bias and percentage of ripple during which the cell is
%active (control for the cell being active in a different number of PFC vs
%CA1 ripples)
clear all;
close all;
%%
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};
day = 1;
epochs = [1:2:17];

%%
blen = 6/1000;

exc_SWRburstingProb = []; %frac of rips with bursting
inh_SWRburstingProb = [];
nonsig_SWRburstingProb = [];
exc_SWRburstingProp = []; %frac of spks in burst over all rips
inh_SWRburstingProp = [];
nonsig_SWRburstingProp = [];
exc_SWRspkLat = [];
inh_SWRspkLat = [];
nonsig_SWRspkLat = [];

for a = 1:length(animalprefixlist)
    animalprefix = char(animalprefixlist{a});
    animdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',animdir,animalprefix,day));% get ripple time
    load(sprintf('%s%sspikes0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sswsALL0%d.mat',animdir,animalprefix,day));
%     load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHighModStrength0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',animdir,animalprefix,day));
    
    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
    
    for e = 1:length(epochs)
        epoch = epochs(e);
        
        ep2 = find(sleeps(:,2) == epoch);
        modcells = epochModulation.cellidx;
        modvals = epochModulation.modVals(:,ep2);
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhvals = modvals(find(epochModulation.modMat(:,ep2) == -1));
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        excvals = modvals(find(epochModulation.modMat(:,ep2) == 1));
        exccells(:,3) = 1;
        noncells = modcells(find(epochModulation.modMat(:,ep2) == 0),:);
        nonsigvals = modvals(find(epochModulation.modMat(:,ep2) == 0));
        noncells(:,3) = 0;
        
        allcells = [inhcells; exccells; noncells];
        modvals = [inhvals; excvals; nonsigvals];
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
        
        riptimes = [ripple{day}{epoch}.starttime ripple{day}{epoch}.endtime];
        swslist = [sws{day}{epoch}.starttime sws{day}{epoch}.endtime];
        swsdur = sws{day}{epoch}.total_duration;
        curreegfile = [animdir,'/EEG/',animalprefix,'eeg', '01','-',epochstring,'-','02']; %use any tetrode
        load(curreegfile);
        time1 = geteegtimes(eeg{day}{epoch}{2}) ; % construct time array
        
        [~,swsvec] = wb_list2vec(swslist,time1);
        
        [~,ripvec] = wb_list2vec(riptimes,time1);
        
        %%
        swsdur = sum(swslist(:,2) - swslist(:,1));
        ripdur = sum(riptimes(:,2) - riptimes(:,1));
        for cellcount = 1:length(allcells(:,1))
            cellmod = allcells(cellcount,3);
            
            if ~isempty(riptimes)
                index = [day,epoch,allcells(cellcount,[1:2])];
                
                if (length(riptimes(:,1)) > 1) && (length(swslist(:,1)) > 1)
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                        fRate = spikes{index(1)}{index(2)}{index(3)}{index(4)}.meanrate;
                        
                        ripSpkCnt = 0;
                        burstSpks = [];
                        spkLat = [];
                        inRipSpks = [];
                        for r = 1:length(riptimes(:,1))
                            spkRipTmp = find((spiketimes > riptimes(r,1)) & (spiketimes < riptimes(r,2)));
                            spkTimeTmp = spiketimes(spkRipTmp);
                            ripSpkCnt = ripSpkCnt + length(spkTimeTmp);
                            if (~isempty(spkTimeTmp)) %Latency to first spk
                                spkLatTmp = spkTimeTmp(1) - riptimes(r,1);
                                if isnan(spkLatTmp)
                                    keyboard
                                end
                                spkLat = [spkLat; spkLatTmp];
                            end
                            if (~isempty(spkTimeTmp)) && (length(spkTimeTmp) > 2)
                                inRipSpks = [inRipSpks; spkTimeTmp];
                                tmpisi = diff(spkTimeTmp);
                                bspikes = find(tmpisi < blen);
                                if ~isempty(bspikes)
                                    burstSpks = [burstSpks; 1];
                                else
                                    burstSpks = [burstSpks; 0];
                                end
                            end
                        end
                        tmpisiAll = diff(inRipSpks);
                        
                        % find the intervals less than the burst length
                        bspikesAll = find(tmpisiAll < blen);
                        
                        nadjacent = 0;
                        if (length(bspikesAll) > 2)
                            nadjacent = length(find(diff(bspikesAll) == 1));
                        end
                        propb = (length(bspikesAll) * 2 - nadjacent) / length(inRipSpks);
                        if length(inRipSpks) > 10
                            if cellmod == 1
                                exc_SWRburstingProb = [exc_SWRburstingProb; [sum(burstSpks)/length(riptimes(:,1)) modvals(cellcount)]];
                                exc_SWRspkLat = [exc_SWRspkLat; [mean(spkLat) modvals(cellcount)]];
                                exc_SWRburstingProp = [exc_SWRburstingProp; [propb modvals(cellcount)]];
                            elseif cellmod == -1
                                inh_SWRburstingProb = [inh_SWRburstingProb; [sum(burstSpks)/length(riptimes(:,1)) modvals(cellcount)]];
                                inh_SWRspkLat = [inh_SWRspkLat; [mean(spkLat) modvals(cellcount)]];
                                inh_SWRburstingProp = [inh_SWRburstingProp; [propb modvals(cellcount)]];
                            elseif cellmod == 0
                                nonsig_SWRburstingProb = [nonsig_SWRburstingProb; [sum(burstSpks)/length(riptimes(:,1)) modvals(cellcount)]];
                                nonsig_SWRspkLat = [nonsig_SWRspkLat; [mean(spkLat) modvals(cellcount)]];
                                nonsig_SWRburstingProp = [nonsig_SWRburstingProp; [propb modvals(cellcount)]];
                            end
                        end
                    end
                end
            end
        end
    end
end
keyboard
datameans = [mean(exc_SWRburstingProb(:,1)) mean(inh_SWRburstingProb(:,1)) ];
datasems = [(std(exc_SWRburstingProb(:,1))/sqrt(length(exc_SWRburstingProb(:,1))))...
    (std(inh_SWRburstingProb(:,1))/sqrt(length(inh_SWRburstingProb(:,1))))];

bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Bursting Probability')
title('Bursting Probability during NC SWRs')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p h] = ranksum(exc_SWRburstingProb(:,1),inh_SWRburstingProb(:,1))

datacombinedBurst = [exc_SWRburstingProb(:,1); inh_SWRburstingProb(:,1)];
g1 = repmat({'EXC'},length(exc_SWRburstingProb(:,1)),1);
g2 = repmat({'INH'},length(inh_SWRburstingProb(:,1)),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedBurst,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.02 0.2])
title('CA1 SWR Bursting')
title(['CA1 SWR Bursting-p = ' num2str(p)])
%%
datameansLat = [mean(exc_SWRspkLat(:,1)) mean(inh_SWRspkLat(:,1))];
datasemsLat = [(std(exc_SWRspkLat(:,1))/sqrt(length(exc_SWRspkLat(:,1))))...
    (std(inh_SWRspkLat(:,1))/sqrt(length(inh_SWRspkLat(:,1))))];

figure
bar([1:2],datameansLat,'k')
hold on
er = errorbar([1:2],datameansLat,datasemsLat);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Latency to first spk')
title('Latency to first spk during NC SWRs')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p1 h1] = ranksum(exc_SWRspkLat(:,1),inh_SWRspkLat(:,1))

datacombinedLatency = [exc_SWRspkLat(:,1); inh_SWRspkLat(:,1)];
g1 = repmat({'EXC'},length(exc_SWRspkLat(:,1)),1);
g2 = repmat({'INH'},length(inh_SWRspkLat(:,1)),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedLatency,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.02 0.2])
title('CA1 SWR Spike Latencuy')
title(['CA1 SWR Spike Latency-p = ' num2str(p)])

[rburstIn pburstIn] = corrcoef(inh_SWRburstingProb)

quartile_sep = floor(length(inh_SWRburstingProb(:,1))/4);
spkburstquar = sortrows(inh_SWRburstingProb,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spkburstquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spkburstquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];

figure; hold on; errorbar(X, data_means(:,2), data_sems(:,2),'k','LineWidth',3);
xlim([0.5 4.5])
ylabel('Modulation Index')
xlabel('Burst Probability Quartile')

[rburstEx pburstEx] = corrcoef(exc_SWRburstingProb)

quartile_sep = floor(length(exc_SWRburstingProb(:,1))/4);
spkburstquar = sortrows(exc_SWRburstingProb,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spkburstquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spkburstquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];

errorbar(X, data_means(:,2), data_sems(:,2),'r','LineWidth',3);
xlim([0.5 4.5])
ylabel('Modulation Index')
xlabel('Burst Probability Quartile')

legend({['INH-p = ' num2str(pburstIn(1,2)) 'r=' num2str(rburstIn(1,2))] ['EXC-p = ' num2str(pburstEx(1,2)) 'r=' num2str(rburstEx(1,2))]})
title('CA1 SWR Bursting - Quartile')

%%
%Latency Quartile

[rlatIn platIn] = corrcoef(inh_SWRspkLat)

quartile_sep = floor(length(inh_SWRspkLat(:,1))/4);
spklatquar = sortrows(inh_SWRspkLat,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spklatquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spklatquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];

figure; hold on; errorbar(X, data_means(:,2), data_sems(:,2),'k','LineWidth',3);
xlim([0.5 4.5])
ylabel('Modulation Index')
xlabel('Spike Latency Quartile')

[rlatEx platEx] = corrcoef(exc_SWRspkLat)

quartile_sep = floor(length(exc_SWRspkLat(:,1))/4);
spklatquar = sortrows(exc_SWRspkLat,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spklatquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spklatquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];

errorbar(X, data_means(:,2), data_sems(:,2),'r','LineWidth',3);
xlim([0.5 4.5])
ylabel('Modulation Index')
xlabel('Spike Latency Quartile')

legend({['INH-p = ' num2str(platIn(1,2)) 'r=' num2str(rlatIn(1,2))] ['EXC-p = ' num2str(platEx(1,2)) 'r=' num2str(rlatEx(1,2))]})
title('CA1 SWR Spike Latency - Quartile')
ylim([-0.7 0.9])

keyboard;