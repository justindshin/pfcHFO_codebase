function jds_ripple_CA1PFCcoreactivationRate(animalprefixlist)

%To do: plot the peak strength within specified ripple events. compare
%ripple types for different assemblies - might need to change binsize for
%raster generation
day = 1;

coReactRates = nan(100,5);
coReactRatesShuf = nan(100,5);
coReactProp = nan(100,5);
for rips = 1:5
    coReactTmp = [];
    coReactTmpShuf = [];
    coReactPropTmp = [];
    for a = 1:length(animalprefixlist)
        
        animalprefix = char(animalprefixlist(a));
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

        %Load reactivation strength file for all assemblies and epochs
        load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
        CA1_R = RtimeStrength; clear RtimeStrength
        load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
        PFC_R = RtimeStrength; clear RtimeStrength

        if rips == 1
            load(sprintf('%s%srippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
        elseif rips == 2
            load(sprintf('%s%srippletime_coordSWS%02d.mat',dir,animalprefix,day));
        elseif rips == 3
            load(sprintf('%s%sripplecoordinationSWS%02d.mat',dir,animalprefix,day));
        elseif rips == 4
            load(sprintf('%s%srippleframes_SWS%02d.mat',dir,animalprefix,day));
        elseif rips == 5
            load(sprintf('%s%srippleframes_SWS%02d.mat',dir,animalprefix,day));
        end

        %Analyze epochs where there are assemblies detected
        epochs = find(~cellfun(@isempty,CA1_R));
        for e = 1:length(epochs)
            ep = epochs(e);
            %Load ripples
            if rips == 1
                riptimes = [ripple{day}{ep}.starttime ripple{day}{ep}.endtime];
            elseif rips == 2
                riptimes = [ripple{day}{ep}.starttime ripple{day}{ep}.endtime];
            elseif rips == 3
                riptimes = [ripplecoordination{day}{ep}.starttime ripplecoordination{day}{ep}.endtime];
            elseif rips == 4
                riptimes = [ripple{day}{ep}.outstarttime ripple{day}{ep}.outendtime];
            elseif rips == 5
                riptimes = [ripple{day}{ep}.instarttime ripple{day}{ep}.inendtime];
            end
            if isempty(riptimes)
                continue
            end
            ripdur = sum(riptimes(:,2) - riptimes(:,1));

            CA1assemblytmp = CA1_R{ep}.reactivationStrength;
            CA1num = length(CA1assemblytmp);
            PFCassemblytmp = PFC_R{ep}.reactivationStrength;
            PFCnum = length(PFCassemblytmp);

            if length(riptimes(:,1)) > 20
                timevec = CA1assemblytmp{1}(:,1);
                react_idx = [];
                for t = 1:length(riptimes)
                    idxst = lookup(riptimes(t,1), timevec);
                    idxend = lookup(riptimes(t,2), timevec);
                    react_idx = [react_idx; [idxst idxend]];
                end

                Rmatrix = zeros(length(riptimes(:,1)),(length(CA1assemblytmp) + length(PFCassemblytmp)));
                pfcMat = zeros(length(riptimes(:,1)),length(PFCassemblytmp));
                ca1Mat = zeros(length(riptimes(:,1)),length(CA1assemblytmp));

                Rmatrixshuf = zeros(length(riptimes(:,1)),(length(CA1assemblytmp) + length(PFCassemblytmp)));
                pfcMatshuf = zeros(length(riptimes(:,1)),length(PFCassemblytmp));
                ca1Matshuf = zeros(length(riptimes(:,1)),length(CA1assemblytmp));
                if (~isempty(CA1assemblytmp)) && (~isempty(PFCassemblytmp))
                    for ii = 1:length(CA1assemblytmp)
                        CA1strengthstmp = CA1assemblytmp{ii}(:,2);
                        CA1strengthstmpshuf = CA1strengthstmp(randperm(length(CA1strengthstmp)));
                        aMn = mean(CA1strengthstmp);
                        aStd = std(CA1strengthstmp);
                        aThresh = aMn + aStd*2;
                        for r = 1:length(react_idx(:,1))
                            tmp = max(CA1strengthstmp(react_idx(r,1):react_idx(r,2))); %get vector of reactivation strenths for specified time period
                            tmpshuf = max(CA1strengthstmpshuf(react_idx(r,1):react_idx(r,2)));
                            if tmp > 5
                                Rmatrix(r,ii) = 1;
                                ca1Mat(r,ii) = 1;
                            end
                            if tmpshuf > 5
                                Rmatrixshuf(r,ii) = 1;
                                ca1Matshuf(r,ii) = 1;
                            end
                        end
                    end

                    for ii = 1:length(PFCassemblytmp)
                        PFCstrengthstmp = PFCassemblytmp{ii}(:,2);
                        PFCstrengthstmpshuf = PFCstrengthstmp(randperm(length(PFCstrengthstmp)));
                        for rr = 1:length(react_idx(:,1))
                            aMn = mean(PFCstrengthstmp);
                            aStd = std(PFCstrengthstmp);
                            aThresh = aMn + aStd*2;
                            tmp = max(PFCstrengthstmp(react_idx(rr,1):react_idx(rr,2))); %get vector of reactivation strenths for specified time period
                            tmpshuf = max(PFCstrengthstmpshuf(react_idx(rr,1):react_idx(rr,2)));
                            if tmp > 5
                                Rmatrix(rr,ii+CA1num) = 1;
                                pfcMat(rr,ii) = 1;
                            end
                            if tmpshuf > 5
                                Rmatrixshuf(rr,ii+CA1num) = 1;
                                pfcMatshuf(rr,ii) = 1;
                            end
                        end
                    end
                    ca1active = sum(ca1Mat,2);
                    pfcactive = sum(pfcMat,2);
                    idx1 = find(ca1active > 0);
                    idx2 = find(pfcactive > 0);
                    idx3 = intersect(idx1,idx2);
                    maxActiveCa1 = ca1active(idx3);
                    maxActivePfc = pfcactive(idx3);
                    mincoactive = min([maxActiveCa1 maxActivePfc]')';
                    rateTmp = sum(mincoactive)/ripdur;
                    coReactTmp = [coReactTmp; rateTmp];
                    propTmp = length(idx3)/length(riptimes(:,1));
                    coReactPropTmp = [coReactPropTmp; propTmp];
                    %% shuffle
                    ca1activeshuf = sum(ca1Matshuf,2);
                    pfcactiveshuf = sum(pfcMatshuf,2);
                    idx1 = find(ca1activeshuf > 0);
                    idx2 = find(pfcactiveshuf > 0);
                    idx3 = intersect(idx1,idx2);
                    maxActiveCa1shuf = ca1activeshuf(idx3);
                    maxActivePfcshuf = pfcactiveshuf(idx3);
                    mincoactive = min([maxActiveCa1shuf maxActivePfcshuf]')';
                    rateTmpshuf = sum(mincoactive)/ripdur;
                    coReactTmpShuf = [coReactTmpShuf; rateTmpshuf];
                end
            end
        end
    end
    coReactRates(1:length(coReactTmp),rips) = coReactTmp; 
    coReactProp(1:length(coReactTmp),rips) = coReactPropTmp;
    coReactRatesShuf(1:length(coReactTmpShuf),rips) = coReactTmpShuf;
end

datameansRate = nanmean(coReactRates);
datasemsRate = [];

datameansProp = nanmean(coReactProp);
datasemsProp = [];
for i = 1:5
    tmp = nanstd(coReactRates(:,i))./sqrt(length(find(~isnan(coReactRates(:,i)))));
    tmp2 = nanstd(coReactProp(:,i))./sqrt(length(find(~isnan(coReactProp(:,i)))));
    datasemsRate(i) = tmp;
    datasemsProp(i) = tmp2;
end

bar(nanmean(coReactRates),'k'); hold on
errorbar(1:5, datameansRate, datasemsRate, '-k','LineStyle','none')
ylabel('CA1-PFC Assembly Coactivation Rate (Events/Sec)')
xticklabels({'IND-CA1','Coord-CA1','Coordinated Ripple','OutFrame','InFrame'})
set(gcf, 'renderer', 'painters')

figure;
bar(nanmean(coReactProp),'k'); hold on
errorbar(1:5, datameansProp, datasemsProp, '-k','LineStyle','none')
ylabel('CA1-PFC Assembly Coactivation (Proportion Ripples Active)')
xticklabels({'IND-CA1','Coord-CA1','Coordinated Ripple','OutFrame','InFrame'})
set(gcf, 'renderer', 'painters')

figure
[p1 h1] = ranksum(coReactRates(:,1),coReactRates(:,2))
datacombinedRate = [coReactRates(:,1); coReactRates(:,2)];
g1 = repmat({'Ind'},length(coReactRates(:,1)),1);
g2 = repmat({'Coord'},length(coReactRates(:,2)),1);
g12 = [g1;g2];

datacombinedRateShuf = [coReactRates(:,1); coReactRates(:,2); coReactRatesShuf(:,1); coReactRatesShuf(:,2)];
g3 = repmat({'Ind-shuf'},length(coReactRatesShuf(:,1)),1);
g4 = repmat({'Coord-shuf'},length(coReactRatesShuf(:,2)),1);
g1234 = [g1;g2;g3;g4];

h = boxplot(datacombinedRate,g12,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');

title(['CA1-PFC coactivation rate-ranksum p = ' num2str(p1)])
ylabel('CA1-PFC assembly coactivation rate') 
ylim([0 4])

figure
h = boxplot(datacombinedRateShuf,g1234,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');

title(['CA1-PFC coactivation rate-ranksum p = ' num2str(p1)])
ylabel('CA1-PFC assembly coactivation rate') 
ylim([0 4])

figure;
datacombinedProp = [coReactProp(:,1); coReactProp(:,2)];
g1 = repmat({'Ind'},length(coReactProp(:,1)),1);
g2 = repmat({'Coord'},length(coReactProp(:,3)),1);
g = [g1;g2];

h = boxplot(datacombinedProp,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title('CA1-PFC coactivation proportion')
ylim([0 0.5])

keyboard;

    