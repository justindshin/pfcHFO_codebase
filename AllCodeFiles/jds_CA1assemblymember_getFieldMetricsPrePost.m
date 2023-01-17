function jds_CA1assemblymember_getFieldMetricsPrePost(animalprefixlist)
%---------------------------------------------------------------%
%This function calculates the number of place fields for cells,
%within/out of field spike ratio, path equivalence,
% behavioral sequence score, and field width
%Compares CA1 exc and inh cells
%---------------------------------------------------------------%

%%
FieldsCnt = [];
Ratio = [];
placeFieldWidth = []; %If multiple fields, average width of fields
pathEquiv = [];
behSeq = [];
cellMod = [];

pre = 1;
post = 0;

day = 1;
epochs = [2:2:16];
numshufs = 10;
% numshufs = 100;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    load(sprintf('%s%sCA1ctxripallmod_epsExcludeHigh0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
    %%
    %-----match neurons across epochs-----%
    for e = 1:length(epochs)
        modcells = epochModulation.cellidx;
        if pre == 1
            eps = [epochs(e) (epochs(e) + 1)];
        elseif post == 1
            eps = [epochs(e) (epochs(e) - 1)];
        end
        epsleep = eps(2);

        sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

        ep2 = find(sleeps(:,2) == epsleep);
        modvals = epochModulation.modVals(:,ep2);
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = modvals(find(epochModulation.modMat(:,ep2) == -1));
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = modvals(find(epochModulation.modMat(:,ep2) == 1));

        allmodcells = [inhcells; exccells];

        %%
        if ~isempty(allmodcells)
            assemblies = RtimeStrength{epsleep};

            for as = 1:length(assemblies.reactivationStrength)
                runchange = assemblies.runchange(as);
                wts = assemblies.weights{as};
                mnWt = mean(wts(:,3));
                stdWt = std(wts(:,3));
                wtThresh = mnWt + stdWt*2;
                highIdx = find(wts(:,3) > wtThresh);
                highCells = wts(highIdx, [1 2]);
                for cellcount = 1:length(highCells(:,1)) %get spikes for each cell
                    [b cellLoc] = ismember(highCells(cellcount,:), allmodcells(:,[1 2]), 'rows', 'legacy');
                    if b
                        cellmod = allmodcells(cellLoc, 3);
                        cellMod = [cellMod; cellmod];
                        ca_run = [];
                        cind = highCells(cellcount,:);
                        run_ep = eps(1);
                        NumSpks = [];
                        LinTmp = [];
                        OccTmp = [];
                        SmOccTmp = [];
                        fields = [];
                        for t = 1:4
                            tmpfield1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,5); %field
                            pos1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,1); %position bin
                            occ1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,2); %not smoothed occupancy
                            smocc1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,6); %smoothed occ
                            NumSpksTmp = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,3); %num spks per bin

                            stdLength1 = linspace(pos1(1),pos1(end),100);
                            tmp1 = lookup(stdLength1,pos1);
                            tmpfield1std = tmpfield1(tmp1);

                            fields{t} = tmpfield1std;

                            OccTmp = [OccTmp; occ1];
                            SmOccTmp = [SmOccTmp; smocc1];
                            NumSpks = [NumSpks; NumSpksTmp];
                            LinTmp = [LinTmp; tmpfield1];
                        end

                        shuf = [];
                        for s = 1:numshufs
                            sampledData = distsample(sum(NumSpks),OccTmp);
                            spksPerBin = histcounts(sampledData,length(SmOccTmp));

                            shufLinField = gaussSmooth(spksPerBin',2)./SmOccTmp;
                            lowocc = find(SmOccTmp < 2*.1); %very low occupancy bins should be excluded
                            shufLinField(lowocc) = nan;

                            shuf = [shuf; shufLinField'];
                        end

                        fieldVec = zeros(1,length(LinTmp));
                        for l = 1:length(LinTmp)
                            shufTmp = shuf(:,l);
                            pvalue = 1 - sum(shufTmp < LinTmp(l))/length(shufTmp);
                            if pvalue <= 0.01
                                fieldVec(l) = 1;
                            end
                        end

                        binaryVector = fieldVec;
                        % Label each region with a label - an "ID" number.
                        [labeledVector, numRegions] = bwlabel(binaryVector);
                        % Measure lengths of each region and the indexes
                        measurements = regionprops(labeledVector, LinTmp', 'Area', 'PixelValues');
                        % Find regions where the area (length) are 5 or greater and
                        % put the values into a cell of a cell array
                        for k = 1:numRegions
                            if measurements(k).Area >= 5
                                % Width is 5 bins or greater, so store the values.
                                ca_run{k} = measurements(k).PixelValues;
                            end
                        end

                        fieldIdx = [];
                        for f = 1:length(ca_run)
                            tmp = ca_run{f};
                            if ~isempty(tmp)
                                idx = find(labeledVector == f);
                                fieldIdx = [fieldIdx; idx'];
                            end
                        end

                        if ~isempty(ca_run)
                            %this is using number of spikes for ratio
                            %calculation
                            %                             SpksInField = sum(NumSpks(fieldIdx));
                            %                             SpksOutField = sum(NumSpks(setdiff(1:length(NumSpks),fieldIdx)));

                            %this is using the mean occupancy normalized
                            %firing rate - better metric than above because
                            %occupancy is taken into consideration
                            NumFields = sum(~cellfun(@isempty, ca_run));
                            FieldWidth = nanmean(cellfun(@length, ca_run))*2;
                            
                            SpksInField = nanmean(LinTmp(fieldIdx));
                            [C notInIdx] = setdiff(1:length(LinTmp),fieldIdx);
                            SpksOutField = nanmean(LinTmp(notInIdx));

                            spkRatio = SpksInField/SpksOutField;

                        else
                            SpksInField = NaN;
                            SpksOutField = NaN;

                            NumFields = NaN;
                            FieldWidth = NaN;

                            spkRatio = NaN;
                        end
                        beh1_4 = corrcoef(fields{1},fields{4},'rows','complete');
                        beh2_3 = corrcoef(fields{2},fields{3},'rows','complete');
                        runBeh = [beh1_4(1,2);beh2_3(1,2)];
                        eq1_2 = corrcoef(fields{1},fields{2},'rows','complete');
                        eq3_4 = corrcoef(fields{3},fields{4},'rows','complete');
                        runEquiv = [eq1_2(1,2); eq3_4(1,2)];

                        FieldsCnt = [FieldsCnt; [NumFields cellmod runchange]];
                        Ratio = [Ratio; [spkRatio cellmod runchange]];
                        placeFieldWidth = [placeFieldWidth; [FieldWidth cellmod runchange]];
                        pathEquiv = [pathEquiv; [max(runEquiv) cellmod runchange]];
                        behSeq = [behSeq; [max(runBeh) cellmod runchange]];
                    end
                end
            end
        end
    end
end

keyboard

%%
%num fields
mnFields = [nanmean(FieldsCnt) nanmean(FieldsCntInh)];
semFields = [nanstd(FieldsCnt)./(sqrt(length(FieldsCnt)))...
    nanstd(FieldsCntInh)./(sqrt(length(FieldsCntInh)))];
figure; hold on
bar(mnFields,'k');
errorbar(1:length(mnFields), mnFields, semFields, '-k', 'LineStyle','none');
title('Number Fields')

[p1 h1] = ranksum(FieldsCnt,FieldsCntInh)

figure;
datacombinedGain = [FieldsCnt; FieldsCntInh];
g1 = repmat({'CA1exc'},length(FieldsCnt),1);
g2 = repmat({'CA1inh'},length(FieldsCntInh),1);
g = [g1;g2];

boxplot(datacombinedGain,g);
title(['NumFields-p = ' num2str(p1)])

%In/Out Spike Ratio
mnRatio = [nanmean(Ratio) nanmean(RatioInh)];
semRatio = [nanstd(Ratio)./(sqrt(length(Ratio)))...
    nanstd(RatioInh)./(sqrt(length(RatioInh)))];
figure; hold on
bar(mnRatio,'k');
errorbar(1:length(mnRatio), mnRatio, semRatio, '-k', 'LineStyle','none');
title('In-Out Spike Ratio')

%Field Width
mnWidth = [nanmean(FieldWidth) nanmean(FieldWidthInh)];
semWidth = [nanstd(FieldWidth)./(sqrt(length(FieldWidth)))...
    nanstd(FieldWidthInh)./(sqrt(length(FieldWidthInh)))];
figure; hold on
bar(mnWidth,'k');
errorbar(1:length(mnWidth), mnWidth, semWidth, '-k', 'LineStyle','none');
title('Field Width')

[p2 h2] = ranksum(FieldWidth,FieldWidthInh)

figure;
datacombinedGain = [FieldWidth; FieldWidthInh];
g1 = repmat({'CA1exc'},length(FieldWidth),1);
g2 = repmat({'CA1inh'},length(FieldWidthInh),1);
g = [g1;g2];

boxplot(datacombinedGain,g);
title(['Field Width-p = ' num2str(p2)])



%%
% %num fields
% mnFields = [nanmean(FieldsCntNM) nanmean(FieldsCntExc) nanmean(FieldsCntInh)];
% semFields = [nanstd(FieldsCntNM)./(sqrt(length(FieldsCntNM)))...
%     nanstd(FieldsCntExc)./(sqrt(length(FieldsCntExc)))...
%     nanstd(FieldsCntInh)./(sqrt(length(FieldsCntInh)))];
% figure; hold on
% bar(mnFields,'k');
% errorbar(1:length(mnFields), mnFields, semFields, '-k', 'LineStyle','none');
% title('Number Fields')
%
% %In/Out Spike Ratio
% mnRatio = [nanmean(RatioNM) nanmean(RatioExc) nanmean(RatioInh)];
% semRatio = [nanstd(RatioNM)./(sqrt(length(RatioNM)))...
%     nanstd(RatioExc)./(sqrt(length(RatioExc)))...
%     nanstd(RatioInh)./(sqrt(length(RatioInh)))];
% figure; hold on
% bar(mnRatio,'k');
% errorbar(1:length(mnRatio), mnRatio, semRatio, '-k', 'LineStyle','none');
% title('In-Out Spike Ratio')
%
% %Field Width
% mnWidth = [nanmean(FieldWidthNM) nanmean(FieldWidthExc) nanmean(FieldWidthInh)];
% semWidth = [nanstd(FieldWidthNM)./(sqrt(length(FieldWidthNM)))...
%     nanstd(FieldWidthExc)./(sqrt(length(FieldWidthExc)))...
%     nanstd(FieldWidthInh)./(sqrt(length(FieldWidthInh)))];
% figure; hold on
% bar(mnWidth,'k');
% errorbar(1:length(mnWidth), mnWidth, semWidth, '-k', 'LineStyle','none');
% title('Field Width')

keyboard
function  out = gaussSmooth(vector, binrange)

paddinglength = round(binrange*2.5);
padding = ones(paddinglength,1);

out = smoothvect([padding*vector(1) ; vector; padding*vector(end)],gaussian(binrange,binrange*5));
out = out(paddinglength+1:end-paddinglength);


