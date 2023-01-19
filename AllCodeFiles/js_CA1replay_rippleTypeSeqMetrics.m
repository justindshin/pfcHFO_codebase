function js_CA1replay_rippleTypeSeqMetrics(animalprefixlist)

day = 1;

bins = 20;

propSig_noncoord = [];
propSig_coord = [];
pval_noncoord = [];
pval_coord = [];
r2_noncoord = [];
r2_coord = [];
jdist_noncoord = [];
jdist_coord = [];
wc_noncoord = [];
wc_coord = [];
proprev_noncoord = [];
proprev_coord = [];

for a = 1:length(animalprefixlist)
    animalprefix = char(animalprefixlist(a));
    epochs = [3:2:17];
    %
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    for ep=1:length(epochs)
        cntrevnoncoord = 0;
        cntrevcoord = 0;
        epoch = epochs(ep);

        load(sprintf('%s%sreplaydecode0%d_%02d.mat',dir,animalprefix,day,epoch));

        sigreplayidx = replaytraj.besttraj;
        riptype = replaytraj.ripType;
        signoncoordIdx1 = find(sigreplayidx ~= 0);
        signoncoordIdx2 = find(riptype == 1);
        signoncoordIdx3 = intersect(signoncoordIdx1,signoncoordIdx2);
        probnoncoordtmp = length(signoncoordIdx3)/length(signoncoordIdx2);

        for i = 1:length(signoncoordIdx3)
            pmat1 = replaytraj.pMat{signoncoordIdx3(i)}{sigreplayidx(signoncoordIdx3(i))}.pMat;
            pmatlength1 = length(replaytraj.pMat{signoncoordIdx3(i)}{sigreplayidx(signoncoordIdx3(i))}.pMat)*2;
            pv1 = replaytraj.pvalue(signoncoordIdx3(i),sigreplayidx(signoncoordIdx3(i)));
            slp1 = replaytraj.slopes(signoncoordIdx3(i),sigreplayidx(signoncoordIdx3(i)));
            if slp1 < 0
                cntrevnoncoord = cntrevnoncoord + 1;
            end
            rs1 = replaytraj.rsquare(signoncoordIdx3(i),sigreplayidx(signoncoordIdx3(i)));
            jDistBins = [];
            for j = 1:size(pmat1,2)-1
                tmp = pmat1(:,j);
                tmp2 = pmat1(:,j+1);
                if (sum(tmp) ~= 0) && (sum(tmp2) ~= 0)
                    [M1,I1] = max(tmp);
                    [M2,I2] = max(tmp2);
                    jDist = abs(I1-I2)*2;
                    jDistBins = [jDistBins; jDist];
                end
            end

            positionVec = 1:pmatlength1/2;
            timevec = 1:size(pmat1,2);
            mloc = sum(sum(pmat1.*positionVec'))/sum(sum(pmat1));
            mt = sum(sum(pmat1.*timevec))/sum(sum(pmat1));
            dloc = positionVec'-mloc;
            dt = timevec - mt;
            cov_loc_t = sum(sum(pmat1.*dloc.*dt))/sum(sum(pmat1));
            cov_loc = sum(sum(pmat1.*(dloc.^2)))/sum(sum(pmat1));
            cov_t = sum(sum(pmat1.*(dt.^2)))/sum(sum(pmat1));
            weighted_corr = cov_loc_t/sqrt(cov_loc*cov_t);

            maxjump = max(jDistBins);
            jd1 = maxjump/pmatlength1;
            jdist_noncoord = [jdist_noncoord; jd1];
            pval_noncoord = [pval_noncoord; pv1];
            wc_noncoord = [wc_noncoord; abs(weighted_corr)];
            r2_noncoord = [r2_noncoord; rs1];
        end

        sigcoordIdx1 = find(sigreplayidx ~= 0);
        sigcoordIdx2 = find(riptype == 2);
        sigcoordIdx3 = intersect(sigcoordIdx1,sigcoordIdx2);
        probcoordtmp = length(sigcoordIdx3)/length(sigcoordIdx2);
        if (length(signoncoordIdx3) >= 10)
            proprev_noncoord = [proprev_noncoord; cntrevnoncoord/length(signoncoordIdx3)];
        end

        for ii = 1:length(sigcoordIdx3)
            pmat2 = replaytraj.pMat{sigcoordIdx3(ii)}{sigreplayidx(sigcoordIdx3(ii))}.pMat;
            pmatlength2 = length(replaytraj.pMat{sigcoordIdx3(ii)}{sigreplayidx(sigcoordIdx3(ii))}.pMat)*2;
            pv2 = replaytraj.pvalue(sigcoordIdx3(ii),sigreplayidx(sigcoordIdx3(ii)));
            slp2 = replaytraj.slopes(sigcoordIdx3(ii),sigreplayidx(sigcoordIdx3(ii)));
            if slp2 < 0
                cntrevcoord = cntrevcoord + 1;
            end
            rs2 = replaytraj.rsquare(sigcoordIdx3(ii),sigreplayidx(sigcoordIdx3(ii)));
            jDistBins = [];
            for j = 1:size(pmat2,2)-1
                tmp = pmat2(:,j);
                tmp2 = pmat2(:,j+1);
                if (sum(tmp) ~= 0) && (sum(tmp2) ~= 0)
                    [M1,I1] = max(tmp);
                    [M2,I2] = max(tmp2);
                    jDist = abs(I1-I2)*2;
                    jDistBins = [jDistBins; jDist];
                end
            end

            %weighted correlation calculation
            positionVec = 1:pmatlength2/2;
            timevec = 1:size(pmat2,2);
            mloc = sum(sum(pmat2.*positionVec'))/sum(sum(pmat2));
            mt = sum(sum(pmat2.*timevec))/sum(sum(pmat2));
            dloc = positionVec'-mloc;
            dt = timevec - mt;
            cov_loc_t = sum(sum(pmat2.*dloc.*dt))/sum(sum(pmat2));
            cov_loc = sum(sum(pmat2.*(dloc.^2)))/sum(sum(pmat2));
            cov_t = sum(sum(pmat2.*(dt.^2)))/sum(sum(pmat2));
            weighted_corr = cov_loc_t/sqrt(cov_loc*cov_t);

            maxjump = max(jDistBins); %maximum jump distance for event
            jd2 = maxjump/pmatlength2; %normalized max jump distance
            jdist_coord = [jdist_coord; jd2];
            pval_coord = [pval_coord; pv2];
            wc_coord = [wc_coord; abs(weighted_corr)]; %abs WC
            r2_coord = [r2_coord; rs2];
        end
        if (length(sigcoordIdx3) >= 10)
            proprev_coord = [proprev_coord; cntrevcoord/length(sigcoordIdx3)];
        end
        if (length(signoncoordIdx2) >= 5) && (length(sigcoordIdx2) >= 5)
            propSig_noncoord = [propSig_noncoord; probnoncoordtmp];
            propSig_coord = [propSig_coord; probcoordtmp];
        end
    end
end


figure;
meanJdist = [mean(jdist_noncoord) mean(jdist_coord)]
medianJdist = [median(jdist_noncoord) median(jdist_coord)]
semJdist = [std(jdist_noncoord)./sqrt(length(jdist_noncoord)) std(jdist_coord)./sqrt(length(jdist_coord))]
datacombinedJdist = [jdist_noncoord; jdist_coord];
g1 = repmat({'Ind'},length(jdist_noncoord),1);
g2 = repmat({'Coord'},length(jdist_coord),1);
g = [g1;g2];

boxplot(datacombinedJdist,g);
title('Norm. Max Jump Distance')

figure;
meanWC = [nanmean(wc_noncoord) nanmean(wc_coord)]
medianWC = [nanmedian(wc_noncoord) nanmedian(wc_coord)]
semWC = [nanstd(wc_noncoord)./sqrt(length(find(~isnan(wc_noncoord))))...
    nanstd(wc_coord)./sqrt(length(find(~isnan(wc_coord))))]
datacombinedWC = [wc_noncoord; wc_coord];
g1 = repmat({'Ind'},length(wc_noncoord),1);
g2 = repmat({'Coord'},length(wc_coord),1);
g = [g1;g2];

boxplot(datacombinedWC,g);
title('Weighted Correlation')

figure;
datacombinedProp = [propSig_noncoord; propSig_coord];
g1 = repmat({'Ind'},length(propSig_noncoord),1);
g2 = repmat({'Coord'},length(propSig_coord),1);
g = [g1;g2];

boxplot(datacombinedProp,g);
title('Proportion Ripple Type Significant/Epoch')
hold on
for e = 1:length(propSig_noncoord)
    x = [1 2];
    plot(x,[propSig_noncoord(e) propSig_coord(e)],'k')
end

figure;
datacombinedr2 = [r2_noncoord; r2_coord];
g1 = repmat({'Ind'},length(r2_noncoord),1);
g2 = repmat({'Coord'},length(r2_coord),1);
g = [g1;g2];
[pR2 hR2] = ranksum(r2_noncoord,r2_coord);
h = boxplot(datacombinedr2,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
boxplot(datacombinedr2,g);
title(['R-squared p=' num2str(pR2)])
hold on

figure;
bar([nanmean(r2_noncoord) nanmean(r2_coord)],'k')
hold on
errorbar([1:2],[nanmean(r2_noncoord) nanmean(r2_coord)],...
    [(nanstd(r2_noncoord)./sqrt(length(r2_noncoord)))...
    (nanstd(r2_coord)./sqrt(length(r2_coord)))],'k','LineStyle','none','LineWidth',2)
title(['R-squared p=' num2str(pR2)])
xticklabels({'Ind','Coord'})

%%

figure;
bar([mean(proprev_noncoord) mean(proprev_coord)],'k')
hold on
errorbar([1:2],[mean(proprev_noncoord) mean(proprev_coord)],...
    [(std(proprev_noncoord)./sqrt(length(proprev_noncoord)))...
    (std(proprev_coord)./sqrt(length(proprev_coord)))],'k','LineStyle','none','LineWidth',2)
[h1 p1] = ttest2(proprev_noncoord,proprev_coord)
title(['Proportion Reverse Replay p=' num2str(p1)])
xticklabels({'Ind','Coord'})

figure;
bar([mean(propSig_noncoord) mean(propSig_coord)],'k')
hold on
errorbar([1:2],[mean(propSig_noncoord) mean(propSig_coord)],...
    [(std(propSig_noncoord)./sqrt(length(propSig_noncoord)))...
    (std(propSig_coord)./sqrt(length(propSig_coord)))],'k','LineStyle','none','LineWidth',2)
[h2 p2] = ttest2(propSig_noncoord,propSig_coord)
title(['Proportion Significant Replay p=' num2str(p2)])
xticklabels({'Ind','Coord'})

% figure
% plot(jd_idx_noncoord./sum(jd_idx_noncoord)); hold on;
% plot(jd_idx_coord./sum(jd_idx_coord))
% ylabel('Jump to location probability')
% xlabel('Position bin')

keyboard