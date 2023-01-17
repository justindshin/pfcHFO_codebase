function jds_sequenceScoreIndCoord(animalprefixlist)

day = 1;

rZcoord = [];
rZind = [];
wccoord = [];
wcind = [];
jdcoord = [];
jdind = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    
    epochs = [3:2:17];
   
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    
    for ep=1:length(epochs)
        epoch = epochs(ep);
        
        load(sprintf('%s%sPCCreplaydecode_alldata_50msRipNoIriCrit_sequenceScores0%d_%02d.mat',dir,animalprefix,day,epoch));
        
        pMats = PCC_events.pMats;
        %Filter for low pvalue
        events_sig = [];
        for event = 1:length(pMats)
            tmp_pmat = pMats{event};
            rtype = tmp_pmat.rippletype;
            traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
            if length(traj) > 1
                continue
            end
            if min(tmp_pmat.traj_pvals) < 0.05
%                 if tmp_pmat.weightedcorr(traj) < 0.5
%                     continue
%                 end
                pmat1 = PCC_events.pMats{event}.alltrajmats{traj};
                pmatlength1 = length(pmat1)*2;
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
                maxjump = max(jDistBins);
                jd1 = maxjump/pmatlength1;
                if rtype == 1
                    wcval = tmp_pmat.weightedcorr(traj);
                    seqScore = tmp_pmat.sequenceScore(traj);
                    wcind = [wcind; abs(wcval)];
                    rZind = [rZind; seqScore];
                    jdind = [jdind; jd1];
                elseif rtype == 2
                    wcval = tmp_pmat.weightedcorr(traj);
                    seqScore = tmp_pmat.sequenceScore(traj);
                    wccoord = [wccoord; abs(wcval)];
                    rZcoord = [rZcoord; seqScore];
                    jdcoord = [jdcoord; jd1];
                end
            end
        end
    end
end

[pWc, ~] = ranksum(wcind,wccoord)
[pSeq, ~] = ranksum(rZind,rZcoord)
[pJd, ~] = ranksum(jdind,jdcoord)

plotcdf(wcind); hold on;
plotcdf(wccoord);
ylabel('Cumulative proportion')
xlabel('Weighted correlation')
title('WC vs shuffle')
set(gcf, 'renderer', 'painters')

figure
plotcdf(jdind); hold on;
plotcdf(jdcoord);
ylabel('Cumulative proportion')
xlabel('Normalized max jump distance')
title('WC vs shuffle')
set(gcf, 'renderer', 'painters')

figure
datacombinedJdist = [jdind; jdcoord];
g1 = repmat({'Ind'},length(jdind),1);
g2 = repmat({'Coord'},length(jdcoord),1);
g = [g1;g2];

boxplot(datacombinedJdist,g);
title('Norm. Max Jump Distance')

figure;
datacombinedWC = [wcind; wccoord];
g1 = repmat({'Ind'},length(wcind),1);
g2 = repmat({'Coord'},length(wccoord),1);
g = [g1;g2];

boxplot(datacombinedWC,g);
title('Weighted Correlation')
keyboard




