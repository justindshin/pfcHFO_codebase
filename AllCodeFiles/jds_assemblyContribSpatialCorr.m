%Get firing rate bias and percentage of ripple during which the cell is
%active (control for the cell being active in a different number of PFC vs
%CA1 ripples)
clear all;
close all;
%%
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};
day = 1;
fieldSimilarityHighContrib = [];
fieldSimilarityLowContrib = [];
for a = 1:length(animalprefixlist)
    animalprefix = char(animalprefixlist{a});
    animdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    load(sprintf('%s%sCA1_icareactivationtimes20SWSSpk0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%slinfields0%d.mat',animdir,animalprefix,day));
    
    for i = 1:length(icareactivationtimes)
        eps = icareactivationtimes{i}.epochs;
        tmpWeights = icareactivationtimes{i}.pc_weights;
        
        for as = 1:length(tmpWeights)
            currWeights = tmpWeights{as};
            currWeights(:,3) = currWeights(:,3);
            meanWeight = mean(currWeights(:,3));
            stdWeight = std(currWeights(:,3));
            highThresh = meanWeight + 2*stdWeight; %Threshold for high contributing cells is 2sd above mean
%             highThresh = graythresh(currWeights(:,3)); %otsu's method for thresholding
            highIdx = (currWeights(:,3) > highThresh);
            highCells = currWeights(highIdx,[1 2]);
            lowCells = currWeights(~highIdx,[1 2]);
            
            for ii = 1:length(highCells(:,1))
                lfield1 = linfields{day}{eps(1)}{highCells(ii,1)}{highCells(ii,2)};
                for j = 1:length(highCells(:,1))
                    lfield2 = linfields{day}{eps(1)}{highCells(j,1)}{highCells(j,2)};
                    if (ii~=j) && (ii<j)
                        fieldstmp1 = [lfield1{1}(:,5); lfield1{2}(:,5); lfield1{3}(:,5); lfield1{4}(:,5)];
                        fieldstmp2 = [lfield2{1}(:,5); lfield2{2}(:,5); lfield2{3}(:,5); lfield2{4}(:,5)];
                        rval = corrcoef(fieldstmp1,fieldstmp2,'rows','complete');
                        rval = rval(1,2);
                        fieldSimilarityHighContrib = [fieldSimilarityHighContrib; rval];
                    end
                end
            end
            
            for ii = 1:length(lowCells(:,1))
                lfield1 = linfields{day}{eps(1)}{lowCells(ii,1)}{lowCells(ii,2)};
                for j = 1:length(lowCells(:,1))
                    lfield2 = linfields{day}{eps(1)}{lowCells(j,1)}{lowCells(j,2)};
                    if (ii~=j) && (ii<j)
                        fieldstmp1 = [lfield1{1}(:,5); lfield1{2}(:,5); lfield1{3}(:,5); lfield1{4}(:,5)];
                        fieldstmp2 = [lfield2{1}(:,5); lfield2{2}(:,5); lfield2{3}(:,5); lfield2{4}(:,5)];
                        rval = corrcoef(fieldstmp1,fieldstmp2,'rows','complete');
                        rval = rval(1,2);
                        fieldSimilarityLowContrib = [fieldSimilarityLowContrib; rval];
                    end
                end
            end
        end
    end
end

figure; bar(1,[nanmean(fieldSimilarityLowContrib)],'w'); hold on
bar(2,[nanmean(fieldSimilarityHighContrib)],'k'); hold on

semLow = nanstd(fieldSimilarityLowContrib)./sqrt(length(find(~isnan(fieldSimilarityLowContrib))))
semHigh = nanstd(fieldSimilarityHighContrib)./sqrt(length(find(~isnan(fieldSimilarityHighContrib))))

er = errorbar([1],[nanmean(fieldSimilarityLowContrib)],[semLow]);
er.Color = 'k';
er.LineStyle = 'none';
er.LineWidth = 1;
er2 = errorbar([2],[nanmean(fieldSimilarityHighContrib)],[semHigh]);
er2.Color = 'k';
er2.LineStyle = 'none';
er2.LineWidth = 2;

ax = gca
ax.FontSize = 16
xticks([1 2])
xlabel('Assembly Contribution')
xticklabels({'Other Pairs','Member Pairs'})
ylabel('Field Similarity (r)')
keyboard