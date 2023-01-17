animalprefixlist = {'ZT2','JS34','JS17','JS21','JS15','JS14','ER1','ER1_2','KL8'};
day = 1;
wavs_p = [];
wavs_i = [];
maxAmp = [];
lRatio = [];
isoDist = [];
pktotr_p = [];
meanFR_p = [];
burstIdx_p = [];
pktotr_i = [];
meanFR_i = [];
burstIdx_i = [];
logIsi_p = [];
logIsi_i = [];
peakVoltage = [];
intervals = -3:0.04:1;
intervals2 = intervals(1:end-1)+.02;

plotisi = 0;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    doneCells = [];
    dir = sprintf('/Volumes/OneTouch/MatclustFiles/%s/',animalprefix);
    if isequal(animalprefix,'ER1_2')
        animalprefix = 'ER1';
    end
    load(sprintf('%s%sspikesInfo%02d.mat',dir,animalprefix,day));
    for i = 1:length(spikesInfo{day})
        for ii = 1:length(spikesInfo{day}{i})
            if ~isempty(spikesInfo{day}{i})
                for l = 1:length(spikesInfo{day}{i}{ii})
                    if ~isempty(spikesInfo{day}{i}{ii})
                        if ~isempty(spikesInfo{day}{i}{ii}{l})
                            isi = {};
                            spks = spikesInfo{day}{i}{ii}{l}.data;
                            isi.log10 = zeros(length(intervals2),size(spks, 2));
                            isi.log10_bins = 10.^intervals2';
                            for j = 1:size(spks  ,2)
                                ISIs = log10(diff(spks));
                                [ISIlog,~] = histcounts(ISIs,intervals);
                                isi.log10(:,j) = ISIlog./(diff(10.^intervals))/length(spks);
                            end

                            figure(1); hold on, set(gca,'xscale','log');
                            ax1 = gca;
                            xlabel(ax1,'Time (s)'), ylabel(ax1,'Occurrence'), title(ax1,'Log ISI distribution')
                            avwv = -1.*normalize(spikesInfo{day}{i}{ii}{l}.avgWav,'range')';
                            %   burst = spikesInfo{day}{i}{ii}{l}.burstindex/spikesInfo{day}{i}{ii}{l}.numSpks;
                            burst = spikesInfo{day}{i}{ii}{l}.burstprobability;
                            rate = spikesInfo{day}{i}{ii}{l}.meanrate;
                            spkwidth = abs(spikesInfo{day}{i}{ii}{l}.p2t);
                            peakVoltage = [peakVoltage; max(spikesInfo{day}{i}{ii}{l}.avgWav)];
                            if ~ismember([ii l],doneCells,'rows','legacy')
                                doneCells = [doneCells; [ii l]];
                                lr = spikesInfo{day}{i}{ii}{l}.Lratio;
                                isoD = spikesInfo{day}{i}{ii}{l}.IsolationDistance;
                                lRatio = [lRatio; lr];
                                isoDist = [isoDist; isoD];
                            end
                            maxAmp = [maxAmp; max(abs(spikesInfo{1, 1}{i}{ii}{l}.avgWav))];
                            if rate < 7
                                logIsi_p = [logIsi_p isi.log10.*(10.^intervals2)'];
                                wavs_p = [wavs_p; avwv];
                                burstIdx_p = [burstIdx_p; burst];
                                meanFR_p = [meanFR_p; rate];
                                pktotr_p = [pktotr_p; spkwidth];
                            elseif rate > 7
                                logIsi_i = [logIsi_i isi.log10.*(10.^intervals2)'];
                                wavs_i = [wavs_i; avwv];
                                burstIdx_i = [burstIdx_i; burst];
                                meanFR_i = [meanFR_i; rate];
                                pktotr_i = [pktotr_i; spkwidth];
                            end
                            if plotisi %example ISI dist
                                if rate > 7
                                    ISI = diff(spks);
                                    ISI = ISI(find(ISI<.5));
                                    ISI = [ISI;-ISI];
                                    edges = [-.3:.0005:.3];
                                    N = histc(ISI,edges);

                                    if (sum(N))
                                        figure
                                        axis([-.3 .3 0 max(N)]);
                                        phandle = bar(edges,N,'histc');
                                        set(phandle,'LineStyle','none');
                                        set(phandle,'FaceColor',[0 0 0]);
                                        xlabel('Time between events (seconds)');
                                        ylabel('Number of events');
                                        if rate > 7
                                            title('CA1Int example ISI')
                                        elseif rate < 7
                                            title('CA1Pyr example ISI')
                                        end
                                    else
                                        axis([-.3 .3 0 1]);
                                    end
                                    keyboard
                                    close(figure(2))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
figure
scatter3(meanFR_p,pktotr_p,burstIdx_p,'ro')
hold on
scatter3(meanFR_i,pktotr_i,burstIdx_i,'bo')
set(gca, 'XScale', 'log')
set(gca, 'ZScale', 'log')
xlabel('Mean firing rate (Hz)')
ylabel('Trough-to-peak (ms)')
zlabel('Burst index')
set(gcf, 'renderer', 'painters')


figure
boundedline(1:40,nanmean(wavs_p),(nanstd(wavs_p)./sqrt(size(wavs_p,1))),'k')
boundedline(1:40,nanmean(wavs_i),(nanstd(wavs_i)./sqrt(size(wavs_i,1))),'r')
xlim([5 40])
xticks([5 14 23 32])
xticklabels({'-0.3','0','0.3','0.6'})
yticks([-1 0])
ylabel('Voltage (normalized)')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

datacombinedisoDist = [isoDist];
g1 = repmat({'Isolation distance'},length(isoDist),1);
g = [g1];

figure;
h = boxplot(datacombinedisoDist,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title('CA1 Isolation distance')
set(gcf, 'renderer', 'painters')
ylim([-1 75])
ylabel('Isolation distance')
set(gcf, 'renderer', 'painters')


datacombinedPeakV = [peakVoltage];
g1 = repmat({'Peak Voltage (uV)'},length(peakVoltage),1);
g = [g1];

figure;
h = boxplot(datacombinedPeakV,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title('CA1 Peak Voltage')
set(gcf, 'renderer', 'painters')
ylim([0 550])
yticks([0:1:4])
ylabel('Voltage (uV)')

logIsi_p = logIsi_p';
logIsi_i = logIsi_i';

figure(1)
plot(ax1,isi.log10_bins,mean(logIsi_p),'r','LineWidth',3)
hold on
plot(ax1,isi.log10_bins,mean(logIsi_i),'b','LineWidth',3)
legend({'Pyramidal','Interneuron'})
set(gcf, 'renderer', 'painters')
%
figure
boundedline(1:100,nanmean(logIsi_p),(nanstd(logIsi_p)./sqrt(size(logIsi_p,1))),'r')
boundedline(1:100,nanmean(logIsi_i),(nanstd(logIsi_i)./sqrt(size(logIsi_i,1))),'b')

