function jds_rippletriggered_compareStrength_M(animalprefixlist)

%To do: plot the peak strength within specified ripple events. compare
%ripple types for different assemblies - might need to change binsize for
%raster generation
day = 1;

b = gaussian(2,6);
bins = 400;
peakbins = find(abs(-bins:bins)<=4);

coordriptrig = [];
noncoordriptrig = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
%     load(sprintf('%s%sCA1PFC_RTimeStrengthSleepNewSpk_50_%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    %Load ripples
%     load(sprintf('%s%srippletime_coordSWS%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%srippletime_leadlag%02d.mat',dir,animalprefix,day));
    %     load(sprintf('%s%srippleframes_SWS%02d.mat',dir,animalprefix,day));
    %     noncoordripple = ripple;
    coordripple = ripplecoupling; clear ripple;
%     load(sprintf('%s%srippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
    %     load(sprintf('%s%srippletime_SWS%02d.mat',dir,animalprefix,day));
    noncoordripple = ripplecoupling;
    %     coordripple = ctxripple;

    %Which epochs to analyze
    epochs = find(~cellfun(@isempty,RtimeStrength));
    %     epochs = [1];
    %         if isequal(animalprefix,'JS34')
    %             epochs = [5 11 13 15 17];
    %         end

    for e = 1:length(epochs)
        ep = epochs(e);
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
        noncoordripstarts = noncoordripple{day}{ep}.ctxleadriptimes;
        coordripstarts = coordripple{day}{ep}.ctxlagriptimes;
        if (length(coordripstarts) >= 10) && (length(noncoordripstarts) >= 10)
            if ~isempty(assemblytmp)
                for ii = 1:length(assemblytmp)
%                     if RtimeStrength{ep}.members{ii}.cross == 1
                        react_idx = [];
                        for t = 1:length(coordripstarts)
                            idxtmp = lookup(coordripstarts(t), assemblytmp{ii}(:,1));
                            react_idx = [react_idx; idxtmp];
                        end
                        atmp = [];
                        for r = 1:length(react_idx)
                            strengthstmp = assemblytmp{ii}(:,2);
                            if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                                tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                                atmp = [atmp; tmp'];
                            end
                        end
                        react_z = smoothvect(zscore(mean(atmp)),b);
%                         react_z = react_z(:,peakbins);
%                         react_z = mean(react_z')';
                        coordriptrig = [coordriptrig;react_z];

                        react_idx = [];
                        for t = 1:length(noncoordripstarts)
                            idxtmp = lookup(noncoordripstarts(t), assemblytmp{ii}(:,1));
                            react_idx = [react_idx; idxtmp];
                        end
                        atmp = [];
                        for r = 1:length(react_idx)
                            strengthstmp = assemblytmp{ii}(:,2);
                            if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                                tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                                atmp = [atmp; tmp'];
                            end
                        end
                        react_z = smoothvect(zscore(mean(atmp)),b);
%                         react_z = react_z(:,peakbins);
%                         react_z = mean(react_z')';
                        noncoordriptrig = [noncoordriptrig;react_z];
%                     end
                end
            end
        end
    end
end

allevents_coordriptrig_bins = coordriptrig(:,peakbins);
allevents_noncoordriptrig_bins = noncoordriptrig(:,peakbins);

allevents_coordriptrig_bins = mean(allevents_coordriptrig_bins')'
allevents_noncoordriptrig_bins = mean(allevents_noncoordriptrig_bins')';

[p h] = signrank(allevents_coordriptrig_bins,allevents_noncoordriptrig_bins)

allevents_coordriptrig_mean = mean(coordriptrig);
allevents_noncoordriptrig_mean = mean(noncoordriptrig);

allevents_coordriptrig_sem = std(coordriptrig)./sqrt(length(coordriptrig(:,1)));
allevents_noncoordriptrig_sem = std(noncoordriptrig)./sqrt(length(noncoordriptrig(:,1)));

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_coordriptrig_mean,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_coordriptrig_mean,allevents_coordriptrig_sem,'-k');
pl2 = plot([-bins:bins],allevents_noncoordriptrig_mean,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_noncoordriptrig_mean,allevents_noncoordriptrig_sem,'-r');
legend([pl1 pl2],{'Lag','Lead'})

title('CA1PFC Peak Reactivation')
ylabel('NonCoordinated CA1 ripples')
xlabel('Coordinated CA1 ripples')
% title(['CA1PFC Peak Reactivation - p' num2str(p(1,2)) '- r' num2str(r(1,2))])
set(gcf, 'renderer', 'painters')

%%
if bins == 400 %using 200ms after ripple for comparison
    [tmp timingindCA1]=sort(mean(coordriptrig(:,401:440),2),'descend');
    figure;
    imagesc([-400:400],1:length(timingindCA1),coordriptrig(timingindCA1,:));
    xlabel('Time (ms)')
    ylabel('#Cell')
    xlim([-1000 1000])
    meanreactcoord = mean(coordriptrig(timingindCA1,401:440),2);

    meanreactnoncoord = mean(noncoordriptrig(timingindCA1,401:440),2);
    figure; scatter(meanreactcoord,meanreactnoncoord,'.k'); hold on; lsline
end

keyboard
