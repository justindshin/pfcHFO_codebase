function jds_noncoordrippletriggered_assemblystrength_M(animalprefixlist,area,state)

%Assembly strength triggered on independent ripples events in CA1 and PFC
day = 1;

% bins = 400; %for 5 ms
% bins = 20; %for 100 ms
bins = 100; %for 20ms - CURRENT ASSEMBLY BINSIZE
peakbins = find(abs(-bins:bins)<=10);

allevents_ctxriptrig = [];
allevents_ca1riptrig = [];
allevents_ctxripstrength = [];
allevents_ca1ripstrength = [];
shuf = 0;

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrength%sNewSpkRips_20_%02d.mat',dir,animalprefix,area,state,day));

    %Load ripples
    load(sprintf('%s%srippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_noncoordSWS%02d.mat',dir,animalprefix,day));

    %Which epochs to analyze
    epochs = find(~cellfun(@isempty,RtimeStrength)); %use epochs where there are assemblies

    for e = 1:length(epochs)
        ep = epochs(e);
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
        ctxripstarts = ctxripple{day}{ep}.starttime;
        ripstarts = ripple{day}{ep}.starttime;

        %Do PFC ripples
        if ~isempty(assemblytmp)
            for ii = 1:length(assemblytmp)
                react_idx = [];
                for t = 1:length(ctxripstarts)
                    idxtmp = lookup(ctxripstarts(t), assemblytmp{ii}(:,1));
                    react_idx = [react_idx; idxtmp];
                end
                atmp = [];
                strengthstmp = assemblytmp{ii}(:,2);
                if shuf == 1
                    strengthstmp = strengthstmp(randperm(length(strengthstmp)));
                end
                for r = 1:length(react_idx)
                    if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                        tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                        atmp = [atmp; tmp'];
                    end
                end
                allevents_ctxriptrig = [allevents_ctxriptrig; zscore(mean(atmp,1))]; %mean and zscore for each assembly
            end
            %Do CA1 ripples
            for ii = 1:length(assemblytmp)
                react_idx = [];
                for t = 1:length(ripstarts)
                    idxtmp = lookup(ripstarts(t), assemblytmp{ii}(:,1));
                    react_idx = [react_idx; idxtmp];
                end
                atmp = [];
                strengthstmp = assemblytmp{ii}(:,2);
                if shuf == 1
                    strengthstmp = strengthstmp(randperm(length(strengthstmp)));
                end
                for r = 1:length(react_idx)
                    if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                        tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                        atmp = [atmp; tmp'];
                    end
                end
                allevents_ca1riptrig = [allevents_ca1riptrig; zscore(mean(atmp,1))];
            end
        end
    end
end

allevents_ctxriptrig_z_mn = mean(allevents_ctxriptrig);
allevents_ca1riptrig_z_mn = mean(allevents_ca1riptrig);

allevents_ctxriptrig_z_sem = (std(allevents_ctxriptrig)./sqrt(length(allevents_ctxriptrig(:,1))));
allevents_ca1riptrig_z_sem = (std(allevents_ca1riptrig)./sqrt(length(allevents_ca1riptrig(:,1))));

mnPeakCtxRip = mean(allevents_ctxriptrig(:,peakbins),2);
mnPeakCa1Rip = mean(allevents_ca1riptrig(:,peakbins),2);

rDiff = mnPeakCtxRip-mnPeakCa1Rip;
stem(rDiff); view(90,90)
g1 = gaussian(3, 10);

%% EXAMPLES
%CA1
subplot(2,1,1)
hold on
plot(smoothvect(allevents_ca1riptrig(102,:),g1))
xlim([50 150])
xticks([50 100 150])
subplot(2,1,2)
hold on
plot(smoothvect(allevents_ca1riptrig(103,:),g1))
xlim([50 150])
xticks([50 100 150])
subplot(2,1,1)
plot(smoothvect(allevents_ctxriptrig(102,:),g1))
subplot(2,1,2)
plot(smoothvect(allevents_ctxriptrig(103,:),g1))

%PFC
% subplot(2,1,1)
% hold on
% plot(smoothvect(allevents_ca1riptrig(9,:),g1))
% xlim([50 150])
% xticks([50 100 150])
% subplot(2,1,2)
% hold on
% plot(smoothvect(allevents_ca1riptrig(14,:),g1))
% xlim([50 150])
% xticks([50 100 150])
% subplot(2,1,1)
% plot(smoothvect(allevents_ctxriptrig(9,:),g1))
% subplot(2,1,2)
% plot(smoothvect(allevents_ctxriptrig(14,:),g1))
% set(gcf, 'renderer', 'painters')
%%

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_ca1riptrig_z_mn,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_ca1riptrig_z_mn,allevents_ca1riptrig_z_sem,'-k');
pl2 = plot([-bins:bins],allevents_ctxriptrig_z_mn,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_ctxriptrig_z_mn,allevents_ctxriptrig_z_sem,'-r');

title(sprintf('NC Ripple Triggered %s Reactivation Strength',area))

ylabel('Reactivation Strength')
xlabel('Time from Ripple Onset (s)')
xlim([(-bins) bins]);
% xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
xticklabels({'-2','-1','0','1','2'});

xx = [0 0];
plot(xx,ax1.YLim,'--k')
ylim = (ax1.YLim);
legend([pl1 pl2],{'CA1 Ripples','PFC Ripples'})
set(gcf, 'renderer', 'painters')

figure;
imagesc(allevents_ctxriptrig(randperm(length(allevents_ctxriptrig(:,1))),:))
title([sprintf('NC Ripple Triggered %s Reactivation Strength',area) '- NonCoord PFC Rips'])
colorbar
set(gcf, 'renderer', 'painters')
clim([-3 3])
xticks([1 101 201])
xticklabels({'-2','0','2'})
xlabel('Time from Ripple onset')

figure
imagesc(allevents_ca1riptrig(randperm(length(allevents_ca1riptrig(:,1))),:))
title([sprintf('NC Ripple Triggered %s Reactivation Strength',area) '- NonCoord CA1 Rips'])
colorbar
set(gcf, 'renderer', 'painters')
clim([-3 3])
xticks([1 101 201])
xticklabels({'-2','0','2'})
xlabel('Time from Ripple onset (s)')

keyboard
