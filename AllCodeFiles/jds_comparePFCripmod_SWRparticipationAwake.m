%Get firing rate bias and percentage of ripple during which the cell is
%active (control for the cell being active in a different number of PFC vs
%CA1 ripples)
clear all;
close all;
%%
savedirX = '/Volumes/JUSTIN/SingleDay/ProcessedData/';

animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};
day = 1;
epochs = [2:2:16];

%%
% epochs = [1:2:17];
exc_SWRfiring = [];
inh_SWRfiring = [];
nonsig_SWRfiring = [];
inh_FRtheta = [];
exc_FRtheta = [];
inh_FRgain = [];
exc_FRgain = [];
nonsig_FRgain = [];
nonsig_FRtheta = [];
inh_FRsws = [];
exc_FRsws = [];
nonsig_FRsws = [];
inh_FRrip = [];
exc_FRrip = [];
nonsig_FRrip = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));% get ripple time
    load(sprintf('%s%sspikes0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sthetatime0%d.mat', dir,animalprefix,day));
    
    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

    for e = 1:length(epochs)
        epoch = epochs(e);
        if epoch == 16
            ep_sleep = epoch + 1;
        else
            ep_sleep = epoch - 1;
        end
        
        ep2 = find(sleeps(:,2) == ep_sleep);
        modcells = epochModulation.cellidx;
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = 1;
        
        allmodcells = [inhcells; exccells];
        
        cellidx = allmodcells(:,[1 2]);
        riptimes = [ripple{day}{epoch}.starttime ripple{day}{epoch}.endtime];
        thetalist = [thetatime{day}{epoch}.starttime thetatime{day}{epoch}.endtime];
        thetadur = thetatime{day}{epoch}.total_duration;
        ripdur = sum(riptimes(:,2) - riptimes(:,1));
        for cell = 1:length(allmodcells(:,1))
            if ~isempty(riptimes)
                index = [day,epoch,cellidx(cell,:)];
                cellmod = allmodcells(cell,3);
                if (length(riptimes(:,1)) > 1) && (length(thetalist(:,1)) > 1)
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                        if length(spiketimes) > 100
                           
                            fRateMean = spikes{index(1)}{index(2)}{index(3)}{index(4)}.meanrate;
                            if fRateMean < 7
                                spikebins = periodAssign(spiketimes, riptimes);
                                
                                thetabins = periodAssign(spiketimes, thetalist);
                                inThetaSpk = length(find(thetabins ~= 0));
                                fRateTheta = inThetaSpk/thetadur;
                                
                                numSpkRip = length(find(spikebins ~= 0));
                                ripFR = numSpkRip/ripdur;
                                
                                inRipSpk = length(find(unique(spikebins) ~= 0));
                                pctRipActivity = (inRipSpk/(length(riptimes(:,1))));
                                
                                if cellmod == 1
                                    exc_SWRfiring = [exc_SWRfiring; pctRipActivity];
                                    exc_FRtheta = [exc_FRtheta; fRateTheta];
                                    exc_FRrip = [exc_FRrip; ripFR];
                                    exc_FRgain = [exc_FRgain; (ripFR/fRateTheta)];
                                elseif cellmod == -1
                                    inh_SWRfiring = [inh_SWRfiring; pctRipActivity];
                                    inh_FRtheta = [inh_FRtheta; fRateTheta];
                                    inh_FRrip = [inh_FRrip; ripFR];
                                    inh_FRgain = [inh_FRgain; (ripFR/fRateTheta)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

datameans = [mean(exc_SWRfiring) mean(inh_SWRfiring)];
datasems = [(std(exc_SWRfiring)/sqrt(length(exc_SWRfiring)))...
    (std(inh_SWRfiring)/sqrt(length(inh_SWRfiring)))];

bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Percent CA1 SWRs Active')
title('PFC Ripple Modulation and Awake CA1 SWR Participation')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p h] = ranksum(exc_SWRfiring,inh_SWRfiring)

datacombinedthetaSWRparticipation = [exc_SWRfiring; inh_SWRfiring];
g1 = repmat({'CA1exc'},length(exc_SWRfiring),1);
g2 = repmat({'CA1inh'},length(inh_SWRfiring),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedthetaSWRparticipation,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Awake SWR participation-p = ' num2str(p)])
ylabel('Proportion ripples active')
set(gcf, 'renderer', 'painters')
%%
datameansFR = [mean(exc_FRtheta) mean(inh_FRtheta)];
datasemsFR = [(std(exc_FRtheta)/sqrt(length(exc_FRtheta)))...
    (std(inh_FRtheta)/sqrt(length(inh_FRtheta)))];

figure
bar([1:2],datameansFR,'k')
hold on
er = errorbar([1:2],datameansFR,datasemsFR);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Baseline Firing Rate (Hz)')
title('CA1 Cell Mean Firing Rate During Run')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p1 h1] = ranksum(exc_FRtheta,inh_FRtheta)

datacombinedthetaFR = [exc_FRtheta; inh_FRtheta];
g1 = repmat({'CA1exc'},length(exc_FRtheta),1);
g2 = repmat({'CA1inh'},length(inh_FRtheta),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedthetaFR,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Theta FR-p = ' num2str(p1)])
ylabel('Firing rate (Hz)')
set(gcf, 'renderer', 'painters')

%%
datameansFRrip = [mean(exc_FRrip) mean(inh_FRrip)];
datasemsFRrip = [(std(exc_FRrip)/sqrt(length(exc_FRrip)))...
    (std(inh_FRrip)/sqrt(length(inh_FRrip)))];

figure
bar([1:2],datameansFRrip,'k')
hold on
er = errorbar([1:2],datameansFRrip,datasemsFRrip);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Firing Rate During SWRs (Hz)')
title('CA1 Cell Firing Rate During Awake SWRs')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p2 h2] = ranksum(exc_FRrip,inh_FRrip)

%%

datameansFRgain = [mean(exc_FRgain) mean(inh_FRgain)];
datasemsFRgain = [(std(exc_FRgain)/sqrt(length(exc_FRgain)))...
    (std(inh_FRgain)/sqrt(length(inh_FRgain)))];

figure
bar([1:2],datameansFRgain,'k')
hold on
er = errorbar([1:2],datameansFRgain,datasemsFRgain);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Awake SWR Firing Rate Gain (Hz)')
title('CA1 Cell Firing Rate Gain')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p3 h3] = ranksum(exc_FRgain,inh_FRgain)

datacombinedswrFRgain = [exc_FRgain; inh_FRgain];
g1 = repmat({'CA1exc'},length(exc_FRgain),1);
g2 = repmat({'CA1inh'},length(inh_FRgain),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedswrFRgain,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['SWR FR gain-p = ' num2str(p3)])
ylabel('SWR FR gain')
set(gcf, 'renderer', 'painters')
ylim([-1 6])

keyboard;