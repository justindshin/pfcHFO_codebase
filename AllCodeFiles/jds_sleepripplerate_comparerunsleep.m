function jds_sleepripplerate_comparerunsleep(animalprefixlist)
%Examines temporal structure of coord and noncoord ripples (e.g. when do
%they occur during SWS (early vs late, intermixed, etc)?
%get rate of coord vs non coord over sessions

%Look at each bout of SWS, concatenated SWS across whole sleep
%epoch, and across all epochs

%compute a crosscorrelogram for coord and non coord ripples and perform
%spectral analysis on result to determine whether there is periodic
%structure

day = 1;
bin = 0.05;

cRipRates_sleep = [];
hRipRates_sleep = [];

cRipRates_sleep_nc = [];
hRipRates_sleep_nc = [];

coordinatedRipRate = [];

cRipRates_run = [];
hRipRates_run = [];
swsAll = [];

for a = 1:length(animalprefixlist)
    hpripcnt = 0;
    ctxripcnt = 0;
    coordripcnt = 0;
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    
    load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));
    ctx_sws = ctxripple;
    load(sprintf('%s%srippletime_SWS0%d.mat',dir,animalprefix,day));
    hp_sws = ripple;
    load(sprintf('%s%sctxrippletime_ALL0%d.mat',dir,animalprefix,day));
    ctx_run = ctxripple;
    load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));
    hp_run = ripple;

    load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    ctx_sws_nc = ctxripple;
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    hp_sws_nc = ripple;
    load(sprintf('%s%sripplecoordinationSWS0%d.mat',dir,animalprefix,day));
    coordRips = ctxripple;
    
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%simmobiletime0%d.mat',dir,animalprefix,day));
    
    epochs = [1:17];
    
    for ep = 1:length(epochs)
        concat_crip = zeros(1000,1)';
        concat_hrip = zeros(1000,1)';
        epoch = epochs(ep);
        
        if mod(epoch,2) == 1
            swsdur = sws{day}{epoch}.total_duration;
            
            if swsdur/60 > 1
                
                ctxriptimestmp = [ctx_sws{day}{epoch}.starttime ctx_sws{day}{epoch}.endtime];
                hpriptimestmp = [hp_sws{day}{epoch}.starttime hp_sws{day}{epoch}.endtime];
                ctxncriptimestmp = [ctx_sws_nc{day}{epoch}.starttime ctx_sws_nc{day}{epoch}.endtime];
                hpncriptimestmp = [hp_sws_nc{day}{epoch}.starttime hp_sws_nc{day}{epoch}.endtime];
                coordriptimestmp = [ripplecoordination{day}{epoch}.starttime ripplecoordination{day}{epoch}.endtime];
                
                c_rate = length(ctxriptimestmp(:,1))/swsdur;
                h_rate = length(hpriptimestmp(:,1))/swsdur;
                c_ncrate = length(ctxncriptimestmp(:,1))/swsdur;
                h_ncrate = length(hpncriptimestmp(:,1))/swsdur;
                coord_rate = length(coordriptimestmp(:,1))/swsdur;
                
                cRipRates_sleep = [cRipRates_sleep c_rate];
                hRipRates_sleep = [hRipRates_sleep h_rate];
                cRipRates_sleep_nc = [cRipRates_sleep_nc c_ncrate];
                hRipRates_sleep_nc = [hRipRates_sleep_nc h_ncrate];
                coordinatedRipRate = [coordinatedRipRate coord_rate];

                swsAll = [swsAll; swsdur];
            else
                cRipRates_sleep = [cRipRates_sleep NaN];
                hRipRates_sleep = [hRipRates_sleep NaN];
                cRipRates_sleep_nc = [cRipRates_sleep_nc NaN];
                hRipRates_sleep_nc = [hRipRates_sleep_nc NaN];
                coordinatedRipRate = [coordinatedRipRate NaN];
            end
            
        elseif mod(epoch,2) == 0
            
            imdur = immobility{day}{epoch}.total_duration;
            
            ctxriptimestmp = [ctx_run{day}{epoch}.starttime ctx_run{day}{epoch}.endtime];
            hpriptimestmp = [hp_run{day}{epoch}.starttime hp_run{day}{epoch}.endtime];
            
            c_rate = length(ctxriptimestmp(:,1))/imdur;
            h_rate = length(hpriptimestmp(:,1))/imdur;
            
            cRipRates_run = [cRipRates_run c_rate];
            hRipRates_run = [hRipRates_run h_rate];
        end
    end
end

meanSleepCtxHp = [nanmean(cRipRates_sleep) nanmean(hRipRates_sleep)]
semSleepCtxHp = [(nanstd(cRipRates_sleep)./sqrt(length(find(~isnan(cRipRates_sleep)))))...
    (nanstd(hRipRates_sleep)./sqrt(length(find(~isnan(hRipRates_sleep)))))]

meanSleepCtxHp_nc = [nanmean(cRipRates_sleep_nc) nanmean(hRipRates_sleep_nc)]
semSleepCtxHp_nc = [(nanstd(cRipRates_sleep_nc)./sqrt(length(find(~isnan(cRipRates_sleep_nc)))))...
    (nanstd(hRipRates_sleep_nc)./sqrt(length(find(~isnan(hRipRates_sleep_nc)))))]

meanSleepCtxHp_c = [nanmean(coordinatedRipRate) nanmean(coordinatedRipRate)]
semSleepCtxHp_c = [(nanstd(coordinatedRipRate)./sqrt(length(find(~isnan(coordinatedRipRate)))))...
    (nanstd(coordinatedRipRate)./sqrt(length(find(~isnan(coordinatedRipRate)))))]

meanRunCtxHp = [nanmean(cRipRates_run) nanmean(hRipRates_run)]
semRunCtxHp = [(nanstd(cRipRates_run)./sqrt(length(find(~isnan(cRipRates_run)))))...
    (nanstd(hRipRates_run)./sqrt(length(find(~isnan(hRipRates_run)))))]

figure; 
errorbar([1 2],meanSleepCtxHp,semSleepCtxHp,'r','LineStyle','none','Marker','o'); hold on
errorbar([1 2],meanRunCtxHp,semRunCtxHp,'k','LineStyle','none','Marker','o'); 
errorbar([1 2],meanSleepCtxHp_nc,semSleepCtxHp_nc,'b','LineStyle','none','Marker','o');
errorbar([1 2],meanSleepCtxHp_c,semSleepCtxHp_c,'m','LineStyle','none','Marker','o'); 


xlim([0.5 2.5])
legend({'Sleep','Run','Independent','Coordinated'})
xticks([1 2])
xticklabels({'PFC','CA1'})
swsAll = swsAll/60;
meanNrem = mean(swsAll)
semNrem = std(swsAll)./sqrt(length(swsAll))

keyboard
