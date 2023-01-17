function [c1vsc2, xcorr_sm_all, Zcrosscov, Zcrosscov_sm_all] = js_CA1_PFCrips_xcorr_singleday(animalprefixlist, epochs)

day = 1;
bin = 0.1;
tmax = 5;
sw1 = bin*3; % for smoothing corrln. Necessary?
shufnum = 1000;
noncoordripscorr = [];
noncoordripscorr_shuf = [];
for a = 1:length(animalprefixlist)
    
    animalprefix = char(animalprefixlist(a));
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));% get ripple time
    
    load(sprintf('%s%srippletime_SWS0%d.mat',dir,animalprefix,day));
%     ripple = ctxripple;
%     load(sprintf('%s%sdeltawavetimes_SWS0%d.mat',dir,animalprefix,day));
%     load(sprintf('%s%sctxspindletime_SWS0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));
%     ctxripple = deltawaves;
%     ctxripple = ctxspindle;
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    
%     ripple = deltawaves;
%     ripple = ctxspindle;
    %     ctxripple = deltawaves;
    
    for ep=1:length(epochs)
        
        epoch = epochs(ep);
        
        ctx = ctxripple{day}{epoch};
        hp = ripple{day}{epoch};
        
        if ~isempty(sws{day}{epoch}.starttime)
            swsdur = sws{day}{epoch}.total_duration;
            
            %         if (~isempty(ctx.starttime)) && (~isempty(hp.starttime))
            if (~isempty(ctx)) && (~isempty(hp))
                
                
                %             ctxmidtimes = ((ctx.endtime - ctx.starttime)./2) + ctx.starttime;
                %             hpmidtimes = ((hp.endtime - hp.starttime)./2) + hp.starttime;
                ctxmidtimes = ctx.starttime;
                hpmidtimes = hp.starttime;
                
                ctxmidtimes_shuf = ctx.starttime;
                iri = diff(ctxmidtimes_shuf); %use distribution of IRIs to jitter ctxriptimes
                iri(find(iri > 5)) = [];
                
                for s = 1:shufnum
                    sizeR = [1 length(ctxmidtimes_shuf)] ;
                    num1 = floor(length(ctxmidtimes_shuf)/2);
                    R = zeros(sizeR);  % set all to zero
                    ix = randperm(numel(R)); % randomize the linear indices
                    ix = ix(1:num1); % select the first
                    R(ix) = 1; % set the corresponding positions to 1
                    R(find(R == 0)) = -1;
                    jittimes = (datasample(iri,length(ctxmidtimes_shuf)).*R');
                    
                    ctxmidtimes_shuf = ctxmidtimes + jittimes;
                    xc_shuf = spikexcorr(hpmidtimes, ctxmidtimes_shuf, bin, tmax);
                    
                    if ~isempty(xc_shuf.c1vsc2)
                        Zcrosscov_shuf = zscore(xc_shuf.c1vsc2);
                        
                        nstd=round(sw1/(xc_shuf.time(2) - xc_shuf.time(1))); % will be 3 std
                        g1 = gaussian(nstd, nstd);
                        timebase = xc_shuf.time;
                        bins_run = find(abs(timebase) <= tmax); % +/- Corrln window
                        
                        %                     xcorr_sm = smoothvect(xc.c1vsc2, g1);
                        Zcrosscov_sm_shuf = smoothvect(Zcrosscov_shuf, g1);% smoothed
                        %                     Zcrosscov_sm = Zcrosscov;
                        
                        noncoordripscorr_shuf = [noncoordripscorr_shuf; Zcrosscov_sm_shuf];
                    end
                end
                
%                 xc = spikexcorr(ctxmidtimes, hpmidtimes, bin, tmax);
                xc = spikexcorr(hpmidtimes, ctxmidtimes, bin, tmax);
                %              xc = spikexcorr(hpmidtimes, ctxmidtimes, bin, tmax);
                %             plot(xc.c1vsc2);
                
                % converted to zscore
                
                %             p1 = xc.nspikes1/swsdur; p2 = xc.nspikes2/swsdur; % Fir rate in Hz
                %             exp_p = p1*p2; % per sec
                %             crosscov = (xc.c1vsc2 ./ (bin*swsdur))-exp_p;
                % Convert to Z-score
                %             factor = sqrt((bin*swsdur) / exp_p);
                %             Zcrosscov = crosscov .* (factor);
                if ~isempty(xc.c1vsc2)
                    Zcrosscov = zscore(xc.c1vsc2);
                    
                    nstd=round(sw1/(xc.time(2) - xc.time(1))); % will be 3 std
                    g1 = gaussian(nstd, nstd);
                    timebase = xc.time;
                    bins_run = find(abs(timebase) <= tmax); % +/- Corrln window
                    
%                     xcorr_sm = smoothvect(xc.c1vsc2, g1);
                    Zcrosscov_sm = smoothvect(Zcrosscov, g1);% smoothed
%                     Zcrosscov_sm = Zcrosscov;
                    
                    noncoordripscorr = [noncoordripscorr; Zcrosscov_sm];
                end
            end
        end
    end
end

% figure
% plot(smooth(nanmean(noncoordripscorr)));
% 
% figure
% boundedline([-(tmax/bin):(tmax/bin)-1],smooth(nanmean(noncoordripscorr)),smooth(nanstd(noncoordripscorr)./sqrt(size(noncoordripscorr,1))),'-k');
% xlabel(['Lag (bins) - ' num2str(bin) 'ms bins'])
% ylabel('Z-score')
% title('Cross Corr - Coordinated PFC and CA1 Ripples')
% xlim([-40 40])
% set(gcf, 'renderer', 'painters')
% x = [0 0];
% y = [-0.5 4];
% hold on
% plot(x,y,'--r')

figure
boundedline([-(tmax/bin):(tmax/bin)-1],nanmean(noncoordripscorr),nanstd(noncoordripscorr)./sqrt(size(noncoordripscorr,1)),'-k');
xlabel(['Lag (bins) - ' num2str(bin) 'ms bins'])
ylabel('Z-score')
xlim([-40 40])
set(gcf, 'renderer', 'painters')

%%
%Compare to shuffled times
if shufnum > 0
    figure; hold on;
    boundedline([-(tmax/bin):(tmax/bin)-1],nanmean(noncoordripscorr),nanstd(noncoordripscorr)./sqrt(size(noncoordripscorr,1)),'-k');
    boundedline([-(tmax/bin):(tmax/bin)-1],nanmean(noncoordripscorr_shuf),nanstd(noncoordripscorr_shuf)./sqrt(size(noncoordripscorr_shuf,1)),'-r');
    xlabel(['Lag (bins) - ' num2str(bin) 'ms bins'])
    ylabel('Z-score')
end
keyboard
