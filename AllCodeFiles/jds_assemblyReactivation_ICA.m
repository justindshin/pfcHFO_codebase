close all;
clear all;
%%
savedata = 1;

PFC = 1;
CA1 = 0;
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8','ER1'};
% animalprefixlist = {'JS14','JS15','ER1','KL8','ER1'};

% animalprefixlist = {'ER1'};
sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17; 10 19];
% sleeps = [1 1];

excludeRips = 0;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    day = 1;
    bin = 20; % ms
    bintime = bin./1000./60; %min
    
    %%
%     load(sprintf('%s%s_spikematrix_ev_allepochallcell100_0%d',dir,animalprefix,day));
    load(sprintf('%s%s_spikematrix_ev_allepochallcell20_0%d',dir,animalprefix,day));
    load(sprintf('%s%ssws0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));
    ripplerun = ripple;
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    
    swsdurs = [];
    for t = 1:2:length(sws{day})
        tmp = sws{day}{t}.total_duration;
        swsdurs = [swsdurs tmp];
    end
    
    swsdurs = swsdurs./60;
    idx = find(swsdurs>=0); 
    
    %compile epochs to analyze
    eps = [];
    for t = 1:length(idx) %need to modify in case of gaps
        if idx(t) > 1
            sleep = sleeps(idx(t),2);
            run = sleep - 1;
            eps = [eps; run sleep];
%         else
%             sleep = 1;
%             run = 2;
%             eps = [eps; run sleep];
        end
    end
    
    ep19idx = find(eps(:,2) <= 17);
    eps = eps(ep19idx,:);
    
    icareactivationtimes = [];
    RtimeStrength = [];
    
    for e = 1:length(eps(:,1)) %need to modify in case of gaps
        ep = eps(e,:);
        
        if PFC == 1
            cellidx = observation_matrix{ep(2)}.ctxidx;
            total_ctx_neuron = observation_matrix{ep(2)}.ncortical;
            %---------Sleep spike matrix-------%
            POST_spike = observation_matrix{ep(2)}.ctxdata;
            %---------W1 and W2 spike matrix-------%
            W_spike = observation_matrix{ep(1)}.ctxdata;
            area = 'PFC';
        elseif CA1 == 1
            cellidx = observation_matrix{ep(2)}.hpidx;
            total_ctx_neuron = observation_matrix{ep(2)}.nhp;
            %---------Sleep spike matrix-------%
            POST_spike = observation_matrix{ep(2)}.hpdata;
            %---------W1 and W2 spike matrix-------%
            W_spike = observation_matrix{ep(1)}.hpdata;
            area = 'CA1';
        end
        %%
        %restrict time to posstart and posend
        st_end_post = [pos{day}{ep(2)}.data(1,1) pos{day}{ep(2)}.data(end,1)];
        %-----POST SWS-----%
        swsep = sws{day}{ep(2)};
        if ~isempty(swsep.starttime)
            swslist_POST = [swsep.starttime swsep.endtime];
            times_filteeg = observation_matrix{ep(2)}.timeeeg;
            idx_st = lookup(st_end_post(1),times_filteeg);
            idx_end = lookup(st_end_post(2),times_filteeg);
            timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
            POST_sws=zeros(size(timevect));
            for sws_seg =1:length(swslist_POST(:,2))
                indtemp = find(timevect >= swslist_POST(sws_seg,1).*1000 & timevect < swslist_POST(sws_seg,2).*1000);
                POST_sws(indtemp) = 1;
            end
            %%
            %--------ripple time------%
            rippost = ripple{day}{ep(2)}; %Use PFC ripple times for the sleep sessions
            riplist_POST(:,1) = rippost.starttime;
            riplist_POST(:,2) = rippost.endtime;
            %         times_filteeg = observation_matrix{ep(2)}.timeeeg;
            %         timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
            POST_rip=zeros(size(timevect));
            for rip_seg =1:length(riplist_POST(:,2))
                indtemp = find(timevect >= riplist_POST(rip_seg,1).*1000 & timevect < riplist_POST(rip_seg,2).*1000);
                POST_rip(indtemp) = 1;
            end
            
            st_end_W = [pos{day}{ep(1)}.data(1,1) pos{day}{ep(1)}.data(end,1)];
            ripW = ripplerun{day}{ep(1)};
            riplist_W(:,1) = ripW.starttime;
            riplist_W(:,2) = ripW.endtime;
            times_filteeg = observation_matrix{ep(1)}.timeeeg;
            idx_st = lookup(st_end_W(1),times_filteeg);
            idx_end = lookup(st_end_W(2),times_filteeg);
            timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
            W_rip=zeros(size(timevect));
            for rip_seg =1:length(riplist_W(:,2))
                indtemp = find(timevect >= riplist_W(rip_seg,1).*1000 & timevect < riplist_W(rip_seg,2).*1000);
                W_rip(indtemp) = 1;
            end
            
            nriplist_W(:,1) = ripW.nstarttime;
            nriplist_W(:,2) = ripW.nendtime;
            times_filteeg = observation_matrix{ep(1)}.timeeeg;
            timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
            W_nrip=zeros(size(timevect));
            for rip_seg =1:length(nriplist_W(:,2))
                indtemp = find(timevect >= nriplist_W(rip_seg,1).*1000 & timevect < nriplist_W(rip_seg,2).*1000);
                W_nrip(indtemp) = 1;
            end
            
            %%
            qPOST_spike = POST_spike;
            qW_spike = W_spike;
            
            % Exclude ripple time for RUN epochs
            if excludeRips == 1
                idx_nrip_w = find(W_rip==0);
                W_spike_raw = W_spike;
                qW_spike = qW_spike(:,idx_nrip_w);
                qW_spike_all = W_spike_raw;
            else
                W_spike_raw = W_spike;
                qW_spike_all = W_spike_raw;
            end
            %%
            
            
            qW_spike_z = zscore(qW_spike')';
            
            qW_spike_z(isnan(qW_spike_z)) = 0; 
            %%
            %-----------correlation matrix----------%
            cW_spike = corr(qW_spike_z'); %Correlation matrix using zscore spike count data
            cW_spike(find(isnan(cW_spike)))=0;
            
            cPOST_spike =corr(qPOST_spike');
            cPOST_spike(find(isnan(cPOST_spike)))=0;
            %%
            %---------------PCA-----------------%
            [u0,s0,v0] = svd(cW_spike);
            sdiag0 = diag(s0);
            %     figure,
            %     plot(sdiag0,'o');
            %     title('eigenvalues')
            
            %-----PCs significant test-----%
            lamdamax = (1+sqrt(length(qW_spike(:,1))/length(qW_spike(1,:)))).^2;%+(length(qW2_spike(:,1))).^(-2/3);
            
            diagind = find(sdiag0 > lamdamax);
            numPCs = length(diagind);
            
            %ICA
            Psign = u0(:,[1:numPCs]); %Restrict ICA to significant eigenvalues
%             Zproj = Psign'*cW_spike;
            Zproj = Psign'*qW_spike_z;
            %     [icasig, A, W] = fastica (Zproj); %Use un-mixing matrix W to get V. V = Psign*W
            [S, H, iter, W] = robustica(Zproj,{});
            V = Psign*W; %Columns of V are the weight vectors of the assembly patterns
            %Scale the weight vectors to unit length and process such that highest
            %absolute value is positive (since sign is arbitrary in ICA)
            
            vtmp = [];
            for t = 1:length(V(1,:))
                w_vectmp = V(:,t);
                w_vec = w_vectmp/norm(w_vectmp);
                min_w = min(w_vec);
                max_w = max(w_vec);
                if abs(min_w) > max_w
                    w_vec = w_vec*(-1);
                end
                vtmp = [vtmp w_vec];
            end
            
            V = vtmp;
            %Calculate projection matrices by taking each column of V and taking
            %the outer product of itself
            
%             figure;
            %Use this section to find high weight cells
            ica_neuron_weights = [];
            cW_spike_ic = cell(1,numPCs);
            for i = 1:numPCs
                tmpmat = V(:,diagind(i))*V(:,diagind(i))';
                tmpmat2 = tmpmat - diag(diag(tmpmat)); %set diagonal to 0
                cW_spike_ic{i} = tmpmat2;
%                 if i<=5
%                     subplot(5,1,i)
%                     stem(V(:,diagind(i)))
%                     title('RUN')
% %                     xlim([0 total_ctx_neuron+1]);
%                 end
                [A B]=sort(V(:,diagind(i)),'descend');
                tmpcellidx = cellidx(B,:);
                tmpcellidx = [tmpcellidx A];
                ica_neuron_weights{i} = tmpcellidx;
                ica_neuron_weights{i} = [cellidx V(:,diagind(i))]; %not sorted
            end
            
            thetaratio = sdiag0(diagind)./lamdamax;
            %%
            %------zscore--------%
            %Finding the reactivation strength over the entire epoch
            qPOST_spike_rip = zscore(POST_spike')';
            qW_spike_rip = zscore(W_spike_raw')';
            qW_spike_rip(isnan(qW_spike_rip)) = 0; %Zscored spike matrix for awake SWRs
            qPOST_spike_rip(isnan(qPOST_spike_rip)) = 0; %Zscored spike matrix for POST SLEEP SWRs
            
            %%
            R_W = zeros(1,length(qW_spike_rip(1,:)));
            R_POST = zeros(1,length(qPOST_spike_rip(1,:)));
            
            for pc = 1:length(diagind)
                PC = cW_spike_ic{pc};
                for t = 1:length(qW_spike_rip(1,:))
                    R1_W_temp = qW_spike_rip(:,t)'*PC*qW_spike_rip(:,t);
                    R_W(pc,t) =  R1_W_temp;
                end
                for t = 1:length(qPOST_spike_rip(1,:))
                    R1_POST_temp = qPOST_spike_rip(:,t)'*PC*qPOST_spike_rip(:,t);
                    R_POST(pc,t) =  R1_POST_temp;
                end
            end
            
            %%
            R_W(find(R_W < 0)) = 0;
            R_POST(find(R_POST <0)) = 0;
            reactivation_strength_run = abs(R_W)';
            trun = (1:length(reactivation_strength_run)).*bin./1000;%s
            reactivation_strength_post = (abs(R_POST)');
            tpost = (1:length(reactivation_strength_post)).*bin./1000;%s
            
            clr = ['k','b','r','m','g','r','k','b','r','m','g','r','k','b','r',...
                'm','g','r','k','b','r','m','g','r','k','b','r','m','g','r','k','b','r','m','g','r'];
            %         figure
            ylimit = max(max(reactivation_strength_post));
%             for i = 1:length(reactivation_strength_run(1,:))
%                 hold on
%                 plot(tpost,reactivation_strength_post(:,i),clr(i) ,'linewidth',1)
%                 xlim([tpost(1), tpost(end)])
%                 ylim([0 ylimit+2])
%                 title('POST')
%                 xlabel('Time (s)')
%             end
            clear riplist_POST riplist_W nriplist_W
            
            times_filteeg_post = observation_matrix{ep(2)}.timeeeg;
            posttimevect = times_filteeg_post(1)*1000:bin:times_filteeg_post(end)*1000;
            [~,swsvec_post] = wb_list2vec(swslist_POST.*1000, posttimevect);
            swsidx_post = find(swsvec_post > 0);
            postripidx = find(POST_rip > 0);
            postidx = intersect(postripidx,swsidx_post);
            y = zeros(length(postidx),1)';
            y(y==0) = 20;
            postidx = postidx.*bin./1000;
            %         figure(2)
            %         hold on
            %         plot(postidx,y,'ro')
            
            %Find the mean and STD of reactivation strength during SWS, Extract
            %times where the signal is +3SD above mean, and record times.
            
            thresh_post = [];
            mean_post = [];
            idx_post = find(POST_sws~=0);
            post_sws_react = reactivation_strength_post(idx_post,:);
            
            for i = 1:length(reactivation_strength_run(1,:))
                posttmp = mean(post_sws_react(:,i));
                mean_post = [mean_post posttmp];
                posttmp = posttmp + (3*std(post_sws_react(:,i)));
                thresh_post = [thresh_post posttmp]; %Find the threshold based on reactivation strengths during SWS
            end
            
            
            times_filteeg_post = observation_matrix{ep(2)}.timeeeg;
            posttimevect = times_filteeg_post(1)*1000:bin:times_filteeg_post(end)*1000;
            
            
            times_filteeg_w = observation_matrix{ep(1)}.timeeeg;
            wtimevect = times_filteeg_w(1)*1000:bin:times_filteeg_w(end)*1000;
            st_end_W = [pos{day}{ep(1)}.data(1,1) pos{day}{ep(1)}.data(end,1)];
            
            stIdx = lookup(st_end_W(1)*1000,wtimevect);
            endIdx = lookup(st_end_W(2)*1000,wtimevect);
            
            wtimevect = wtimevect(stIdx:endIdx);
            %truncate to get rid of empty bins
            reactivation_strength_run = reactivation_strength_run([stIdx:endIdx],:);
            
            %Bin run reactivation strength and determine whether
            %gradually increasing, unchanging, or decreasing
            
            windowWidth = 10; %for 1 second using 100ms bins
            
            rChange = [];
            rvals = [];
            pvals = [];
            for runas = 1:length(reactivation_strength_run(1,:))
                runReactTmp = reactivation_strength_run(:,runas);
                averagedData = zeros(1, floor(length(runReactTmp)/windowWidth));
                binidx = 1:length(averagedData);
                winst = 1;
                for k = 1:floor(length(runReactTmp)/windowWidth)
                    averagedData(k) = mean(runReactTmp(winst:winst+windowWidth-1));
                    winst = winst + windowWidth;
                end
                [r p] = corrcoef(binidx,averagedData);
                rvals = [rvals; r(1,2)];
                pvals = [pvals; p(1,2)];
                if p(1,2) >= 0.05
                    rChange = [rChange; 0];
                elseif (p(1,2) < 0.05) && (r(1,2) > 0)
                    rChange = [rChange; 1];
                elseif (p(1,2) < 0.05) && (r(1,2) < 0)
                    rChange = [rChange; -1];
                end
            end
            
            react_time_strength_all = [];
            react_time_strength_sleep = [];
            
            for ii = 1:length(reactivation_strength_run(1,:))
                
                [~,swsvec_post] = wb_list2vec(swslist_POST.*1000, posttimevect);
                wpctmp = reactivation_strength_run(:,ii);
                postpctmp = reactivation_strength_post(:,ii);
                
                react_time_strength_sleep{ii} = [(wtimevect./1000)' wpctmp];
                
                react_time_strength_all{ii} = [(posttimevect./1000)' postpctmp];
                
                icareactivationtimes{e}.post_strengths{ii} = [(posttimevect./1000)' postpctmp]; %All event strengths, no filter
                icareactivationtimes{e}.run_strengths{ii} = [(wtimevect./1000)' wpctmp];
                icareactivationtimes{e}.runchange = rChange; %increasing, unchanged, decreasing
                icareactivationtimes{e}.runrval = rvals; %correlation rvals
                icareactivationtimes{e}.runpval = pvals; %correlation pvals
                icareactivationtimes{e}.post_means = mean_post;
                icareactivationtimes{e}.postthresh = thresh_post;
                icareactivationtimes{e}.thetaratio = thetaratio;
                icareactivationtimes{e}.pc_weights = ica_neuron_weights;
                icareactivationtimes{e}.cellidx = cellidx;
                icareactivationtimes{e}.timebinsize = bin;
            end
            icareactivationtimes{e}.epochs = ep;
            
            RtimeStrength{ep(2)}.reactivationStrength = react_time_strength_all;
            RtimeStrength{ep(2)}.reactivationStrengthRun = react_time_strength_sleep;
            RtimeStrength{ep(2)}.weights = ica_neuron_weights;
            RtimeStrength{ep(2)}.epochs = ep;
            RtimeStrength{ep(2)}.cellidx = cellidx;
            RtimeStrength{ep(2)}.runchange = rChange;
            clear swslist_POST swslist
        end
    end
    
    if savedata == 1
        if excludeRips == 1
            save(sprintf('%s%s%s_icareactivationtimes20SWSSpk%02d.mat', dir,animalprefix,area,day), 'icareactivationtimes');
            save(sprintf('%s%s%s_RTimeStrengthSleepNewSpk_%02d_%02d.mat', dir,animalprefix,area,bin,day), 'RtimeStrength');
        else
            save(sprintf('%s%s%s_icareactivationtimes20SWSSpkRips%02d.mat', dir,animalprefix,area,day), 'icareactivationtimes');
            save(sprintf('%s%s%s_RTimeStrengthSleepNewSpkRips_%02d_%02d.mat', dir,animalprefix,area,bin,day), 'RtimeStrength');
        end
    end
end
