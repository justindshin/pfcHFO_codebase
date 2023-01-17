function [gain, s_gain] = jds_GLM_CA1spikePrediction_allep_yu(animalprefixlist)


%%
%Things to incorporate
%-Need to restrict to cells that pass specific criterion (Active in
%more than 10 SWR events (Rothschild 2017)
%%
savedata = 1;


allanim_pg_log_steps_iter = [];
allanim_pg_log_steps_iter_s = [];
p_vals_steps_iter = [];

% timesteps = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
% timesteps = [0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2];
timesteps = [1; 2; 3];
% timesteps = [-0.2 0.2];

datadir = '/Volumes/JUSTIN/FigureWorking/CA1ReactivationSuppression/DataMats/';

for steps = 1:length(timesteps(:,1))
    allanim_pg_log_steps = [];
    allanim_pg_log_steps_s = [];
    p_vals_steps = [];
    for runs = 1:10
        allanim_pg = [];
        allanim_pg_log = [];
        allanim_pg_s = [];
        allanim_pg_log_s = [];
        allanim_CA1resp = [];
        allanim_PFCmatrix = [];
        for a = 1:length(animalprefixlist)
            animalprefix = animalprefixlist{a};
            
            if strcmp(animalprefix,'ZT2')
                savedir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
            elseif strcmp(animalprefix,'JS34')
                savedir = ('/Volumes/JUSTIN/SingleDay/JS34_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/JS34_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/JS34_direct/');
            elseif strcmp(animalprefix,'BG1')
                savedir = ('/Volumes/JUSTIN/SingleDay/BG1_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/BG1_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/BG1_direct/');
            elseif strcmp(animalprefix,'JS17')
                savedir = ('/Volumes/JUSTIN/SingleDay/JS17_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/JS17_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/JS17_direct/');
            elseif strcmp(animalprefix,'JS21')
                savedir = ('/Volumes/JUSTIN/SingleDay/JS21_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/JS21_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/JS21_direct/');
            elseif strcmp(animalprefix,'JS14')
                savedir = ('/Volumes/JUSTIN/SingleDay/JS14_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/JS14_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/JS14_direct/');
            elseif strcmp(animalprefix,'JS15')
                savedir = ('/Volumes/JUSTIN/SingleDay/JS15_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/JS15_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/JS15_direct/');
            elseif strcmp(animalprefix,'ER1')
                savedir = ('/Volumes/JUSTIN/SingleDay/ER1_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/ER1_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/ER1_direct/');
            elseif strcmp(animalprefix,'KL8')
                savedir = ('/Volumes/JUSTIN/SingleDay/KL8_direct/');
                ripdir = ('/Volumes/JUSTIN/SingleDay/KL8_direct/');
                dir = ('/Volumes/JUSTIN/SingleDay/KL8_direct/');
            end
            
            
            %%
            eps = [1:2:17];
            
            %-----match neurons across epochs-----%
            day = 1; %always single day expt
            %match across all epochs since we're concatenating analysis over all eps
            [ctxidx, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
            ctxnum = length(ctxidx(:,1));
            hpnum = length(hpidx(:,1));
            
            %-----create the event matrix during SWRs-----%
            spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
            
            % get ripple time
%             load(sprintf('%s%sripplecoordinationSWS0%d.mat',ripdir,animalprefix,day));
            load(sprintf('%s%srippleframes_SWS0%d.mat',ripdir,animalprefix,day));
            
            %Get matrix of spikes per cell per SWR event (colums = diff cells, rows =
            %diff SWRs) - Same matrix for each PFC cello
            PFCmatrix = [];
            CA1_resp_tmp = [];
            for ep = eps
                clear pfc_riptimes hpc_riptimes riptimes
                cell2use_idx = zeros(ctxnum,1);
                disp([animalprefix '-epoch-' num2str(ep)])
                
%                 rip = ripplecoordination{day}{ep};
                rip = ripple{day}{ep};
                if length(rip.outstarttime) > 1
                    %             randlength = randi(2000,length(rip.starttime),1)/1000;
                    %             halfriptime = (rip.endtime - rip.starttime)./2;
                    %             ripmids = rip.starttime + halfriptime;
                    %
                    %             riptimes(:,1) = ripmids - 0.2;
                    %             riptimes(:,2) = ripmids + 0.2;

                    %                 riptimes(:,1) = rip.starttime - 0.2;
                    %                 riptimes(:,2) = rip.endtime + 0.2;

                    %             p_riptimes(:,1) = rip.starttime - 0.2;
                    %             p_riptimes(:,2) = rip.starttime;
                    %

                    if steps == 1
                        riptimes(:,1) = rip.outstarttime - 0.6;
                        riptimes(:,2) = rip.outstarttime;
                    elseif steps == 2
                        riptimes(:,1) = rip.outstarttime;
                        riptimes(:,2) = rip.outendtime;
                    elseif steps == 3
                        riptimes(:,1) = rip.outendtime;
                        riptimes(:,2) = rip.outendtime + 0.6;
                    end
                    
                    %             if ~isempty(rip.starttime)
                    %                 rip_starttime = riptimes(:,1)*1000;  % in ms
                    % %
                    % %                 rip_endtime = riptimes(:,2)*1000;  % in ms
                    %
                    % %             Find ripples separated by at least 500ms -- sj_getpopulationevents2
                    %                 iri = diff(rip_starttime);
                    %                 keepidx = [1;find(iri>=500)+1];
                    %
                    %             %     riplength = rip_endtime - rip_starttime;
                    %             %     keepidx = find(riplength >= 50);% use the ripple last for more than 50 ms
                    %             %     keepidx = intersect(keepidx,keepidx2);
                    %             %
                    %                 riptimes = riptimes(keepidx,:);
                    %             end
                    %800ms before rip to 600ms before
                    %     ripplelengths = riptimes(:,2) - riptimes(:,1);
                    
                    %What peri-SWR ripple to analyze
                    pfc_riptimes = riptimes(:,1);
                    pfc_riptimes(:,2) = riptimes(:,2);
                    
                    hpc_riptimes = riptimes(:,1);
                    hpc_riptimes(:,2) = riptimes(:,2);
                    
                    %     crit_num = floor(length(pfc2_riptimes(:,1))*0.20);
                    %     if crit_num < 10
                    %         crit_num = 10;
                    %     end
                    pfc_matrix_tmp = [];
                    if length(pfc_riptimes(:,1)) > 1
                        for cellcount = 1:ctxnum %get spikes for each cell
                            index = [day,ep,ctxidx(cellcount,:)] ;
                            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                            else
                                spiketimes = [];
                            end
                            spikebins = periodAssign(spiketimes, pfc_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
                            if ~isempty(spiketimes)
                                validspikes = find(spikebins);
                                spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                                spikebins = spikebins(validspikes);
                            end
                            spikecount = zeros(1,size(pfc_riptimes,1));
                            for i = 1:length(spikebins)
                                spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                            end
                            pfc_matrix_tmp = [pfc_matrix_tmp spikecount'];
                        end
                        PFCmatrix = [PFCmatrix; pfc_matrix_tmp]; %concatenating num spikes per cell, per event
                    end
                    
                    %GET CA1 CELL DATA
                    if length(hpc_riptimes(:,1)) > 1
                        for cellcount = 1:hpnum %get spikes for each cell
                            index = [day,ep,hpidx(cellcount,:)] ;
                            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                            else
                                spiketimes = [];
                            end
                            spikebins = periodAssign(spiketimes, hpc_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
                            if ~isempty(spiketimes)
                                validspikes = find(spikebins);
                                spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                                spikebins = spikebins(validspikes);
                            end
                            spikecount = zeros(1,size(hpc_riptimes,1));
                            for i = 1:length(spikebins)
                                spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                            end
                            CA1_resp_tmp{ep}{cellcount}.data = spikecount';
                            CA1_resp_tmp{ep}{cellcount}.CA1cellidx = hpidx(cellcount,:);
                        end
                    end
                end
            end
            
            CA1_resp = [];
            for hp = 1:hpnum
                thiscell_all = [];
                for ep = eps
                    if ~isempty(CA1_resp_tmp{ep})
                        tmp = CA1_resp_tmp{ep}{hp}.data;
                        thiscell_all = [thiscell_all; tmp];
                    end
                end
                CA1_resp{hp} = thiscell_all;
            end
            
            allanim_CA1resp{a} = CA1_resp;
            allanim_PFCmatrix{a} = PFCmatrix;
            
            %GLM with n-fold cross validation
            mse_CA1 = [];
            s_mse_CA1 = [];
            gain = [];
            s_gain = [];
            for i = 1:hpnum
                spk_cnt_str = num2str(CA1_resp{i});
                K = 10;
                cv = cvpartition(spk_cnt_str, 'kfold',K);
                disp(['Cell number ' num2str(i) ' out of ' num2str(hpnum)])% ' - step ' num2str(steps)])
                mse = zeros(K,1);
                shuf_mse = zeros(K,1);
                allshufs = [];
                yhats = [];
                for k=1:K
                    % training/testing indices for this fold
                    trainIdx = cv.training(k);
                    testIdx = cv.test(k);
                    
                    % train GLM model
                    CA1mat = CA1_resp{i};
                    CA1mat2 = CA1mat(trainIdx);
                    warning('off','all');
                    mdl = fitglm(PFCmatrix(trainIdx,:), CA1mat2,'linear','Distribution', 'poisson'); %Shuffle PFCmat2 to determine shuffled data to get error bars
                    
                    constarray = table2array(mdl.Coefficients(:,1));
                    % predict regression output
                    %             Y_hat = predict(mdl, PFCmatrix(testIdx,:));
                    Y_hat = glmval(constarray, PFCmatrix(testIdx,:),'log');
                    
                    %Do shuffling
                    for s = 1:5000
                        shuf = Y_hat(randperm(length(Y_hat)));
                        shuf_err(s) = mean(abs(CA1mat(testIdx) - shuf));
                    end
                    
                    % compute mean squared error
                    mse(k) = mean(abs(CA1mat(testIdx) - Y_hat));
                    shuf_mse(k) = mean(shuf_err);
                    allshufs = [allshufs; shuf_err'];
                    yhats = [yhats; (CA1mat(testIdx) - Y_hat)];
                end
                mse_CA1{i}.mse = mse;
                mse_CA1{i}.mean_mse = mean(mse);
                mse_CA1{i}.CA1idx = hpidx(i,:);
                mse_CA1{i}.mean_shuf_mse = mean(shuf_mse);
                mse_CA1{i}.shuf_mse = shuf_mse;
                mse_CA1{i}.allshufs = allshufs;
                mse_CA1{i}.yhats = yhats;
                %SHUFFLE ORIG DATA AND GET PREDICTION GAIN
                cv2 = cvpartition(spk_cnt_str, 'kfold',K);
                disp('Generating shuffled dataset...')
                allshufs_s = [];
                yhats_s = [];
                for kk=1:K
                    disp([num2str(K) ' fold cross validation - fold number ' num2str(kk)])
                    trainIdx2 = cv2.training(kk);
                    testIdx2 = cv2.test(kk);
                    %             for p = 1:1 %number of shuf models
                    %                 % train GLM model
                    %                 shuf_CA1mat = CA1_resp{i};
                    %                 shuf_CA1mat = shuf_CA1mat(randperm(length(shuf_CA1mat)));
                    %                 CA1mat2_shuf = shuf_CA1mat(trainIdx2);
                    %                 warning('off','all');
                    %                 mdl2 = fitglm(PFCmatrix(trainIdx2,:), CA1mat2_shuf,'Distribution', 'poisson'); %Shuffle PFCmat2 to determine shuffled data to get error bars
                    %                 constarray_s = table2array(mdl2.Coefficients(:,1));
                    %
                    %                 % predict regression output
                    %                 Y_hat_s = glmval(constarray_s, PFCmatrix(testIdx2,:),'log');
                    %                 %                 Y_hat_s = predict(mdl2, PFCmatrix(testIdx2,:));
                    %
                    %                 %Do shuffling
                    %                 for s = 1:5000
                    %                     shuf_shuf = Y_hat_s(randperm(length(Y_hat_s)));
                    %                     shuf_shuf_err(s) = mean(abs(shuf_CA1mat(testIdx2) - shuf_shuf));
                    %                 end
                    %
                    %                 mse_s(kk) = mean(abs(shuf_CA1mat(testIdx2) - Y_hat_s));
                    %                 shuf_mse_s(kk) = mean(shuf_shuf_err);
                    %                 allshufs_s = [allshufs_s; shuf_shuf_err'];
                    %                 yhats_s = [yhats_s; (shuf_CA1mat(testIdx2) - Y_hat_s)];
                    %
                    %                 %                 s_gain{i}{kk}.data = (mean(shuf_shuf_err)/mean(abs(shuf_CA1mat(testIdx2) - Y_hat_s)));
                    %                 %                 s_gain{i}{kk}.foldnumber = kk;
                    %                 %                 s_gain{i}{kk}.cellidx = hpidx(i,:);
                    %                 %                 shuf_shuf_err = [];
                    %             end
                    for p = 1:1 %number of shuf models
                        % train GLM model
                        PFCmatrixshuf = PFCmatrix(randperm(length(PFCmatrix(:,1))),:);
                        CA1mat = CA1_resp{i};
                        CA1mat2 = CA1mat(trainIdx2);
                        warning('off','all');
                        mdl2 = fitglm(PFCmatrixshuf(trainIdx2,:), CA1mat2,'Distribution', 'poisson'); %Shuffle PFCmat2 to determine shuffled data to get error bars
                        constarray_s = table2array(mdl2.Coefficients(:,1));
                        
                        % predict regression output
                        Y_hat_s = glmval(constarray_s, PFCmatrixshuf(testIdx2,:),'log');
                        %                 Y_hat_s = predict(mdl2, PFCmatrix(testIdx2,:));
                        
                        %Do shuffling
                        for s = 1:5000
                            shuf_shuf = Y_hat_s(randperm(length(Y_hat_s)));
                            shuf_shuf_err(s) = mean(abs(CA1mat(testIdx2) - shuf_shuf));
                        end
                        
                        mse_s(kk) = mean(abs(CA1mat(testIdx2) - Y_hat_s));
                        shuf_mse_s(kk) = mean(shuf_shuf_err);
                        allshufs_s = [allshufs_s; shuf_shuf_err'];
                        yhats_s = [yhats_s; (CA1mat(testIdx2) - Y_hat_s)];
                        
                        %                 s_gain{i}{kk}.data = (mean(shuf_shuf_err)/mean(abs(shuf_CA1mat(testIdx2) - Y_hat_s)));
                        %                 s_gain{i}{kk}.foldnumber = kk;
                        %                 s_gain{i}{kk}.cellidx = hpidx(i,:);
                        %                 shuf_shuf_err = [];
                    end
                end
                s_mse_CA1{i}.mse_s = mse_s;
                s_mse_CA1{i}.mean_mse_s = mean(mse);
                s_mse_CA1{i}.CA1idx = hpidx(i,:);
                s_mse_CA1{i}.mean_shuf_mse_s = mean(shuf_mse_s);
                s_mse_CA1{i}.shuf_mse_s = shuf_mse_s;
                s_mse_CA1{i}.allshufs_s = allshufs_s;
                s_mse_CA1{i}.yhats_s = yhats_s;
            end
            
            pred_gain = [];
            pred_gain_log = [];
            for i = 1:length(mse_CA1)
                mse = mean(mse_CA1{i}.mse);
                shuf_mse = mean(mse_CA1{i}.shuf_mse);
                pg = (shuf_mse/mse);
                pg_log = log10(shuf_mse/mse);
                pred_gain = [pred_gain pg];
                pred_gain_log = [pred_gain_log pg_log];
            end
            
            allanim_pg = [allanim_pg; pred_gain'];
            allanim_pg_log = [allanim_pg_log; pred_gain_log'];
            
            pred_gain_s = [];
            pred_gain_log_s = [];
            for i = 1:length(s_mse_CA1)
                mse_s = mean(s_mse_CA1{i}.mse_s);
                shuf_mse_s = mean(s_mse_CA1{i}.shuf_mse_s);
                pg_s = (shuf_mse_s/mse_s);
                pg_log_s = log10(shuf_mse_s/mse_s);
                pred_gain_s = [pred_gain_s pg_s];
                pred_gain_log_s = [pred_gain_log_s pg_log_s];
            end
            
            allanim_pg_s = [allanim_pg_s; pred_gain_s'];
            allanim_pg_log_s = [allanim_pg_log_s; pred_gain_log_s'];
            
            %     m_pred_gain = mean(pred_gain);
            %     gain.meanGain = mean(pred_gain);
            %     gain.allCellGain = pred_gain;
            %     gain.meanGain_log = mean(pred_gain_log);
            %     gain.allCellGain_log = pred_gain_log;
            %     gain.animal = animalprefix;
            %     gain.epoch = 'all epochs combined - matched cells';
            %     gain.mse = mse_CA1;
            %
            
            
            disp([animalprefix ' processing complete'])
        end
        
        allanim_pg_log_steps{runs} = allanim_pg_log;
        allanim_pg_log_steps_s{runs} = allanim_pg_log_s;
        
        %     p_vals_steps(steps) = permutationTest(allanim_pg_log, allanim_pg_log_s, 10000, 'plotresult', 1,'sidedness','larger');
        p_vals_steps(runs) = ranksum(allanim_pg_log,allanim_pg_log_s);
    end
    allanim_pg_log_steps_iter{steps} = allanim_pg_log_steps;
    allanim_pg_log_steps_iter_s{steps} = allanim_pg_log_steps_s;
    p_vals_steps_iter{steps} = p_vals_steps;
    
    pgData.data = allanim_pg_log_steps_iter;
    pgData.shuffled = allanim_pg_log_steps_iter_s;
    pgData.pvals = p_vals_steps_iter;
    
%     if savedata == 1
%         save(sprintf('%sGLMripcoord_timesteps_PFCpredCA1.mat', datadir), 'pgData');
%     end
%     if savedata == 1
%         save(sprintf('%sGLMinframeripcoord_timesteps_PFCpredCA1.mat', datadir), 'pgData');
%     end
    if savedata == 1
        save(sprintf('%sGLMoutframeripcoord_timesteps_PFCpredCA1.mat', datadir), 'pgData');
    end
end
figure; plot(p_vals_steps)
permutationTest(allanim_pg_log, allanim_pg_log_s, 10000, 'plotresult', 1,'sidedness','larger');
mean_data = mean(allanim_pg_log);
mean_shuf = mean(allanim_pg_log_s);
data_sem = std(allanim_pg_log)./sqrt(size(allanim_pg_log,1))
shuf_sem = std(allanim_pg_log_s)./sqrt(size(allanim_pg_log_s,1))
figure
combdata = [mean_data mean_shuf];
combsem = [data_sem shuf_sem];
bar(combdata,'k'); hold on;
er = errorbar([1 2],combdata,combsem);
er.Color = 'k';
er.LineStyle = 'none';
er.LineWidth = 2;

figure; hold on
allDat = [pgData.data{1}{1} pgData.data{2}{1} pgData.data{3}{1}];
meanDat = mean(allDat);
semDat = std(allDat)./sqrt(length(allDat(:,1)));
bar(meanDat,'k')
errorbar(1:3,meanDat, semDat,'-k','LineStyle', 'none')
for s = 1:2
    x = [s s+1];
    tmp = pgData.data{s}{1};
    tmp2 = pgData.data{s+1}{1};
    plot(x,[tmp tmp2],'-r')
end
xlim([0.5 3.5])


keyboard

