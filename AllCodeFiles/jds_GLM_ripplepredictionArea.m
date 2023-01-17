%SVM Ripple Type decoding batch
clear all
close all
saveDir = '/Volumes/JUSTIN/SingleDay/ProcessedDataNew/';
savedata = 0;
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS15','JS14','ER1','KL8'};
epochs = [1:2:17];
day = 1;
rtype = 'CA1';
animData = [];
animDataRipSplit = [];
splits = 10;
precision_data = [];
recall_data = [];
fscore_data = [];
precision_s = [];
recall_s = [];
fscore_s = [];
for spl = 1:splits
    for a = 1:length(animalprefixlist)
        prediction_pvals = [];
        prediction_pvals_s = [];
        allPctCorrect = [];
        ripplespkmat = [];
        rippletype = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

        load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
        if isequal(rtype,'PFC')
            [ctxidx, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day, epochs, []); %(tet, cell)
            ctxnum = length(ctxidx(:,1));
            numcells = ctxnum;
            cellidx = ctxidx;
            load(sprintf('%s/%sctxrippletime_noncoordSWS%02d.mat', dir, animalprefix, day));
            nc_rip = ctxripple; clear ctxripple
            load(sprintf('%s/%sctxrippletime_coordSWS%02d.mat', dir, animalprefix, day));
            c_rip = ctxripple; clear ctxripple

        elseif isequal(rtype,'CA1')
            [ctxidx, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day, epochs, []); %(tet, cell)
            hpnum = length(hpidx(:,1));
            numcells = hpnum;
            cellidx = hpidx;
            load(sprintf('%s/%srippletime_noncoordSWS%02d.mat', dir, animalprefix, day));
            nc_rip = ripple; clear ripple
            load(sprintf('%s/%srippletime_coordSWS%02d.mat', dir, animalprefix, day));
            c_rip = ripple; clear ripple
        end

        dat = [];        

        for e = 1:length(epochs)
            epoch = epochs(e);

            nc_riptimes = [nc_rip{day}{epoch}.starttime nc_rip{day}{epoch}.endtime];
            nc_riptimes(:,3) = 0;

            c_riptimes = [c_rip{day}{epoch}.starttime c_rip{day}{epoch}.endtime];
            c_riptimes(:,3) = 1;

            %combine riptimes
            riptimes = sortrows([nc_riptimes; c_riptimes],1);
            ripnum = length(riptimes(:,1));

            if ripnum > 1
                celldata = [];
                spikecounts = [];
                for cellcount = 1:numcells %get spikes for each cell
                    index = [day,epoch,cellidx(cellcount,:)] ;
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                    else
                        spiketimes = [];
                    end
                    spikebins = periodAssign(spiketimes, riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
                    if ~isempty(spiketimes)
                        validspikes = find(spikebins);
                        spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                        spikebins = spikebins(validspikes);
                        tmpcelldata = [spiketimes spikebins];
                    end
                    if ~isempty(spiketimes)
                        tmpcelldata(:,3) = cellcount; %keep count of how many cells active during rip event
                    else
                        tmpcelldata = [0 0 cellcount];
                    end
                    celldata = [celldata; tmpcelldata];
                    spikecount = zeros(1,size(riptimes,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    %                 spikecount(find(spikecount>0)) = 1; %If that cell is active or not
                    spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
                end
                ripplespkmat = [ripplespkmat; spikecounts]; %get feature matrix here
                rippletype = [rippletype; riptimes(:,3)];
            end
            clear nc_riptimes c_riptimes
        end
        %even out ripple numbers
        coordNum = length(find(rippletype == 1));
        indNum = length(find(rippletype ==0));
        numDiff = abs(indNum-coordNum);
        if coordNum < indNum
            delIdx = find(rippletype ==0);
            randidx = delIdx(randperm(length(delIdx)));
            rippletype(randidx(1:numDiff),:) = [];
            ripplespkmat(randidx(1:numDiff),:) = [];
        end
        %GLM with n-fold cross validation
        rip_type_str = num2str(rippletype);
        K = 10;
        cv = cvpartition(rip_type_str, 'kfold',K);
        mse = zeros(K,1);
        shuf_mse = zeros(K,1);
        allshufs = [];
        yhats = [];
        for k=1:K
            % training/testing indices for this fold
            trainIdx = cv.training(k);
            testIdx = cv.test(k);

            % train GLM model
            rtypemat = rippletype(trainIdx);
            warning('off','all');
            mdl = fitclinear(ripplespkmat(trainIdx,:), rtypemat);

            % predict regression output
            Y_hat = predict(mdl, ripplespkmat(testIdx,:));
            
            %Do shuffling

%             for s = 1:5000
%                 shuf = Y_hat(randperm(length(Y_hat)));
%                 corrPred_shuf = rippletype(testIdx) == shuf;
%                 pctCorr_shuf = sum(corrPred_shuf)/length(corrPred_shuf);
%                 shufPct(s) = pctCorr_shuf;
%             end

            % compute p_value
            corrPred = rippletype(testIdx) == Y_hat;

            %% Precision Recall and F score
            [confMat,order] = confusionmat(rippletype(testIdx),Y_hat);
            for i =1:size(confMat,1)
                recall(i)=confMat(i,i)/sum(confMat(i,:));
            end
            recall(isnan(recall))=[];
            Recall = sum(recall)/size(confMat,1);
            
            for i =1:size(confMat,1)
                precision(i)=confMat(i,i)/sum(confMat(:,i));
            end
            Precision=sum(precision)/size(confMat,1);
            F_score=2*Recall*Precision/(Precision+Recall); 
            precision_data = [precision_data; Precision];
            recall_data = [recall_data; Recall];
            fscore_data = [fscore_data; F_score];

            %%
            pctCorr = sum(corrPred)/length(corrPred);
%             p_value = mean(pctCorr <= shufPct);
%             prediction_pvals = [prediction_pvals; p_value];
            allPctCorrect = [allPctCorrect; pctCorr]; %should have 5 values per animal
        end

        %SHUFFLE ORIG DATA AND GET PREDICTION GAIN
        cv2 = cvpartition(rip_type_str, 'kfold',K);
        disp('Generating shuffled dataset...')
        allshufs_s = [];
        yhats_s = [];
        allPctCorrect_shuf = [];
        for p = 1:1000
            p
            shufPctCorr = [];
            for kk=1:K
                %             disp([num2str(K) ' fold cross validation - fold number ' num2str(kk)])
                trainIdx2 = cv2.training(kk);
                testIdx2 = cv2.test(kk);
                %number of shuf models
                % train GLM model
                ripplespkmatshuf = ripplespkmat(randperm(length(ripplespkmat(:,1))),:);
                rtypemat_s = rippletype(trainIdx2);
                warning('off','all');
                mdl2 = fitclinear(ripplespkmatshuf(trainIdx2,:), rtypemat_s);

                % predict regression output
                Y_hat_s = predict(mdl2, ripplespkmat(testIdx2,:));

                %Do shuffling
%                 for s = 1:5000
%                     shuf_shuf = Y_hat_s(randperm(length(Y_hat_s)));
%                     corrPred_shuf = rippletype(testIdx2) == shuf_shuf;
%                     pctCorr_shuf_shuf = sum(corrPred_shuf)/length(corrPred_shuf);
%                     shuf_shufPct(s) = pctCorr_shuf_shuf;
%                 end
                %% Precision recall fscore
                [confMat_s,order] = confusionmat(rippletype(testIdx2),Y_hat_s);
                for i =1:size(confMat_s,1)
                    recall(i)=confMat_s(i,i)/sum(confMat_s(i,:));
                end
                recall(isnan(recall))=[];
                Recall_s = sum(recall)/size(confMat_s,1);

                for i =1:size(confMat_s,1)
                    precision(i)=confMat_s(i,i)/sum(confMat_s(:,i));
                end
                Precision_s = sum(precision)/size(confMat_s,1);
                F_score_s = 2*Recall_s*Precision_s/(Precision_s+Recall_s);
                
                precision_s = [precision_s; Precision_s];
                recall_s = [recall_s; Recall_s];
                fscore_s = [fscore_s; F_score_s];
                % compute p_value for shuf model
                corrPred_s = rippletype(testIdx2) == Y_hat_s;
                pctCorr_s = sum(corrPred_s)/length(corrPred_s);
%                 p_value_s = mean(pctCorr_s <= shuf_shufPct);
%                 prediction_pvals_s = [prediction_pvals_s; p_value_s];
                shufPctCorr = [shufPctCorr; pctCorr_s];
            end
            allPctCorrect_shuf = [allPctCorrect_shuf; mean(shufPctCorr)];
        end
        p_val_anim = mean(mean(allPctCorrect) <= allPctCorrect_shuf);
        animData{a}.meanPctCorrect = mean(allPctCorrect);
        animData{a}.AllPctCorrectCrossVal = allPctCorrect;
        animData{a}.PctCorrectCrossValShuf = allPctCorrect_shuf;
        animData{a}.pval = p_val_anim;
    end
    animDataRipSplit{spl} = animData;
end
if savedata
    save(sprintf('%s%sRipplePrediction_%s.mat',saveDir,rtype,rtype), 'animDataRipSplit');
end

keyboard