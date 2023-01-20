%Get modulation state of cells for each epoch and save. Code is very slow...
clear all;
close all;
%%
savedirX = '/Volumes/JUSTIN/SingleDay/ProcessedDataNew/';

animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};

day = 1;
%animal idx: ZT2 - 1, BG1 - 2, JS34 - 3, JS17 - 4, JS21 - 5
%%

area = 'PFC';
filename = sprintf('Allanim_noncoord250ca1ripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',area);
% filename = sprintf('Allanim_noncoordctxripplemod_run_%s_alldata_sepeps_gather_X6.mat',area);

load([savedirX filename],'allripplemod','allripplemod_idx');

%%
exc_contrib = [];
inh_contrib = [];
allsig = vertcat(allripplemod.rasterShufP2) < 0.05;
flipmatrix_anim = [];
cellidx_mod_anim = [];
modvals_anim = [];
FRvals_anim = [];
FRvals_sleep_anim = [];
posmod = [];
negmod = [];
for a = 1:length(animalprefixlist)
    flipmatrix = [];
    cellidx_mod = [];
    modvals = [];
    FRvals = [];
    FRvals_sleep = [];
    animalprefix = char(animalprefixlist{a});
    animdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    load(sprintf('%s%sspikes%02d.mat',animdir,animalprefix,day));
    
    animidx = find(allripplemod_idx(:,1) == a);
    animdata_idx  = allripplemod_idx(animidx,:);
    %get mod val
    animmod = vertcat(allripplemod.Dm);
    animmod = animmod(animidx);
    animdata_idx = [animdata_idx animmod];
    
    %     animsig = vertcat(allripplemod.sig_ttest);
    animsig = vertcat(allripplemod.rasterShufP2) < 0.05;
    animsig = animsig(animidx);
    
    dm_sig = [animmod animsig];
    
    inh_idx = find(dm_sig(:,1) < 0);
    exc_idx = find(dm_sig(:,1) > 0);
    cellidx = sort([inh_idx; exc_idx]);
    sig_idx = find(dm_sig(:,2) == 1);
    
    allsig_cells = animdata_idx(intersect(cellidx, sig_idx),[3:6]);
    animidx(:,1) = a;
    animidx = animidx(1:length(allsig_cells(:,1)));
    allsig_cells = [animidx allsig_cells];
    
    dummyindex = [animidx allsig_cells(:,[2 3 4])];
    cellidx = dummyindex;
    matchedCells = [];
    doneCells = [0 0];
    for i=1:length(cellidx(:,1))
        if ~ismember(cellidx(i,[3 4]),doneCells,'rows','legacy')
            animtetcell=cellidx(i,[1 3 4]);
            ind=[];
            while rowfind(animtetcell([1 2 3]),dummyindex(:,[1 3 4]))~=0          % collect all rows (epochs)
                ind = [ind rowfind(animtetcell([1 2 3]),dummyindex(:,[1 3 4]))];        % finds the first matching row
                dummyindex(rowfind(animtetcell([1 2 3]),dummyindex(:,[1 3 4])),:)=[0 0 0 0]; % after adding index, remove the corresponding row
                % so you could find the next one
            end
            thisCell = allsig_cells(ind,:);
            %             cellidx_mod = [cellidx_mod; cellidx(i,[3 4])];
            epochs = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
            filterepochs = [2:2:16];
            
            activeEps = [];
            runMeanFr = nan(1,8);
            sleepMeanFr = nan(1,9);
            for f = 1:length(filterepochs)
                currEp = filterepochs(f);
                [ctxidx, hpidx] =  jds_getallepcells_includeall(animdir, animalprefix, day, currEp, []);
                if area == 'PFC'
                    areacellidx = ctxidx;
                elseif area == 'CA1'
                    areacellidx = hpidx;
                end
                if ~isempty(thisCell)
                    [b cellidxtmp] = ismember(thisCell(1,[3 4]), areacellidx, 'rows', 'legacy');
                    if cellidxtmp ~= 0
                        activeEps = [activeEps; currEp];
                        meanFR = spikes{day}{currEp}{thisCell(1,3)}{thisCell(1,4)}.meanrate;
                        runMeanFr(f) = meanFR;
                    end
                end
            end
            %         epochs = [1 2; 2 4; 3 6; 4 8; 5 10; 6 12; 7 14; 8 16];
            if ~isempty(runMeanFr)
%                 if max(runMeanFr) < 7
                if max(runMeanFr) ~= 0 %Include all cells for PFC
                    cellidx_mod = [cellidx_mod; cellidx(i,[3 4])];
                    for e = 1:length(epochs(:,1))
                        modvec = [];
                        ep = epochs(e,2);
                        idx = find(thisCell(:,2) == ep);
                        if ~isempty(idx)
                            cellmod = thisCell(idx,5);
                            if cellmod > 0
                                posmod = [posmod; thisCell(idx,5)];
                                cellmod = 1;
                            elseif cellmod < 0
                                negmod = [negmod; thisCell(idx,5)];
                                cellmod = -1;
                            end
                            modvalvec(e) = thisCell(idx,5);
                        else
                            cellmod = 0;
                            modvalvec(e) = 0;
                        end
                        modVec(e) = cellmod;
                        if ~isempty(spikes{day}{ep}{thisCell(1,3)}{thisCell(1,4)})
                            if isfield(spikes{day}{ep}{thisCell(1,3)}{thisCell(1,4)},'meanrate')
                                sleepMeanFr(e) = spikes{day}{ep}{thisCell(1,3)}{thisCell(1,4)}.meanrate;
                            end
                        end
                    end
                                        
                    flipmatrix = [flipmatrix; modVec];
                    doneCells = [doneCells; cellidx(i,[3 4])];
                    modvals = [modvals; modvalvec];
                    FRvals = [FRvals; runMeanFr];
                    FRvals_sleep = [FRvals_sleep; sleepMeanFr];
                end
            end
        end
    end
    flipmatrix_anim{a} = flipmatrix;
    cellidx_mod_anim{a} = cellidx_mod;
    modvals_anim{a} = modvals;
    FRvals_anim{a} = FRvals;
    FRvals_sleep_anim{a} = FRvals_sleep;
end
anims = length(animalprefixlist);
subplot(anims,1,1)
allAnimMat = [];
savedata = 1;
for an  = 1:length(flipmatrix_anim)
    animalprefix = animalprefixlist{an};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    epochModulation = [];
    subplot(anims, 1, an);
    mattmp = flipmatrix_anim{an};
    rowsum = sum(abs(mattmp),2);
    zeridx = find(rowsum ~= 0);
    newmat = mattmp(zeridx,:);
    cellidx_mod_anim{an} = cellidx_mod_anim{an}(zeridx,:);
    imagesc(newmat);
    colormap(gray)
    allAnimMat = [allAnimMat; newmat];
    epochModulation.modMat = newmat;
    epochModulation.cellidx = cellidx_mod_anim{an};
    epochModulation.modVals = modvals_anim{an};
    epochModulation.runFRvals = FRvals_anim{an};
    epochModulation.sleepFRvals = FRvals_sleep_anim{an};
    
    if savedata == 1
        save(sprintf('%s%s%sca1ripmodsig_epsExcludeHigh%02d.mat', dir,animalprefix,area,day), 'epochModulation');
    end
end
figure(1);
subplot(8,1,1)
title('Significant CA1 Ripple Modulation in PFC')
ylabel('Cell #')
subplot(8,1,8)
xlabel('Sleep Epoch')

for i = 1:length(allAnimMat(:,1))
    tmp = allAnimMat(i,:);
    modidx = find(tmp ~= 0);
    tmp2 = tmp(modidx);
    eachcell{i} = tmp2;
end

numCol = max(cellfun(@length, eachcell));
numRow = length(allAnimMat(:,1));

consistencyMatrix = zeros(numRow,numCol);
transitions = 0;
allFlips = 0;
for i = 1:length(eachcell)
    tmp = eachcell{i};
    transTmp = length(tmp) - 1;
    consistencyMatrix(i,1:length(tmp)) = tmp;
    transitions = transitions + transTmp;
    modFlips = length(find(diff(tmp) ~= 0));
    allFlips = allFlips + modFlips;
end

figure; imagesc(consistencyMatrix); colormap(gray)

modLengths = cellfun(@length, eachcell)
[b, idx] = sort(modLengths,'descend')
figure; imagesc(consistencyMatrix(idx,:)); colormap(gray)

% inhmean = mean(abs(negmod));
% excmean = mean(posmod);
%
% datameans = [excmean inhmean];
% datasems = [(std(posmod)/sqrt(length(posmod))) (std(abs(negmod))/sqrt(length(negmod)))];
%
% figure;
% bar([1:2],datameans,'k')
% hold on
% er = errorbar([1:2],datameans,datasems);
% er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
% xticklabels({'PFCexc','PFCinh'}); xtickangle(45)
keyboard


