function jds_PCCreplay_decoding_CA1_allcells(animalprefix,day,ep,cellcountthresh)

savedata = 1;

savedir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
ripdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

if mod(ep,2) == 0
    eprun = ep;
elseif ep == 1
    eprun = 2;
elseif mod(ep,2) ~= 0
    eprun = (ep-1);
end

%%
%-----match neurons across epochs-----%
[ctxidx, hpidx] = matchidx_acrossrunsleep_singleday(savedir, animalprefix, day, ep, []); %(tet, cell)
ctxnum = length(ctxidx(:,1));
hpnum = length(hpidx(:,1));
tBinSz = 10; %default temporal bin in ms [typically around 15ms]
wellcutoff = 0; %cm; make 0 for sleep
%%
%-----create the ratemaps [nPosBin x nHPCells]-----%
rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
cellidxm = [];
load(sprintf('%s%slinfields0%d.mat',savedir,animalprefix,day)); % get linearized place fields
for i = 1:hpnum
    cind = hpidx(i,:);
    if (length(linfields{day}{eprun})>= cind(1))
        if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
            linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
        else
            linfield1 =[];
        end
    else
        linfield1=[];
    end
    
    if ~isempty(linfield1)
        linfield_hp = [];
        lintrack_hp = [];
        pos_hp = [];
        for track = 1:4 %concatenating all linearized firing rates
            temp1 = linfield1{track};
            pos1 = temp1(:,1);
            lintrack1 = ones(size(pos1))*track; %get length of lin traj and convert to vector of traj number
            occnormrate1 = temp1(:,5);
            linfield_hp = [linfield_hp;occnormrate1];
            pos_hp = [pos_hp;pos1];
            lintrack_hp = [lintrack_hp;lintrack1];
        end
        if (max(linfield_hp) >= 3) % peak firing rate max larger than 3 Hz
            a = find(isnan(linfield_hp));
            %pad nan
            if ~isempty(a)
                [lo,hi]= findcontiguous(a);  %find contiguous NaNs
                for ii = 1:length(lo)
                    if lo(ii) > 1 & hi(ii) < length(linfield_hp)
                        fill = linspace(linfield_hp(lo(ii)-1), ...
                            linfield_hp(hi(ii)+1), hi(ii)-lo(ii)+1);
                        linfield_hp(lo(ii):hi(ii)) = fill;
                    end
                end
            end
            rm = [rm;linfield_hp']; %horizontal vector of firing rates (each row, diff cell)
            pm = [pm;pos_hp']; %horizontal vector of pos
            tm = [tm;lintrack_hp'];
            cellidxm = [cellidxm; cind];
        end
    end
end
rm = rm'; %[nPosBin x nHPCells]
pm = pm';
tm = tm';
for i = 1:4 %Cutt off reward well positions
    pm_traj = pm(find(tm == i));
    maxpos = max(max(pm_traj));
    rm(find(tm == i & pm <= wellcutoff)) = 0;
    rm(find(tm == i & pm >= maxpos-wellcutoff)) = 0;
end
rm = rm+ (eps.^8); %Add a small number so there are no zeros
expecSpk =rm.*tBinSz./1000; %[nPos x nCells] Expected number of spikes per bin
hpnum = length(rm(1,:)); % update (only cells with >3hz FR)

%%
%-----create the event matrix during SWRs-----%
spikes = loaddatastruct(savedir, animalprefix, 'spikes', day); % get spikes

%%DESIGNATE DECODING TIMES
% get ripple time
% load(sprintf('%s%srippletime_ALL0%d.mat',ripdir,animalprefix,day));
% rip = ripple{day}{ep};
% riptimes(:,1) = rip.starttime;
% riptimes(:,2) = rip.endtime;
%
% rip_starttime = 1000*riptimes(:,1);  % in ms
% rip_endtime = 1000*riptimes(:,2);  % in ms
%
% % Find ripples separated by at least a second -- sj_getpopulationevents2
%
% iri = diff(rip_starttime);
% keepidx = [1;find(iri>=1000)+1];
%
% riplength = rip_endtime - rip_starttime;
% keepidx2 = find(riplength >= 50);% use the ripple last for more than 50 ms
% keepidx = intersect(keepidx,keepidx2);
% %
% rip_starttime = rip_starttime(keepidx);
% riptimes = riptimes(keepidx,:);

% load(sprintf('%s%srippletime_SWS0%d.mat',savedir,animalprefix,day));
load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',savedir,animalprefix,day));

sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
ep2 = find(sleeps(:,2) == ep);

modcells = epochModulation.cellidx;

inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
inhcells(:,3) = -1;
exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
exccells(:,3) = 1;

allmodcells = [inhcells; exccells];

% rip = ripple{day}{ep};
% 
% riptimes = [rip.starttime rip.endtime];
load(sprintf('%s/%srippletime_noncoordSWS%02d.mat', savedir, animalprefix, day));
nc_rip = ripple; clear ripple
load(sprintf('%s/%srippletime_coordSWS%02d.mat', savedir, animalprefix, day));
load(sprintf('%s/%sspikes%02d.mat', savedir, animalprefix, day));
c_rip = ripple; clear ripple

nc_riptimes = [nc_rip{day}{ep}.starttime nc_rip{day}{ep}.endtime];
nc_riptimes(:,3) = 1;

c_riptimes = [c_rip{day}{ep}.starttime c_rip{day}{ep}.endtime];
c_riptimes(:,3) = 2;

%combine riptimes
riptimes = sortrows([nc_riptimes; c_riptimes],1);

%
rip_starttime = 1000*riptimes(:,1);  % in ms
rip_endtime = 1000*riptimes(:,2);  % in ms
%
% % Find ripples separated by at least a second -- sj_getpopulationevents2
%
% iri = diff(rip_starttime);
% keepidx = [1;find(iri>=100)+1];
%
riplength = rip_endtime - rip_starttime;
keepidx = find(riplength >= 50);% use the ripple last for more than 100 ms
% keepidx = intersect(keepidx,keepidx2);
% %
riptimes = riptimes(keepidx,:);

ripnum = size(riptimes,1);

%%
PCCmat = [];
PCCmatall = [];
pMats = [];
if ripnum > 1
    celldata = [];
    spikecounts = [];
    for cellcount = 1:hpnum %get spikes for each cell
        index = [day,ep,cellidxm(cellcount,:)] ;
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
        spikecounts = [spikecounts; spikecount]; %concatenating num spikes per cell, per event
    end
    cellcounts = sum((spikecounts > 0));
    eventindex = find(cellcounts >= cellcountthresh); %Find event indices that have more than cellthresh # of cells
    if ~isempty(eventindex)
        for event = 1:length(eventindex)
            event
            rippletype = riptimes(eventindex(event),3);
            PCCmattmp = [];
            tmpind = find(celldata(:,2) == eventindex(event));
            spiketimes = celldata(tmpind,1); %spike times during ripple event
            cellindex = celldata(tmpind,3);
            cellsi = celldata(find(celldata(:,2)==eventindex(event)),3);
            
            [cellsi,ia] = unique(cellsi,'first');
            [~,sortorder] = sort(ia);
            event_cellSeq = cellsi(sortorder); %sort in order of spikes
            cellsactive = cellidxm(cellsi,:);
            %-----create the event matrix during SWRs (spkT{cells}.spiketimes) -----%
            for cell = event_cellSeq'
                validspikeidx = find(cellindex == cell);
                spkT{cell} = spiketimes(validspikeidx).*1000;
            end
            %--------calculate the posterior probability on an event-by-event basis--------%
            startevent = riptimes(eventindex(event),1).*1000;
            endevent = riptimes(eventindex(event),2).*1000;
            timebins = startevent:tBinSz:endevent; % timebins are the binedges
            nTBin = length(timebins)-1;
            nCell = hpnum;
            spkPerBin = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
            for nn  = 1:hpnum
                cellInd = nn; %current cell
                if length(spkT) >= cellInd
                    if ~isempty(spkT{cellInd})
                        %                     spkPerBin(1,:,cellInd) = histcounts(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
                        temp = histc(spkT{cellInd}, timebins); %[1 x nTBin x nCell] Number of spikes per timebing per cell
                        spkPerBin(1,:,cellInd) = temp(1:end-1);
                    end
                end
            end
            
            nSpkPerTBin = squeeze(sum(spkPerBin,3)); %[nTBin x 1] number of spikes in tBin
            nPBin = size(rm,1); %N positional bin
            expecSpk  = reshape(expecSpk, [nPBin,1, nCell]); %[nPos x 1 x nCell]
            expon = exp(-expecSpk); %Exponent of equation.
            factSpkPerBin = factorial(spkPerBin); %Factorial to divide by
            
            wrking = bsxfun(@power, expecSpk, spkPerBin); %[nPos x nTbin x nCell]
            wrking = bsxfun(@rdivide, wrking, factSpkPerBin); %[nPos x nTbin x nCell]
            wrking = bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
            post = prod(wrking,3); %Non normalised prob [nPos x Tbin]
            post(:,nSpkPerTBin==0) = 0; % FO changed this 27/5/15 so the posterior matrix can be smoothed.
            post(isnan(post)) = 0;
            trajinfo = mean(tm,2);% trajactory number
            
            numShuf = 1000;
            shufDecode = [];
            for s = 1:numShuf
                alltrajshuf = [];
                for tt = 1:4
                    celltrajshuf = [];
                    trajidx = find(trajinfo == tt);
                    thisTraj = rm(trajidx,:);
                    for c = 1:length(thisTraj(1,:))
                        shiftval = randi(length(thisTraj(:,1)));
                        tmp = circshift(thisTraj(:,c),shiftval);
                        celltrajshuf = [celltrajshuf tmp];
                    end
                    alltrajshuf = [alltrajshuf; celltrajshuf];
                end
                
                rm_s = alltrajshuf; %update rm_s with shuffled fields
                
                expecSpk_s = rm_s.*tBinSz./1000;
                expecSpk_s  = reshape(expecSpk_s, [nPBin,1, nCell]); %[nPos x 1 x nCell]
                expon_s = exp(-expecSpk_s); %Exponent of equation.
                factSpkPerBin = factorial(spkPerBin); %Factorial to divide by
                
                wrking_s = bsxfun(@power, expecSpk_s, spkPerBin); %[nPos x nTbin x nCell]
                wrking_s = bsxfun(@rdivide, wrking_s, factSpkPerBin); %[nPos x nTbin x nCell]
                wrking_s = bsxfun(@times,wrking_s, expon); %[nPos x nTbin x nCell]
                post_s = prod(wrking_s,3); %Non normalised prob [nPos x Tbin]
                post_s(:,nSpkPerTBin==0)  =0; % FO changed this 27/5/15 so the posterior matrix can be smoothed.
                post_s(isnan(post_s)) = 0;
                trajinfo = mean(tm,2);% trajactory number
                shufDecode{s} = post_s;
            end
            
            %shuffle the linfield of cells active during this ripple
            for ii = 1:length(cellsactive(:,1))
                shufcell = cellsactive(ii,[1 2]);
                
                [b cellLoc] = ismember(shufcell, cellidxm, 'rows', 'legacy');
                %Shuffling each trajectory independently
                cellshufmat = [];
                cellLin = rm(:,cellLoc);
                for shiftNum = 1:1000
                    sepTrajShuf = [];
                    for tr = 1:4
                        trajidx = find(trajinfo == tr);
                        cellLin_tr = cellLin(trajidx);
                        shiftval = randi(length(trajidx),1);
                        dirshift = randi(2,1);
                        if dirshift == 1
                            dirshift = 1;
                        else
                            dirshift = -1;
                        end
                        shuftr = circshift(cellLin_tr,dirshift*shiftval);
                        sepTrajShuf = [sepTrajShuf; shuftr];
                    end
                    cellshufmat(:,shiftNum) = sepTrajShuf;
                end
                
                cellExpecSpk_shuf = cellshufmat;
                
                shufDecodeCell = [];
                for s = 1:length(cellExpecSpk_shuf(1,:))
                    expecSpk_ss = expecSpk;
                    expecSpk_ss(:,cellLoc) = cellExpecSpk_shuf(:,s);
                    expecSpk_ss  = reshape(expecSpk_ss, [nPBin,1, nCell]); %[nPos x 1 x nCell]
                    expon_s = exp(-expecSpk_ss); %Exponent of equation.
                    factSpkPerBin = factorial(spkPerBin); %Factorial to divide by
                    
                    wrking_s = bsxfun(@power, expecSpk_ss, spkPerBin); %[nPos x nTbin x nCell]
                    wrking_s = bsxfun(@rdivide, wrking_s, factSpkPerBin); %[nPos x nTbin x nCell]
                    wrking_s = bsxfun(@times,wrking_s, expon); %[nPos x nTbin x nCell]
                    post_s = prod(wrking_s,3); %Non normalised prob [nPos x Tbin]
                    post_s(:,nSpkPerTBin==0)  =0; % FO changed this 27/5/15 so the posterior matrix can be smoothed.
                    post_s(isnan(post_s)) = 0;
                    trajinfo = mean(tm,2);% trajactory number
                    shufDecodeCell{s} = post_s;
                end
                trajtmp = [];
                pMat_traj = [];
                wc_traj = [];
                replay_dir = [];
                for traj = 1:4 % four trajactory
                    trajidx = find(trajinfo == traj);
                    pMat = post(trajidx,:);% create a posterior matrix for each traj (split post by traj)
                    
                    distvector = pm(trajidx,1)';
                    szPM1 = size(pMat,1);
                    szPM2 = size(pMat,2);
                    
                    for i = 1:szPM2
                        if (sum(pMat(:,i))>0)
                            pMat(:,i) = pMat(:,i)./sum(pMat(:,i)); % normalized across positions to 1 for each time bin
                        end
                    end
                    
                    pMat_traj{traj} = pMat;
                    %% Weighted correlation calculation
                    positionVec = distvector;
                    timevec = 1:szPM2;
                    mloc = sum(sum(pMat.*positionVec'))/sum(sum(pMat));
                    mt = sum(sum(pMat.*timevec))/sum(sum(pMat));
                    dloc = positionVec'-mloc;
                    dt = timevec - mt;
                    cov_loc_t = sum(sum(pMat.*dloc.*dt))/sum(sum(pMat));
                    cov_loc = sum(sum(pMat.*(dloc.^2)))/sum(sum(pMat));
                    cov_t = sum(sum(pMat.*(dt.^2)))/sum(sum(pMat));
                    weighted_corr = cov_loc_t/sqrt(cov_loc*cov_t);
                    wc_traj(traj) = weighted_corr;
                    
                    rvalues_s = [];
                    for shufMat = 1:length(shufDecode)
                        pMat_s = shufDecode{shufMat};
                        pMat_s = pMat_s(trajidx,:);
                        
                        for i = 1:szPM2
                            if (sum(pMat_s(:,i))>0)
                                pMat_s(:,i) = pMat_s(:,i)./sum(pMat_s(:,i)); % normalized across positions to 1 for each time bin
                            end
                        end
                        
                        mloc_s = sum(sum(pMat_s.*positionVec'))/sum(sum(pMat_s));
                        mt_s = sum(sum(pMat_s.*timevec))/sum(sum(pMat_s));
                        dloc_s = positionVec'-mloc_s;
                        dt_s = timevec - mt_s;
                        cov_loc_t_s = sum(sum(pMat_s.*dloc_s.*dt))/sum(sum(pMat_s));
                        cov_loc_s = sum(sum(pMat_s.*(dloc_s.^2)))/sum(sum(pMat_s));
                        cov_t_s = sum(sum(pMat_s.*(dt_s.^2)))/sum(sum(pMat_s));
                        weighted_corr_s = cov_loc_t_s/sqrt(cov_loc_s*cov_t_s);
                        rvalues_s = [rvalues_s; weighted_corr_s];
                    end
                    
                    rvalues_ss = [];
                    for shufMat = 1:length(shufDecodeCell)
                        pMat_s2 = shufDecodeCell{shufMat};
                        pMat_s2 = pMat_s2(trajidx,:);
                        
                        for i = 1:szPM2
                            if (sum(pMat_s2(:,i))>0)
                                pMat_s2(:,i) = pMat_s2(:,i)./sum(pMat_s2(:,i)); % normalized across positions to 1 for each time bin
                            end
                        end
                        
                        mloc_s = sum(sum(pMat_s2.*positionVec'))/sum(sum(pMat_s2));
                        mt_s = sum(sum(pMat_s2.*timevec))/sum(sum(pMat_s2));
                        dloc_s = positionVec'-mloc_s;
                        dt_s = timevec - mt_s;
                        cov_loc_t_s = sum(sum(pMat_s2.*dloc_s.*dt))/sum(sum(pMat_s2));
                        cov_loc_s = sum(sum(pMat_s2.*(dloc_s.^2)))/sum(sum(pMat_s2));
                        cov_t_s = sum(sum(pMat_s2.*(dt_s.^2)))/sum(sum(pMat_s2));
                        weighted_corr_s = cov_loc_t_s/sqrt(cov_loc_s*cov_t_s);
                        rvalues_ss = [rvalues_ss; weighted_corr_s];
                    end
                    
                    rZ = (abs(weighted_corr) - abs(nanmean(rvalues_s)))/nanstd(abs(rvalues_s));
                    rZ_cellshuf = (abs(weighted_corr) - abs(nanmean(rvalues_ss)))/nanstd(abs(rvalues_ss));
                    
                    numcellsactive = length(cellsactive(:,1));
                    
                    PCC = (rZ - rZ_cellshuf)*numcellsactive;

                    trajRzs(traj) = rZ;
                    
                    if weighted_corr < 0
                        pvalue = sum(weighted_corr > rvalues_s)/length(rvalues_s);
                        dir = -1;
                        replaydir(traj) = dir;
                    elseif weighted_corr > 0
                        pvalue = sum(weighted_corr < rvalues_s)/length(rvalues_s);
                        dir = 1;
                        replaydir(traj) = dir;
                    end
                    
                    trajpvalues(traj) = pvalue;
                    
                    tmp = [PCC traj shufcell rZ rZ_cellshuf];
                    
                    eventtime = riptimes(eventindex(event),[1 2]);
                    tmp = [tmp event];
                    trajtmp = [trajtmp; tmp];
                end
                PCCmatall = [PCCmatall; trajtmp];
            end
            pMats{event}.alltrajmats = pMat_traj;
            pMats{event}.eventnumber = event;
            pMats{event}.eventindex = eventindex(event);
            pMats{event}.traj_pvals = trajpvalues;
            pMats{event}.dir = replaydir;
            pMats{event}.weightedcorr = wc_traj;
            pMats{event}.cellsactive = cellsactive;
            pMats{event}.sequenceScore = trajRzs;
            pMats{event}.spksPerCell = spkT;%cellidxm(cellsi,:); %CHANGE THIS
            pMats{event}.eventtime = eventtime;
            pMats{event}.rippletype = rippletype;
        end
        PCC_events.pMats = pMats;
        PCC_events.alleventsAllcell = PCCmatall;
        PCC_events.cellidx = cellidxm;
        PCC_events.modcells = allmodcells;
        PCC_events.allriptimes = riptimes;
        PCC_events.descript = 'PCC, trajectory, shuff cell, rZ, rZshuf, event';
        %%
        if savedata
            save(sprintf('%s%sPCCreplaydecode_alldata_50msRipNoIriCritTbin%02d_sequenceScores%02d_%02d.mat', savedir,animalprefix,tBinSz,day,ep), 'PCC_events');
        end
    end
end