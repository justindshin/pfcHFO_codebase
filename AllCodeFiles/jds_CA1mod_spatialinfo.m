function jds_CA1mod_spatialinfo(animalprefixlist)
%---------------------------------------------------------------%
%This function calculates the spatial infor for CA1 cells that
%are modulated during PFC ripples - Compares CA1 exc and inh cells
%---------------------------------------------------------------%

%%
spatialInfoSpk_exc = [];
spatialInfoSpk_inh = [];
spatialInfoSec_exc = [];
spatialInfoSec_inh = [];
meanFR_exc = [];
meanFR_inh = [];
meanFRep_exc = [];
meanFRep_inh = [];
day = 1;
epochs = [2:2:16];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    %%
    %-----match neurons across epochs-----%
    for e = 1:length(epochs)
        ep = epochs(e);
        [ctxidx, hpidx] =  jds_getallepcells(dir, animalprefix, day, ep, []); %(tet, cell)
        hpnum = length(hpidx(:,1));
        
        load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',dir,animalprefix,day));
        
        sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

        if pre == 1
            ep_sleep = ep + 1;
        elseif post == 1
            ep_sleep = ep - 1;
        end
        
        ep2 = find(sleeps(:,2) == ep_sleep);
        modcells = epochModulation.cellidx;
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = 1;
        
        allmodcells = [inhcells; exccells];
        
        %%
        if ~isempty(hpidx)
            hpidx = hpidx(:,[1 2]);
            hpnum = length(hpidx(:,1));
            %-----create the ratemaps [nPosBin x nHPCells]-----%
            rm = []; % ratemap matrix
            oc = [];
            cellidxm = [];
            
            load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
            load(sprintf('%s%sspikes0%d.mat',dir,animalprefix,day));
            load(sprintf('%s%sthetatime0%d.mat',dir,animalprefix,day));
            for i = 1:hpnum
                cind = hpidx(i,:);
                if (length(linfields{day}{ep})>= cind(1))
                    if  (length(linfields{day}{ep}{cind(1)})>= cind(2))
                        linfield1 = linfields{day}{ep}{cind(1)}{cind(2)};
                    else
                        linfield1 =[];
                    end
                else
                    linfield1=[];
                end
                
                if ~isempty(linfield1)
                    linfield = [];
                    occupancy = [];
                    for track = 1:4
                        temp1 = linfield1{track};
                        occnormrate1 = temp1(:,5);
                        occTmp = temp1(:,2);
                        occupancy = [occupancy; occTmp];
                        linfield = [linfield;occnormrate1];
                    end
                    if (max(linfield) >= 3) % peak firing rate max larger than 3 Hz
                        aa = find(isnan(linfield));
                        %pad nan
                        if ~isempty(aa)
                            [lo,hi]= findcontiguous(aa);  %find contiguous NaNs
                            for ii = 1:length(lo)
                                if lo(ii) > 1 & hi(ii) < length(linfield)
                                    fill = linspace(linfield(lo(ii)-1), ...
                                        linfield(hi(ii)+1), hi(ii)-lo(ii)+1);
                                    linfield(lo(ii):hi(ii)) = fill;
                                end
                            end
                        end
                        rm = [rm;linfield'];
                        oc = [oc; occupancy'];
                        cellidxm = [cellidxm; cind];
                    end
                end
            end
            rm = rm'; %[nPosBin x nHPCells]
            oc = oc';
            rm = rm+ (eps.^8); %Add a small number so there are no zeros
            hpnum = length(rm(1,:)); % update
            
            if ep <10
                epochstring = ['0',num2str(ep)];
            else
                epochstring = num2str(ep);
            end
            
            thetalist = [thetatime{day}{ep}.starttime thetatime{day}{ep}.endtime];
            thetadur = thetatime{day}{ep}.total_duration;
            runTime = sum(oc(:,1));
            posProb = oc(:,1)./runTime;
            
            for modcell = 1:length(allmodcells(:,1))
                thisCell = allmodcells(modcell,[1 2]);
                cellmod = allmodcells(modcell,3);
                [b cellidx] = ismember(thisCell, cellidxm, 'rows', 'legacy');
                spkInfo = [];
                secInfo = [];
                if cellidx ~= 0
                    cellSpks = spikes{day}{ep}{thisCell(1)}{thisCell(2)}.data(:,1);
                    thetabins = periodAssign(cellSpks,thetalist);
                    thetaSpks = length(find(thetabins ~= 0));
                    meanFR = thetaSpks/thetadur;
                    meanFRep = spikes{day}{ep}{thisCell(1)}{thisCell(2)}.meanrate;
%                     meanFR = meanFRep;
                    linFR = rm(:,cellidx);
                    for bin = 1:length(linFR)
                        spkInfoTmp = posProb(bin)*(linFR(bin)/meanFR)*log2(linFR(bin)/meanFR);
                        
                        secInfoTmp = posProb(bin)*linFR(bin)*log2(linFR(bin)/meanFR);
                        
                        spkInfo = [spkInfo; spkInfoTmp];
                        secInfo = [secInfo; secInfoTmp];
                    end
                    if cellmod == 1
                        spatialInfoSpk_exc = [spatialInfoSpk_exc; nansum(spkInfo)];
                        spatialInfoSec_exc = [spatialInfoSec_exc; nansum(secInfo)];
                        meanFR_exc = [meanFR_exc; meanFR];
                        meanFRep_exc = [meanFRep_exc; meanFRep];
                    elseif cellmod == -1
                        spatialInfoSpk_inh = [spatialInfoSpk_inh; nansum(spkInfo)];
                        spatialInfoSec_inh = [spatialInfoSec_inh; nansum(secInfo)];
                        meanFR_inh = [meanFR_inh; meanFR];
                        meanFRep_inh = [meanFRep_inh; meanFRep];
                    end
                end
            end
        end
    end
end

datameans = [mean(spatialInfoSpk_exc) mean(spatialInfoSpk_inh)];
datasems = [(std(spatialInfoSpk_exc)/sqrt(length(spatialInfoSpk_exc)))...
    (std(spatialInfoSpk_inh)/sqrt(length(spatialInfoSpk_inh)))];

figure
bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Spatial Info (Bits/Spike)')

xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p1 h1] = ranksum(spatialInfoSpk_exc,spatialInfoSpk_inh)

datacombinedInfoSpk = [spatialInfoSpk_exc; spatialInfoSpk_inh];
g1 = repmat({'EXC'},length(spatialInfoSpk_exc),1);
g2 = repmat({'INH'},length(spatialInfoSpk_inh),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedInfoSpk,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.02 0.2])
title(['CA1 Spatial Info/spk-p = ' num2str(p1)])

datameans2 = [mean(spatialInfoSec_exc) mean(spatialInfoSec_inh)];
datasems2 = [(std(spatialInfoSec_exc)/sqrt(length(spatialInfoSec_exc)))...
    (std(spatialInfoSec_inh)/sqrt(length(spatialInfoSec_inh)))];

figure
bar([1:2],datameans2,'k')
hold on
er = errorbar([1:2],datameans2,datasems2);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Spatial Info (Bits/Second)')

xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p2 h2] = ranksum(spatialInfoSec_exc,spatialInfoSec_inh)

datacombinedInfoSec = [spatialInfoSec_exc; spatialInfoSec_inh];
g1 = repmat({'EXC'},length(spatialInfoSec_exc),1);
g2 = repmat({'INH'},length(spatialInfoSec_inh),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedInfoSec,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-2 12])
title(['CA1 Spatial Info/sec-p = ' num2str(p2)])

keyboard
