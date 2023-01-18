function jds_PPC_thetalockingstrength(animalprefixlist)

savedata = 1;
day = 1;

daystring = '01';

%%
epochs = [2:2:16];

PPC_Exc = [];
PPC_Inh = [];
PPC_ExcKappa = [];
PPC_InhKappa = [];

pre = 1;
post = 0;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    %%
    %-----create the event matrix during SWRs-----%
    load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day)); % get spikes
    load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sthetatime%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    modcells = epochModulation.cellidx;
    
    tets = tetinfo{1}{epochs(1)}; 
    
    reftet = []; 
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'CA1Ref')
                reftet = [reftet; t];
            end
        end
    end
    
    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
    
    for ep = 1:length(epochs)

        if pre == 1
            eps = [epochs(ep) (epochs(ep) + 1)];
        elseif post == 1
            eps = [epochs(ep) (epochs(ep) - 1)];
        end
        epsleep = eps(2);
        eprun = eps(1);

        ep2 = find(sleeps(:,2) == epsleep);
        modvals = epochModulation.modVals(:,ep2);
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = modvals(find(epochModulation.modMat(:,ep2) == -1));
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = modvals(find(epochModulation.modMat(:,ep2) == 1));

        allmodcells = [inhcells; exccells];

        cellnum = length(allmodcells(:,1));
                
        if (eprun <10) && (isequal(animalprefix, 'ZT2'))
            epochstring = ['0',num2str(eprun)];
        else
            epochstring = num2str(eprun);
        end
        
        thetalist = thetatime{day}{eprun};
        thetalist = [thetalist.starttime thetalist.endtime];
        
        for cellcount = 1:cellnum %get spikes for each cell
            cellmod = allmodcells(cellcount, 3);
            index = [day,eprun,allmodcells(cellcount,[1 2])];
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            
            goodspikes = isExcluded(spiketimes, thetalist);
                        
            if (reftet<10)
                reftetstring = ['0',num2str(reftet)];
            else
                reftetstring = num2str(reftet);
            end
            
            curreegfile = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstring,'-',reftetstring];
            load(curreegfile);
            thetagnd_run = thetagnd; clear thetagnd;
            
            phasedata = thetagnd_run{day}{eprun}{reftet}.data(:,2);
            
            t = geteegtimes(thetagnd_run{day}{eprun}{reftet});
            
            if ~isempty(goodspikes)
                sph = phasedata(lookup(spiketimes, t));
                sph = double(sph(logical(goodspikes))) / 10000;  % If no spikes, this will be empty
            else
                sph = [];
            end
            
            if length(sph)>1 
                
                % Rayleigh and Modulation: Originally in lorenlab Functions folder
                stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                [m, ph] = modulation(sph);
                phdeg = ph*(180/pi);
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                thetahat_deg = thetahat*(180/pi);
                
                allphases = sph*(180/pi);
                
                [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                cellPPC = [];
                if prayl < 0.05
                    for p = 1:length(allphases)
                        for p2 = 1:length(allphases)
                            if (p ~= p2) && (p2 > p)
                                phaseDiff = allphases(p2) - allphases(p);
                                %calculates the pairwise phase consistency
                                %(PPC)
                                PPCtmp = cosd(phaseDiff);
                                cellPPC = [cellPPC; PPCtmp];
                            end
                        end
                    end
                    if cellmod == 1
                        PPC_Exc = [PPC_Exc; mean(cellPPC)];
                        PPC_ExcKappa = [PPC_ExcKappa; kappa];
                    elseif cellmod == -1
                        PPC_Inh = [PPC_Inh; mean(cellPPC)];
                        PPC_InhKappa = [PPC_InhKappa; kappa];
                    end
                end
            end
        end
    end
end

meanPCCexc = mean(PPC_Exc)
meanPCCinh = mean(PPC_Inh)

% keyboard;

datameansPPC = [mean(PPC_Exc) mean(PPC_Inh)];
datasemsPPC = [(std(PPC_Exc)/sqrt(length(PPC_Exc)))...
    (std(PPC_Inh)/sqrt(length(PPC_Inh)))];

[p h] = ranksum(PPC_Exc,PPC_Inh)

datacombinedPPC = [PPC_Exc; PPC_Inh];
g1 = repmat({'EXC'},length(PPC_Exc),1);
g2 = repmat({'INH'},length(PPC_Inh),1);
g = [g1;g2];

% boxplot(datacombinedCoact,g);
figure
h = boxplot(datacombinedPPC,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.02 0.2])
title(['Pairwise phase consistency-p = ' num2str(p)])

[p2 h2] = ranksum(PPC_ExcKappa,PPC_InhKappa)

datacombinedKappa = [PPC_ExcKappa; PPC_InhKappa];
g1 = repmat({'EXC'},length(PPC_ExcKappa),1);
g2 = repmat({'INH'},length(PPC_InhKappa),1);
g = [g1;g2];

% boxplot(datacombinedCoact,g);
figure
h = boxplot(datacombinedKappa,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.2 1.2])
title(['Kappa (post)-p = ' num2str(p2)])

kappaMn = [mean(PPC_ExcKappa) mean(PPC_InhKappa)];
kappaSems =  [(std(PPC_ExcKappa)/sqrt(length(PPC_ExcKappa)))...
    (std(PPC_InhKappa)/sqrt(length(PPC_InhKappa)))];

figure
bar([1:2],datameansPPC,'k')
hold on
er = errorbar([1:2],datameansPPC,datasemsPPC);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Pairwise Phase Consistency (PPC)')
title('CA1 Cell Phase Locking Strength (PPC)')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

figure
bar([1:2],kappaMn,'k')
hold on
er = errorbar([1:2],kappaMn,kappaSems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Kappa')
title('CA1 Cell Phase Locking Strength (Kappa)')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

figure;
datacombined = [PPC_Exc; PPC_Inh];
g1 = repmat({'CA1exc'},length(PPC_Exc),1);
g2 = repmat({'CA1inh'},length(PPC_Inh),1);
g = [g1;g2];

boxplot(datacombined,g,'PlotStyle','compact');

keyboard

