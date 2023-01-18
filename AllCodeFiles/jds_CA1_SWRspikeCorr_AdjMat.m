function out = jds_CA1_SWRspikeCorr_AdjMat(animalprefixlist)

day = 1; %always single day expt
%%
epochs = [1:2:17];
inh_con = [];
exc_con = [];
nm_con = [];
inh_pctcells = [];
exc_pctcells = [];
inh_spks = [];
exc_spks = [];
compileModmat = [];
cellcnt = [];
cellModValues = [];
for a = 1:length(animalprefixlist)
    
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
    % get ripple time
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',dir,animalprefix,day));
    cellperm = epochModulation.modMat(randperm(length(epochModulation.modMat(:,1))),:);
    compileModmat = [compileModmat; cellperm];
    cellcnt(a) = length(epochModulation.modMat(:,1));
    for e = 1:length(epochs)
        ep = epochs(e);
        rip = ripple{day}{ep};
        riptimes(:,1) = rip.starttime;
        riptimes(:,2) = rip.endtime;

        [ctxidx, hpidx] = jds_getallepcells(dir, animalprefix, day, ep, [])
        ctxnum = length(ctxidx(:,1));
        hpnum = length(hpidx(:,1));
        
        if length(riptimes(:,1)) > 5
            celldata = [];
            CA1matrix = [];
            for cellcount = 1:hpnum %get spikes for each cell
                index = [day,ep,hpidx(cellcount,:)] ;
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
                CA1matrix = [CA1matrix spikecount']; %concatenating num spikes per cell, per event
            end
        end
        if (isequal(animalprefix,'ZT2')) && (ep == 11)
            figure; imagesc(CA1matrix); colormap(inferno)
            colorbar
            title(['Example SWR matrix-' animalprefix '-ep' num2str(ep)])
        end
        
        %%
        %Get correlations between cells during specified times
        corridxs = zeros((hpnum*hpnum),6); rcorr = []; pcorr = []; rip_co_occur = [];
        corlnMat = zeros(hpnum,hpnum);
        allCellSWRrespmat = CA1matrix;
        cnt = 1;
        for n=1:hpnum
            CA1_cell1 = CA1matrix(:,n);
            %             rip_co_occur{n}.descipt = 'num_SWRco_occur, hptet1, hpcell1, hptet2, hpcell2';
            %         allCellSWRrespmat = [allCellSWRrespmat PFC_cell];
            for nn=1:hpnum
                if n ~= nn
                    CA1_cell2 = CA1matrix(:,nn);
                    [r, p] = corrcoef(CA1_cell1,CA1_cell2);
                    rcorr = [rcorr; r(1,2)]; %partial correlation
                    pcorr = [pcorr; p(1,2)];
                    corridxs(cnt,:) = [day ep hpidx(n,1) hpidx(n,2) hpidx(nn,1) hpidx(nn,2)];
                    if p(1,2) < 0.025
                        corlnMat(n,nn) = r(1,2);
                    end
                    cnt = cnt + 1;
                end
            end
        end
        corlnMat(find(corlnMat < 0)) = 0;
        corlnMat(find(isnan(corlnMat))) = 0;

        if (isequal(animalprefix,'ZT2')) && (ep == 11)
            figure; imagesc(corlnMat); colormap(inferno)
            colorbar
            title(['Example SWR Corln matrix-' animalprefix '-ep' num2str(ep)])
        end
        
        co_occur = [];
        for l = 1:length(allCellSWRrespmat(1,:))
            c1Mat = allCellSWRrespmat(:,l);
            cnt = 1;
            for ll = 1:length(allCellSWRrespmat(1,:))
                c2Mat = allCellSWRrespmat(:,ll);
                co_occur{l}.data(cnt,1) = length(find((c1Mat~=0) & (c2Mat~=0)));
                totalnum = length(find((c1Mat~=0) & (c2Mat~=0)));
                nospk1idx = c1Mat~=0;
                nospk2idx = c2Mat~=0;
                both_nospk = nospk1idx + nospk2idx;
                nospk_overlap = length(find(both_nospk == 0));
                percent_cofire = (totalnum/(length(riptimes) - nospk_overlap))*100;
                co_occur{l}.data(cnt,2) = percent_cofire;
                cnt = cnt + 1;
            end
        end
        
%         adjMat = [];
%         for t = 1:hpnum
%             for tt = 1:length(co_occur)
%                 if (co_occur{t}.data(tt,2) >= 20) && (co_occur{t}.data(tt,2) ~= 100); %floor(length(riptimes)*0.1)
%                     %             adjMat(t,tt) = co_occur{t}.data(tt,2)*.01;
%                     adjMat(t,tt) = 1;
%                 else
%                     adjMat(t,tt) = 0;
%                 end
%             end
%         end
        adjMat = corlnMat;
        adjMat(find(adjMat > 0)) = 1;
        
        %calculate the degree matrix degMat
        degMat = zeros(hpnum);
        for g = 1:(hpnum)
            degMat(g,g) = sum(adjMat(g,:));
        end
        myColorMap = lines(length(corlnMat));
        
        myLabel = cell(length(corlnMat),1);
        indices = [];
        for i = 1:hpnum
            cellidx = hpidx(i,:);
            [b,cellLoc] = ismember(cellidx, epochModulation.cellidx,'rows','legacy');
            if cellLoc ~= 0
                cellmod = epochModulation.modMat(cellLoc,e);
                cellmodval = epochModulation.modVals(cellLoc,e);
                cellModValues = [cellModValues; cellmodval];
                if cellmod == 1
                    myLabel{i} = 'EXC';
                    indices(i) = 1;
                elseif cellmod == -1
                    myLabel{i} = 'INH';
                    indices(i) = -1;
                elseif cellmod == 0
                    myLabel{i} = 'NM';
                    indices(i) = 0;
                end
            else
                myLabel{i} = 'NM';
                indices(i) = 0;
            end
        end
        excidx = find(indices == 1);
        inhidx = find(indices == -1);
        nmidx = find(indices == 0);
        numCon = sum(degMat);
        
        if ~isempty(excidx)
            tmp = numCon(excidx);
            tmp = tmp./hpnum;
            exc_con = [exc_con; tmp'];
            for thiscell = 1:length(excidx)
                tmpidx = find(CA1matrix(:,excidx(thiscell)) ~= 0);
                tmpmatrix = CA1matrix(tmpidx,:);
                spkvec = CA1matrix(tmpidx,excidx(thiscell));
                tmpmatrix2 = tmpmatrix;
                tmpmatrix2 = sum(tmpmatrix2,2);
                tmpmatrix(find(tmpmatrix > 0)) = 1;
                tmpmatrix = sum(tmpmatrix,2);
                pctactive = mean(tmpmatrix./hpnum);
%                 normspks = mean(tmpmatrix2./tmpmatrix);
                nummspks = mean(spkvec);
                exc_pctcells = [exc_pctcells; pctactive];
                exc_spks = [exc_spks; nummspks];
            end
        end
        
        if ~isempty(inhidx)
            tmp = numCon(inhidx);
            tmp = tmp./hpnum;
            inh_con = [inh_con; tmp'];
            for thiscell = 1:length(inhidx)
                tmpidx = find(CA1matrix(:,inhidx(thiscell)) ~= 0);
                tmpmatrix = CA1matrix(tmpidx,:); %all active cells
                spkvec = CA1matrix(tmpidx,inhidx(thiscell));
                tmpmatrix2 = tmpmatrix;
                tmpmatrix2 = sum(tmpmatrix2,2); %number of spikes in each ripple summed
                tmpmatrix(find(tmpmatrix > 0)) = 1; 
                tmpmatrix = sum(tmpmatrix,2); %number cells active in each ripple
                pctactive = mean(tmpmatrix./hpnum); %mean % active cells active when this cell is active
%                 normspks = mean(tmpmatrix2./tmpmatrix); %mean number of spikes per cell EXCLUDE THIS
                numspks = mean(spkvec);
                inh_pctcells = [inh_pctcells; pctactive];
                inh_spks = [inh_spks; numspks];
            end
        end
        
        if ~isempty(nmidx)
            tmp = numCon(nmidx);
            tmp = tmp./hpnum;
            nm_con = [nm_con; tmp'];
        end
        
        figure;
%         circularGraph(corlnMat,'Colormap',myColorMap,'Label',myLabel)
        circularGraph(adjMat,'Colormap',myColorMap,'Label',myLabel)
        keyboard
        close

        clear riptimes
    end
end
keyboard
inh_mean = mean(inh_con);
exc_mean = mean(exc_con);

datameans = [exc_mean inh_mean];
datasems = [(nanstd(exc_con)/sqrt(length(find(~isnan(exc_con)))))...
    (nanstd(inh_con)/sqrt(length(find(~isnan(inh_con)))))];

figure
bar([1:2],datameans(1:2),'k')
hold on
er = errorbar([1:2],datameans(1:2),datasems(1:2));
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('% Cells Sig Corr')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)


[p h] = ranksum(exc_con,inh_con)

figure;
datacombinedSWRSigCorr = [exc_con; inh_con];
g1 = repmat({'CA1exc'},length(exc_con),1);
g2 = repmat({'CA1inh'},length(inh_con),1);
g = [g1;g2];

h = boxplot(datacombinedSWRSigCorr,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-6 10])
title(['%Cells significantly corr-p = ' num2str(p)])

datameans2 = [mean(exc_pctcells) mean(inh_pctcells)];
datasems2 = [(nanstd(exc_pctcells)/sqrt(length(find(~isnan(exc_pctcells)))))...
    (nanstd(inh_pctcells)/sqrt(length(find(~isnan(inh_pctcells)))))];

figure
bar([1:2],datameans2(1:2),'k')
hold on
er = errorbar([1:2],datameans2(1:2),datasems2(1:2));
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('% Cells Active')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

datameans3 = [mean(exc_spks) mean(inh_spks)];
datasems3 = [(nanstd(exc_spks)/sqrt(length(find(~isnan(exc_spks)))))...
    (nanstd(inh_spks)/sqrt(length(find(~isnan(inh_spks)))))];

figure
bar([1:2],datameans3(1:2),'k')
hold on
er = errorbar([1:2],datameans3(1:2),datasems3(1:2));
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Avg. Spikes/Cell')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

figure; imagesc(compileModmat); colormap(gray)
title(num2str(cellcnt))

figure;
cellModValues = cellModValues(find(cellModValues ~= 0));
cellModValues = sort(cellModValues,'descend')
bar(cellModValues); view(90,90);

posmodvals = cellModValues(find(cellModValues > 0));
negmodvals = abs(cellModValues(find(cellModValues < 0)));

[p2 h2] = ranksum(posmodvals,negmodvals)
datacombinedmodVals = [posmodvals; negmodvals];
g1 = repmat({'EXC'},length(posmodvals),1);
g2 = repmat({'INH'},length(negmodvals),1);
g = [g1;g2];

% boxplot(datacombinedCoact,g);
figure
h = boxplot(datacombinedmodVals,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.2 1.5])
title(['Abs. modulation index-p = ' num2str(p2)])
keyboard;





   



