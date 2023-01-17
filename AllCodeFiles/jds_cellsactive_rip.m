function jds_cellsactive_rip(animalprefixlist)


%%
%Things to incorporate
%-Need to restrict to cells that pass specific criterion (Active in
%more than 10 SWR events (Rothschild 2017)
%%
% savedata = 0;

PFCmatrix_anim = [];
CA1matrix_anim = [];
PFC_pct_active_coord = [];
PFC_pct_active_noncoord = [];
CA1_pct_active_coord = [];
CA1_pct_active_noncoord = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
 
    %%
    
    eps = [1:2:17];
    
    
    %-----match neurons across epochs-----%
    day = 1; %always single day expt

    %-----create the event matrix during SWRs-----%
    spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
    
    % get ripple time
    load(sprintf('%s/%srippletime_noncoordSWS%02d.mat', dir, animalprefix, day));
    nc_rip = ripple; clear ripple
    load(sprintf('%s/%srippletime_coordSWS%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
    c_rip = ripple; clear ripple
    
    for ep = 1:length(eps)
        PFCmatrix = [];
        CA1matrix = [];
        e = eps(ep);
        disp([animalprefix '-epoch-' num2str(e)])

        [ctxidx, hpidx] = jds_getallepcells(dir, animalprefix, day, e, []); %(tet, cell)
        ctxnum = length(ctxidx(:,1));
        hpnum = length(hpidx(:,1));

        nc_riptimes = [nc_rip{day}{e}.starttime nc_rip{day}{e}.endtime];
        nc_riptimes(:,3) = 1;

        c_riptimes = [c_rip{day}{e}.starttime c_rip{day}{e}.endtime];
        c_riptimes(:,3) = 2;

        %combine riptimes
        riptimes = sortrows([nc_riptimes; c_riptimes],1);
        
        if ~isempty(riptimes)
            
            pfc_matrix_tmp = [];
            if length(riptimes(:,1)) > 1
                for cellcount = 1:ctxnum %get spikes for each cell
                    index = [day,e,ctxidx(cellcount,:)] ;
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
                    end
                    spikecount = zeros(1,size(riptimes,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    pfc_matrix_tmp = [pfc_matrix_tmp spikecount'];
                end
                matrixtmp = pfc_matrix_tmp;
                numcells = length(matrixtmp(1,:));
                for rip = 1:length(matrixtmp(:,1))
                    thisripple = matrixtmp(rip,:);
                    rtype = riptimes(rip,3);
                    numactive = length(find(thisripple ~= 0));
                    pct_active = (numactive/numcells)*100;
                    if rtype == 2
                        PFC_pct_active_coord = [PFC_pct_active_coord; pct_active];
                    elseif rtype == 1
                        PFC_pct_active_noncoord = [PFC_pct_active_noncoord; pct_active];
                    end
                end
            end

            %GET CA1 CELL DATA
            ca1_matrix_tmp = [];
            if length(riptimes(:,1)) > 1
                for cellcount = 1:hpnum %get spikes for each cell
                    index = [day,e,hpidx(cellcount,:)] ;
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
                    end
                    spikecount = zeros(1,size(riptimes,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    ca1_matrix_tmp = [ca1_matrix_tmp spikecount'];
                end
                matrixtmp = ca1_matrix_tmp;
                numcells = length(matrixtmp(1,:));
                for rip = 1:length(matrixtmp(:,1))
                    thisripple = matrixtmp(rip,:);
                    rtype = riptimes(rip,3);
                    numactive = length(find(thisripple ~= 0));
                    pct_active = (numactive/numcells)*100;
                    if rtype == 2
                        CA1_pct_active_coord = [CA1_pct_active_coord; pct_active];
                    elseif rtype == 1
                        CA1_pct_active_noncoord = [CA1_pct_active_noncoord; pct_active];
                    end
                end            
            end
        end
    end
end

[p1 h1] = ranksum(PFC_pct_active_noncoord,PFC_pct_active_coord)
datacombinedPctActivepfc = [PFC_pct_active_noncoord; PFC_pct_active_coord];
g1 = repmat({'Ind'},length(PFC_pct_active_noncoord),1);
g2 = repmat({'Coord'},length(PFC_pct_active_coord),1);
g = [g1;g2];

figure; 
h = boxplot(datacombinedPctActivepfc,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['PFC - %cells active-p = ' num2str(p1)])
ylabel('%cells active')
set(gcf, 'renderer', 'painters')
ylim([-5 60])

[p2 h2] = ranksum(CA1_pct_active_noncoord,CA1_pct_active_coord)
datacombinedPctActiveca1 = [CA1_pct_active_noncoord; CA1_pct_active_coord];
g1 = repmat({'Ind'},length(CA1_pct_active_noncoord),1);
g2 = repmat({'Coord'},length(CA1_pct_active_coord),1);
g = [g1;g2];

figure; 
h = boxplot(datacombinedPctActiveca1,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['CA1 - %cells active-p = ' num2str(p2)])
ylabel('%cells active')
set(gcf, 'renderer', 'painters')
ylim([-5 50])
keyboard

