%Plots a raster during specified epoch - Includes CA1 and PFC cells, CA1
%and PFC EEG, Ripple Band
clear; close all;
animalprefix='ZT2';
day=1; ep = 13;
prepostRip = 0.05;

% --------------- Parameters ---------------
% Lot of this is not needed if not doing EEG
%-------------------------------------------
Fspos = 30; %Hz
respos = 1/30; % sec
Fseeg = 1500; %Hz
reseeg = 1/1500; % sec
Fsspikes = 10000; %Hz
resspikes = 1/10000; %sec

%%
if strcmp(animalprefix,'ZT2')
    savedir = ('/Volumes/JUSTIN/SingleDay/ZT2_direct/');
    riptetlist = [10 11 12 14 16 17 18 19 29 24 25 27 32 36];  % No Need if no ripples
    maineegtet = 24;  % CA1 tet % No need if no EEG
    peegtet = 30; % PFCtet % No need if no EEG
end

eegtets = maineegtet;
% Also get a PFC eeg
eegtets = [eegtets, peegtet];
peegidx= find(eegtets==peegtet);
maineegidx = find(eegtets==maineegtet);
saveg = 0;


%%
%-----match neurons across epochs-----%
[ctxidx, hpidx] = jds_getallepcells(savedir, animalprefix, day, ep, []); %(tet, cell)
ctxnum = length(ctxidx(:,1));
hpnum = length(hpidx(:,1));

%%
%-----create the event matrix during SWRs-----%
spikes = loaddatastruct(savedir, animalprefix, 'spikes', day); % get spikes

load(sprintf('%s%srippletime_noncoordSWS0%d.mat',savedir,animalprefix,day));
load(sprintf('%s%sctxrippletime_SWS0%d.mat',savedir,animalprefix,day));
load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',savedir,animalprefix,day));

rip = ripple{day}{ep};

riptimes = [(rip.starttime-prepostRip) (rip.endtime+prepostRip)];

triggers = (riptimes(:,1)*1000);
triggers_end = (riptimes(:,2)*1000);
pt = triggers; pt_sec = pt./1000;
pt_end = triggers_end; pt_endsec = pt_end./1000;

ripnum = size(riptimes,1);

EEGfile = sprintf('%s/EEG/%seegref%02d-%02d-%02d.mat', savedir, animalprefix, day,ep,maineegtet);
load(EEGfile);
eeg=eegref;
e = eeg{day}{ep}{maineegtet};

teeg = geteegtimes(e);
eind = lookup(pt, teeg);
e.samprate=round(e.samprate);
eegstart = eeg{day}{ep}{maineegtet}.starttime; % in secs - Epoch start
eegdata = eeg{day}{ep}{maineegtet}.data;
eegend = teeg(end); % in secs - Epoch end

ripfile = sprintf('%s/EEG/%sripple%02d-%02d-%02d.mat', savedir, animalprefix, day,ep,maineegtet);
load(ripfile);
ripamp = ripple{day}{ep}{maineegtet}.data(:,1);

sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
ep2 = find(sleeps(:,2) == ep);

modcells = epochModulation.cellidx;

inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
inhcells(:,3) = -1;
exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
exccells(:,3) = 1;

allmodcells = [inhcells; exccells];
cellcountthresh = 10
%%
if ripnum > 1
    celldata = [];
    spikecounts = [];
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
        spikecounts = [spikecounts; spikecount]; %concatenating num spikes per cell, per event
    end
    cellcounts = sum((spikecounts > 0));
    eventindex = find(cellcounts >= cellcountthresh); %Find event indices that have more than cellthresh # of cells
    if ~isempty(eventindex)
        for event = 1:length(eventindex)
            event
            PCCmattmp = [];
            tmpind = find(celldata(:,2) == eventindex(event));
            spiketimes = celldata(tmpind,1); %spike times during ripple event
            cellindex = celldata(tmpind,3);
            uniqueId = unique(cellindex);
            cellSpksCom = [];
            for c = 1:length(uniqueId)
                idx = find(uniqueId(c) == cellindex);
                com = mean(spiketimes(idx));
                cellSpksCom = [cellSpksCom; [uniqueId(c) com]];
            end
            cellSpksCom = sortrows(cellSpksCom,2);
            cellSpks = [];
            for c = 1:length(cellSpksCom(:,1))
                idx = find(cellSpksCom(c,1) == cellindex);
                spks = spiketimes(idx);
                idvec(1:length(spks)) = cellSpksCom(c,1);
                cellSpks = [cellSpks; [idvec' spks]];
                idvec = [];
            end

%             cellSpksCom = sortrows(cellSpksCom,2);
%             cellSpks = [cellindex spiketimes];
%             cellSpks = sortrows(cellSpks,2);
            spiketimes = cellSpks(:,2);
            cellindex = cellSpks(:,1);
            event_cellSeq = unique(cellSpks(:,1),'stable');
            cellsi = celldata(find(celldata(:,2)==eventindex(event)),3);
            
%             [cellsi,ia] = unique(cellsi,'first');
%             [cellsi,ia] = unique(cellSpksOrdered,'first');
%             [~,sortorder] = sort(ia);
%             event_cellSeq = cellsi(sortorder); %sort in order of spikes
            cellsactive = hpidx(event_cellSeq,:);
            %-----create the event matrix during SWRs (spkT{cells}.spiketimes) -----%
            for cell = event_cellSeq'
                validspikeidx = find(cellindex == cell);
                spkT{cell} = spiketimes(validspikeidx).*1000;
            end
            winst = pt_sec(eventindex(event));
            winend = pt_endsec(eventindex(event));
            win = winend-winst;
            
            % Get eeg times axis - DONT NEED THIS IF NO EEG
            eind1 = lookup(winst, teeg);
            eind2 = lookup(winend, teeg);
            taxisEEG = teeg(eind1:eind2);
            taxisEEG = taxisEEG - winst;
            
            taxis = winst:winend; % IF winst is in secs, this will create an axis in seconds
            taxis = taxis - winst; % Start axis from 0
            
            winst_ms = winst*1000;
            winend_ms = winend*1000;
            

            baseline = 0;
            
            % CA1 spikes
            cnt = 0;
            activeca1cnt = 0;
            for c=1:length(event_cellSeq)
                thiscell = cellsactive(c,:);
                %sig mod?
                idx1 = find(thiscell(:,1) == allmodcells(:,1));
                idx2 = find(thiscell(:,2) == allmodcells(:,2));
                idx3 = intersect(idx1, idx2);
                
                if ~isempty(idx3)
                    cellmod = allmodcells(idx3,3);
                else
                    cellmod = 0;
                end
                
                cellorderidx = event_cellSeq(c);
                currspkt = spkT{cellorderidx}/1000;
                
                figure(1); hold on;
                if ~isempty(currspkt)
                    activeca1cnt = activeca1cnt+1;
                    cnt=cnt+2;
                    currspkt = currspkt - winst;
                    if cellmod == 0
                    plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','k','LineWidth',1);
                    elseif cellmod == 1
                        plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','r','LineWidth',1);
                    elseif cellmod == -1
                        plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','b','LineWidth',1);
                    end
                end
            end
            
            % Update baseline to give a little vertical gap
            baseline = baseline + (activeca1cnt)*4;
            baseline = baseline+1;
            
            % EEG On main tet
            n = maineegidx;
            eegtet = eegdata;
            curreeg = eegtet(eind1:eind2);
            
            % Plot
            eegscale = max(curreeg)-min(curreeg);
            downeeg = baseline; upeeg = downeeg+8;
            plotscale = 8;
            curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
            plot(taxisEEG,curreeg,'k-','LineWidth',1); hold on
            
            % Update baseline
            baseline = baseline + 8;
            
            % EEG in ripple band
            curreeg = double(ripamp(eind1:eind2));
            
            % Plot
            eegscale = max(curreeg)-min(curreeg);
            downeeg = baseline; upeeg = downeeg+8;
            plotscale = 8;
            curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
            plot(taxisEEG,curreeg,'k-','LineWidth',1);
            
            %Update baseline
            baseline = baseline + 1;
            
            set(gca, 'YTick',[]);
            
            winsecs = [0:10:winend-winst];
            secs = [winst:10:winend];
            secs = round(secs);
            msecs = [winst_ms:10000:winend_ms];
            xlabel('Time (secs)','FontSize',18,'Fontweight','normal');
            title(['Window Time: ' num2str(roundn(winst,-1)) ' - ' num2str(roundn(winend,-1)) ' - CellsActive ' num2str(length(cellsactive(:,1)))],'FontSize',18,'Fontweight','normal');
            
            % Draw Lines
            ylim = get(gca,'YLim');
            ypts = ylim(1):ylim(2);
            x1 = [0+prepostRip 0+prepostRip];
            x2 = [(winend-winst)-prepostRip (winend-winst)-prepostRip];
            y1 = [0 baseline];
            plot(x1,y1,'--k'); plot(x2,y1,'--k');
            xlim([0 (winend-winst)])
            set(gcf,'renderer','Painters')
            keyboard
            close
        end
    end
end
keyboard;


