clear; close all;
prefix='ZT2';
day=1;
% epochs=[1:2:17];
epochs=[15];

win = 10;
overlap = 0;

lowsp_thrs = 4; %cm/sec
highsp_thrs = lowsp_thrs;
dospeed = 1;
% passemnum = 2;
% cassemnum = 10;
passemnum = 1;
cassemnum = 11;
sd=3;
%ZT2 ep 15 pfc #2 ca1 #10 (cell cnts 67p 85c)
%ZT2 ep 13 pfc #6 ca1 #7(7 is good)
% Get Data Location
% -------------------------------
switch prefix
    case 'ZT2'
        directoryname = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        animdirect = directoryname;
        riptetlist = [10 11 12 14 16 17 18 19 29 24 25 27 32 36];  % No Need if no ripples
        maineegtet = 24;  % CA1 tet % No need if no EEG
        peegtet = 30; % PFCtet % No need if no EEG
        
    case 'KL8'
        directoryname = '/Volumes/Seagate Portable Drive/SingleDay/KL8_direct/';
        dire = '/Volumes/Seagate Portable Drive/SingleDay/KL8_direct/';
        animdirect = directoryname;
        riptetlist = [10];  % No Need if no ripples
        maineegtet = 10;  % CA1 tet % No need if no EEG
        peegtet = 5; % PFCtet % No need if no EEG
        
    case 'ER1'
        directoryname = '/Volumes/Seagate Portable Drive/SingleDay/ER1_direct/';
        dire = '/Volumes/Seagate Portable Drive/SingleDay/ER1_direct/';
        animdirect = directoryname;
        riptetlist = [13];  % No Need if no ripples
        maineegtet = 13;  % CA1 tet % No need if no EEG
        peegtet = 27; % PFCtet % No need if no EEG
end

currdir = pwd;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
if (dire(end) == '/')
    dire = dire(1:end-1);
end

if (day < 10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

% --------------- Parameters ---------------
% Lot of this is not needed if not doing EEG
%-------------------------------------------
Fspos = 30; %Hz
respos = 1/30; % sec
Fseeg = 1500; %Hz
reseeg = 1/1500; % sec
Fsspikes = 10000; %Hz
resspikes = 1/10000; %sec

%% -----------------------------------------
% SET DATA
% -------------------------------------------
eegtets = riptetlist;
% Also get a PFC eeg
eegtets = [eegtets, peegtet];
peegidx= find(eegtets==peegtet);
maineegidx = find(eegtets==maineegtet);
saveg = 0;

for e = 1:length(epochs)
    epoch = epochs(e);
    % Get data files
    % Spike data
    %-----------
    spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
    load(spikefile);
    tetinfofile = sprintf('%s/%stetinfo.mat', directoryname, prefix);
    load(tetinfofile);
    cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
    load(cellinfofile);
    % GET Rid of ripples and position info
    %-------------------------------
    pfcreactfile = sprintf('%s/%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat', directoryname, prefix, day);
    load(pfcreactfile);
    p_Rtime = RtimeStrength;
    pfcmainreactfile = sprintf('%s/%sPFC_icareactivationtimes20SWSSpk%02d.mat', directoryname, prefix, day);
    load(pfcmainreactfile);
    p_icareactivation = icareactivationtimes;
    
    ca1reactfile = sprintf('%s/%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat', directoryname, prefix, day);
    load(ca1reactfile);
    c_Rtime = RtimeStrength;
    ca1mainreactfile = sprintf('%s/%sCA1_icareactivationtimes20SWSSpk%02d.mat', directoryname, prefix, day);
    load(ca1mainreactfile);
    c_icareactivation = icareactivationtimes;
    posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
    load(posfile);
    linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
    load(linposfile);
    load(sprintf('%s/%ssws%02d.mat', animdirect, prefix, day));
    
    % Pos - Linpos
    % -------------
    filterString = 'strcmp($tag2, ''CA1Pyr'') && ($numspikes > 100)';
    
    cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsi = [repmat([day epoch], size(cellindices,1),1 ), cellindices]; % day-epoch-tet-cell for CA1 cells
    usecellsi = 1:size(cellsi,1);
    
    % Change below to put all PFC cells together, instead of ripmod or ripunmod
    % ---------------------------------------------------------
    filterString = 'strcmp($area, ''PFC'') && ($numspikes > 100)';
    pcellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsp = [repmat([day epoch], size(pcellindices,1),1 ), pcellindices]; % day-epoch-tet-cell for PFC cells
    usecellsp = 1:size(cellsp,1);
    
    % GET Spike data
    %-----------
    % CA1
    for i=1:size(cellsi,1)
        eval(['spiketimei{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,7);']);
    end
    
    % PFC
    for i=1:size(cellsp,1)
        eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,7);']);
    end
    
    % ------------------------------
    % Figure Parametersand Font Sizes
    % ------------------------------
    forppr = 0;
    
    set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
    
    if forppr==1
        set(0,'defaultaxesfontsize',16);
        tfont = 18; % title font
        xfont = 16;
        yfont = 16;
    else
        set(0,'defaultaxesfontsize',16);
        tfont = 20;
        xfont = 14;
        yfont = 14;
    end
    clr = {'b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y'};
    
    clr1='k';
    clr2='r';
    
    figdir = '/Volumes/JUSTIN/SingleDay/ReactivationRasters/';
    % ---------------------------------------
    %Set window time
    winst = 24298.6;
    winend = 24308.6; % secs
    epochend = 24308.6;
    
%     winst = sws{day}{epoch}.starttime(1,1);
%     winend = winst + win; % secs
%     epochend = pos{day}{epoch}.data(end,1);
    
    ii = 1;
    rastnum = 1;
    while winend <= epochend
        
        rastnum
        
        
        %         figure(1); xlim([0 win]); hold on;
        %         redimscreen;
        
        
        % WE NEED TO MAKE A taxis from WINST to WINEND (start to end of epoch)
        taxis = winst:winend; % IF winst is in secs, this will create an axis in seconds
        taxis = taxis - winst; % Start axis from 0
        
        winst_ms = winst*1000;
        winend_ms = winend*1000;
        
        % Look up position index at this time in stateatrix
        % posidx_rip = lookup(pt(i),statematrix.time); % All in secs
        
        baseline = 0;
        
        % First PFC Spikes on bottom. Each tick has space of height2, and using 1.8 of it for line
        % -----------------------
        if rastnum == 1
            preactivationTmp = p_icareactivation{7};
            pReactivation = preactivationTmp.post_strengths{passemnum}; 
            
            %update cellindices
            
            pWeights = preactivationTmp.pc_weights{passemnum};
            mnWeight = mean(pWeights(:,3));
            stdWeight = std(pWeights(:,3));
            highIdx = find(pWeights(:,3) > (mnWeight + 2*stdWeight));
            highCells = pWeights(highIdx,1:3);
            pWeights(highIdx,:) = [];
            updateWeightsP = [highCells; pWeights];
            phighCellIdx = highCells(:,[1 2]);
            
            pdeleteIdx = [];
            for pcell = 1:length(phighCellIdx(:,1))
                idx1 = find(cellsp(:,3) == phighCellIdx(pcell,1));
                idx2 = find(cellsp(:,4) == phighCellIdx(pcell,2));
                idx3 = intersect(idx1,idx2);
                pdeleteIdx = [pdeleteIdx; idx3];
            end
            usecellsp(pdeleteIdx) = [];
        end
        cnt = 0;
        activepfc = 0;
        for c=usecellsp %For member cells, plot last
            cellid = pcellindices(c,:);
            idx1 = find(updateWeightsP(:,1) == cellid(1));
            idx2 = find(updateWeightsP(:,2) == cellid(2));
            idx3 = intersect(idx1,idx2);
            if ~isempty(idx3)
                
                currspkt = spiketimep{c};
                currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
                
                if ~isempty(currspkt)
                    currspkt = currspkt - winst;
                end
                subplot(4,1,1); xlim([0 win]); hold on;
                if ~isempty(currspkt)
                    cnt=cnt+1;
                    activepfc = activepfc + 1;
                    
                    if size(currspkt,2)~=1, currspkt=currspkt'; end
                    plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','k','LineWidth',1);
                end
            end
        end
        for c=pdeleteIdx' %For member cells, plot last
            
            currspkt = spiketimep{c};
            currspkt = currspkt;
            currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
            
            if ~isempty(currspkt)
                currspkt = currspkt - winst;
            end
            subplot(4,1,1); xlim([0 win]); hold on;
            if ~isempty(currspkt)
                cnt=cnt+1;
                activepfc = activepfc + 1;
                
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','b','LineWidth',1);
            end
        end
        ylim([0 cnt*2]);
        reactidx = lookup([winst winend],pReactivation(:,1));
        
        title([prefix ' - Day',num2str(day) 'Ep' num2str(epoch) ' Window Time: ' num2str(roundn(winst,-1))...
            ' - ' num2str(roundn(winend,-1)) ' - CA1' num2str(cassemnum) ' - PFC' num2str(passemnum)]);
        
        subplot(4,1,2);
        plot(pReactivation(reactidx(1):reactidx(2),2),'b');
        xlim([0 win*50]);
        
        %plot the reactivation strength here
        
        baseline = baseline + (activepfc)*2;
        
        % Update baseline to give a little vertical gap
        baseline = 0;
        
        % Now, CA1 spikes
        % ---------------
        if rastnum == 1
            creactivationTmp = c_icareactivation{7};
            cReactivation = creactivationTmp.post_strengths{cassemnum};
            cWeights = creactivationTmp.pc_weights{cassemnum};
            mnWeight = mean(cWeights(:,3));
            stdWeight = std(cWeights(:,3));
            highIdx = find(cWeights(:,3) > (mnWeight + 2*stdWeight));
            highCells = cWeights(highIdx,1:3);
            cWeights(highIdx,:) = [];
            updateWeightsC = [highCells; cWeights];
            chighCellIdx = highCells(:,[1 2]);
            subplot(4,1,3); xlim([0 win]); hold on;
            cdeleteIdx = [];
            for ccell = 1:length(chighCellIdx(:,1))
                idx1 = find(cellsi(:,3) == chighCellIdx(ccell,1));
                idx2 = find(cellsi(:,4) == chighCellIdx(ccell,2));
                idx3 = intersect(idx1,idx2);
                cdeleteIdx = [cdeleteIdx; idx3];
            end
            usecellsi(cdeleteIdx) = [];
        end
        cnt = 0;
        activeca1cnt = 0;
        for c=usecellsi
            cellid = cellindices(c,:);
            idx1 = find(updateWeightsC(:,1) == cellid(1));
            idx2 = find(updateWeightsC(:,2) == cellid(2));
            idx3 = intersect(idx1,idx2);
            if ~isempty(idx3)
                
                currspkt = spiketimei{c};
                currspkt = currspkt;
                currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
                
                % If spikes, subtract from subtract from start time and bin
                if ~isempty(currspkt)
                    currspkt = currspkt - winst;
                end
                
                subplot(4,1,3); hold on;
                if ~isempty(currspkt)
                    cnt=cnt+1;
                    activeca1cnt = activeca1cnt+1;
                    
                    if size(currspkt,2)~=1, currspkt=currspkt'; end
                    plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','k','LineWidth',1);
                end
            end
        end
        for c=cdeleteIdx' %For member cells, plot last
            
            currspkt = spiketimei{c};
            currspkt = currspkt;
            currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
            
            if ~isempty(currspkt)
                currspkt = currspkt - winst;
            end
            subplot(4,1,3); xlim([0 win]); hold on;
            if ~isempty(currspkt)
                cnt=cnt+1;
                activeca1cnt = activeca1cnt+1;
                
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','r','LineWidth',1);
            end
        end
        ylim([0 cnt*2]);
        %Plot reactivation strength here
        reactidx = lookup([winst winend],cReactivation(:,1));
        subplot(4,1,4);
        plot(cReactivation(reactidx(1):reactidx(2),2),'r');
        xlim([0 win*50]);
        
        if rastnum == 1
            figure; subplot(2,1,1); stem(updateWeightsP(:,3),'k','filled')
            hold on
            stem(updateWeightsP(1:length(phighCellIdx(:,1)),3),'b','filled')
            view(90,90)
            ylim([-0.6 0.6])
            yticks([-0.5:0.1:0.5])
            xlim([-5 length(updateWeightsP(:,1))+5])
            xticks([1 length(updateWeightsP(:,1))])
            title('PFC Weights')
            
            subplot(2,1,2); stem(updateWeightsC(:,3),'k','filled')
            hold on
            stem(updateWeightsC(1:length(chighCellIdx(:,1)),3),'r','filled')
            view(90,90)
            ylim([-0.6 0.6])
            yticks([-0.5:0.1:0.5])
            xlim([-5 length(updateWeightsC(:,1))+5])
            xticks([1 length(updateWeightsC(:,1))])
            title('CA1 Weights')
        end
        
        % Update baseline to give a little vertical gap
        baseline = baseline + (activeca1cnt)*2;
        baseline = baseline+1;
        
                keyboard; % Pause after each plot
%         pause(3);
        
        saveg=0;
        if saveg==1
            figfile = [figdir,prefix,'Day',num2str(day),'Ep',num2str(epoch),'RectivationEgNo',num2str(ii),'Window',num2str(win)];
            print('-djpeg', figfile);
        end
        
        ii = ii+1;
        close all
        
        % Update winst and winend
        winst = winst + win - overlap;
        winend = winst + win;
        rastnum = rastnum + 1;
    end
end

keyboard;

