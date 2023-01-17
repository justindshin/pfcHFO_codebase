clear; close all;
prefix='ZT2';
day=1;
% epochs=[1:2:17];
epochs=[11];

% -------------------------------
win = 10;
overlap = 5;

% Get Data Location
% -------------------------------
switch prefix
    case 'ZT2'
        directoryname = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        animdirect = directoryname;
        riptetlist = [10 11 12 14 16 17 18 19 29 24 25 27 32 36];  % No Need if no ripples
        maineegtet = 14;  % CA1 tet % No need if no EEG
        peegtet = 40; % PFCtet % No need if no EEG

    case 'KL8'
        directoryname = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        animdirect = directoryname;
        riptetlist = [10];  % No Need if no ripples
        maineegtet = 10;  % CA1 tet % No need if no EEG
        peegtet = 5; % PFCtet % No need if no EEG

    case 'ER1'
        directoryname = '/Volumes/JUSTIN/SingleDay/ER1_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/ER1_direct/';
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
    % ripfile = sprintf('%s/%srippletime_noncoordSWS%02d.mat', directoryname, prefix, day);
    ripfile = sprintf('%s/%srippletime_SWS%02d.mat', directoryname, prefix, day);
    load(ripfile);
    nc_ripple = ripple;
    c_ripple = ripple;
    pripfile = sprintf('%s/%sctxrippletime_SWS%02d.mat', directoryname, prefix, day);
    load(pripfile);

    load(sprintf('%s/%sswsALL%02d.mat', directoryname, prefix, day));

    swstimes = sws{day}{epoch}.starttime;

    % Get cells
    % ---------
    % CA1 cells (black)
    filterString = 'strcmp($tag2, ''CA1Pyr'') && ($numspikes > 100)';

    cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsi = [repmat([day epoch], size(cellindices,1),1 ), cellindices]; % day-epoch-tet-cell for CA1 cells
    usecellsi = 1:size(cellsi,1);

    % Change below to put all PFC cells together, instead of ripmod or ripunmod
    % ---------------------------------------------------------
    % PFC cells
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

    for i=1:size(cellsp,1)
        eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,7);']);
    end

    riptimes = [nc_ripple{day}{epoch}.starttime nc_ripple{day}{epoch}.endtime]; %non coordinated ca1 ripples
    rip_starttime = riptimes(:,1);
    rip_endtime = riptimes(:,2);

    p_riptimes = [ctxripple{day}{epoch}.starttime ctxripple{day}{epoch}.endtime]; %Noncoordinated ctx ripples
    p_rip_starttime = p_riptimes(:,1);
    p_rip_endtime = p_riptimes(:,2);

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
        set(0,'defaultaxesfontsize',24);
        tfont = 28;
        xfont = 20;
        yfont = 20;
    end
    clr = {'b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y'};

    clr1='k';
    clr2='r';

    figdir = '/Volumes/JUSTIN/Rasters/';

    %     winst = 16608.5; %EXAMPLE PLOT
    %     winend = 16618.5; % secs
    %%
    % epochend = 6840;

    %5second window
    for rast = 2:length(swstimes(:,1))
        winst = swstimes(rast,1) - 2; %start 2 seconds before start of SWS
        winend = swstimes(rast,1) + 8; % full 10 seconds
        epochend = sws{day}{epoch}.endtime(end);
        epochend = winend;
        %while winend <= eegend
        ii = 1;
        rastnum = 1;


        rastnum
        rastnum = rastnum + 1;

        figure(1); xlim([0 win]); hold on;
        redimscreen;

        % WE NEED TO MAKE A taxis from WINST to WINEND (start to end of epoch)
        taxis = winst:winend; % IF winst is in secs, this will create an axis in seconds
        taxis = taxis - winst; % Start axis from 0

        winst_ms = winst*1000;
        winend_ms = winend*1000;

        % Find ripples within this window

        baseline = 0;

        % First PFC Spikes on bottom. Each tick has space of height2, and using 1.8 of it for line
        % -----------------------
        % First PFC spikes, ripunmod, then ripmod
        % Can change this to just have all PFC cells
        % ---------------------------------------------
        cnt = 0;
        activepfc = 0;
        [B, I] = sort(cellfun(@length,spiketimep),'descend');
        for c=usecellsp
            cc = I(c);
            eval(['currspkt = spiketimep{',num2str(cc),'};']);
            currspkt = currspkt;
            currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));

            % If spikes, subtract from subtract from start time and bin
            % NOTE THAT HERE WE ARE FORCING SPIKETIMES TO START FROM 0 IN OUR WINDOW,
            % BY SUBTRACTING WINST FROM SPIKETIMES
            % -------------------------------------------------------------
            if ~isempty(currspkt)
                currspkt = currspkt - winst;
            end
            figure(1); xlim([0 win]); hold on;
            if ~isempty(currspkt)
                activepfc = activepfc+1;
                cnt=cnt+1;
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','r','LineWidth',1);
            end
        end

        baseline = baseline + (activepfc)*2;
        baseline = baseline+1;

        % Now, CA1 spikes
        % ---------------
        cnt = 0;
        activeca1cnt = 0;
        [B, I] = sort(cellfun(@length,spiketimei),'descend');
        for c=usecellsi
            cc = I(c);
            eval(['currspkt = spiketimei{',num2str(cc),'};']);
            currspkt = currspkt;
            currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));

            % If spikes, subtract from subtract from start time and bin
            if ~isempty(currspkt)
                currspkt = currspkt - winst;
            end

            figure(1); hold on;
            if ~isempty(currspkt)
                activeca1cnt = activeca1cnt+1;
                cnt=cnt+1;
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','k','LineWidth',1);
            end
        end
        cnt = 0;
        baseline = baseline + (activeca1cnt)*2;
        baseline = baseline+4;

        ripsinwin = rip_starttime(find(rip_starttime>=winst & rip_starttime<=winend ));
        ripsinwinend = rip_endtime(find(rip_endtime>=winst & rip_endtime<=winend ));

        if ~isempty(ripsinwin)
            ripsinwin = ripsinwin - winst;
            ripsinwinend = ripsinwinend - winst;
        end
        plotraster(ripsinwin,(baseline+2*(cnt-1))*ones(size(ripsinwin)),1.8,[],'Color','b','LineWidth',5);
        plotraster(ripsinwinend,(baseline+2*(cnt-1))*ones(size(ripsinwinend)),1.8,[],'Color','g','LineWidth',5);

        cripsinwin = p_rip_starttime(find(p_rip_starttime>=winst & p_rip_starttime<=winend ));
        cripsinwinend = p_rip_endtime(find(p_rip_endtime>=winst & p_rip_endtime<=winend ));

        if ~isempty(cripsinwin)
            cripsinwin = cripsinwin - winst;
            cripsinwinend = cripsinwinend - winst;
        end

        plotraster(cripsinwin,(baseline+2*(cnt-1))*ones(size(cripsinwin)),1.8,[],'Color','m','LineWidth',5);
        plotraster(cripsinwinend,(baseline+2*(cnt-1))*ones(size(cripsinwinend)),1.8,[],'Color','k','LineWidth',5);

        % ------------------------------------------------------------------------
        % NOTE: Can change axis tick mark spacing. Skip abels ideally if plotting
        % entire epoch. You just need a scale bar
        % ------------------------------------------------------------------------
        winsecs = [0:10:winend-winst]; %eg. 0:10:60
        secs = [winst:10:winend];
        secs = round(secs); %secs = roundn(secs,-1);
        msecs = [winst_ms:10000:winend_ms];
        %set(gca,'XTick',winsecs,'XTickLabel',num2str(secs'));
        xlabel('Time (secs)','FontSize',18,'Fontweight','normal');
        title([prefix ' - Day',num2str(day) 'Ep' num2str(epoch) ' Window Time: ' num2str(roundn(winst,-1))...
            ' - ' num2str(roundn(winend,-1)) '  hpctet - ' num2str(maineegtet)...
            ' pfctet - ' num2str(peegtet)],'FontSize',18,'Fontweight','normal');
        %title(['SWR Time: ',num2str(roundn(pt(i),-1)), 'secs']);

        baseline = baseline+1;
        set(gca,'XLim',[0 winend-winst]);
        set(gca,'YLim',[0 baseline+10]);

        keyboard; % Pause after each plot

        saveg=0;
        if saveg==1
            figfile = [figdir,prefix,'Day',num2str(day),'Ep',num2str(epoch),'RasteregNo',num2str(ii),'Window',num2str(win),'Overlap',num2str(overlap)];
            print('-djpeg', figfile);
        end

        ii = ii+1;
        close all

    end
end

keyboard;

