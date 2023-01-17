clear all
close all
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS15','JS14','ER1','KL8'};
epochs = [3:2:17];
day = 1;
area = 'CA1';
mean_coordCorr = [];
mean_noncoordCorr = [];
mean_diff = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    load(sprintf('%s/%srippletime_noncoordSWS%02d.mat', dir, animalprefix, day));
    load(sprintf('%s%s%s_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,area,day));
        nc_rip = ripple; clear ripple
%     nc_rip = ripple; clear ripple
    load(sprintf('%s/%srippletime_coordSWS%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
        c_rip = ripple; clear ripple
%     c_rip = ripple; clear ripple

    dat = [];
    for e = 1:length(epochs)
        epoch = epochs(e);

        nc_riptimes = [nc_rip{day}{epoch}.starttime nc_rip{day}{epoch}.endtime];

        c_riptimes = [c_rip{day}{epoch}.starttime c_rip{day}{epoch}.endtime];

        if (length(nc_riptimes(:,1)) < 20) || (length(c_riptimes(:,1)) < 20)
            continue
        end

        rtmp = RtimeStrength{epoch};
        
        for r = 1:length(rtmp.reactivationStrength)
            assembly = rtmp.weights{r};
            thresh = mean(assembly(:,3)) + 2*std(assembly(:,3));
            cellidx = assembly(find(assembly(:,3) > thresh),[1 2]);
            numcells = length(cellidx(:,1));
            if numcells > 1
                if ((length(nc_riptimes) > 10) && (length(c_riptimes) > 10))
                    celldata = [];
                    spikecounts = [];
                    for cellcount = 1:numcells %get spikes for each cell
                        index = [day,epoch,cellidx(cellcount,:)] ;
                        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                        else
                            spiketimes = [];
                        end
                        spikebins = periodAssign(spiketimes, nc_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                        spikecount = zeros(1,size(nc_riptimes,1));
                        for i = 1:length(spikebins)
                            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                        end
                        spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
                    end
                    for i = 1:numcells
                        for ii = 1:numcells
                            if (i ~= ii) && (ii > i)
                                n1 = spikecounts(:,i);
                                n2 = spikecounts(:,ii);
                                coactiveZ = coactivezscore(n1, n2);
                                mean_noncoordCorr = [mean_noncoordCorr; coactiveZ];
                            end
                        end
                    end
                    %do coord rips
                    celldata = [];
                    spikecounts = [];
                    for cellcount = 1:numcells %get spikes for each cell
                        index = [day,epoch,cellidx(cellcount,:)] ;
                        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                        else
                            spiketimes = [];
                        end
                        spikebins = periodAssign(spiketimes, c_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                        spikecount = zeros(1,size(c_riptimes,1));
                        for i = 1:length(spikebins)
                            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                        end
                        spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
                    end
                    for i = 1:numcells
                        for ii = 1:numcells
                            if (i ~= ii) && (ii > i)
                                n1 = spikecounts(:,i);
                                n2 = spikecounts(:,ii);
                                coactiveZ = coactivezscore(n1, n2);
                                mean_coordCorr = [mean_coordCorr; coactiveZ];
                            end
                        end
                    end
                end
            end
        end
    end
end

[p h] = signrank(mean_noncoordCorr,mean_coordCorr);
datacombinedCofiring = [mean_noncoordCorr; mean_coordCorr];
g1 = repmat({'Independent'},length(mean_noncoordCorr),1);
g2 = repmat({'Coordinated'},length(mean_coordCorr),1);
g = [g1;g2];

figure
histogram(mean_noncoordCorr-mean_coordCorr)
x = [nanmedian(mean_noncoordCorr-mean_coordCorr) nanmedian(mean_noncoordCorr-mean_coordCorr)];
y = [0 60];
hold on
plot(x,y,'-r','LineWidth',4)
[h2 p2] = ttest(mean_noncoordCorr-mean_coordCorr)
title(['RippleCofiringDifference(ind-coord)-p = ' num2str(p2)])

figure;
h = boxplot(datacombinedCofiring,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['RippleCofiring-p = ' num2str(p)])
ylabel('Cofiring (z)')
set(gcf, 'renderer', 'painters')
keyboard

