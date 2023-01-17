function jds_getSpikeWidth_meanFR(dataDir,animID,sessionNum)
%createNQSpikesFiles(rawDir,dataDir,animID,sessionNum)
%
%This function extracts sorted spike information from matclust files
%and saves data in the FF format. It is assumed that there is
%one (and only one) folder named '*.matclust' in the current working
%directory that contains matclust ('matclust_*.mat') files (result of createAllMatclustFiles.m and saving clustered
%data in matclust).
%
%The function also assumes that there is a '*.time' folder in the current directory
%conaining time information about the session (from extractTimeBinaryFile.m).
%
% It also assumes position information has already  been extracted, and is located in the dataDir.
%
%rawDir -- the directory where the raw dat folders are located
%dataDir -- the directory where the processed files should be saved
%animID -- a string identifying the animal's id (appended to the
%beginning of the files).
%sessionNum -- the session number (in chronological order for the animal)
cd(dataDir);
load('times.mat');
epochList = ranges(2:end,:); 

spikesInfo = [];
matclustfiles = dir('matclust_*');
disp('Processing MatClust files...');

for matclustfileInd = 8:length(matclustfiles)

    [filePath,prefixName,ext] = fileparts(matclustfiles(matclustfileInd).name);
    nt = strfind(prefixName,'nt');
    nTrodeNum = str2num(prefixName(nt(end)+2:end));
    channelString = num2str(nTrodeNum);

    clear clustattrib clustdata
    load(matclustfiles(matclustfileInd).name);

    relWaves= load(['waves_nt' channelString]); 

    load('times.mat', 'ranges');
    ranges= ranges(2:18,:);

    load(sprintf('%s%spos0%d.mat',dataDir,animID,sessionNum));
    load(sprintf('%s%scellinfo.mat',dataDir,animID));

    for e = 1:size(ranges,1)

        epochString = getTwoDigitNumber(e);
        currentSession = sessionNum;
        sessionString = getTwoDigitNumber(sessionNum);
        currentTimeRange = ranges(e,:);

        starttime = pos{sessionNum}{e}.data(1,1);
        endtime = pos{sessionNum}{e}.data(end,1);
        if ~isempty(clustattrib.clusters)
            for clustnum = 1:length(clustattrib.clusters)
                if (~isempty(clustattrib.clusters{clustnum}))
                    if (is_cluster_defined_in_epoch(clustnum,e))
                        if ~isempty(cellinfo{sessionNum}{e}{nTrodeNum})
                            if length(cellinfo{sessionNum}{e}{nTrodeNum}) >= clustnum
                                if isfield(cellinfo{sessionNum}{e}{nTrodeNum}{clustnum},'area')
                                    if isequal(cellinfo{sessionNum}{e}{nTrodeNum}{clustnum}.area,'CA1')
                                        tag = cellinfo{sessionNum}{e}{nTrodeNum}{clustnum}.tag2;
                                        spikesInfo{currentSession}{e}{nTrodeNum}{clustnum} = [];

                                        timepoints = clustdata.params(clustattrib.clusters{clustnum}.index,1);
                                        timepoints = timepoints(find((timepoints >= starttime) & (timepoints <= endtime)));
                                        if length(timepoints) > 50
                                            if (size(clustdata.params,2) > 4)
                                                amps = clustdata.params(clustattrib.clusters{clustnum}.index,2:5);
                                            else
                                                amps = nan(size(clustdata.params,1),1);
                                            end
                                            if (isequal(channelString,'25')) && (clustnum == 1) && (isequal(animID,'JS14'))
                                                continue
                                            end
                                            amps = amps(find((timepoints >= starttime) & (timepoints <= endtime)),:);
                                            ampvar = var(amps);
                                            [maxVal maxChan] = max(mean(amps));
                                            [trash, maxvar] = max(ampvar);
                                            amps = amps(:,maxvar);
                                            timepoints = timepoints(:);

                                            [csi propbursts] = computecsi(timepoints, amps, 6); %using 6ms Harris et al 2001

                                            %                             [burstIdx] = computeburstindex(timepoints, 6, 200, 300)
                                            xc = spikexcorr(timepoints, timepoints, 0.001, 0.500);
                                            burstIdx = mean(xc.c1vsc2(502:506))/mean(xc.c1vsc2(700:800));
                                            
                                            Lratio = ReportLRatio(clustnum);
                                            IsolDist = ReportIsolationDistance(clustnum);

                                            wavesIdx= clustattrib.clusters{clustnum}.index; % get idx of spikes in clust
                                            curWaves= relWaves.waves(:,:,wavesIdx);  % get waveforms of spikes
                                            chanWavs = curWaves(:,maxChan,:);
                                            chanWavs = reshape(chanWavs,40,length(chanWavs(1,1,:)));
                                            [peakIdx Imax] = max(chanWavs); % get idx of spike peak on each channel
                                            [trouIdx Imin] = min(chanWavs); % get idx of spike trough on each channel

                                            width = (double(Imin-Imax)*(1/30000))';    % get width in idx in each channel, convert to time in s
                                            width = width(find((timepoints >= starttime) & (timepoints <= endtime)),:); % winnow spike widths by current epoch

                                            currWavs = chanWavs(:,find((timepoints >= starttime) & (timepoints <= endtime)));
                                            av_waveform = mean(currWavs,2);
                                            waveform_sem = std(currWavs')./sqrt(length(currWavs(1,:)));

                                            spikedur= mean(width)*1000; %convert to ms
                                            [peakIdx Imaxavg] = max(av_waveform); % get idx of spike peak on each channel
                                            [trouIdx Iminavg] = min(av_waveform);
                                            spikedurAvg = double(abs(Imaxavg-Iminavg))*(1/30000);
                                            spikedurAvg = spikedurAvg*1000; %spike width from mean waveform

                                            if ~isempty(timepoints)
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.data(:,1) = timepoints;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.numSpks = length(timepoints);
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.meanrate = length(timepoints)/(endtime-starttime);
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.p2t = spikedur;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.p2tAvg = spikedurAvg;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.tag = tag;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.Lratio = Lratio;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.IsolationDistance = IsolDist;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.burstprobability = propbursts;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.burstindex = burstIdx;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.complexspkidx = csi;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.avgWav = av_waveform;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.wavSEM = waveform_sem;
                                                spikesInfo{currentSession}{e}{nTrodeNum}{clustnum}.timerange = currentTimeRange;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    else
                        continue
                    end
                end
            end
        end
    end
end

cd(dataDir)
save([animID,'spikesInfo',sessionString], 'spikesInfo');

%----------------------------------------------------------
%Helper functions
%-----------------------------------------------------------

function numString = getTwoDigitNumber(input)

if (input < 10)
    numString = ['0',num2str(input)];
else
    numString = num2str(input);
end
