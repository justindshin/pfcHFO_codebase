function out = jds_triggered_wavelet(animalprefixlist, epochs)
% out = calcripspectrum(index, excludeperiods, eeg,ripples,cellinfo, options)
%  Computes the spectrogram around the middle of each decoded event.
%  Options:
%       appendindex-- Determines whether index is included in output vector
%           Default: 1
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.2]

%  out is a structure with the following fields
%       S-- This is a MxNxT matrix of the spectrogram for each tetrode
%           M is time relative to triggering event, N is frequency, T is event
%       F-- Frequency vector
%       T-- time relative to triggering event
%       fit-- This is the fit based on the spectrum computed for the entire
%           epoch to normalize S. To reconstruct S without normalization,
%           add log10(frequency)*fit(2)+fit(1) to S.
%       index-- Only if appendindex is set to 1 (default)


day = 1;

daystring = '01';

%parse the options
Fs = 1500;
win = [0.5 0.5];
S0 = 2*(1/Fs);
DJ = 0.15;
DT = 1/Fs;
N = sum(win)*Fs;
J1 = fix(log(N*DT/S0)/log(2))/DJ;
allAnimSpecs = [];
for a = 1:length(animalprefixlist)
    allEpSpecs = [];
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    
    for i = 1:length(epochs)
        epoch = epochs(i);
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
        
        load(sprintf('%s%srippletime_coordSWS0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%stetinfo.mat',dir,animalprefix));
        
        rip = ripple{day}{epoch};
        riptimes = [rip.starttime rip.endtime];
        if isempty(riptimes)
            continue
        end
        
        iri = diff(riptimes(:,1));
        keepidx = [1;find(iri>=0.5)+1];
        
        if isempty(keepidx)
            continue
        end
        riptimes = riptimes(keepidx,:);
        
        load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));
        rTets = find(~cellfun(@isempty,ctxripples{day}{epoch}));
        
        tetsNumRips = [];
        for scan = 1:length(rTets)
            t = rTets(scan);
            numR = length(ctxripples{day}{epoch}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
        end
        [ripcnt idx] = max(tetsNumRips);
        tetlist = rTets(idx);
        if ~isempty(riptimes)
            tets = tetinfo{1}{epoch};
            
            %to use all tetrodes
%             tetlist = [];
%             for t = 1:length(tets)
%                 tmp = tets{t};
%                 if isfield(tmp, 'descrip')
%                     if isequal(tmp.descrip, 'ctxriptet')
%                         tetlist = [tetlist; t];
%                     end
%                 end
%             end
            
            
            allTetSpecs = [];
            
            for ii = 1:length(tetlist)
                
                tet = tetlist(ii);
                
                if (tet<10)
                    tetstring = ['0',num2str(tet)];
                else
                    tetstring = num2str(tet);
                end
                
                curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-',tetstring];
                load(curreegfile);
                
                % Define EEG
                e = eegref{day}{epoch}{tet}.data';
                starttime = eegref{day}{epoch}{tet}.starttime;
                endtime = (length(e)-1) * (1 / Fs);
                
                % Define triggering events as the start of each ripple
                triggers = riptimes(:,1)-starttime;
                
                %Remove triggering events that are too close to the beginning or end
                
                while triggers(1)<win(1)
                    triggers(1) = [];
                end
                while triggers(end)> endtime-win(2)
                    triggers(end) = [];
                end
                
                % Compute a z-scored spectrogram using the mean and std for the entire session
                meanbase = [];
                idxcnt = 0;
                for s = 1:length(e)/Fs
                    basetmp = e((idxcnt+1):(Fs+idxcnt));
                    ptmp = wavelet(basetmp,1/Fs,0,DJ,S0,J1, 'MORLET');
                    ptmp = abs(ptmp);
                    idxcnt = idxcnt + Fs;
                    meanbase = cat(3, meanbase, ptmp);
                end
                
                allmeanbase = mean(meanbase,3);
                
                meanP = allmeanbase;
                stdP = std(meanbase,[],3);
                
                % Calculate the event triggered spectrogram
                windoweddata = jds_createdatamatc_wavelet(e,triggers,Fs,win);
                allwavsp = [];
                for s = 1:length(windoweddata(:,1))
                    [temp,period] = wavelet(windoweddata(s, 1:end),1/Fs,0,DJ,S0,J1,'MORLET');
                    pow=abs(temp);
                    allwavsp(:,:,s)=(pow-meanP)./stdP;
                end
                
                allTetSpecs = cat(3, allTetSpecs, mean(allwavsp,3));
                clear allwavsp
            end
        end
        allEpSpecs = cat(3, allEpSpecs, mean(allTetSpecs,3));
    end
    clear triggers riptimes
    meanspectAnim = mean(allEpSpecs,3);

    imagesc(meanspectAnim)
    ax = gca;
    ax.FontSize = 14
    frequency=round(arrayfun(@(x) 1/x, period));
    yLabels = arrayfun(@num2str,frequency,'UniformOutput',false);
    yticklabels({yLabels{10},yLabels{20},yLabels{30},yLabels{40},yLabels{50},yLabels{60}})
    ylabel('Frequency')
    xlim([0 1500])
    xticks([0:750:1500])
    xticklabels({'-500','0','500'})
    xlabel('Time From Ripple Onset (ms)')
    colorbar
    set(gcf, 'renderer', 'painters')
    saveg=0;
    if saveg==1
        figfile = ['/Volumes/JUSTIN/FigureWorking/CA1ReactivationSuppression/EPS/Wavelet/PFC/','Isolated',animalprefix,'Day',num2str(day),'Tetrode',num2str(tet)];
        print('-djpeg', figfile);
        print('-depsc', figfile);
    end
    close
    allAnimSpecs = cat(3, allAnimSpecs, mean(allEpSpecs,3));
end
meanspect = mean(allAnimSpecs,3);

imagesc(meanspect)
ax = gca;
ax.FontSize = 14
frequency=round(arrayfun(@(x) 1/x, period));
yLabels = arrayfun(@num2str,frequency,'UniformOutput',false);
yticklabels({yLabels{10},yLabels{20},yLabels{30},yLabels{40},yLabels{50},yLabels{60}})
ylabel('Frequency')
xlim([0 1500])
xticks([0:750:1500])
xticklabels({'-500','0','500'})
xlabel('Time From Ripple Onset (ms)')
colorbar
set(gcf, 'renderer', 'painters')
saveg = 1;
if saveg==1
    figfile = ['/Volumes/JUSTIN/FigureWorking/CA1ReactivationSuppression/EPS/Wavelet/PFC/','AllAnimWaveletCA1tetIndPFCripTrig'];
    print('-djpeg', figfile);
    print('-depsc', figfile);
end
keyboard;
end

