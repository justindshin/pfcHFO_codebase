function jds_extractdeltawaves_zugaro(animalprefixlist)
%extracts delta waves using the method/parameters in Todorova/Zugaro

day = 1;
epochs = [1:2:17];
daystring = '01';
savedata = 1;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    
    tets = tetinfo{1}{epochs(1)};
    
    ctxtets = []; %get all ctxriptet tetrodes
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'ctxriptet')
                ctxtets = [ctxtets; t];
            end
        end
    end
    
    for e = 1:length(epochs)
        ampdataall = [];
        epoch = epochs(e);
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
        
        for i = 1:length(ctxtets)
            ctxtet = ctxtets(i);
            
            if (ctxtet<10)
                ctxtetstring = ['0',num2str(ctxtet)];
            else
                ctxtetstring = num2str(ctxtet);
            end
            
            curreegfile = [dir,'/EEG/',animalprefix,'delta', daystring,'-',epochstring,'-',ctxtetstring];
            load(curreegfile);
            
            ampdatatmp = delta{day}{epoch}{ctxtet}.data(:,1);
            ampdataall = [ampdataall; ampdatatmp'];
        end
        
        ampdata = mean(ampdataall); %mean amplitude across all tetrodes
        
        %Zscore amplitude data and find epochs where peak is >=-2 (LFP is
        %flipped)
        z_amp = zscore(ampdata);
        zci = find(diff(sign(z_amp))); %find the zero crossing indices
        startidx = zci(1:end-1); %look at pretty much every interval
        endidx = zci(2:end);
        indices = [startidx' endidx'];
        
        overthresh_idx = [];
        peak_idx = [];
        for l = 1:length(indices(:,1))
            idx = indices(l,:);
            ampvector = z_amp(idx(1):idx(2));
            [min_z minidx] = min(ampvector);
            %             downwardMono = ~any(diff(ampvector(1:minidx))>0);
            %             upwardMono = ~any(diff(ampvector(minidx+1:end))<0);
            minidx = (idx(1) + minidx(1)) - 1;
            if min_z < -2
                overthresh_idx = [overthresh_idx; idx];
                peak_idx = [peak_idx; minidx];
            end
        end
        
        times = geteegtimes(delta{day}{epoch}{ctxtet}) ; % construct time array
        
        deltawavetimes = times(overthresh_idx);
        peaktimes = times(peak_idx)';
        
        %Filter by length (>150ms but <500ms, Todorova)
        deltalengths = (deltawavetimes(:,2) - deltawavetimes(:,1))*1000;
        keepidx = find(deltalengths > 150 & deltalengths < 500);
        deltawavetimes = deltawavetimes(keepidx,:);
        peaktimes = peaktimes(keepidx,:);
        
        %Constrain by SWS
        load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time
        
        swstime = sws{day}{epoch};
        if ~isempty(swstime.starttime)
            swslist = [swstime.starttime swstime.endtime];
            [~,swsvec] = wb_list2vec(swslist,times);
            [~,wavevec] = wb_list2vec(deltawavetimes,times);
            
            waves_new = swsvec & wavevec;
            
            wavetimes = vec2list(waves_new,times);
            
            if length(wavetimes(:,1)) > 1
                peaktimes_bins = periodAssign(peaktimes,wavetimes);
                peaktimes_keep = peaktimes(find(peaktimes_bins ~= 0));
                peaktimes_keep2 = peaktimes_bins(find(peaktimes_bins ~= 0));
                
                waveswithpeak = wavetimes(peaktimes_keep2,:);
                finalpeaktimes = peaktimes_keep;
                
                %error test
                wavebins = periodAssign(finalpeaktimes,waveswithpeak);
                dif = diff(wavebins);
                mismatch = find(dif > 1);
                if ~isempty(mismatch)
                    keyboard;
                end
                
                deltawaves{day}{epoch}.starttime = waveswithpeak(:,1);
                deltawaves{day}{epoch}.endtime = waveswithpeak(:,2);
                deltawaves{day}{epoch}.peaktime = finalpeaktimes;
                deltawaves{day}{epoch}.descrip = 'delta waves extracted from mean of ctxriptets, 150-500ms';
            else
                deltawaves{day}{epoch}.starttime = [];
                deltawaves{day}{epoch}.endtime = [];
                deltawaves{day}{epoch}.peaktime = [];
            end
        else
            deltawaves{day}{epoch}.starttime = [];
            deltawaves{day}{epoch}.endtime = [];
            deltawaves{day}{epoch}.peaktime = [];
        end
        clear swslist
    end
    
    if savedata == 1;
        save(sprintf('%s%sdeltawavetimes_SWS%02d.mat', dir,animalprefix,day), 'deltawaves');
    end
    clear deltawaves
end