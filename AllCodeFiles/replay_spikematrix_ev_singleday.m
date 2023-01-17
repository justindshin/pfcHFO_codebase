function replay_spikematrix_ev_singleday(animalprefixlist,day,ep)


%%
bin = 20; % ms
%%
%----save the results?-----%
savedata = 1;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    %%
    %-----------define params and animal dir----------%
    if strcmp(animalprefix,'ZT2')
        animaldir = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/ER1_NEW_direct2/');
        exclude_list = [];
    elseif strcmp(animalprefix,'BG1')
        animaldir = '/Volumes/JUSTIN/SingleDay/BG1_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    elseif strcmp(animalprefix,'JS34')
        animaldir = '/Volumes/JUSTIN/SingleDay/JS34_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    elseif strcmp(animalprefix,'JS17')
        animaldir = '/Volumes/JUSTIN/SingleDay/JS17_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    elseif strcmp(animalprefix,'JS21')
        animaldir = '/Volumes/JUSTIN/SingleDay/JS21_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    elseif strcmp(animalprefix,'JS14')
        animaldir = '/Volumes/JUSTIN/SingleDay/JS14_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    elseif strcmp(animalprefix,'JS15')
        animaldir = '/Volumes/JUSTIN/SingleDay/JS15_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    elseif strcmp(animalprefix,'ER1')
        animaldir = '/Volumes/JUSTIN/SingleDay/ER1_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
        %         exclude_list = [1,2; 2, 2; 8,2; 9, 4; 9, 5; 21, 3; 23, 1;23, 3; 24,3; 24, 4; 24, 6; 24, 10; 26, 2]; %Feb, 16, 2017, for ER1;--Wenbo
        
    elseif strcmp(animalprefix,'KL8')
        animaldir = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
        %         exclude_list = [10,7; 15,7; 18,2; 18,5; 23,2; 23,4; 23,3; 23,5; 23,5; 23,8; 24,4]; %Apr, 21, 2017, for KL8;--Wenbo
    elseif strcmp(animalprefix,'AM2')
        animaldir = '/Volumes/JUSTIN/NovelFamiliarNovel/AM2_direct/';
        %         animaldir = ('/Volumes/Seagate Backup Plus Drive/SingledayExp/KL8_direct/');
        exclude_list = [];
    end
    eegdir = [animaldir,'EEG/'];
    
    %%
    %---------------------spike matrix----------------------%
    for e = 1:length(ep)
        
        epoch = ep(e);
        
%         if (mod(epoch,2) == 0 || epoch == 1)
%             eps = [epoch epoch+1];
%         else
%             eps = [epoch epoch-1];
%         end

        if epoch == 1
            eps = [epoch epoch+1];
        elseif epoch == 2
            eps = [epoch epoch-1];
        end

%         eps = [1:17]; %For matched
        
        [ctxidx, hpidx] = matchidx_acrossep_singleday(animaldir, animalprefix, day, eps, exclude_list); %(tet, cell)
        ctxnum = length(ctxidx(:,1));
        hpnum = length(hpidx(:,1));
        %     if ep <= 15
        %         eps = [ep ep+2]; %use two contiguous epochs of the same type
        %     elseif ep == 16
        %         eps = [ep ep-2];
        %     elseif ep ==17
        %         eps = [ep ep-2];
        %     end
        
        
        tmpflist1 = sprintf('%s%seeg%02d-%02d-%02d.mat', eegdir,animalprefix, day, epoch, ctxidx(1,1));
        
        load(tmpflist1);
        times_filteeg = geteegtimes(eeg{day}{epoch}{ctxidx(1,1)}) ;
        times_filteeg = times_filteeg(:)';% reference time is eeg time
        %     duration = (times_filteeg(end) - times_filteeg(1)).*1000; % ms
        %     nbins = ceil(duration/bin);
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
        nbins = length(timevect);
        spikes = loaddatastruct(animaldir, animalprefix, 'spikes', day); % get spikes
        for i = 1:ctxnum
            % get cortical neurons
            index = [day,epoch,ctxidx(i,:)];
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            spkstime = ( double(spiketimes)-times_filteeg(1) ).*1000; % ms
            
            if bin == 1
                spkstime = round((double(spiketimes)-times_filteeg(1) ).*1000); % ms
                for ntime = 1:nbins
                    nspks = length(find(spkstime == ntime*bin));
                    spike_matrix_ctx(i,ntime) = nspks;
                end
            else
                for ntime = 1:nbins
                    nspks = length(find(spkstime >= (ntime-1)*bin + 1 & spkstime < ntime*bin));
                    spike_matrix_ctx(i,ntime) = nspks;
                end
            end
        end
        
        for i = 1:hpnum
            index = [day,epoch,hpidx(i,:)] ;
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            spkstime = ( double(spiketimes)-times_filteeg(1) ).*1000; % ms
            
            if bin == 1
                spkstime = round((double(spiketimes)-times_filteeg(1) ).*1000); % ms
                for ntime = 1:nbins
                    nspks = length(find(spkstime == ntime*bin));
                    spike_matrix_hp(i,ntime) = nspks;
                end
            else
                for ntime = 1:nbins
                    nspks = length(find(spkstime >= (ntime-1)*bin + 1 & spkstime < ntime*bin));
                    spike_matrix_hp(i,ntime) = nspks;
                end
            end
        end
        
        observation_matrix{epoch}.ctxdata =  spike_matrix_ctx;
        observation_matrix{epoch}.hpdata =  spike_matrix_hp;
        observation_matrix{epoch}.hpidx =  hpidx;
        observation_matrix{epoch}.ctxidx =  ctxidx;
        observation_matrix{epoch}.ncortical =  ctxnum;
        observation_matrix{epoch}.nhp =  hpnum;
        observation_matrix{epoch}.timeeeg =  times_filteeg;
        clear spike_matrix_hp spike_matrix_ctx
    end
    %%
    if savedata
%         save(sprintf('%s%s_spikematrix_ev_allepochallcell%d_%02d.mat', animaldir,animalprefix,bin,day), 'observation_matrix', '-v7.3');
%         save(sprintf('%s%s_spikematrix_ev_trackedallepochallcell%d_%02d.mat', animaldir,animalprefix,bin,day), 'observation_matrix', '-v7.3');
        save(sprintf('%s%s_spikematrix_ev_ep1allcell%d_%02d.mat', animaldir,animalprefix,bin,day), 'observation_matrix', '-v7.3');

    end
end