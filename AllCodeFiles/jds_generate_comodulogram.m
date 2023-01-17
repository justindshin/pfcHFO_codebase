%% Constructing an example
%
clear
day = 1;
epochs = [1:2:17];
allComod = [];
allComodSm = [];
allComodSmAnim = [];
comodMean = [];
cnt = 0;
animalprefixlist = {'ER1','KL8','ZT2','JS17','JS21','JS34','JS14','JS15'};

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    figdir = '/Volumes/JUSTIN/SingleDay/Comodulograms/';
    load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sripples0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    animComod = [];
    for e = 1:length(epochs)
        epoch = epochs(e);
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
        
        rTetsctx = find(~cellfun(@isempty,ctxripples{day}{epoch}));
        rTetsca1 = find(~cellfun(@isempty,ripples{day}{epoch}));
        
        tetsNumRips = [];
        for scan = 1:length(rTetsctx)
            t = rTetsctx(scan);
            numR = length(ctxripples{day}{epoch}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
        end
        [ripcnt idx] = max(tetsNumRips);
        pfctet = rTetsctx(idx);
        
        if (pfctet<10)
            pfctetstring = ['0',num2str(pfctet)];
        else
            pfctetstring = num2str(pfctet);
        end
        
        tetsNumRips = [];
        for scan = 1:length(rTetsca1)
            t = rTetsca1(scan);
            numR = length(ripples{day}{epoch}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
        end
        [ripcnt idx] = max(tetsNumRips);
        ca1tet = rTetsca1(idx);
        
        if (ca1tet<10)
            ca1tetstring = ['0',num2str(ca1tet)];
        else
            ca1tetstring = num2str(ca1tet);
        end
        
        load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
        swslist = [sws{1}{epoch}.starttime sws{1}{epoch}.endtime];
        if ~isempty(swslist)
            %         for s = 1:length(swslist(:,1))
            %             if (swslist(s,2) - swslist(s,1)) >= 60
            if sum((swslist(:,2) - swslist(:,1))) >= 60
                cnt = cnt+1;
                curreegfile = [dir,'/EEG/',animalprefix,'eegref', '01' ,'-',epochstring,'-',ca1tetstring];
                load(curreegfile);
                lfp1 = eegref{1}{epoch}{ca1tet}.data';
                clear eegref
                
                curreegfile = [dir,'/EEG/',animalprefix,'eegref', '01' ,'-',epochstring,'-',pfctetstring];
                load(curreegfile);
                lfp2 = eegref{1}{epoch}{pfctet}.data';
                
                time1 = geteegtimes(eegref{1}{epoch}{pfctet}) ; % construct time array
                
                %                 [~,swsvec] = wb_list2vec(swslist(s,:),time1);
                [~,swsvec] = wb_list2vec(swslist,time1);
                lfp1 = lfp1(logical(swsvec));
                lfp2 = lfp2(logical(swsvec));
                
                data_length = length(lfp1);
                srate=1500;
                dt = 1/srate;
                t=dt*(1:data_length);
                
                %% Define the amplitude frequencies
                AmpFreqVector1=0:10:250;
                AmpFreqVector2=0:10:250;
                
                AmpFreq_BandWidth1=10;
                AmpFreq_BandWidth2=10;
                
                %% Filtering and Hilbert transform
                
                'filtering'
                tic
                Comodulogram=single(zeros(length(AmpFreqVector1),length(AmpFreqVector2)));
                AmpFreqTransformed1 = zeros(length(AmpFreqVector1), data_length);
                AmpFreqTransformed2 = zeros(length(AmpFreqVector2), data_length);
                
                for ii=1:length(AmpFreqVector1)
                    Af1 = AmpFreqVector1(ii);
                    Af2=Af1+AmpFreq_BandWidth1;
                    AmpFreq=eegfilt(lfp1,srate,Af1,Af2); % filtering
                    AmpFreqTransformed1(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
                end
                mn = mean(AmpFreqTransformed1,2);
                stddev = std(AmpFreqTransformed1');
                AmpFreqTransformed1 = (AmpFreqTransformed1 - mn)./stddev';
                
                for jj=1:length(AmpFreqVector2)
                    Pf1 = AmpFreqVector2(jj);
                    Pf2 = Pf1 + AmpFreq_BandWidth2;
                    AmpFreq=eegfilt(lfp2,srate,Pf1,Pf2); % filtering
                    AmpFreqTransformed2(jj, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope of LFP2
                end
                mn2 = mean(AmpFreqTransformed2,2);
                stddev2 = std(AmpFreqTransformed2');
                AmpFreqTransformed2 = (AmpFreqTransformed2 - mn2)./stddev2';
                toc
                
                %% Compute MI and comodulogram
                
                'Comodulation loop'
                
                counter1=0;
                for ii=1:length(AmpFreqVector1)
                    counter1=counter1+1;
                    
                    Pf1 = AmpFreqVector1(ii);
                    Pf2 = Pf1+AmpFreq_BandWidth1;
                    
                    counter2=0;
                    for jj=1:length(AmpFreqVector2)
                        counter2=counter2+1;
                        
                        Af1 = AmpFreqVector2(jj);
                        Af2 = Af1+AmpFreq_BandWidth2;
                        [MI]=ModIndex_Amp(AmpFreqTransformed1(ii, :), AmpFreqTransformed2(jj, :));
                        Comodulogram(counter1,counter2)=MI;
                    end
                end
                if sum(sum(Comodulogram)) == 0
                    keyboard
                end
                K = (1/9)*ones(3);
                comodSmooth = conv2(Comodulogram,K,'same');
                %             contourf(AmpFreqVector1+AmpFreq_BandWidth1/2,AmpFreqVector2+AmpFreq_BandWidth2/2,comodSmooth',30)
                %             keyboard
                %                 allComod(:,:,cnt) = Comodulogram;
                %                 allComodSm(:,:,cnt) = comodSmooth;
                allComod(:,:,cnt) = Comodulogram;
                allComodSm(:,:,cnt) = comodSmooth;
                animComod(:,:,cnt) = comodSmooth;
                %             end
                
                saveg=1;
                if saveg==1
                    figure; 
                    contourf(AmpFreqVector1+AmpFreq_BandWidth1/2,AmpFreqVector2+AmpFreq_BandWidth2/2,comodSmooth',20)
                    colormap(inferno); set(gcf,'renderer','Painters'); colorbar
                    figfile = [figdir,animalprefix,'Day',num2str(day),'Ep',num2str(epoch),'Comodulogram',num2str(epoch),'CA1tet',num2str(ca1tet),'PFCtet',num2str(pfctet)];
                    print('-depsc', figfile);
                    print('-djpeg', figfile);
                    close
                end
            end
        end
    end
    allComodSmAnim{a}.comodulograms = animComod;
    allComodSmAnim{a}.animalprefix = animalprefix;
    allComodSmAnim{a}.descrip = 'CA1-PFC comodulogram';
end
%% Plot comodulogram
keyboard
clf
contourf(AmpFreqVector1+AmpFreq_BandWidth1/2,AmpFreqVector2+AmpFreq_BandWidth2/2,mean(allComodSm,3)',20)
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar
