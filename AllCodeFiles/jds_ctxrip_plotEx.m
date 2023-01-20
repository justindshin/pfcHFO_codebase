function jds_ctxrip_plotEx(animalprefixlist, day, epochs)

% SET DATA
%--------------------add your animals' directories--------------------%
allanimamp = [];
alltetavgenv = [];
pre = 750;
post = 750;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    if (day<10)
        daystring = ['0',num2str(day)];
    else
        daystring = num2str(day);
    end

    plotindiv = 1;
    smoothing_width = 0.004;
    samprate = 1500;
    for ep=1:length(epochs)
        epoch = epochs(ep);
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end

        epochstring2 = num2str(epoch);

        load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));% get ripple time
        load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get immmobility time
        load(sprintf('%s%stetinfo.mat',dir,animalprefix));

        tets = tetinfo{1}{1};

        %to use all tetrodes
        ctxtets = [];
        for t = 1:length(tets)
            tmp = tets{t};
            if isfield(tmp, 'descrip')
                if isequal(tmp.descrip, 'ctxriptet')
                    ctxtets = [ctxtets; t];
                end
            end
        end

        ctxtimetet = ctxtets(1);
        sws = sws{day}{epoch};
        swslist(:,1) = sws.starttime;
        swslist(:,2) = sws.endtime;

        rip_ctx = ctxripple{day}{epoch};
        riplist_ctx(:,1) = rip_ctx.starttime;%.starttime;
        riplist_ctx(:,2) = rip_ctx.endtime;%.endtime;

        if (ctxtimetet<10)
            ctxtimetetstring = ['0',num2str(ctxtimetet)];
        else
            ctxtimetetstring = num2str(ctxtimetet);
        end

        curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-',ctxtimetetstring];
        load(curreegfile);

        time1 = geteegtimes(eegref{day}{epoch}{ctxtimetet}) ; % construct time array

        [~,ctxripvec] = wb_list2vec(riplist_ctx(:,[1,2]),time1);
        [~,swsvec] = wb_list2vec(swslist,time1);

        ripvec_new = swsvec & ctxripvec;

        %standardize orientation
        vector = ripvec_new(:);
        timevec = time1(:);

        paddedvector = [ 0 ; vector ; 0];

        starttimes = find(diff(paddedvector) == 1);
        endtimes = find(diff(paddedvector)==-1);

        ctxripidxlist = [starttimes , endtimes];

        alltetavg = [];
        subplot(5,1,1);
        hold on

        for i = 1210%1:length(ctxripidxlist) 
            for t = 1:5%length(ctxtets)
                tet = ctxtets(t);
                if (tet<10)
                    ctxtetstring = ['0',num2str(tet)];
                else
                    ctxtetstring = num2str(tet);
                end
                currteteegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-',ctxtetstring];
                if isequal(animalprefix,'ER1')
                    ripeegfile = [dir,'/EEG/',animalprefix,'ripple', daystring,'-',epochstring2,'-',ctxtetstring];
                else
                    ripeegfile = [dir,'/EEG/',animalprefix,'ripple', daystring,'-',epochstring,'-',ctxtetstring];
                end
                load(currteteegfile);
                load(ripeegfile);
                eegdata = eegref{day}{epoch}{tet}.data;
                eegscaling = eegref{day}{epoch}{tet}.voltage_scaling  
                kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
                temprenv = ripple{day}{epoch}{tet}.data(:,3);
                ramp = ripple{day}{epoch}{tet}.data(:,1);
                ripenv = smoothvect(temprenv, kernel);
                meanEnv = mean(ripenv);
                thresh = meanEnv + (std(double(ripenv)))*3; %3 SD above mean
                allevents = [];
                alleventsenv = [];
                if ctxripidxlist(i,1)-pre > 0 && ctxripidxlist(i,2)+post < ctxripidxlist(length(ctxripidxlist),2)
                    datatoplot = eegdata(ctxripidxlist(i,1)-pre:ctxripidxlist(i,2)+post);
                    if length(datatoplot) >= 150 && max(datatoplot) < 500 && min(datatoplot) > -500
                        if plotindiv == 1
                            subplot(5,1,t)
                            plot(eegdata(ctxripidxlist(i,1)-pre:ctxripidxlist(i,2)+post)/eegscaling);
                            hold on;
                            lengthofrip = ctxripidxlist(i,2) - ctxripidxlist(i,1);
                            plot(pre+1:pre+1+lengthofrip,eegdata(ctxripidxlist(i,1):ctxripidxlist(i,2))/eegscaling,'r');
                            %                             X = [1500 1500 post+lengthofrip post+lengthofrip];
                            %                             Y = [min(datatoplot) max(datatoplot) min(datatoplot) max(datatoplot)];
                            %                             plot(X(1:2),Y(1:2), '-r'); plot(X(3:4),Y(3:4), '-r');
                            %                             close
                            y = [400 400];
                            x = [750 900];
                            y2 = [0 1000];
                            x2 = [750 750];
                            plot(x,y,'-b')
                            plot(x2,y2,'-m')
                            set(gcf, 'renderer', 'painters')
                        end
                        tmpmid = ceil(ctxripidxlist(i,1) + (ctxripidxlist(i,2) - ctxripidxlist(i,1))/2);
                        tmpamp = ramp(tmpmid-750:tmpmid+750)';
                        allevents = [allevents; tmpamp];
                        allenvtmp = ripenv(ctxripidxlist(i,1)-1500:ctxripidxlist(i,1)+1500)';
                        alleventsenv = [alleventsenv; allenvtmp];
                        %                     figure; subplot(2,1,1); plot(alleventstmp,'-k'); xticklabels([-1500:500:1500]); xlim([0 3000]); subplot(2,1,2); plot(smoothdata(allenvtmp),'-k'); xlim([0 3000]);
                        %                     xticklabels([-1500:500:1500]);
                        %                     xlabel('Samples from Ripple Onset');
                        %                     subplot(2,1,1); ylabel('EEG Trace');
                        %                     subplot(2,1,2); ylabel('Ripple power'); hold on
                        %                     yy = [thresh thresh];
                        %                     y2 = [meanEnv meanEnv];
                        %                     xx = [0 3000];
                        %                     plot(xx,yy,'-r')
                        %                     plot(xx,y2,'-b')
                        %                     keyboard;
                    end
                end
            end
            %             figure(tet); plot(mean(allevents));
            alltetavg = [alltetavg; mean(allevents,1)];
            alltetavgenv = [alltetavgenv; alleventsenv];
        end
        allanimamp = [allanimamp; mean(alltetavg,1)];
        clear swslist riplist_ctx
    end
end
figure(100);
subplot(2,1,1)
plot(mean(allanimamp)); hold on;
subplot(2,1,2)
plot(mean(alltetavgenv));
keyboard;

