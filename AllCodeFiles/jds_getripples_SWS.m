function jds_getripples_SWS(animalprefixlist, eps)
%Using UMAP (Uniform Manifold Approximation and Projection) to find
%clusters of coordinated SWR events based on CA1 or PFC population firing rate
%vectors

%Compile population firing rate vectors for all events (concatenate across
%all epochs?)
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    savedata = 1;
    dir = sprintf('/Volumes/JUSTIN/NovelFamiliarNovel/%s_direct/',animalprefix);;

    day = 1; %always single day expt

    % get ripple time
    load(sprintf('%s%sctxrippletime_ALL0%d.mat',dir,animalprefix,day));
%     load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));

    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time

    for ep = eps
        if ep <10
            epochstring = ['0',num2str(ep)];
        else
            epochstring = num2str(ep);
        end
        %     rip = ripplecoordination{ep}.coordinated;
        rip = ctxripple{day}{ep};

        %     riptimes(:,1) = rip(:,1);
        %     riptimes(:,2) = rip(:,2);
        riptimes(:,1) = rip.starttime;
        riptimes(:,2) = rip.endtime;

        %     rip_starttime = riptimes(:,1)*1000;  % in ms
        %
        %     rip_endtime = riptimes(:,2)*1000;  % in ms
        %
        % %     Find ripples separated by at least 500ms -- sj_getpopulationevents2
        %     iri = diff(rip_starttime);
        %     keepidx = [1;find(iri>=500)+1];
        %
        %     riplength = rip_endtime - rip_starttime;
        %     keepidx2 = find(riplength >= 50);% use the ripple last for more than 50 ms
        %     keepidx = intersect(keepidx,keepidx2);
        %
        %     rip_starttime = rip_starttime(keepidx);
        %     riptimes = riptimes(keepidx,:);

        swsep = sws{day}{ep};
        if ~isempty(swsep.starttime)
            swslist(:,1) = swsep.starttime;
            swslist(:,2) = swsep.endtime;

            curreegfile = [dir,'/EEG/',animalprefix,'eeg', '01','-',epochstring,'-','02']; %use any tetrode
            load(curreegfile);
            time1 = geteegtimes(eeg{day}{ep}{2}) ; % construct time array

            [~,swsvec] = wb_list2vec(swslist,time1);

            [~,ripvec] = wb_list2vec(riptimes(:,[1,2]),time1);

            ripvec_new = swsvec & ripvec;

            riptimes = vec2list(ripvec_new,time1);
            ctxripplenew{day}{ep}.starttime = riptimes(:,1);
            ctxripplenew{day}{ep}.endtime = riptimes(:,2);
%                     ripplenew{day}{ep}.starttime = riptimes(:,1);
%                     ripplenew{day}{ep}.endtime = riptimes(:,2);
        else
%                     ripplenew{day}{ep}.starttime = [];
%                     ripplenew{day}{ep}.endtime = [];
            ctxripplenew{day}{ep}.starttime = [];
            ctxripplenew{day}{ep}.endtime = [];
        end
        clear riptimes swslist
    end
    %CHANGE THIS TO OVERWRITE
    clear ctxripple
    ctxripple = ctxripplenew;
%     clear ripple
%     ripple = ripplenew;

    if savedata == 1
        save(sprintf('%s/%sctxrippletime_SWS%02d.mat', dir, animalprefix, day), 'ctxripple');
%             save(sprintf('%s/%srippletime_SWS%02d.mat', dir, animalprefix, day), 'ripple');
    end
end
