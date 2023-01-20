function jds_getripples_SWS(animalprefixlist, eps, area)
%Get ripples contstrained by SWS

day = 1; %always single day expt

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    savedata = 1;
    dir = sprintf('/Volumes/JUSTIN/NovelFamiliarNovel/%s_direct/',animalprefix);

    % get ripple time
    if isequal(area,'PFC')
        load(sprintf('%s%sctxrippletime_ALL0%d.mat',dir,animalprefix,day));
        rip = ctxripple;
    elseif isequal(area,'CA1')
        load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));
        rip = ripple;
    end

    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time

    for ep = eps
        if ep <10
            epochstring = ['0',num2str(ep)];
        else
            epochstring = num2str(ep);
        end

        riptimes(:,1) = rip{day}{ep}.starttime;
        riptimes(:,2) = rip{day}{ep}.endtime;

        if ~isempty(swsep.starttime)
            swslist(:,1) = sws{day}{ep}.starttime;
            swslist(:,2) = sws{day}{ep}.endtime;

            curreegfile = [dir,'/EEG/',animalprefix,'eeg', '01','-',epochstring,'-','02']; %use any tetrode
            load(curreegfile);
            time1 = geteegtimes(eeg{day}{ep}{2}) ; %construct time array

            [~,swsvec] = wb_list2vec(swslist,time1);

            [~,ripvec] = wb_list2vec(riptimes(:,[1,2]),time1);

            ripvec_new = swsvec & ripvec;

            riptimes = vec2list(ripvec_new,time1);
            ripplenew{day}{ep}.starttime = riptimes(:,1);
            ripplenew{day}{ep}.endtime = riptimes(:,2);
        else
            ripplenew{day}{ep}.starttime = [];
            ripplenew{day}{ep}.endtime = [];
        end
        clear riptimes swslist
    end

    clear riptimes
    if isequal(area,'PFC')
        ctxripple = ripplenew;
        if savedata == 1
            save(sprintf('%s/%sctxrippletime_SWS%02d.mat', dir, animalprefix, day), 'ctxripple');
        end
    elseif isequal(area,'CA1')
        ripple = ripplenew;
        if savedata == 1
            save(sprintf('%s/%srippletime_SWS%02d.mat', dir, animalprefix, day), 'ripple');
        end
    end
end
