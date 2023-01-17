function jds_getnoncoordripples_SWS(animalprefixlist, eps, area)
savedata = 1;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/NovelFamiliarNovel/%s_direct/',animalprefix);;

    day = 1; %always single day expt

    % get ripple time
    load(sprintf('%s%sripplecoordination%srips_ALL0%d.mat',dir,animalprefix,area,day));
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time
    ripplenew = [];
    for ep = eps

        if ep <10
            epochstring = ['0',num2str(ep)];
        else
            epochstring = num2str(ep);
        end
        if isequal(area,'CA1')
            rip = ripplecoordination_ca1{ep}.noncoordinated;
        elseif isequal(area,'Ctx')
            rip = ripplecoordination_ctx{ep}.noncoordinated;
        end

        riptimes(:,1) = rip(:,1);
        riptimes(:,2) = rip(:,2);

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

            ripplenew{day}{ep}.starttime = riptimes(:,1);
            ripplenew{day}{ep}.endtime = riptimes(:,2);
        else
            ripplenew{day}{ep}.starttime = [];
            ripplenew{day}{ep}.endtime = [];
        end
        clear riptimes swslist
    end
    clear ripple

    if (savedata == 1) && (isequal(area,'CA1'))
        ripple = ripplenew;
        save(sprintf('%s/%srippletime_noncoordSWS%02d.mat', dir, animalprefix, day), 'ripple');
    elseif (savedata == 1) && (isequal(area,'Ctx'))
        ctxripple = ripplenew;
        save(sprintf('%s/%sctxrippletime_noncoordSWS%02d.mat', dir, animalprefix, day), 'ctxripple');
    end
end
