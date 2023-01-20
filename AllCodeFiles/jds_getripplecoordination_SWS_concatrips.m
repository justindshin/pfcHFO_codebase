function jds_getripplecoordination_SWS_concatrips(animalprefixlist, eps)
%Get concatenated coordinated ripples across CA1 and PFC

savedata = 1;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    day = 1; %always single day expt
    
    % get ripple time
    load(sprintf('%s%sripplecoordinationCA1rips_ALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sripplecoordinationCtxrips_ALL0%d.mat',dir,animalprefix,day));
    
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time 
    
    ctxripplenew = [];
    ripplenew = [];
    
    for ep = eps
        if ep <10
            epochstring = ['0',num2str(ep)];
        else
            epochstring = num2str(ep);
        end
        ca1rip_coord = ripplecoordination_ca1{ep}.coordinated;
        ctxrip_coord = ripplecoordination_ctx{ep}.coordinated;
        
        ca1rip_noncoord = ripplecoordination_ca1{ep}.noncoordinated;
        ctxrip_noncoord = ripplecoordination_ctx{ep}.noncoordinated;
        
        swsep = sws{day}{ep};
        if ~isempty(swsep.starttime)
        swslist(:,1) = swsep.starttime;
        swslist(:,2) = swsep.endtime;
        
        curreegfile = [dir,'/EEG/',animalprefix,'eeg', '01','-',epochstring,'-','02']; %use any tetrode
        load(curreegfile);
        time1 = geteegtimes(eeg{day}{ep}{2}); % construct time array
        
        [~,swsvec] = wb_list2vec(swslist,time1);
        
        [~,ctxripvec_c] = wb_list2vec(ctxrip_coord(:,[1,2]),time1);
        [~,ctxripvec_nc] = wb_list2vec(ctxrip_noncoord(:,[1,2]),time1);
        
        [~,ca1ripvec_c] = wb_list2vec(ca1rip_coord(:,[1,2]),time1);
        [~,ca1ripvec_nc] = wb_list2vec(ca1rip_noncoord(:,[1,2]),time1);
        
        ctxripvec_new_c = swsvec & ctxripvec_c;
        ctxripvec_new_nc = swsvec & ctxripvec_nc;
        ca1ripvec_new_c = swsvec & ca1ripvec_c;
        ca1ripvec_new_nc = swsvec & ca1ripvec_nc;
        
        %Concatenate ctx and ca1 rippls
        concat_c = ctxripvec_new_c | ca1ripvec_new_c;
        concat_times = vec2list(concat_c,time1);
        
        ripplecoordination{day}{ep}.starttime = concat_times(:,1);
        ripplecoordination{day}{ep}.endtime = concat_times(:,2);
        else
        ripplecoordination{day}{ep}.starttime = [];
        ripplecoordination{day}{ep}.endtime = [];
        end
        clear ca1rip_coord ca1rip_noncoord ctxrip_coord ctxrip_noncoord swslist
    end
    
    if savedata == 1
        save(sprintf('%s/%sripplecoordinationSWS%02d.mat', dir, animalprefix, day), 'ripplecoordination');
    end
end
