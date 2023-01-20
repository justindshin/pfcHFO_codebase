function js_ctxhprip_coord_HPrips_leadlag(animalprefixlist)

day = 1;
ctxrip_lead = 0;
hprip_lead = 0;
total_coordrip_num = 0;
number_coordrips = 0;
separation = [0];
percentage_data = zeros(3,length(separation));
onset_sep_ctx = [];
onset_sep_hp = [];
hprip_lead_eps = nan(length(animalprefixlist),9);
ctxrip_lead_eps = nan(length(animalprefixlist),9);
hprip_lead_all = [];
ctxrip_lead_all = [];

savedata = 0;
for a = 1:length(animalprefixlist)
    ripplecoupling = [];
    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    epochs = [1:2:17];

    for ep=1:length(epochs)
        epoch = epochs(ep);

        load(sprintf('%s%srippletime_coordSWS0%d.mat',dir,animalprefix,day));% get ripple time
        load(sprintf('%s%sctxrippletime_coordSWS0%d.mat',dir,animalprefix,day));% get ripple time

        ctx_rippletimes = ctxripple{day}{epoch}.starttime;
        hp_rippletimes = ripple{day}{epoch}.starttime;

        if (~isempty(ctx_rippletimes)) && (~isempty(hp_rippletimes))

            ctx_rippletimes_all = [ctx_rippletimes ctxripple{day}{epoch}.endtime];
            hp_rippletimes_all = [hp_rippletimes ripple{day}{epoch}.endtime];

            iri_ctx = diff(ctx_rippletimes(:,1));
            iri_hp = diff(hp_rippletimes(:,1));

            for r = 1:length(separation)
                ripsep = separation(r);
                keepidx_ctx = [1;find(iri_ctx>=ripsep)+1];
                keepidx_hp = [1;find(iri_hp>=ripsep)+1];

                ctx_rippletimes = ctx_rippletimes_all(keepidx_ctx,:);
                hp_rippletimes = hp_rippletimes_all(keepidx_hp,:);

                if (length(ctx_rippletimes(:,1)) > 20) && (length(hp_rippletimes(:,1)) > 20)
                    hprip_starttime = hp_rippletimes(:,1);
                    hpripbins_start = periodAssign(hprip_starttime, ctx_rippletimes(:,[1 2])); %lag
                    hprip_endtime = hp_rippletimes(:,2);
                    hpripbins_end = periodAssign(hprip_endtime, ctx_rippletimes(:,[1 2])); %lead

                    lagidx = hpripbins_start > 0;
                    leadidx = hpripbins_end > 0;

                    embed_find = sum([lagidx leadidx],2);
                    %Get rid of ripples that are embedded (st and end)
                    embed_idx = find(embed_find == 2);

                    lagidx(embed_idx) = 0;
                    leadidx(embed_idx) = 0;

                    if r == 1
                        %starttime differences for pfc and ca1 ripples
                        %where hp leads
                        leadhpripidx = find(leadidx == 1);
                        lagctxripidx = hpripbins_end(leadhpripidx);

                        ctxlagtimes_tmp = ctx_rippletimes(lagctxripidx,:);
                        hpleadtimes_tmp = hp_rippletimes(leadhpripidx,:);

                        timediff = ctx_rippletimes(lagctxripidx,1) - hprip_starttime(leadhpripidx);
                        onset_sep_hp = [onset_sep_hp; timediff];

                        %starttime differences for pfc and ca1 ripples
                        %where ctx leads
                        laghpripidx = find(lagidx == 1);
                        leadctxripidx = hpripbins_start(laghpripidx);

                        ctxleadtimes_tmp = ctx_rippletimes(leadctxripidx,:);
                        hplagtimes_tmp = hp_rippletimes(laghpripidx,:);

                        timediff2 = hprip_starttime(laghpripidx) - ctx_rippletimes(leadctxripidx,1);
                        onset_sep_ctx = [onset_sep_ctx; timediff2];
                    end

                    ctxrip_lead = sum(lagidx);
                    hprip_lead = sum(leadidx);
                    total_coordrip_num = sum(lagidx) + sum(leadidx);

                    hprip_lead_eps(a,ep) = hprip_lead/total_coordrip_num;
                    ctxrip_lead_eps(a,ep) = ctxrip_lead/total_coordrip_num;

                    hprip_lead_all = [hprip_lead_all;hprip_lead/total_coordrip_num];
                    ctxrip_lead_all = [ctxrip_lead_all;ctxrip_lead/total_coordrip_num];

                    number_coordrips = number_coordrips + sum(lagidx) + sum(leadidx);
                    percentage_data(1,r) = percentage_data(1,r) + total_coordrip_num;
                    percentage_data(2,r) = percentage_data(2,r) + ctxrip_lead;
                    percentage_data(3,r) = percentage_data(3,r) + hprip_lead;
                    ripplecoupling{day}{epoch}.ctxleadriptimes = ctxleadtimes_tmp;
                    ripplecoupling{day}{epoch}.hplagriptimes = hplagtimes_tmp;
                    ripplecoupling{day}{epoch}.ctxlagriptimes = ctxlagtimes_tmp;
                    ripplecoupling{day}{epoch}.hpleadriptimes = hpleadtimes_tmp;
                else
                    ripplecoupling{day}{epoch}.ctxleadriptimes = [];
                    ripplecoupling{day}{epoch}.hplagriptimes = [];
                    ripplecoupling{day}{epoch}.ctxlagriptimes = [];
                    ripplecoupling{day}{epoch}.hpleadriptimes = [];
                end
            end

        end
    end
    if savedata == 1
        save(sprintf('%s%srippletime_leadlag%02d.mat', dir,animalprefix,day), 'ripplecoupling');
    end
end

[p r] = ranksum(hprip_lead_all,ctxrip_lead_all)
dataMn = [nanmean(hprip_lead_all) nanmean(ctxrip_lead_all)]
dataSem = [(nanstd(hprip_lead_all)./sqrt(length(find(~isnan(hprip_lead_all)))))...
    (nanstd(ctxrip_lead_all)./sqrt(length(find(~isnan(ctxrip_lead_all)))))]

bar([mean(hprip_lead_all) mean(ctxrip_lead_all)],'k'); hold on
errorbar(1:2,dataMn,dataSem,'k','LineStyle','none')

datacombinedRipLead = [hprip_lead_all; ctxrip_lead_all];
g1 = repmat({'HP-Leads'},length(hprip_lead_all),1);
g2 = repmat({'PFC-Leads'},length(ctxrip_lead_all),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedRipLead,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['CA1PFC - Ripples Leading-p = ' num2str(p)])
ylabel('Proportion')
set(gcf, 'renderer', 'painters')

ctxrip_lead_all = percentage_data(2,:)./percentage_data(1,:);
hprip_lead_all = percentage_data(3,:)./percentage_data(1,:);
plot(ctxrip_lead_all, 'LineWidth',2)
hold on
plot(hprip_lead_all, 'LineWidth',2)

figure
histogram(onset_sep_ctx,20)
hold on
histogram(onset_sep_hp,20)

figure
plotcdf(onset_sep_ctx)
hold on
plotcdf(onset_sep_hp)

mnCtxEps = nanmean(ctxrip_lead_eps);
mnHpEps = nanmean(hprip_lead_eps);
semCtx = [];
semHp = [];
for i = 1:length(ctxrip_lead_eps(1,:))
    semCtx(i) = nanstd(ctxrip_lead_eps(:,i))./sqrt(length(find(~isnan(ctxrip_lead_eps(:,i)))));
    semHp(i) = nanstd(hprip_lead_eps(:,i))./sqrt(length(find(~isnan(hprip_lead_eps(:,i)))));
end

figure;
plot(mnCtxEps,'ro')
hold on
plot(mnHpEps,'ko')
errorbar(1:length(mnCtxEps), mnCtxEps, semCtx,'r','LineStyle','none')
errorbar(1:length(mnHpEps), mnHpEps, semHp,'k','LineStyle','none')
xlim([0.5 9.5])
x = [0.5 9.5];
y = [0.5 0.5];
plot(x,y,'--k')

keyboard
