function jds_ripplecoordinationtriggered_assemblystrength_M(animalprefixlist,state)

%To do: plot the peak strength within specified ripple events. compare
%ripple types for different assemblies - might need to change binsize for
%raster generation
day = 1;

bins = 100;

areas = {'CA1','PFC'}
numCA1 = 0;
numPFC = 0;

for ar = 1:length(areas)
    area = areas{ar};

    allevents_coordriptrig = [];

    for a = 1:length(animalprefixlist)

        animalprefix = char(animalprefixlist(a));
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

        %Load reactivation strength file for all assemblies and epochs
        load(sprintf('%s%s%s_RTimeStrength%sNewSpk_20_%02d.mat',dir,animalprefix,area,state,day));

        %Load ripples
        load(sprintf('%s%sripplecoordinationSWS%02d.mat',dir,animalprefix,day));
        ripple = ripplecoordination;

        %Analyze epochs where there are assemblies detected
        epochs = find(~cellfun(@isempty,RtimeStrength));

        for e = 1:length(epochs)
            ep = epochs(e);
            assemblytmp = RtimeStrength{ep}.reactivationStrength;
            ripstarts = ripple{day}{ep}.starttime;

            if ar == 1
                numCA1 = numCA1 + length(assemblytmp);
            elseif ar == 2
                numPFC = numPFC + length(assemblytmp);
            end

            %Do CA1 ripples
            for ii = 1:length(assemblytmp)
                react_idx = [];
                for t = 1:length(ripstarts)
                    idxtmp = lookup(ripstarts(t), assemblytmp{ii}(:,1));
                    react_idx = [react_idx; idxtmp];
                end
                atmp = [];
                for r = 1:length(react_idx)
                    strengthstmp = assemblytmp{ii}(:,2);
                    if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                        tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                        atmp = [atmp; tmp'];
                    end
                end
                allevents_coordriptrig = [allevents_coordriptrig; zscore(mean(atmp,1))];
            end
        end
    end

    allevents_coordriptrig_z_mn = mean(allevents_coordriptrig);

    allevents_coordriptrig_z_sem = (std(allevents_coordriptrig)./sqrt(length(allevents_coordriptrig(:,1))));

    figure(1); hold on
    ax1 = gca;
    ax1.FontSize = 14;
    if ar == 1
        pl1 = plot([-bins:bins],allevents_coordriptrig_z_mn,'-k','LineWidth',1)
        boundedline([-bins:bins],allevents_coordriptrig_z_mn,allevents_coordriptrig_z_sem,'-k');
    end

    if ar == 2
        pl2 = plot([-bins:bins],allevents_coordriptrig_z_mn,'-r','LineWidth',1)
        boundedline([-bins:bins],allevents_coordriptrig_z_mn,allevents_coordriptrig_z_sem,'-r');
        title('Coord Ripple Triggered Reactivation Strength')

        ylabel('Reactivation Strength')
        xlabel('Time from Ripple Onset (s)')
        xlim([(-bins) bins]);
%         xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
        xticklabels({'-2','-1','0','1','2'});
        xx = [0 0];
        plot(xx,ax1.YLim,'--k')
        ylim = (ax1.YLim);
        legend([pl1 pl2],areas)
        set(gcf, 'renderer', 'painters')
    end

    figure
    imagesc(allevents_coordriptrig(randperm(length(allevents_coordriptrig(:,1))),:))
    xticks(101)
    title([sprintf('Ripple Triggered %s Reactivation Strength',area) '- Coord Rips'])
    colorbar
    set(gcf, 'renderer', 'painters')
end
keyboard

