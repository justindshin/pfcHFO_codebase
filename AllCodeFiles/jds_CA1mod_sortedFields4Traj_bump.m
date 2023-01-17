function jds_CA1mod_sortedFields4Traj_bump(animalprefixlist,epochs)

preFieldsExc = [];
preFieldsInh = [];
postFieldsExc = [];
postFieldsInh = [];

day = 1;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %%
    %-----match neurons across epochs-----%
    for e = 1:length(epochs)
        if epochs(e) ~= 16
            eps = epochs(e):(epochs(e) + 2);
            epsleep = eps(2);
            [ctxidx, hpidx] =  matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
            hpnum = length(hpidx(:,1));

            load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',dir,animalprefix,day));

            sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

            ep2 = find(sleeps(:,2) == epsleep);
            modcells = epochModulation.cellidx;
            inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
            inhcells(:,3) = -1;
            exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
            exccells(:,3) = 1;

            allmodcells = [inhcells; exccells];

            load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
            for t = 1:4
                preFieldsTmp = [];
                postFieldsTmp = [];
                for i = 1:hpnum
                    [b, cellLoc] = ismember(hpidx(i,[1 2]),hpidx,'rows','legacy');
                    cind = hpidx(i,:);
                    pre_ep = eps(1);
                    post_ep = eps(3);
                    preLinTmp = [];
                    postLinTmp = [];

                    tmpfield1 = linfields{day}{pre_ep}{cind(1)}{cind(2)}{t}(:,5);
                    tmpfield1(isnan(tmpfield1)) = 0;
                    tmpfield2 = linfields{day}{post_ep}{cind(1)}{cind(2)}{t}(:,5);
                    tmpfield2(isnan(tmpfield2)) = 0;
                    pos1 = linfields{day}{pre_ep}{cind(1)}{cind(2)}{t}(:,1);
                    pos2 = linfields{day}{post_ep}{cind(1)}{cind(2)}{t}(:,1);
                    stdLength1 = linspace(pos1(1),pos1(end),100);
                    stdLength2 = linspace(pos2(1),pos2(end),100);
                    tmp1 = lookup(stdLength1,pos1);
                    tmp2 = lookup(stdLength2,pos2);
                    tmp1 = tmpfield1(tmp1);
                    tmp2 = tmpfield2(tmp2);
                    preFieldsTmp = [preFieldsTmp; normalize(tmp1','range')];
                    postFieldsTmp = [postFieldsTmp; normalize(tmp2','range')];
                end
                figure; subplot(1,2,1); hold on
                spacing = 1;
                preFieldsTmp(find(nansum(preFieldsTmp') == 0),:) = [];
                hpidxTmp = hpidx;
                hpidxTmp(find(nansum(preFieldsTmp') == 0),:) = [];
                [tmp timingindCA1pre]=max(preFieldsTmp');
                [x sortedIdxPre] = sort(timingindCA1pre);
                cellSorted = hpidxTmp(sortedIdxPre,:);
                preFieldsSorted = preFieldsTmp(sortedIdxPre,:);
                modcellcnt = 1;
                for i = 1:length(cellSorted)
                    [b, cellLoc] = ismember(cellSorted(i,[1 2]),allmodcells(:,[1 2]),'rows','legacy');
                    xAx = 1:length(preFieldsSorted(i,:));
                    if b
                        if allmodcells(modcellcnt,3) == 1
                            com = ceil(mean(find(preFieldsSorted(i,:) == max(preFieldsSorted(i,:)))));
                            plot(xAx, preFieldsSorted(i,:)+spacing, 'r')
                            plot([com com], [preFieldsSorted(i,com)+spacing+1 preFieldsSorted(i,com)+spacing-1],'m')
                        elseif allmodcells(modcellcnt,3) == -1
                            com = ceil(mean(find(preFieldsSorted(i,:) == max(preFieldsSorted(i,:)))));
                            plot(xAx, preFieldsSorted(i,:)+spacing, 'b')
                            plot([com com], [preFieldsSorted(i,com)+spacing+1 preFieldsSorted(i,com)+spacing-1],'m')
                        end
                        modcellcnt = modcellcnt + 1;
                    else
                        com = ceil(mean(find(preFieldsSorted(i,:) == max(preFieldsSorted(i,:)))));
                        plot(xAx, preFieldsSorted(i,:)+spacing, 'k')
                        plot([com com], [preFieldsSorted(i,com)+spacing+1 preFieldsSorted(i,com)+spacing-1],'m')
                    end
                    spacing = spacing + 2;
                end
                ylim([0 max(preFieldsSorted(i,:))+spacing])
                set(gcf, 'renderer', 'painters')
                subplot(1,2,2)
                imagesc(flipud(preFieldsSorted)); colormap(inferno)
                title([animalprefix 'Run-' num2str(eps(1)) ' Sleep-' num2str(eps(2)) ' Traj-' num2str(t)])

                figure; subplot(1,2,1); hold on
                spacing = 1;
                postFieldsTmp(find(nansum(postFieldsTmp') == 0),:) = [];
                hpidxTmp = hpidx;
                hpidxTmp(find(nansum(postFieldsTmp') == 0),:) = [];
                [tmp timingindCA1post]=max(postFieldsTmp');
                [x sortedIdxPost] = sort(timingindCA1post);
                cellSorted = hpidxTmp(sortedIdxPost,:);
                postFieldsSorted = postFieldsTmp(sortedIdxPost,:);
                modcellcnt = 1;
                for i = 1:length(cellSorted)
                    [b, cellLoc] = ismember(cellSorted(i,[1 2]),allmodcells(:,[1 2]),'rows','legacy');
                    xAx = 1:length(postFieldsSorted(i,:));
                    if b
                        if allmodcells(modcellcnt,3) == 1
                            com = ceil(mean(find(postFieldsSorted(i,:) == max(postFieldsSorted(i,:)))));
                            plot(xAx, postFieldsSorted(i,:)+spacing, 'r')
                            plot([com com], [postFieldsSorted(i,com)+spacing+1 postFieldsSorted(i,com)+spacing-1],'m')
                        elseif allmodcells(modcellcnt,3) == -1
                            com = ceil(mean(find(postFieldsSorted(i,:) == max(postFieldsSorted(i,:)))));
                            plot(xAx, postFieldsSorted(i,:)+spacing, 'b')
                            plot([com com], [postFieldsSorted(i,com)+spacing+1 postFieldsSorted(i,com)+spacing-1],'m')
                        end
                        modcellcnt = modcellcnt + 1;
                    else
                        com = ceil(mean(find(postFieldsSorted(i,:) == max(postFieldsSorted(i,:)))));
                        plot(xAx, postFieldsSorted(i,:)+spacing, 'k')
                        plot([com com], [postFieldsSorted(i,com)+spacing+1 postFieldsSorted(i,com)+spacing-1],'m')
                    end
                    spacing = spacing + 2;
                end
                ylim([0 max(postFieldsSorted(i,:))+spacing])
                set(gcf, 'renderer', 'painters')
                subplot(1,2,2)
                imagesc(flipud(postFieldsSorted)); colormap(inferno)
                title([animalprefix ' Run-' num2str(eps(3)) ' Sleep-' num2str(eps(2)) ' Traj-' num2str(t)])
                keyboard 
                close all
            end
        else
            eps = [epochs(e) 17];
            epsleep = eps(2);
            [ctxidx, hpidx] =  matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
            hpnum = length(hpidx(:,1));

            load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',dir,animalprefix,day));

            sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

            ep2 = find(sleeps(:,2) == epsleep);
            modcells = epochModulation.cellidx;
            inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
            inhcells(:,3) = -1;
            exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
            exccells(:,3) = 1;

            allmodcells = [inhcells; exccells];

            %%
            load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
            for i = 1:length(hpnum(:,1))
                [b, cellLoc] = ismember(allmodcells(i,[1 2]),hpidx,'rows','legacy');
                cind = hpidx(i,:);
                pre_ep = eps(1);
                preLinTmp = [];
                preMax = [];
                preFieldsTmp = [];
                for t = 1:4
                    tmpfield1 = linfields{day}{pre_ep}{cind(1)}{cind(2)}{t}(:,5);
                    pos1 = linfields{day}{pre_ep}{cind(1)}{cind(2)}{t}(:,1);
                    stdLength1 = linspace(pos1(1),pos1(end),100);
                    tmp1 = lookup(stdLength1,pos1);
                    tmp1 = tmpfield1(tmp1);
                    preMax = [preMax; max(tmp1)];
                    preLinTmp = [preLinTmp; tmp1];
                    preFieldsTmp = [preFieldsTmp; tmp1'];
                end
                preField = normalize(preLinTmp(:,1),'range');
                preField(:,2) = preLinTmp(:,2);
                preTraj = preField(find(preField(:,2) == traj),1);
            end
            figure; hold on
            [tmp timingindCA1]=max(preFieldsExc{t}');
            [x sortedIdx] = sort(timingindCA1);
            imagesc(preFieldsExc{t}(sortedIdx,:)); colormap(inferno)
            title(['PreMod traj-' num2str(t) ' EXC'])
            [tmp timingindCA1]=max(preFieldsInh{t}');
            [x sortedIdx2] = sort(timingindCA1);
            imagesc(preFieldsInh{t}(sortedIdx2,:)); colormap(inferno)
            title(['PreMod traj-' num2str(t) ' INH'])
        end
    end
end

keyboard
