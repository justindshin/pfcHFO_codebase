% based on ripTrigFiringOnsetAnalysis3
clear all
close all

scrsz = get(0,'ScreenSize');
savedirX = '/Volumes/EnvoyPro/ProcessedDataNew/';
area='PFC';
riptype = {'noncoord','Coord','All'};

savefig1=0;
set(0,'defaultaxesfontsize',16);

allallripmodhists = {};
allallripmodhistssig = {};
allallripmodhistssigInd=  {};
allallripmodhiststype = {};
allallripmodpval = {};
allallripmodMI = {};

% how to smooth psths
b=gaussian(20,61);
totalCells = [];
for ripInd=1:length(riptype)
    rtype=riptype{ripInd};

    load([savedirX sprintf('Allanim_%s250ca1ripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',rtype,area)])
    allripplemod_idx=[];
    for w=1:length(allripplemod)
        allripplemod_idx=[allripplemod_idx;allripplemod(w).index];
    end

    allripmodhists=[];
    allripmodhistssig=[];
    allripmodhistssigInd=[];
    allripmodonset3=[];
    allripmodhiststype=[];
    allripmodpvals = [];
    allripmodMI = [];
    totalCells(ripInd) = length(allripplemod);
    for i=1:length(allripplemod)
        if allripplemod(i).rasterShufP2<0.05
            allripmodhists=[allripmodhists; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];

            curranim = allripplemod(i).index(1);
            allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
            allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
            % 1 for swr-exc, 0 for swr-inh
            allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
            allripmodpvals = [allripmodpvals; allripplemod(i).rasterShufP2];
            allripmodMI = [allripmodMI; allripplemod(i).Dm];
        end
    end

    allallripmodhists{ripInd}=allripmodhists;
    allallripmodhistssig{ripInd}=allripmodhistssig;
    allallripmodhistssigInd{ripInd}=allripmodhistssigInd;
    allallripmodhiststype{ripInd}=allripmodhiststype;
    allallripmodpval{ripInd} = allripmodpvals;
    allallripmodMI{ripInd} = allripmodMI;

    cntcellsRip=length(allripplemod);
end
ncpsths=allallripmodhistssig{1};
ncinds=allallripmodhistssigInd{1};
ncpvals = allallripmodpval{1};
ncMI = allallripmodMI{1};
allpsths=allallripmodhistssig{2};
allinds=allallripmodhistssigInd{2};
allpvals = allallripmodpval{2};
allMI = allallripmodMI{2};

figure;
propData = [length(allallripmodhiststype{1})/totalCells(1) length(allallripmodhiststype{2})/totalCells(2) length(allallripmodhiststype{3})/totalCells(3); ...
    length(find(allallripmodhiststype{1} == 0))/length(allallripmodhiststype{1}) length(find(allallripmodhiststype{2} == 0))/length(allallripmodhiststype{2}) length(find(allallripmodhiststype{3} == 0))/length(allallripmodhiststype{3}); ...
    length(find(allallripmodhiststype{1} == 1))/length(allallripmodhiststype{1}) length(find(allallripmodhiststype{2} == 1))/length(allallripmodhiststype{2}) length(find(allallripmodhiststype{3} == 1))/length(allallripmodhiststype{3})];

bar(propData)
xticklabels({'Prop. Modulated', 'Prop. Inh', 'Prop. Exc'})
title(sprintf('%s - NC vs Coord vs All mod',area))
ylabel('Proportion')
set(gcf, 'renderer', 'painters')
pvs = [];
chiStats = [];

%get sig prop diff for two pops of overall mod prop
for m = 1:3
    for n = 1:3
        if (m ~= n) && (n > m)
            n1 = length(allallripmodhiststype{m}); N1 = totalCells(m);
            n2 = length(allallripmodhiststype{n}); N2 = totalCells(n);
            x1 = [repmat('a',N1,1); repmat('b',N2,1)];
            x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
            [tbl,chi2stat,pval] = crosstab(x1,x2)
            pvs = [pvs; pval];
            chiStats = [chiStats; chi2stat];
        end
    end
end

for i = 1:2
    mod = i-1;
    for m = 1:3
        for n = 1:3
            if (m ~= n) && (n > m)
                n1 = length(find(allallripmodhiststype{m} == mod)); N1 = length(allallripmodhiststype{m});
                n2 = length(find(allallripmodhiststype{n} == mod)); N2 = length(allallripmodhiststype{n});
                x1 = [repmat('a',N1,1); repmat('b',N2,1)];
                x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
                [tbl,chi2stat,pval] = crosstab(x1,x2)
                pvs = [pvs; pval];
                chiStats = [chiStats; chi2stat];
            end
        end
    end
end


%% unsorted
zoomWindow=[-500 500];
figure('Position',[900 100 scrsz(3)/5 scrsz(4)])
xaxis=-1049:1050;
subplot(2,1,1);
imagesc(xaxis,1:size(allpsths,1),allpsths);
xlim(zoomWindow)
title('Ripple-triggered PSTHs, all rip-mod cells, CA1');xlabel('Time (ms)');ylabel('Cells')
caxis([-5 5])
subplot(2,1,2);
imagesc(xaxis,1:size(ncpsths,1),ncpsths);
caxis([-4 4])
xlim(zoomWindow)
title('Ripple-triggered PSTHs, all rip-mod cells, PFC');xlabel('Time (ms)');ylabel('Cells')

%% Match cells
NcnonModcells = [];
NcinhModcells = [];
NcexcModcells = [];
NcMI = [];
AllnonModcells = [];
AllinhModcells = [];
AllexcModcells = [];
changeInSig = [];
AllMI = [];
if length(allpsths(:,1)) > length(ncpsths(:,1))
    for c = 1:length(allpsths(:,1))
        indTmp = allinds(c,:);
        [B I] = ismember(indTmp, ncinds, 'rows','legacy');
        allVec = allpsths(c,:);
        allMItmp = allMI(c);
        allPtmp = allpvals(c);
        if B
            ncVec = ncpsths(I,:);
            ncMItmp = ncMI(I);
            ncPtmp = ncpvals(I);
            if ((ncPtmp < 0.05) || (allPtmp < 0.05))
                modChange = diff([ncMItmp allMItmp]);
                changeInSig = [changeInSig; modChange];
            end
            if ncPtmp < 0.05
                if ncMItmp > 0
                    NcexcModcells = [NcexcModcells; ncVec];
                    AllexcModcells = [AllexcModcells; allVec];
                elseif ncMItmp < 0
                    NcinhModcells = [NcinhModcells; ncVec];
                    AllinhModcells = [AllinhModcells; allVec];
                end
            else
                NcnonModcells = [NcnonModcells; ncVec];
                AllnonModcells = [AllnonModcells; allVec];
            end
        else %Nan vector, but check if this cell is now significant
            ncVec = nan(1,length(allpsths(c,:)));
            %             if allPtmp < 0.05
            %                 changeInSig = [changeInSig; allMItmp];
            %             end
            NcnonModcells = [NcnonModcells; ncVec];
            AllnonModcells = [AllnonModcells; allVec];
        end
    end

elseif length(ncpsths(:,1)) > length(allpsths(:,1))
    for c = 1:length(ncpsths(:,1))
        indTmp = ncinds(c,:);
        [B I] = ismember(indTmp, allinds, 'rows','legacy');
        ncVec = ncpsths(c,:);
        ncMItmp = ncMI(c);
        ncPtmp = ncpvals(c);
        if B
            allVec = allpsths(I,:);
            allMItmp = allMI(I);
            allPtmp = allpvals(I);
            if ((ncPtmp < 0.05) || (allPtmp < 0.05))
                modChange = diff([ncMItmp - allMItmp]);
                changeInSig = [changeInSig; modChange];
            end

            if ncPtmp < 0.05
                if ncMItmp > 0
                    NcexcModcells = [NcexcModcells; ncVec];
                    AllexcModcells = [AllexcModcells; allVec];
                elseif ncMItmp < 0
                    NcinhModcells = [NcinhModcells; ncVec];
                    AllinhModcells = [AllinhModcells; allVec];
                else
                    NcnonModcells = [NcnonModcells; ncVec];
                    AllnonModcells = [AllnonModcells; allVec];
                end
            end

        else %Nan vector, but check if this cell is now significant
            %             if ncPtmp < 0.05
            %                 changeInSig = [changeInSig; allMItmp];
            allVec = nan(1,length(ncpsths(c,:)));
            NcnonModcells = [NcnonModcells; ncVec];
            AllnonModcells = [AllnonModcells; allVec];
            %             end
        end
    end
end

%%  exc/inh

NCexcinh=allallripmodhiststype{2};
Allexcinh=allallripmodhiststype{1};

NCexcpsths=ncpsths(NCexcinh>0,:);
NCinhpsths=ncpsths(NCexcinh==0,:);

Allexcpsths=allpsths(Allexcinh>0,:);
Allinhpsths=allpsths(Allexcinh==0,:);

%% another look at the modulated PFC cells, and CA1 cells
% CA1 all ripmod cells
[tmp timingindCA1a]=min(Allinhpsths(:,801:1300)');
[tmp timingindCA1b]=max(Allexcpsths(:,801:1300)');
[A1 B1]=sort(timingindCA1a);
[A2 B2]=sort(timingindCA1b);
figure;
subplot(2,1,1)
imagesc(xaxis,1:length([B1';B2']), [Allinhpsths(B1,:); Allexcpsths(B2,:)]);caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title(sprintf('CA1 all CA1rip-mod cells - %dINH %dEXC',length(Allinhpsths(:,1)), length(Allexcpsths(:,1))))
xlim([-500 500])
hold on
subplot(2,1,2); hold on
boundedline(xaxis,mean(Allexcpsths,1),std(Allexcpsths)./sqrt(size(Allexcpsths,1)),'-r');
boundedline(xaxis,mean(Allinhpsths,1),std(Allinhpsths)./sqrt(size(Allinhpsths,1)),'-b');
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

figure;
subplot(2,1,1)
[tmp timingindPFCa]=min(NCinhpsths(:,801:1300)');
[tmp timingindPFCb]=max(NCexcpsths(:,801:1300)');
[A1 B1]=sort(timingindPFCa);
[A2 B2]=sort(timingindPFCb);
subplot(2,1,1)
imagesc(xaxis,1:length([B1';B2']), [NCinhpsths(B1,:); NCexcpsths(B2,:)]);caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title(sprintf('PFC all CA1rip-mod cells - %dINH %dEXC',length(NCinhpsths(:,1)), length(NCexcpsths(:,1))))
xlim([-500 500])
hold on
subplot(2,1,2); hold on
boundedline(xaxis,mean(NCexcpsths,1),std(NCexcpsths)./sqrt(size(NCexcpsths,1)),'-r');
boundedline(xaxis,mean(NCinhpsths,1),std(NCinhpsths)./sqrt(size(NCinhpsths,1)),'-b');
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

% With CA1 (also looking only at SWR-excited)
Allexcinh=allallripmodhiststype{1};
Allexcpsths=allpsths(Allexcinh>0,:);
Allinhpsths=allpsths(Allexcinh==0,:);
[tmp timingindCA1]=max(Allexcpsths(:,801:1300)');
[A2 B2]=sort(timingindCA1);
figure;
subplot(4,1,1)
imagesc(xaxis,1:length(B2),Allexcpsths(B2,:));caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title('CA1 NC-ctxRip-excited cells')
xlim([-1000 1000])

%look at CA1 inh
[tmp timingindCA1inh]=min(Allinhpsths(:,801:1300)');
[A2 B2]=sort(timingindCA1inh);
subplot(4,1,2)
imagesc(xaxis,1:length(B2),Allinhpsths(B2,:));caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title('CA1 NC-ctxRip-inhibited cells')
xlim([-1000 1000])


% another look at the SWR-excited PFC cells
[tmp timingind]=max(NCexcpsths(:,801:1300)');
[A B]=sort(timingind);
subplot(4,1,3);imagesc(xaxis,1:length(B),NCexcpsths(B,:));caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title('PFC NC-ctxRip-excited cells')

% another look at the SWR-inhibited PFC cells
[tmp timingind2]=min(NCinhpsths(:,801:1300)');
[A2 B2]=sort(timingind2);
xlim([-1000 1000])
subplot(4,1,4);imagesc(xaxis,1:length(B2),NCinhpsths(B2,:));caxis([-3 3])
xlim([-1000 1000])
xlabel('Time (ms)')
ylabel('#Cell')
title('PFC NC-ctxRip-inhibited cells')

%% STACKED PLOTS

[tmp timingindCA1inh]=min(Allinhpsths(:,801:1300)');
[A2 B2]=sort(timingindCA1inh);

figure; yyaxis left; imagesc(xaxis,1:length(B2),Allinhpsths(B2,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(Allinhpsths,1),std(Allinhpsths)./sqrt(size(Allinhpsths,1)),'-b')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')
%%

[tmp timingindCA1exc]=max(Allexcpsths(:,801:1300)');
[A2 B2]=sort(timingindCA1exc);

figure; yyaxis left; imagesc(xaxis,1:length(B2),Allexcpsths(B2,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(Allexcpsths,1),std(Allexcpsths)./sqrt(size(Allexcpsths,1)),'-r')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')

[tmp timingind]=max(NCexcpsths(:,801:1300)');
[A B]=sort(timingind);

figure; yyaxis left; imagesc(xaxis,1:length(B),NCexcpsths(B,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(NCexcpsths,1),std(NCexcpsths)./sqrt(size(NCexcpsths,1)),'-r')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')

[tmp timingind]=min(NCinhpsths(:,801:1300)');
[A B]=sort(timingind);

figure; yyaxis left; imagesc(xaxis,1:length(B),NCinhpsths(B,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(NCinhpsths,1),std(NCinhpsths)./sqrt(size(NCinhpsths,1)),'-b')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')

figure
shadedErrorBar(xaxis,mean(Allinhpsths,1),std(Allinhpsths)./sqrt(size(Allinhpsths,1)),'-b',1);
hold on
shadedErrorBar(xaxis,mean(Allexcpsths,1),std(Allexcpsths)./sqrt(size(Allexcpsths,1)),'-k',0);
xlim([-1000 1000])
title('CA1 mod');

figure
shadedErrorBar(xaxis,mean(NCinhpsths,1),std(NCinhpsths)./sqrt(size(NCinhpsths,1)),'-k',1);
hold on
shadedErrorBar(xaxis,mean(NCexcpsths,1),std(NCexcpsths)./sqrt(size(NCexcpsths,1)),'-r',0);
xlim([-1000 1000])
title('PFC mod');

figure
boundedline(xaxis,mean(ncpsths,1),std(ncpsths)./sqrt(size(ncpsths,1)),'-r');
hold on
boundedline(xaxis,mean(allpsths,1),std(allpsths)./sqrt(size(allpsths,1)),'-b');

xlim([-1000 1000])
title('Overall Mod');

keyboard
