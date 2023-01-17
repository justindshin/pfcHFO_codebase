% based on ripTrigFiringOnsetAnalysis3
clear all
close all

scrsz = get(0,'ScreenSize');
savedirX = '/Volumes/EnvoyPro/ProcessedDataNew/';
area='PFC';
riptype = {'noncoord','Coord'};

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

    if (ripInd == 2) && (isequal(area,'CA1'))
        load([savedirX sprintf('Allanim_%s250ca1ripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',rtype,area)])
    elseif (ripInd == 1) && (isequal(area,'CA1'))
        load([savedirX sprintf('Allanim_%s250ctxripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',rtype,area)])
    end

    if (ripInd == 2) && (isequal(area,'PFC'))
        load([savedirX sprintf('Allanim_%s250ca1ripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',rtype,area)])
    elseif (ripInd == 1) && (isequal(area,'PFC'))
        load([savedirX sprintf('Allanim_%s250ctxripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',rtype,area)])
    end
    
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


%% Match cells
inhModcells = [];
excModcells = [];

for c = 1:length(allpsths(:,1))
    indTmp = allinds(c,:);
    [B I] = ismember(indTmp, ncinds, 'rows','legacy');
    if B
        allMItmp = allMI(c);
        ncMItmp = ncMI(I);
%         if isequal(area,'CA1')
            if ncMItmp > 0
                excModcells = [excModcells; [ncMItmp allMItmp]];
            elseif ncMItmp < 0
                inhModcells = [inhModcells; [ncMItmp allMItmp]];
            end
%         elseif isequal(area,'PFC')
%             if allMItmp > 0
%                 excModcells = [excModcells; [ncMItmp allMItmp]];
%             elseif allMItmp < 0
%                 inhModcells = [inhModcells; [ncMItmp allMItmp]];
%             end
%         end
    end
end

[r p] = corrcoef(inhModcells);
[r2 p2] = corrcoef(excModcells);
figure; hold on
scatter(inhModcells(:,1),inhModcells(:,2),'b.')
lsline
scatter(excModcells(:,1),excModcells(:,2),'r.')
lsline
fontSize = 9; 
title(sprintf(['%s-CoordCA1rip IndPFCrip ripmod corr - pEXC-' num2str(p2(1,2)) ' rEXC-' num2str(r2(1,2))...
    ' pINH-' num2str(p(1,2)) ' rINH-' num2str(r(1,2))],area), 'FontSize', fontSize)
xlabel('IND PFCripMod Index')
ylabel('Coord CA1ripMod Index')
% xlim([-4 5])
% ylim([-4 10])
set(gcf, 'renderer', 'painters')

za  = 0.5*log((1+abs(r(1,2)))/(1-abs(r(1,2))));
zb  = 0.5*log((1+abs(r2(1,2)))/(1-abs(r2(1,2))));
szab = sqrt(1/(length(inhModcells(:,1))-3) + 1/(length(excModcells(:,1))-3));
z = abs(za-zb)/szab;
oneTl = 1-normcdf(z, 0, 1);
twoTl = 2*oneTl;

excSE = sqrt((1-(abs(r2(1,2))*abs(r2(1,2))))/(length(excModcells(:,1))-2));
inhSE = sqrt((1-(abs(r(1,2))*abs(r(1,2))))/(length(inhModcells(:,1))-2));
figure;
bar([abs(r(1,2)) abs(r2(1,2))],'k')
hold on
errorbar(1:2, [abs(r(1,2)) abs(r2(1,2))], [inhSE excSE],'k','LineStyle','none');
title(sprintf(['%s-CoordCA1rip IndPFCrip ripmod corrcoeff Diff - zpval-' num2str(twoTl)],area), 'FontSize', fontSize)
set(gcf, 'renderer', 'painters')
xticklabels({'INH','EXC'})
ylabel('Abs. correlation coefficient')
keyboard
