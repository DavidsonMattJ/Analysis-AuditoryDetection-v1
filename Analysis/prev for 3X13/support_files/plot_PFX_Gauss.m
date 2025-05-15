% plot_PFX_Gauss

% useppants= [1:16,18:19]; % 17 and 20 have poor gaussian fits.
cd(savedatadir);
%
figure(1);clf
set(gcf,'units', 'normalized', 'Position', [.1 .1 .8 .6], 'color', 'w');
load('GFX_SOAdata.mat');
load('GFX_GaussData.mat');

nppants = size(GFX_Data_SOAs.SOAs,1);

for ippant = 1:nppants

figure(1);clf
set(gcf,'units', 'normalized', 'Position', [.1 .1 .8 .6], 'color', 'w');

% plot slow and normal speedS:
useCols= {'b', 'r', 'k'};
legh=[];
SOAs= GFX_Data_SOAs.SOAs(1,:);
plotData= {GFX_Data_SOAs.propSame_slow, GFX_Data_SOAs.propSame_fast,  GFX_Data_SOAs.propSame_all}; 

nTrialsall= {GFX_Data_SOAs.nTrialsSOA_slow, GFX_Data_SOAs.nTrialsSOA_fast,  GFX_Data_SOAs.nTrialsSOA_all}; 
 
for iplot = 1:2
    pltd= plotData{iplot};
subplot(131);
hold on
 
% pM= squeeze(mean(pltd,1));

pM = pltd(ippant,:);

legh(iplot)=plot(SOAs, pM, 'o-', 'color', useCols{iplot}, ...
    'MarkerFaceColor', useCols{iplot}, 'MarkerSize',2,...
    'LineWidth',2);

hold on;
for itext= 1:length(SOAs)
    
    text(SOAs(itext), pM(itext)+.01, num2str(nTrialsall{iplot}(ippant,itext)), 'HorizontalAlignment', 'center');
end
end
%%
xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title(['participant ' num2str(ippant) ' results']);
set(gca,'fontsize',12);
%%
% now plot the gauss fits,

subplot(132); cla
useCols= {'k', 'b', 'r'};
% include group data for comparison?
[nsubs,nconds, nsamps]= size(GFX_gauss_Cond_x_Yfit);
%

for icond= 2:3%:3 % first cond is 'all'
    plotme = squeeze(GFX_gauss_Cond_x_Yfit(:,icond,:));

%     pM = mean(plotme,1);
    pM = plotme(ippant,:);
    
plot(gauss_xvec, pM, '-','color', useCols{icond}, 'linew',2)
    hold on;
end
shg

%%

xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title('Participant gauss fits');
set(gca,'fontsize',12);

%%quick s tats check
% [nsubs, nconds, param(mean,sd)]= size(GFX_gauss_Cond_xparam)

%%  now for the gausian parameters:
usesubspots= [3,6];%
titlesAre = {'Mean (mu)', 'STD (sig)', 'amp'};
ylimsat= [0 .1; 0 0.5 ];


for iparam=1:2 % mean and std of gaussian:

    paramData = squeeze(GFX_gauss_Cond_xparam(:, :,iparam)); % all 3 walk speeds:
subplot(2,3, usesubspots(iparam));

% barM = mean(paramData,1);
barM = paramData(ippant,:);

bh=bar(1:3, [barM(1), nan, nan]); hold on; % 'all'
bh.FaceColor= 'k';
bh=bar(1:3, [nan, barM(2), nan]); hold on; % 'slow'
bh.FaceColor= 'b';
bh=bar(1:3, [nan,nan, barM(3)]); hold on; % 'normal'
bh.FaceColor= 'r';

ylim(ylimsat(iparam,:))
%%
ylabel(titlesAre{iparam});
title(['Gaussian ' titlesAre{iparam}]);
hold on;

xlim([0 4])
xlabel('walk speed');
set(gca,'XTickLabel', {'all','slow', 'fast'})
end
%%
cd(figdir)
cd('PFX Gaussian fits');
print('-dpng', ['Participant ' num2str(ippant) '_compareSlowFast_Gaussparams']);
clf;
%% Now for the gait split:

figure(2);clf
set(gcf,'units', 'normalized', 'Position', [.1 .1 .8 .8], 'color', 'w');
% load("GFX_SOAdata.mat");

% 

% change colormap
slowCols = cbrewer('seq', 'Blues', 10);
slowCols = slowCols(5:10,:);

fastCols = cbrewer('seq', 'Reds', 10);
fastCols = fastCols(5:10,:);

allCols = cbrewer('seq', 'Greys', 10);
allCols = allCols(5:10,:);

% note that we have 3 types of alignment.
% Can base on gait quartile on first stim, (1)
% Can base on vis stim gait quartile, (2);
% Can base on aud stim gait quartile, (3);
% Can base on resp gait quartile, (4);

% cycle through to plot each:
alignmentWas= {'First stim', 'Vis stim', 'Aud stim', 'Last stim', 'Response'};
alignFlds= {'_FirstAligned', '_VisAligned', '_AudAligned', '_LastAligned', '_RespAligned'};
%

nGaitQ = size(GFX_Data_SOAs.propSameByGait_slow_FirstAligned,2);
% GFX_gauss_Cond_x_Gait_x_Yfit

for iAlign= 1:length(alignFlds)

figure(2);clf
% plot slow and normal speedS:
% plotData= {squeeze(GFX_propSamebyGaitandSOA_slow(:,iAlign,:,:)),...
%     squeeze(GFX_propSamebyGaitandSOA_fast(:,iAlign,:,:))};

plotData= { GFX_Data_SOAs.(['propSameByGait_all' alignFlds{iAlign}]),...
    GFX_Data_SOAs.(['propSameByGait_slow' alignFlds{iAlign}]),...
    GFX_Data_SOAs.(['propSameByGait_fast' alignFlds{iAlign}])};

nTrialsall= {GFX_Data_SOAs.nThisSOAthisGaitPhase_all,...
    GFX_Data_SOAs.nThisSOAthisGaitPhase_slow,...}; 
    GFX_Data_SOAs.nThisSOAthisGaitPhase_fast};

useCols= {allCols,slowCols, fastCols};

ispots = [1, 4, 7];
for iplot = 1:3
    pltd= plotData{iplot};
subplot(3,3,ispots(iplot));

%  pM= squeeze(nanmean(pltd,1));
pM = squeeze(pltd(ippant,:,:));

 % result is a nGaitQ x nSOA

legh=[];
legS={};
 for igait= 1:nGaitQ

legh(igait)=plot(SOAs, pM(igait,:), 'o-', 'color', [useCols{iplot}(igait,:)], ...
    'MarkerFaceColor', [useCols{iplot}(igait,:)], 'MarkerSize',2,...
    'LineWidth',2);
hold on;


legS{igait} = ['q' num2str(igait)];



 end
 
 % also show min to max targ count per q.
 nTrials = squeeze(nTrialsall{iplot}(ippant,:,:));
 
text(0.2, .1, ['[min, max] = ']);
 text(0.2, .05, [num2str(min(nTrials(:))) ', ' num2str(max(nTrials(:)))]); 

    

 legend(legh, legS, 'location', 'south')
text(-0.4, .9, ['gait based on ']);
text(-0.4, .8, ['\bf' alignmentWas{iAlign}]);
%%
xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title(['Participant ' num2str(ippant) ' results']);
set(gca,'fontsize',12);
end
%%
% now plot the gauss fits,
% load('GFX_GaussData.mat');
% include group data for comparison?
% [nsubs,nconds, ngaits, nsamps]= size(GFX_gauss_Cond_x_Gait_Yfit);
%%
% [nsubs,nspeeds,naligns,ngaits,SOAS]=size(GFX_gauss_Cond_x_Gait_x_Yfit);

for icond= 1:3%:3 % first cond is 'all', slow , fast.
subplot(3,3,ispots(icond)+1);
cla
% GFX_gauss_Cond_x_Gait_x_Yfit
    plotme = squeeze(GFX_gauss_Cond_x_Gait_x_Yfit(:,icond,iAlign,:,:));

    for igait=1:nGaitQ
        
%     pM = squeeze(mean(plotme(:,igait,:),1));    
    pM = squeeze(plotme(ippant,igait,:));
    
    plot(gauss_xvec, pM, '-','color',[useCols{icond}(igait,:)], 'linew',2)
    hold on;
    end
    xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title('Individual gauss fits');
set(gca,'fontsize',12);
end
%%
%%quick stats check
% [nsubs, nconds, param(mean,sd)]= size(GFX_gauss_Cond_xparam)

% % now for the gait splits:
usesubspots= [3,6,9];
titlesAre = {'Mean (mu)', 'STD (sig)'};
ylimsat= [0 .2; 0 0.6 ];

% cla
% useppants= [1:21];
%%

%prespecify bar locations:
ngroups = 3;
nbars = 5;
xall=[];
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
xall(:,i)= x;
end
%%
for iparam=1:2
paramData = squeeze(GFX_gauss_Cond_x_Gait_param(:, :, iAlign,:,iparam)); % speeds: all, slow, fast

subplot(2,3, usesubspots(iparam));

% barM = squeeze(mean(paramData,1));
barM = squeeze(paramData(ippant,:,:));
% stE = CousineauSEM(paramData);
% each separately:


% all combined:
bh=bar([barM(1,:); nan(1,nGaitQ); nan(1,nGaitQ)], 'FaceColor', 'flat');
hold on;
for ibar=1:nGaitQ
bh(ibar).FaceColor= allCols(ibar,:);
end


% slow
bh=bar([nan(1, nGaitQ); barM(2,:);nan(1, nGaitQ) ], 'FaceColor', 'flat');
hold on;
for ibar=1:nGaitQ
bh(ibar).FaceColor= slowCols(ibar,:);
end


% fast
bh=bar([nan(1, nGaitQ); nan(1, nGaitQ); barM(3,:) ], 'FaceColor', 'flat');
hold on;
for ibar=1:nGaitQ
bh(ibar).FaceColor= fastCols(ibar,:);
end

%%
ylim(ylimsat(iparam,:))
ylabel(titlesAre{iparam});
title(['Gaussian ' titlesAre{iparam}]);

xlim([0 4])
ylim
xlabel('walk speed');
set(gca,'xtick', 1:3, 'XTickLabel', {'all','slow', 'fast'})
% axis tight

my_scaleYrange(2, barM);

end
%%
cd(figdir)
cd('PFX Gaussian fits');
print('-dpng', ['Participant ' num2str(ippant) '_compareSlowFast_XGait_alignedto_' alignmentWas{iAlign}]);

end
end % ppant

