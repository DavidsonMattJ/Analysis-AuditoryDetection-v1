% plot GFX gauss data:

% plot_GFX_Gauss

job.plot_walkSpeedComparison =1; % focuses on the gaussian params comparing 2 walk speeds
job.plot_walkQuantileComparison=0; % focuses oon g params over the gait quantiles.

job.plot_walkQuantileComparison_SOAview=0; %focuses on changes within an SOA.

job.plot_PropSame_perPcnt_perAlign=0;

% useppants= [1:16,18:19]; % 17 and 20 have poor gaussian fits.
cd(savedatadir);
%
figure(1);clf
set(gcf,'units', 'normalized', 'Position', [.1 .1 .8 .6], 'color', 'w');
load('GFX_SOAdata.mat');

% plot slow and normal speedS:
useCols= {'b', 'r', 'k'};
legh=[];
SOAs= GFX_Data_SOAs.SOAs(1,:);

if job.plot_walkSpeedComparison==1
plotData= {GFX_Data_SOAs.propSame_slow, GFX_Data_SOAs.propSame_fast,  GFX_Data_SOAs.propSame_all}; 
for iplot = 1:2
    pltd= plotData{iplot};
subplot(131);
 pM= squeeze(mean(pltd,1));
legh(iplot)=plot(SOAs, pM, 'o-', 'color', useCols{iplot}, ...
    'MarkerFaceColor', useCols{iplot}, 'MarkerSize',2,...
    'LineWidth',2);
hold on;
stE= CousineauSEM(pltd);
errorbar(SOAs, pM, stE, 'Color', useCols{iplot}, 'Linestyle', 'none', 'LineWidth',2);
end
xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title('Group level results');
set(gca,'fontsize',12);
%%
% now plot the gauss fits,
load('GFX_GaussData.mat');
subplot(132); cla
useCols= {'k', 'b', 'r'};
% include group data for comparison?
[nsubs,nconds, nsamps]= size(GFX_gauss_Cond_x_Yfit);
%

for icond= 2:3%:3 % first cond is 'all'
    plotme = squeeze(GFX_gauss_Cond_x_Yfit(:,icond,:));

    pM = mean(plotme,1);
    stE= CousineauSEM(plotme);

%     shadedErrorBar(gauss_xvec, pM, stE, {'-','color',useCols{icond}, 'linew',1}, 1)
    plot(gauss_xvec, pM, '-','color', useCols{icond}, 'linew',2)
    hold on;
end
shg
%%
% pvals=[];
% for isamp= 1:length(gauss_xvec)
% 
%     [h, pvals(isamp)] = ttest(GFX_gauss_Cond_x_Yfit(:,2,isamp),GFX_gauss_Cond_x_Yfit(:,3,isamp));
% 
%     if pvals(isamp) < .05;
%         text(gauss_xvec(isamp), 0.1, '*', 'color', 'k', 'FontSize',20)
%     end
% end
%%

text(SOAs(1), 0.05, 'p < .05')
xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title('Average of individual gauss fits');
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
barM = mean(paramData,1);
stE = CousineauSEM(paramData);

bh=bar(1:3, [barM(1), nan, nan]); hold on; % 'all'
bh.FaceColor= 'k';
bh=bar(1:3, [nan, barM(2), nan]); hold on; % 'slow'
bh.FaceColor= 'b';
bh=bar(1:3, [nan,nan, barM(3)]); hold on; % 'normal'
bh.FaceColor= 'r';


errorbar(1:3, barM, stE, 'LineWidth',2, 'Color','k', 'LineStyle','none');
ylim(ylimsat(iparam,:))
%%
[h,pval, ci,stats]= ttest(paramData(:,2), paramData(:,3)); % compare slow and fast speeds

txtmsg = ['\itt\rm(' num2str(stats.df) ')=' sprintf('%.2f',stats.tstat) ', \itp \rm = ' sprintf('%.3f', pval) ];
text(2.5, diff(ylimsat(iparam,:))*.9, txtmsg, 'HorizontalAlignment','center');
ylabel(titlesAre{iparam});
title(['Gaussian ' titlesAre{iparam}]);
hold on;

sc=scatter(repmat(2, [1,size(paramData,1)]), paramData(:,2),12);
sc.MarkerEdgeColor= 'k';
hold on;
sc=scatter(repmat(3, [1,size(paramData,1)]), paramData(:,3),12);
sc.MarkerEdgeColor= 'k';
xlim([0 4])
xlabel('walk speed');
set(gca,'XTickLabel', {'all','slow', 'fast'})
end
%%
cd(figdir)
cd('GFX Gaussian fits');
print('-dpng', ['GFX_compareSlowFast_Gaussparams']);
end


if job.plot_walkQuantileComparison==1
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
alignFlds= {'_FirstAligned', '_VisAligned', '_AudAligned', '_LastAligned','_RespAligned'};
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

useCols= {allCols,slowCols, fastCols};

ispots = [1, 4, 7];
for iplot = 1:3
    pltd= plotData{iplot};
subplot(3,3,ispots(iplot));
 pM= squeeze(nanmean(pltd,1));
 % result is a 4 x 7

legh=[];
legS={};
 for igait= 1:nGaitQ

legh(igait)=plot(SOAs, pM(igait,:), 'o-', 'color', [useCols{iplot}(igait,:)], ...
    'MarkerFaceColor', [useCols{iplot}(igait,:)], 'MarkerSize',2,...
    'LineWidth',2);
hold on;

stE= CousineauSEM(squeeze(pltd(:, igait,:)));
errorbar(SOAs, pM(igait,:), stE, 'Color', [useCols{iplot}(igait,:)], 'Linestyle', 'none', 'LineWidth',2);

legS{igait} = ['q' num2str(igait)];
 end
 legend(legh, legS, 'location', 'south')
text(-0.4, .9, ['gait based on ']);
text(-0.4, .8, ['\bf' alignmentWas{iAlign}]);
%%
xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title('Group level results');
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
    pM = squeeze(mean(plotme(:,igait,:),1));
    stE= CousineauSEM(squeeze(plotme(:,igait,:)));

%     shadedErrorBar(gauss_xvec, pM, stE, {'-','color',[useCols{icond-1}(igait,:)]}, 1)
    plot(gauss_xvec, pM, '-','color',[useCols{icond}(igait,:)], 'linew',2)
    hold on;
    end
    xlim([-.45 .45]);
ylim([0 1]);
ylabel('proportion "Same" ')
xlabel('SOA {AV <--> VA}')
title('Average of individual gauss fits');
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
nbars = nGaitQ;
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
barM = squeeze(mean(paramData,1));
% stE = CousineauSEM(paramData);
% each separately:


% all combined:
bh=bar([barM(1,:); nan(1,nGaitQ); nan(1,nGaitQ)], 'FaceColor', 'flat');
hold on;
for ibar=1:nGaitQ
bh(ibar).FaceColor= allCols(ibar,:);
end

% test LME: % includes param for x location (incase plotting).
statsOUT_all = plotGoF_gauss(squeeze(paramData(:,1,:)), xall(1,:));

% slow
bh=bar([nan(1, nGaitQ); barM(2,:);nan(1, nGaitQ) ], 'FaceColor', 'flat');
hold on;
for ibar=1:nGaitQ
bh(ibar).FaceColor= slowCols(ibar,:);
end

statsOUT_slow = plotGoF_gauss(squeeze(paramData(:,2,:)), xall(2,:));

% fast
bh=bar([nan(1, nGaitQ); nan(1, nGaitQ); barM(3,:) ], 'FaceColor', 'flat');
hold on;
for ibar=1:nGaitQ
bh(ibar).FaceColor= fastCols(ibar,:);
end

statsOUT_fast = plotGoF_gauss(squeeze(paramData(:,3,:)), xall(3,:));

%%
stE= [CousineauSEM(squeeze(paramData(:,1,:)));...
    CousineauSEM(squeeze(paramData(:,2,:)));...
    CousineauSEM(squeeze(paramData(:,3,:)))]; % calculate within speeds.
eh=errorbar_groupedfit(barM, stE);
% errorbar(1:2, barM, stE, 'LineWidth',2, 'Color','k', 'LineStyle','none');
ylim(ylimsat(iparam,:))
% 
% [h,pval, ci,stats]= ttest(paramData1, paramData2);
% txtmsg = ['\itt\rm(' num2str(stats.df) ')=' sprintf('%.2f',stats.tstat) ', \itp \rm = ' sprintf('%.2f', pval) ];
% text(1.5, diff(ylimsat(iparam,:))*.9, txtmsg, 'HorizontalAlignment','center');
ylabel(titlesAre{iparam});
title(['Gaussian ' titlesAre{iparam}]);
% hold on;
% 
% sc=scatter(repmat(1, [1,nsubs]), paramData1,12);
% sc.MarkerEdgeColor= 'k';
% hold on;
% sc=scatter(repmat(2, [1,nsubs]), paramData2,12);
% sc.MarkerEdgeColor= 'k';
xlim([0 4])
ylim
xlabel('walk speed');
set(gca,'xtick', 1:3, 'XTickLabel', {'all','slow', 'fast'})
% axis tight

my_scaleYrange(2, barM);

end
%%
cd(figdir)
cd('GFX Gaussian fits');

print('-dpng', ['GFX_compareSlowFast_XGait_alignedto_' alignmentWas{iAlign}]);

% export for JASP:
%%
jtable =table();
speedsAre= {'all','slow', 'normal'};
psare= {'_mean', '_sd'};

for ispeed= 1:3


    for imusig= 1:2
        for igait = 1:nGaitQ

            tmpData = GFX_gauss_Cond_x_Gait_param(:, ispeed, iAlign, igait, imusig);

            jtable.([speedsAre{ispeed} psare{imusig} '_g' num2str(igait)]) = tmpData;
        end
    end
end
%%
writetable(jtable,['Gauss_GaitParams_' alignmentWas{iAlign}   '.csv'])
end


end
%%

if job.plot_walkQuantileComparison_SOAview==1; %focuses on changes within an SOA.


    % similar to the above, except now instead of plotting Gaussians, we
    % plot change in propSame over the gait.

figure(2);clf
set(gcf,'units', 'normalized', 'Position', [.1 .1 .8 .8], 'color', 'w');
% load("GFX_SOAdata.mat");

% 

% change colormap
slowCols = cbrewer('seq', 'Blues', 10);
slowCols = slowCols(4:10,:);

fastCols = cbrewer('seq', 'Reds', 10);
fastCols = fastCols(4:10,:);

allCols = cbrewer('seq', 'Greys', 10);
allCols = allCols(4:10,:);

useCols= {allCols,slowCols, fastCols};

spdCols={'k', 'b', 'r'};
% note that we have 3 types of alignment.
% Can base on gait quartile on first stim, (1)
% Can base on vis stim gait quartile, (2);
% Can base on aud stim gait quartile, (3);
% Can base on resp gait quartile, (4);

% cycle through to plot each:
alignmentWas= {'First stim', 'Vis stim', 'Aud stim', 'Last stim', 'Response'};
alignFlds= {'_FirstAligned', '_VisAligned', '_AudAligned', '_LastAligned','_RespAligned'};
%

nGaitQ = size(GFX_Data_SOAs.propSameByGait_slow_FirstAligned,2);
% GFX_gauss_Cond_x_Gait_x_Yfit
for iAlign= 1:length(alignFlds)
% plot slow and normal speedS:
% plotData= {squeeze(GFX_propSamebyGaitandSOA_slow(:,iAlign,:,:)),...
%     squeeze(GFX_propSamebyGaitandSOA_fast(:,iAlign,:,:))};

plotData= { GFX_Data_SOAs.(['propSameByGait_all' alignFlds{iAlign}]),...
    GFX_Data_SOAs.(['propSameByGait_slow' alignFlds{iAlign}]),...
    GFX_Data_SOAs.(['propSameByGait_fast' alignFlds{iAlign}])};

icounter=1;

for avSOAs=1:8

    if avSOAs==8
        averageSOAs= 1:7; % all
    else
        % a single index.
        averageSOAs= avSOAs;
    end
    figure(2);clf
for i=1:nGaitQ
xstring{i} = ['Q' num2str(i)];
end
    for iplot = 1:3

    
    
    pltd= plotData{iplot};
 pM= squeeze(nanmean(pltd,1));
 % result is a 4 x 7

legh=[];
legS={};
 for iSOA= 1:length(SOAs)
     subplot(2,7,iSOA); hold on;
    legh(iSOA)= plot(1:nGaitQ,pM(:,iSOA)' , 'o:', 'color', [spdCols{iplot}], ...
    'MarkerFaceColor', [spdCols{iplot}], 'MarkerSize',1,...
    'LineWidth',2);
hold on;

stE= CousineauSEM(squeeze(pltd(:, :,iSOA)));
errorbar(1:nGaitQ, pM(:,iSOA), stE, 'Color', [spdCols{iplot}], 'Linestyle', 'none', 'LineWidth',2);

legS{iSOA} = [ num2str(SOAs(iSOA))];
 
axis tight
ylim([0 1])
xlim([0 nGaitQ+1])
set(gca,'xtick', [1:nGaitQ], 'xticklabels', xstring);
xlabel('Gait quartile');
title({['SOA {AV <--> VA}'];[ num2str(SOAs(iSOA)) ' sec']});
% my_scaleYrange(2, pM(:,iSOA));
set(gca,'fontsize',12);

ylabel('proportion "Same" ')

 end
 shg
%  legend(legh, legS, 'location', 'south')
text(0.5, .85, ['Gait alignment based on '], 'fontsize', 12);
text(0.5, .8, ['\bf' alignmentWas{iAlign}], 'fontsize', 12);
%%
subplot(2,7,[8:9]+ 2*(iplot-1));
hold on;
% plot the change in propSame overall. 
pData =squeeze(nanmean(pltd(:,:,averageSOAs),3)); %mean over SOAs.
% norm per ppant.
pData= pData  - repmat(mean(pData,2), [1,nGaitQ]);

plot(1:nGaitQ, mean(pData,1), '-o','color',spdCols{iplot})
stE= CousineauSEM(pData);
errorbar(1:nGaitQ, mean(pData,1), stE, 'LineStyle','none', 'color',spdCols{iplot}, 'linew',2)
% test GoF
statsOUT = plotGoF_gauss(pData, 1:nGaitQ);

xlim([0 nGaitQ+1])
shg
xstring={};

set(gca,'xtick', [1:4], 'xticklabels', xstring);
xlabel('Gait quartile');
hold on; plot(xlim , [0 0], 'k:');
ylabel('Normalised change to "same"')
% ylim([-.1 .1])

if length(averageSOAs)==1
title(['Average of SOAs: ' num2str(avSOAs)])
else
title(['Average of SOAs: all'])
end
% test gof?

set(gca,'xtick', [1:4], 'xticklabels', xstring);
xlabel('Gait quartile');
hold on; plot(xlim , [0 0], 'k:');
ylabel('Normalised change to "same"')
% axis tight
end
%%

%%
cd(figdir)
cd('GFX Gaussian fits');
% shg
if length(averageSOAs)==1
print('-dpng', ['GFX_compareEachSOA_XGait_alignedto_' alignmentWas{iAlign} '_SOA'  num2str(avSOAs)]);
else
print('-dpng', ['GFX_compareEachSOA_XGait_alignedto_' alignmentWas{iAlign} '_SOAsAll']);
end
end % avSOAs

%%
% 
% jtable =table();
% speedsAre= {'all','slow', 'normal'};
% 
% for ispeed= 1:3
% tmpData= plotData{ispeed};    
%         for igait = 1:nGaitQ
%             for iSOA=1:length(SOAs)
% 
%             jtable.([speedsAre{ispeed} '_g' num2str(igait) '_SOA' num2str(iSOA)  alignFlds{iAlign}]) = tmpData(:,igait,iSOA);
%             end
%         end
%    
% end
% %%
% writetable(jtable,['PropSame_byGaitQ by ' alignmentWas{iAlign}   '.csv'])
% % export for JASP:

end% ialign
%%
end


if job.plot_PropSame_perPcnt_perAlign==1;

%%
figure(2);clf
set(gcf,'units', 'normalized', 'Position', [.1 .1 .8 .8], 'color', 'w');
% load("GFX_SOAdata.mat");

%% 

% change colormap
slowCols = cbrewer('seq', 'Blues', 10);
slowCols = slowCols(4:10,:);

fastCols = cbrewer('seq', 'Reds', 10);
fastCols = fastCols(4:10,:);

allCols = cbrewer('seq', 'Greys', 10);
allCols = allCols(4:10,:);

useCols= {allCols,slowCols, fastCols};

spdCols={'k', 'b', 'r'};
% note that we have 3 types of alignment.
% Can base on gait quartile on first stim, (1)
% Can base on vis stim gait quartile, (2);
% Can base on aud stim gait quartile, (3);
% Can base on resp gait quartile, (4);

% cycle through to plot each:
alignmentWas= {'First stim', 'Vis stim', 'Aud stim', 'Last stim', 'Response'};
alignFlds= {'_FirstAligned', '_VisAligned', '_AudAligned', '_LastAligned','_RespAligned'};
%%
clf
for iAlign= 2%1:length(alignFlds);
 
  hold on
  % all speeds:
    plotData= { GFX_Data_SOAs.(['propSameByGaitPcnt_perSOA__all' alignFlds{iAlign}]),...
        GFX_Data_SOAs.(['propSameByGaitPcnt_perSOA__slow' alignFlds{iAlign}]),...
        GFX_Data_SOAs.(['propSameByGaitPcnt_perSOA__fast' alignFlds{iAlign}])};
for iSOA=1:7
    for iplot = 1:3
         subplot(1,7,iSOA); hold on
    pltd= squeeze(plotData{iplot}(:,iSOA,:));
 pM= squeeze(nanmean(pltd,1));
 stE= CousineauSEM(pltd);
%  shadedErrorBar(1:length(pM), smooth(pM,20), smooth(stE,20), {'color', spdCols{iplot}},1)
plot(1:length(pM), pM, 'color', spdCols{iplot})
ylim([0 1])
    end %plotspd
end %iSOA
end % align
shg

end % job
    
    
    
    
    