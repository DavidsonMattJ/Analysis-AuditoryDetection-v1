%% plots the fourier series fits, for observed and shuffled data (PFX and GFX)
%note that previous data wrangling must be performed:
% j5_create null GFX (for testing) - creates null for both PFX and GFX
% j8_testfourier_Obs_Null_GFX(precompute for plots)
% j
% j...Gshuffled data is prepared in j5, and fits to shuffled data in
%j8.




%% Choose an options:
%1 = RTs relative to target onset
%2 = RTs relative to response onset
%3 = Accuracy relative to target onset
%4 = Likelihood relative to response onset


% load data
cd([procdatadir filesep 'GFX'])
%load obs and shuffld data for plots:
load('GFX_Data_inGaits.mat')
load('GFX_Data_inGaits_FourierFits.mat')
%%

job.plotGFX_summaryforDVs =1; % used in posters/slides (1 x 3, with group fit at end)

job.plotGFX_fitsOverlayed=0 ; %extra figure to compare timing of GFX fits (overlayed)/

job.plotGFX_groupFits_MSver1 =0; % MS version. trial Data, targOns etc.

job.plotGFX_groupFits_MSver2 =1; % MS version.restricted to dual gait

% testtype=1;
usebin=1; %to do: adapt script to also fit fourier to raw data.


%% figure set up:
% rts = red
% acc = blue
% resps = purp

barCols= {'r',  'b', [.7 0 .7], [.2 .2 .2]};  
fntsize= 18;

fitsF=[]; % for summary plot (after for loop).
clf; 
pc=1; %plotcounter
%
% figure(1); clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0.1 0.1 1 .5]);
% figure(2); clf;
legHZ=[]; % store the figure handles for FITS of our observed data.

%% cycle through plot jobs.
if job.plotGFX_summaryforDVs
%%
   %%PSYC  3X13 ver
testData= {GFX_TargPosData,GFX_TargPosData};
    typeOnsets={'Target', 'Target'};
    typeDVs = { 'Accuracy', 'RT'};

    cfg=[];
    
    cfg.subjIDs = subjIDs;
%     cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel = 'GFX'; % plot separate figures per participant
    cfg.norm=0;
    cfg.ylims = [-.15 .15]; % if norm =0;
    cfg.normtype= 'relative';

    cfg.nGaits_toPlot=1;

    cfg.usebin=1;
    cfg.Hzspace= Hzspace;
    cfg.fntsize= 15;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.plotShuff=1;

%     plot_FourierFitGFX_summary
    plot_GFX_groupFits_vAUD_MSver2(cfg,testData, GFX_FourierNull);


end
 

%% plot summary (overlayed GFX fits).
if job.plotGFX_fitsOverlayed
%%
       testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = {'Accuracy', 'RTs', 'Counts'};

    cfg=[];
    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel = 'GFX'; % plot separate figures per participant
    cfg.norm=0;
    cfg.ylims = [-.15 .15]; % if norm =0;
    cfg.normtype= 'relative';
    cfg.nGaits_toPlot=2;
    cfg.usebin=1;
    for iDV=1:3
        cfg.type = typeOnsets{iDV};
        cfg.DV = typeDVs{iDV};
        cfg.iDV=iDV;
        % UNFINISHED! ! ! ! ! !
        plot_FourierFitGFX_overlayed(cfg,testData{iDV});
            % ! ! ! ! ! ! 
    end


end

if job.plotGFX_groupFits_MSver1 ==1; % MS version. trial Data, targOns etc.
%%
% plot a pretty summary of the group level fourifer fit data.
% similar to the first plot job above, but including an example trial, \
% and trial onsets per gait bin. 
%%
testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};

    cfg=[];
    
    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel = 'GFX'; % plot separate figures per participant
    cfg.norm=0;
    cfg.ylims = [-.15 .15]; % if norm =0;
    cfg.normtype= 'relative';

    cfg.nGaits_toPlot=1;

    cfg.usebin=1;
    cfg.Hzspace= Hzspace;
    cfg.fntsize= 14;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    
%     plot_FourierFitGFX_summary
    plot_GFX_groupFits_MSver1(cfg,testData, GFX_FourierNull);


end

if job.plotGFX_groupFits_MSver2 ==1; % MS version. trial Data, targOns etc.
%%
% plot a pretty summary of the group level fourifer fit data.
% similar to the first plot job above, but including an example trial, \
% and trial onsets per gait bin. 
%%
testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};

    cfg=[];
    
    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel = 'GFX'; % plot separate figures per participant
    cfg.norm=0;
    cfg.ylims = [-.15 .15]; % if norm =0;
    cfg.normtype= 'relative';

    cfg.nGaits_toPlot=1;

    cfg.usebin=1;
    cfg.Hzspace= Hzspace;
    cfg.fntsize= 14;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    
%     plot_FourierFitGFX_summary
    plot_GFX_groupFits_MSver2(cfg,testData, GFX_FourierNull);


end

