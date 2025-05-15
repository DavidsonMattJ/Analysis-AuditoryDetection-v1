% plotjobs_Data_withinGait_basic -

% loads the ppant data collated in J4

%various plot jobs (with descriptions below)

% adapted from the visual detection experiments,
% this version for Auditory
% detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AUDITORY DETECTION V1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%

cd(savedatadir);
pfols= dir([pwd  filesep '*_data.mat']);
nsubs= length(pfols);
procdatadir=savedatadir;
%%
%show ppant list:
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)

%%
% Participant level effects (PFX) 
% plots of raw distributions (no binning).
job.plot_TargOnspos =0; % target onset relative to gait.
job.plot_RespOnspos =0; % response onset relative to gait (includes FAs)

job.plot_RT_perTargpos=0; % RT per targ onset position
job.plot_RT_perResppos=0; % as above, but response .

job.plot_Acc_perTargpos=0;


%%%%
%plot group level effects (GFX):
job.plotGFX_Targcount_perGaitpos=0; %  distribution of targ onsets.
job.plotGFX_Respcount_perGaitpos=0; % distribution of resp onsets.



job.plotGFX_TargContrastPos=0; % single gc only.
job.plotGFX_TargRT_perGaitpos_binned=0; % binned versions
job.plotGFX_RespRT_perGaitpos_binned=0; %
job.plotGFX_TargAcc_perGaitpos_binned=1; %

%% 
% job.plotGFX_eyetimeseries =0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data wrangling first: j4 and j5 need to be completed:
%%

    %%
    disp(['loading GFX data into workspace...'])
    cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat');



%% plotting jobs
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plots at Participant level:
%colours for left / right foot examples.


%PFX
if job.plot_TargOnspos % target onset relative to gait.
 %% set up figure:


 %for each ppant, plot the distribution of targ onset positions:
 % pass in some details needed for accurate plots:
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.usebin=1;
 cfg.figdir=figdir;
 cfg.plotlevel='PFX';
 % cycles through ppants, plots with correct labels.
 plot_onsetDistribution(GFX_TargPosData, cfg);
 
end

%% PFX
if job.plot_RespOnspos % Response onset relative to gait.
 %%
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Response';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.usebin=0; % raw distribution, or binned versions.
 

 cfg.figdir=figdir;
 cfg.plotlevel='PFX';
 % cycles through ppants, plots with correct labels.
 plot_onsetDistribution(GFX_RespPosData, cfg);
 
end
%%


if job.plot_RT_perResppos==1
 % plot the RT in each bin, calculated in concatGFX above.
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Response';
 cfg.DV = 'RT';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.usebin=1; 
 cfg.figdir=figdir;
 cfg.plotlevel = 'PFX'; % plot separate figures per participant
 
    plot_GaitresultsBinned(GFX_RespPosData, cfg);
   %% 
end%
    %
% PFX plotting:
if job.plot_RT_perTargpos==1
 % plot the RT in each bin, calculated in concatGFX above.
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.DV = 'RT';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel = 'PFX'; % plot separate figures per participant
    cfg.usebin=1;
    cfg.figdir=figdir;
    plot_GaitresultsBinned(GFX_TargPosData, cfg);
   %% 
end%

% PFX plotting:
if job.plot_Acc_perTargpos==1
 % plot the RT in each bin, calculated in concatGFX above.
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.DV = 'Accuracy';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.plotShuff=0;
 cfg.usebin=1;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel = 'PFX'; % plot separate figures per participant
        
    cfg.figdir=figdir;
    plot_GaitresultsBinned(GFX_TargPosData, cfg);
   %% 
end%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GFX : group effects plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the position of all target onsets (regardless of outcome), as  % of gait
% cycle:


if job.plotGFX_TargContrastPos==1
    %%
    cfg=[];
    cfg.subjIDs = subjIDs;
    cfg.type = 'Target';
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel='GFX';
    cfg.figdir=figdir;
    plot_contrastDistribution(GFX_TargPosData, cfg);
end
%GFX

if job.plotGFX_Targcount_perGaitpos==1   
    % helper function, shared script for both targ onset, and resp onset data
    
  %% plot the distribution of targ onset positions:
 % pass in some details needed for accurate plots:
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Auditory stim'; % for y axis
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.usebin=1;
 cfg.figdir=figdir;
 cfg.plotlevel='GFX';
 % cycles through ppants, plots with correct labels.
%  plot_onsetDistribution(GFX_TargPosData, cfg); 
 plot_AUDonsetDistribution(GFX_TargPosData, cfg); 
end
%% GFX
% plot the position of all responses (regardless of outcome), as  % of gait
% cycle:
if job.plotGFX_Respcount_perGaitpos==1   
    %%
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Response';
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.usebin=1; % raw distribution, or binned versions.
 
 cfg.figdir=figdir;
 cfg.plotlevel='GFX';
 % cycles through ppants, plots with correct labels.
 plot_AUDonsetDistribution(GFX_RespPosData, cfg);
end 
%% plot accuracy based on target onset
if job.plotGFX_TargAcc_perGaitpos_binned==1 %
    %%
  cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.DV = 'dprime';
 cfg.HeadData= GFX_headY;
 cfg.usebin=1;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel = 'GFX'; % plot separate figures per participant
 cfg.norm=0;
 cfg.plotShuff=0;
 cfg.normtype= 'relative'; %['absolute', 'relative', ''relchange', 'normchange', db']
 cfg.ylims = [-.10 0.10]; % if normalized
    cfg.figdir=figdir;
%  cfg.plotppants=[2:22];
 
 plot_AUDGaitresultsBinned(GFX_TargPosData, cfg);
    
end
% if job.plotGFX_RespAcc_perGaitpos_binned==1 %
%     %%
%   cfg=[];
%  cfg.subjIDs = subjIDs;
%  cfg.type = 'Response';
%  cfg.DV = 'Accuracy';
%  cfg.datadir= datadir; % for orienting to figures folder
%  cfg.HeadData= GFX_headY;
%  cfg.pidx1= pidx1;
%  cfg.pidx2= pidx2;
%  cfg.plotlevel = 'GFX'; % plot separate figures per participant
%  cfg.norm=1;
%  cfg.normtype='relchange';
%  cfg.plotShuff=0;
%  cfg.ylims= [-.05 .05];
%     plot_GaitresultsBinned(GFX_RespPosData, cfg);
%     
% end


if job.plotGFX_TargRT_perGaitpos_binned==1 %
    %%
  cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.DV = 'RT';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel = 'GFX'; % plot separate figures per participant
 cfg.norm=0; % already z scored, so don't tweak.
 cfg.ylims = [-.15 .15]; % if norm =0;
 cfg.plotShuff=0;
 cfg.usebin=1;
 cfg.normtype= 'absolute';   
    plot_AUDGaitresultsBinned(GFX_TargPosData, cfg);
    
end
% %%
% if job.plotGFX_RespRT_perGaitpos_binned==1 %
%     %%
%   cfg=[];
%  cfg.subjIDs = subjIDs;
%  cfg.type = 'Response';
%  cfg.DV = 'RT';
%  cfg.datadir= datadir; % for orienting to figures folder
%  cfg.HeadData= GFX_headY;
%  cfg.pidx1= pidx1;
%  cfg.pidx2= pidx2;
%  cfg.plotlevel = 'GFX'; % plot separate figures per participant
%  cfg.norm=1; % zscored, so don't tweak.
%  cfg.normtype='relative';
%  cfg.plotShuff=0;
%  cfg.ylims = [-.1 .1]; % if norm =0;
%  cfg.usebin=1;
%  cfg.figdir=figdir;
% 
%  plot_GaitresultsBinned(GFX_RespPosData, cfg);
%     
% end

%%




if job.plotGFX_eyetimeseries==1

    plot_Gaittimeseries_Eye(GFX_headY, GFX_EyeD)



end