% plot individual fourier fits (Participant level).
%% plots the fourier series fits, for observed and shuffled data.
%note that shuffled data is prepared in j5, and fits to shuffled data in
%j8.
% j11__Plot_Data_PFX_FourierFits_vGabor
%% Choose an options:
%1 = RTs relative to target onset
%2 = RTs relative to response onset
%3 = Accuracy relative to target onset
%4 = Likelihood relative to response onset


% * NEW AUDITORY Version.

% testtype=1;
usebin=1; %to do: adapt script to also fit fourier to raw data.

set_myOnlineDirectories_AUD;

%
cd([procdatadir filesep 'GFX'])
%load obs and shuffld data for plots:
load('GFX_Data_inGaits.mat')
load('GFX_Data_inGaits_FourierFits.mat')


Hzspace = [0.01:.2:10];

job=[];
% job.plot_PFX_groupsummary =0; % messy figure overlaying all fits, and all forced fits. (one DV at a time).
job.plot_PFX_vAud=1; % plots the forced fits per DV, with histogram (or imagesc) of peak F

% job.plot_PFX_prevalence= 0; % first attempt at summarising prevalence info. needs the individual fits and perms to be performed.
% job.plot_PFX_prevalence_MSver1=0; % 
% 
% job.plot_PFX_prevalence_MSver2=0; %alternate (simplified) version. 
% job.plot_PFX_prevalence_MSver3=0; %alternate (simplified) version, now includes saccade onsets. 
% 
% 
job.plot_PFX_phaseResults= 0; % first attempt at summarising phase data for sig oscillations in pop.
% job.plot_PFX_phaseResults_MSver= 0; % first attempt at summarising phase data for sig oscillations in pop.

%%
if job.plot_PFX_groupsummary


    % % CHOOSE A DV to summarise:
    % testype=1; % 1:4
    % 
    % typeOnsets= {'Target,', 'Target', 'Reponse', 'Target'};
    % typeDVs= {'RT', 'Accuracy', 'Counts', 'Counts'};
    % 
    % if itype==3
    %     dataIN= GFX_RespPosData;
    % else
    %     dataIN= GFX_TargPosData;
    % end
    % 
    % 
    % %%
    % 
    % 
    % cfg=[];
    % cfg.subjIDs = subjIDs;
    % cfg.type = typeOnset;
    % cfg.datadir= datadir; % for orienting to figures folder
    % 
    % cfg.testtype= testtype; % index of DV
    % cfg.DV = typeDV;
    % cfg.HeadData= GFX_headY;
    % cfg.pidx1= pidx1;
    % cfg.pidx2= pidx2;
    % cfg.plotlevel = 'GFX'; % plot separate figures per participant
    % cfg.norm=0; % already z scored, so don't tweak.
    % cfg.ylims = [-.15 .15]; % if norm =0;
    % cfg.normtype= 'relative';
    % cfg.reOrderbyMean= 0; % reorder the data in ascending accuracy?
    % cfg.usebin=1;
    % cfg.fntsize=12;
    % cfg.Hzspace= Hzspace;
    % 
    % % %
    % plot_PFX_FourierFits_overlayed(cfg, dataIN, PFX_FourierNull);
    % 
    % %%
end % JOB


if job.plot_PFX_vAud==1
    %%
    testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Target', 'Target'};
    typeDVs = { 'Accuracy', 'RT', 'dprime', 'crit'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;

    cfg.rankPpants =1 ; %[0,1], 1 to sort imagesc by strength in 2cps

    plot_PFX_fourierFits_verGabor(cfg, testData, PFX_FourierNull);

end

%%

if job.plot_PFX_prevalence==1
    %%
    testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    %%
    plot_PFX_prevalence(cfg, testData, PFX_FourierNull);

end


if job.plot_PFX_prevalence_MSver1==1; % here we plot te

  testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.figdir= figdir;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.fntsize= 12;
    %%
    plot_PFX_prevalence_MSver1(cfg, testData, PFX_FourierNull);


end
%% 

if job.plot_PFX_prevalence_MSver2==1; % here we plot te
%%
  testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.figdir= figdir;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.fntsize= 12;
    %
    plot_PFX_prevalence_MSver2(cfg, testData, PFX_FourierNull);


end
%%

if job.plot_PFX_prevalence_MSver3==1; % here we plot the DVs and saccade 
%%
  testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData, GFX_SaccadePosData};
    typeOnsets={'Target', 'Target', 'Response','Saccade'};
    typeDVs = { 'Accuracy', 'RT', 'Counts','Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.figdir= figdir;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.fntsize= 12;
    %
    plot_PFX_prevalence_MSver3(cfg, testData, PFX_FourierNull);


end


%%








%%


if job.plot_PFX_phaseResults==1
  %%

  % At the moment, plots the individual sine fits for ALL subjects at a
  % fixed freuency. 1 column for 2cps, 1 column for 4cps rows are DVs.
    testData= {GFX_TargPosData,GFX_TargPosData,GFX_TargPosData,GFX_TargPosData};

    typeOnsets={'Target', 'Target', 'Target', 'Target'};
    typeDVs = { 'Accuracy', 'RT', 'dprime','crit'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.figdir= figdir;
    %
    cfg.useShuffdata=0;
    plot_PFX_phaseResults_gabor(cfg, testData, PFX_FourierNull);



end


if job.plot_PFX_phaseResults_MSver==1
    %%
     testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData };
    typeOnsets={'Target', 'Target', 'Response'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.figdir= figdir;
    cfg.fntsize=12;
    cfg.plotShuff=0;
    %
    plot_PFX_phaseResults_MSver1(cfg, testData, PFX_FourierNull);



    
end

if job.plot_PFX_phaseResults_MSver2==1
    %%
     testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData};% ,GFX_SaccadePosData};
    typeOnsets={'Target', 'Target', 'Response'};%,'Saccade'};
    typeDVs = { 'Accuracy', 'RT', 'Counts'};%'Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.figdir= figdir;
    cfg.fntsize=12;
    cfg.plotShuff=0;
    %
    plot_PFX_phaseResults_MSver2(cfg, testData, PFX_FourierNull);



    
end


%% new job. compre saccade frequency with DV frequency:
if job.plot_PFX_saccadeFreqCorr ==1
testData= {GFX_TargPosData,GFX_TargPosData,GFX_RespPosData ,GFX_SaccadePosData};
    typeOnsets={'Target', 'Target', 'Response','Saccade'};
    typeDVs = { 'Accuracy', 'RT', 'Counts','Counts'};

    cfg=[];

    cfg.subjIDs = subjIDs;
    cfg.datadir= datadir; % for orienting to figures folder
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.Hzspace= Hzspace;
    cfg.type = typeOnsets;
    cfg.DV = typeDVs;
    cfg.figdir= figdir;
    cfg.fntsize=12;
    %%
 plot_PFX_saccFreqCorrResults_MSver1(cfg, testData, PFX_FourierNull)

end
