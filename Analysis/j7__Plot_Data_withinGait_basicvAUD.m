% j6__Plot_Data_withinGait_basicvAUD-

% loads the ppant data collated in j2_binData_bycycle./
% j3_binDatabylinkedcycles.

%plots the average RT, relative to target position in gait cycle.
%pariticpant, as a position of the gait cycle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% AUDITORY version %%%%%% 
%   2 walking speeds, extra DVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
jobs=[];

set_myOnlineDirectories_AUD;
cd(procdatadir);
pfols= dir([pwd  filesep '*summary_data.mat']);
nsubs= length(pfols);
%
%show ppant list:
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)


    disp(['loading GFX data into workspace...'])
    cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat');

    load('GFX_Data_inGaits_FourierFits.mat');

%%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Note it is all done from the same function. details below.

% !Output is a 3 x 3, separate figures for single step or stride (pfx).


% all foot (L,R,combined), speed (slow,natural,combined) combinations.

%select PFX (or GFX)
%select Target onset in gait (or Response onset).
%select DV at each position :
   % - Accuracy, RT, HR, FA, dprime, crit, counts.
 
%see Figures/ for output.
   
clf
   cfg=[];
   cfg.plotlevel = 'PFX';                   % PFX, GFX plot separate figures per participant
   cfg.usebin    =  1;                      % average within bins (or 1-100% separately)
   cfg.type      = 'Target';                % Target / Response (not supported)
   cfg.DV        = 'Acc';                    % Accuracy/ RT/ dprime / crit/ counts/ FA / HR
   cfg.plotCOL = 'b';
   cfg.fitCOL= 'b';
   cfg.yyaxis = 'left';
   cfg.ispeed=1:3; % limit to single speed?
   cfg.nGaitstoplot = 2;%:2; %[1,2] step or stride.
   cfg.datadir= datadir;              % for orienting to figures folder
   cfg.HeadData= GFX_headY;
   cfg.pidx1= pidx1;
   cfg.pidx2= pidx2; 
   cfg.subjIDs = subjIDs;
   cfg.figdir=figdir;
   cfg.plotShuff=0; % requires job    
   cfg.normON=1;
   cfg.normtype='relchange'; %relative, relchange, normchange, absolute
   cfg.plotHead = 0; % overlay the head trace? 
   cfg.plotFIT = 1;  % apply best Fourier fit? 

   % Note that the PFX style plots in a 3 x 3 (columns for speeds, rows for
   % L/R support). 
   % Addiing a cfg to plot in this way (grid) or focused.
   cfg.plotGaitGrid=0;
   cfg.plotnewLayout=1; % rather than all single or doub gc, just a selection (LR, RL, LRL).

   % plot helper:
   
       plot_GaitresultsBinned_2speed(GFX_TargPosData, GFX_FourierNull,cfg);
   
   % legend off
   % title('')

   %% different version for ACNS (the above version separates by speed, this version has new DVs per panel (1 speed).
if jobs.plotforACNS==1
% for the ACNS version, plot 4 in a row:
clf
  cfg=[];
   cfg.plotlevel = 'GFX';                   % PFX, GFX plot separate figures per participant
   cfg.usebin    =  1;                      % average within bins (or 1-100% separately)
   cfg.type      = 'Target';                % Targer / Response
   cfg.DV_all        = {'Accuracy','RT','dprime','crit'};                    % Accuracy/ RT/ dprime / crit/ counts/ FA / HR
   cfg.plotCOL_all = {'b','r','m',[.5 .5 .9]};
    cfg.fitCOL_all= {'b','r','m',[.5 .5 .9]};
    cfg.yyaxis_all = {'left', 'right', 'left','right'};
    cfg.ispeed=2; % limit to single speed?
   cfg.datadir= datadir;              % for orienting to figures folder
   cfg.HeadData= GFX_headY;
   cfg.pidx1= pidx1;
   cfg.pidx2= pidx2; 
   cfg.subjIDs = subjIDs;
   cfg.figdir=figdir;
    cfg.plotShuff=0;
    
        cfg.add_cpsFit=1;
    cfg.plotHead = 0; % overlay the head trace? 
    cfg.plotFIT = 1;  % apply best Fourier fit? 

%     subsplts
    for iDV=1:4

        cfg.DV = cfg.DV_all{iDV};
        cfg.plotCOL= cfg.plotCOL_all{iDV};
        cfg.fitCOL= cfg.fitCOL_all{iDV};
        cfg.yyaxis= cfg.yyaxis_all{iDV};
        cfg.iDV=iDV;


        if strcmp(cfg.type,'Target')
            plot_GaitresultsBinned_2speed_ACNS(GFX_TargPosData, GFX_FourierNull,cfg);
        else
            plot_GaitresultsBinned_2speed_ACNS(GFX_RespPosData, cfg);
        end
    end

    %% final plot, this one will overlay Vis vs Audio sensitivity, and Vis vs Audio criterion:
    % for the ACNS version, plot 4 in a row:

setmydirs_detectvGabor;
    cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat','GFX_TargPosData','pidx2');

    % we need to rename the group effects as we load two versions:
    GFX_TargPos_Data_gabor = GFX_TargPosData;
    
    GFXNull_gabor = GFX_FourierNull;
    pidx2_gab = pidx2;
%%

    % load the AUD:
     set_myOnlineDirectories_AUD;
   
     cd(procdatadir); cd('GFX');
    load('GFX_Data_inGaits.mat','GFX_TargPosData','pidx1','pidx2', 'subjIDs');
    GFX_TargPos_Data_aud= GFX_TargPosData;

    GFXNull_aud= GFX_FourierNull;
    pidx2_aud=pidx2;
clf
  cfg=[];
   cfg.plotlevel = 'GFX';                   % PFX, GFX plot separate figures per participant
   cfg.usebin    =  1;                      % average within bins (or 1-100% separately)
   cfg.type      = 'Target';                % Targer / Response
   cfg.DV_all        = {'dprime','dprime','crit','crit'};                    % Accuracy/ RT/ dprime / crit/ counts/ FA / HR
    cfg.DV_modalityall= {'visual', 'auditory', 'visual','auditory'};
   cfg.plotCOL_all = {'m', 'm', [.5 .5 .9],[.5 .5 .9]};
    cfg.fitCOL_all= {'m', 'm', [.5 .5 .9],[.5 .5 .9]};
    % cfg.yyaxis_all = {'left', 'right', 'left','right'};
    cfg.yyaxis_all = {'left', 'left', 'left','left'};
  
    cfg.FITstyle= {'-',':','-',':'};

    cfg.ispeed=3; % limit to single speed?
   cfg.datadir= datadir;              % for orienting to figures folder
   cfg.HeadData= GFX_headY;
   cfg.pidx1= pidx1;
   cfg.pidx2= pidx2; 
   cfg.subjIDs = subjIDs;
   cfg.figdir=figdir;
    cfg.plotShuff=0;
    
        cfg.add_cpsFit=0;
    cfg.plotHead = 0; % overlay the head trace? 
    cfg.plotFIT = 1;  % apply best Fourier fit? 

%     subsplts
    for iDV=1:4

        cfg.DV = cfg.DV_all{iDV};
        cfg.plotCOL= cfg.plotCOL_all{iDV};
        cfg.fitCOL= cfg.fitCOL_all{iDV};
        cfg.yyaxis= cfg.yyaxis_all{iDV};
          cfg.DV_modality = cfg.DV_modalityall{iDV};
        cfg.iDV=iDV;
        if iDV==1 || iDV==3 % visual data
            GFX_T = GFX_TargPos_Data_gabor;
            GFX_N = GFXNull_gabor;
            cfg.pidx2= pidx2_gab;

        elseif iDV==2 || iDV==4

            GFX_T = GFX_TargPos_Data_aud;
            GFX_N = GFXNull_aud;
            cfg.pidx2= pidx2_aud;
        end

      
            plot_GaitresultsBinned_2speed_ACNS_overlaid(GFX_T, GFX_N,cfg);
       
    end
end