% j6__Plot_Data_withinGait_basicvAUDDISCRIM-

% loads the ppant data collated in j2_binData_bycycle./
% j3_binDatabylinkedcycles.

%plots the average RT, relative to target position in gait cycle.
%pariticpant, as a position of the gait cycle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% AUDITORY DISCRIMINATION
%   2 walking speeds, extra DVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % now the auditory task:
    set_myOnlineDirectories_AUD;
   

%%
% load data wrangling first: j4 and j5 need to be completed:
%%

    %%
     cd(savedatadir); cd('GFX');
    load('GFX_Data_inGaits.mat','GFX_TargPosData','pidx1','pidx2', 'subjIDs', 'GFX_headY');
    load('GFX_Data_inGaits_FourierFits.mat')
    disp(['loading GFX data into workspace...'])

    cd(savedatadir)
pfols= dir([pwd  filesep 'p_.mat']);
nsubs= length(pfols);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Note it is all done from the same function. details below.
% output is a 3 x 3, separate figures for single step or stride.
% all foot (L,R,combined), speed (slow,natural,combined) combinations.

%select PFX (or GFX)
%select Target onset in gait (or Response onset).
%select DV at each position :
   % - Accuracy, RT, HR, FA, dprime, crit, counts.
 
%see Figures/ for output.
%    clf
   cfg=[];
   cfg.plotlevel = 'GFX';                   % PFX, GFX plot separate figures per participant
   cfg.usebin    =  1;                      % average within bins (or 1-100% separately)
   cfg.type      = 'Target';                % Targer / Response
   cfg.DV        = 'RT';                    % Accuracy/ RT/ dprime / crit/ counts/ FA / HR
   cfg.plotCOL = 'r';
    cfg.fitCOL= 'r'%[.5 .5 .9];
    cfg.yyaxis = 'right';
    cfg.ispeed=3; % limit to single speed?
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

   if strcmp(cfg.type,'Target')
       plot_GaitresultsBinned_2speed(GFX_TargPosData, GFX_FourierNull,cfg);
   else
       plot_GaitresultsBinned_2speed(GFX_RespPosData, cfg);
   end
   %
   legend off
   title('');
   
%% different version for ACNS (the above version separates by speed, this version has new DVs per panel (1 speed).


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
    cfg.ispeed=3; % limit to single speed?
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
