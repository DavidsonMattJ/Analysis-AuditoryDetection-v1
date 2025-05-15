%% %%%%%% In this script, we will import the raw data, convert to friendly matlab format....etc

% - import 
% - convert to table
% - save in pariticipant folder.

% laptop (mac) directory:
% homedir = '/Users/matthewdavidson/Documents/GitHub/AV-Synchrony-during-locomotion';
%GT directory

%% UPDATED:
% (MD) 
% - renaming scripts to show job order.
% - have changed first import job to repair trial index
%   -also now creates/saves HeadPos struct.



set_myOnlineDirectories_AUD;

jobs=[];

cd(savedatadir);
%%

% 
% Job 1: import
jobs.importandsave_ptable = 1; 

% Preprocess / gait analysis:
%Job 2: Preprocess, extract step start/ends
jobs.SplitHead_bycycle= 1; % Adapted from other analyses (MD)
jobs.AssignWalkPhases = 1;% appends gait-percentile information to trial events.
jobs.epochHead_timeseries=1;
jobs.crunchPFX_concatGFX=1; % main analysis to loop over ppants for eithin stride effects

jobs.plotMain_strideFX =1; 
% jobs.


%%  Perform Calculations of thresholds and effects

%(MD updated for extra analyses)
jobs.calculate_participantlevelPsycho = 1; % psychometric fit per ppant. 

jobs.checkppantvalues = 1; % some plotting (Gabe ver 1)

jobs.plot_summaryData_PSYC3X13= 1; %basic plots to include in assignment




%%%% BEGIN ANALYSIS:

%% import csv and save as matlab table per ppant

if jobs.importandsave_ptable==1
% saves summary and frame x frame files in matlab format.
J1_import; 
end

%% % Perform gait extraction. (run this only once - and not in MATLAB online)
if jobs.SplitHead_bycycle==1
    
    J2_split_bycycle;
    
end

if jobs.AssignWalkPhases == 1
% using the trough times, allocate targ and response onset to %step cycle
%appends to the summary table (new columns)
    J2A_AssignWalkPhases; 
  
end

if jobs.Assign_footsupport  
    J2B_Assign_stepSupport;
    
end
if jobs.AssignStridephases
    % do the same but for 2 linked steps (a full stride).
    % NB: compared to the AV project, we now have the data density for this
    
    J2C_AssignStridephases;
end


if jobs.epochHead_timeseries        
    J3_epoch_gait_timeseries;
end


if jobs.crunchPFX_concatGFX==1
   
J4_crunchPFX_concatGFX; % cycle through and cycle the acc, HR, etc, per gait percent (and binned versions).    
end
    
if jobs.plotMain_strideFX
% see within for separate plot jobs:
    J6__PlotAudData_winGait
end

%% Psychometric fits (no gait split)

if  jobs.calculate_participantlevelPsycho == 1

jobs_calculatePpantPsychometric;

end

if jobs.checkppantvalues == 1
    jobs_ppantcheck;
end

%% basic data analysis - for psyc3x13 assignment
% calculate mean Acc, RT, sensitivity, criterion for both speeds (and
% overall)
if jobs.plot_summaryData_PSYC3X13
    
    JX_3X13_output;
end
% also use (above script), to show psychometric fit.

%
