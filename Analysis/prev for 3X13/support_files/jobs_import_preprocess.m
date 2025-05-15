%% %%%%%% In this script, we will import the raw data, convert to friendly matlab format....etc
% - import 
% - convert to table
% - save in pariticipant folder.


jobs=[];

addpath('support_files');




%specify save directory:
savedatadir='/MATLAB Drive/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

jobs.importandsave_ptable= 1;

jobs.calculate_participantlevelSOA =1;
%print the result of the above:
jobs.print_participantlevel_SOA =1;

jobs.concatenate_groupSOA=1;
jobs.print_grouplevel_SOA=1;




% jobs.loadppanttable_calculategaitq= 0;

% ....

%%%% BEGIN ANALYSIS:

%% import csv and save as matlab table per ppant
if jobs.importandsave_ptable==1

jobs_import;

end

if jobs.calculate_participantlevelSOA ==1

jobs_calculateParticipantlevelSynchfn;

end

%
if jobs.jobs.print_participantlevel_SOA ==1

jobs_plot_PFX_SOAs;



end

%% 
% jobs.loadppanttable_calculategaitq= 1;


