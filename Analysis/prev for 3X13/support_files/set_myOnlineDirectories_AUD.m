% set_myOnlineDirectories


% homedir = '/MATLAB Drive/Auditory Detection Task 1/FROM MATT';
% homedir='C:\Users\mdav0285\Documents\GitHub\Analysis-AuditoryDetection-v1';


%UTS mac 
cd('/Users/164376/Documents/GitHub/Analysis-AuditoryDetection-v1/');
homedir = pwd;
addpath(genpath([homedir filesep 'Analysis']));
addpath(genpath([homedir filesep 'Raw_Data']));
addpath(genpath([homedir filesep 'Processed_Data']));
addpath(genpath([homedir filesep 'Figures']));

rawdatadir=  [homedir filesep 'Raw_Data'];
savedatadir= [homedir filesep 'Processed_Data'];
figdir = [homedir filesep 'Figures'];