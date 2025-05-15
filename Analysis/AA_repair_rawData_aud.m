%Repair_rawData_vAUD


set_myOnlineDirectories_AUD;

% datadir='C:\Users\mdav0285\Documents\GitHub\Analysis-detection-walking-gabor\data_Raw\exclusions'; 
cd(datadir)

% cd("from mobilab")
pfols = dir([pwd filesep '*trialsummary.csv']);
nsubs= length(pfols);
tr= table([1:length(pfols)]',{pfols(:).name}' );
disp(tr)
%% %%
% first job is to load the (huge) csv files, and repair column headers.
pfolsframe= dir([pwd filesep '*framebyframe.csv']);
tr= table([1:length(pfolsframe)]',{pfolsframe(:).name}' );
disp(tr)
%%
 filename = pfolsframe(1).name;
 opts = detectImportOptions(filename,'NumHeaderLines',0);
 oldvarnames= opts.VariableNames;
%remove 'isPrac' if it exists.
 rem = find(contains(oldvarnames, 'isPrac'));
    newvarnames = oldvarnames;
    newvarnames(rem)=[];

    %% ! UNFINISHED! ! ! ! ! ! (below is from gabor, different columns in AUD
for ippant= 1%6:12%:length(pfolsframe)
    % prepare new column headers:
    filename= pfolsframe(ippant).name;
    opts = detectImportOptions(filename,'NumHeaderLines',0);

    oldvarnames= opts.VariableNames;
    
    rem = find(contains(oldvarnames, 'isPrac'));


    fixtable = readtable(filename, opts);
    newvarnames = oldvarnames;
    newvarnames(rem)=[];
    if size(fixtable,2) ~= length(newvarnames)
        fixtable= fixtable(:,1: length(newvarnames));
    end
    fixtable.Properties.VariableNames = newvarnames;
    
    writetable(fixtable,filename);

end
%% subject specific fixes:
if fixIBCP; % this participant had datasets split in two... 
subjID = 'IBCP' % was IGCP
lframe = dir([pwd filesep subjID '*framebyframe.csv']);
lsummary = dir([pwd filesep subjID '*trialsummary.csv']);

t1 = readtable(lsummary(1).name, 'NumHeaderLines',0);
ntrials = max(t1.trial);
nblocks= max(t1.block);
t2 = readtable(lsummary(2).name, 'NumHeaderLines',0);
%adjust trial count and block count.
t2.trial = t2.trial + ntrials + 1; % 0 indexing in Unity.
t2.block = t2.block + nblocks + 1;

% append repaired.
t3= [t1;t2];
%% overwrite.
writetable(t3, lsummary(1).name);

%% now frame data

t1 = readtable(lframe(1).name, 'NumHeaderLines',0);
ntrials = max(t1.trial);
%
t2 = readtable(lframe(2).name, 'NumHeaderLines',0);
%adjust trial count and block count.
t2.trial = t2.trial + ntrials + 1; % 0 indexing in Unity.
% t2.block = t2.block + nblocks + 1;

% % append repaired.
t3= [t1;t2];
%%
% %overwrite.
writetable(t3, lframe(1).name);
end