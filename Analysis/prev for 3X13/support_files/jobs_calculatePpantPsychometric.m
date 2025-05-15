% jobs_calculateParticipantlevelSynchfn
% -- called from JOBS_import_preprocess;

%% THIS SECTION UPDATES FILE PATH FOR DESKTOP OR ONLINE VERSION

% % Set root path to update directory paths
%
% MatlabOnline = 0
%
% if MatlabOnline == 0
%     rootPath = 'C:/Users/mobil/Gabriel/MATLAB_Drive';
% else
    %rootPath = '/MATLAB Drive/Auditory Detection Task 1/FROM MATT/Analysis/support_files';
% end

%% CALCULATE PARTICIPANT LEVEL SOA DATA BEGINS HERE

% this script loads the presaved participant level data (matlab tables).
% calculates proportion 'same' responses per SOA.
% saves at the same location.

% - load data
% - extract trials per SOA type
% - calculate proportion same
% save
% savedatadir= string(rootPath) + '/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

allppants = dir([pwd filesep 'p_*.mat']);



%% Add paths for Support Scripts
% Add path for MLHGaussFit Script, which is called later in this script

%BFPpath = "/MATLAB Drive/Auditory Detection Task 1/Support_Scripts/";


BFPpath = "/MATLAB Drive/Auditory Detection Task 1/FROM MATT/Analysis/support_files/";
addpath(BFPpath);

%% Set Up table to hold ppant data

% colNames = {'ppant', 'SubjID', 'alphaAll', 'betaAll', 'RsqAll', 'alphaSlow', 'betaSlow', 'RsqSlow', 'alphaNorm', 'betaNorm', 'RsqNorm', 'ppantfile'};
% nCols = length(colNames);

% allppantsData = cell2table(cell(0, numel(colNames)), 'VariableNames', colNames);

ppantData= table();
allppantsData=table();
allppantsData_psychPlot= [];

savefile = 'allppantsTable.mat';
savefilePath = fullfile(savedatadir,savefile);

% Create file to wtite all participants data to
save(savefilePath, 'allppantsData');

%% Set up container to hold individual ppant data - written to file above at end of each loop.


% ppantData = NaN(1,nCols);
% 
% ppantData = array2table(ppantData, "VariableNames",colNames);


%% Number of slices of Gait Cycle

nGaitPhases = 4; % MD updated, more flexible this way. (fits fail with nQuants >4)

% have a list of percentages for allocation:
qbounds = [1, ceil(quantile(1:100, nGaitPhases-1)), 100];
gQindx=[]; % stash the percentage index for each quantile (might be different sizes).

for iq= 1:nGaitPhases
    if iq < nGaitPhases
        gQindx{iq} = qbounds(iq):(qbounds(iq+1)-1);
    else
        gQindx{iq} = qbounds(iq):(qbounds(iq+1));
    end
end

%% Begin Participant Loop

%for ippant= 13:17  testing loop contents

for ippant= 1:length(allppants)


    
    cd(savedatadir);
    load(allppants(ippant).name, 'summary_table');
    
    disp("Processing Participant Level Data for Participant "+ippant);

    subjID = summary_table.participant{1};

   
    
    mytab = summary_table;
    
%% Update ppantData

ppantData.ppant = ippant;
ppantData.SubjID = string(subjID);
ppantData.ppantfile = string(allppants(ippant).name);

% also store the fit per ppant.
ppantData_psychPlot=[]; 



    %% change no responses to nans (to avoid stuffing up averages).
    missedRTs = mytab.targRT==-10;
    missedrts_indx = find(missedRTs);
    mytab.targRT(missedrts_indx)=nan;

% first calculate averages (no gait split).

allexpT = find(mytab.isPrac==0);
slowT = find(mytab.walkSpeed==1);
naturalT = find(mytab.walkSpeed==2);

%subselect trials:
allD = mytab(allexpT,:);
slowD = mytab(intersect(allexpT, slowT),:);
natD = mytab(intersect(allexpT, naturalT),:);

% for each, calculate various summary stats:
trialLoop= {allD,slowD,natD};
trialtypePrint= {'All', 'Slow', 'Natural'};

for itrialtype=1:length(trialLoop)
    
    useD = trialLoop{itrialtype};
    % first calc av signal intensity.
    signalInt = useD.toneAmp;
    %remove absent cases (-100)
    signalInt=signalInt(signalInt>-100);
    
    p_signalInt = nanmean(signalInt);
    
    
    % calculate Acc, HR, FAR, d', criterion.
    presTrials = find(useD.tonePresent==1);
    absentTrials = find(useD.tonePresent==0);
    posrespTrials = find(useD.targResponse==1);
    negrespTrials = find(useD.targResponse==0);
    
    % overall acc
    p_Acc = nanmean(useD.correctResponse);
    
    
    %overall RT
    p_RT = nanmean(useD.targRT-useD.targOnset);
    
    %Hit rate
    nHits = length(intersect(presTrials, posrespTrials));
    
    p_HR = nHits/ length(presTrials);
    
    % False alarm rate
    nFA =  length(intersect(absentTrials, posrespTrials));
    
    p_FAR = nFA/ length(absentTrials);
    
    
    % dPrime: 
       
        p_dp = norminv(p_HR)- norminv(p_FAR);
    
    %criterion:
    p_crit = -0.5*(norminv(p_HR)+ norminv(p_FAR));
    
    %store this data:
    ppantData.(['Accuracy_' trialtypePrint{itrialtype}]) = p_Acc;
    ppantData.(['RT_' trialtypePrint{itrialtype}]) = p_RT;
    ppantData.(['HitRate_' trialtypePrint{itrialtype}]) = p_HR;
    ppantData.(['FAR_' trialtypePrint{itrialtype}]) = p_FAR;
    ppantData.(['dprime_' trialtypePrint{itrialtype}]) = p_dp;
    ppantData.(['criterion_' trialtypePrint{itrialtype}]) = p_crit;
    
end



%%  for psychometric fits, we want to restrict to signal present trials only.



    % Filter to only data where tone was present
    
    presentOnly = mytab(mytab.tonePresent==1, :);

    %Remove non-response trials
    presentOnly = presentOnly(presentOnly.targResponse ~= -10, :);

    %Filter Data to only walking trials
    walkingIdx = find(presentOnly.walkSpeed~=0);
    walkingData = presentOnly(walkingIdx,:);
    
    % create slow data and fast data.
    slowtrials = walkingData.walkSpeed==1;
    normtrials = walkingData.walkSpeed==2;
    
    slowDat = walkingData(slowtrials,:);
    normDat = walkingData(normtrials,:);

    allDat = presentOnly;
 

%%

% slowDat = presentOnly(presentOnly.walkSpeed == 1, :);
% 
% normDat = presentOnly(presentOnly.walkSpeed == 2, :);


useDat = { allDat, slowDat, normDat } ;
speedNames = { 'All', 'Slow', 'Normal'};

for ispeed = 1:length(useDat)

%Create 2 column matrix to hold data
dataMat = NaN(height(useDat{ispeed}), 2);

thisSpeedDat = useDat{ispeed};
thisSpeed = string(speedNames{ispeed});

%Fill dataMat, with Intensities in Column 1 and Responses in Column 2
dataMat(:,1) = thisSpeedDat.toneAmp(:);
dataMat(:,2) = thisSpeedDat.targResponse(:);



%% Define Variables

gamma = 0;
figureFlag = 0;
nBins = 9;

%% Call binFitPlot -     function [ alpha , beta, Rsq ] = binFitPlot( dataMat, gamma, figureFlag , nBins);

dataMat  = dataMat(21:end,:);

[ alpha , beta, Rsq , binnedData,plotOUT] = binFitPlot(dataMat, gamma, figureFlag , nBins);

% titleText4 = strcat(subjID, " - ", thisSpeed, " Walk");   
% title(titleText4);

% xlim([-35 -15]);
% 
% cd(figdir);
% try cd('participant_psychoFits');
% catch
%     mkdir('participant_psychoFits')
%     cd('participant_psychoFits')
% end
% 
% print('-dpng', titleText4);

if ispeed == 1
    ppantData.alphaAll = alpha; 
    ppantData.betaAll = beta; 
    ppantData.RsqAll = Rsq; 
elseif ispeed == 2
    ppantData.alphaSlow = alpha; 
    ppantData.betaSlow = beta; 
    ppantData.RsqSlow = Rsq; 
elseif ispeed == 3
    ppantData.alphaNorm = alpha; 
    ppantData.betaNorm = beta; 
    ppantData.RsqNorm = Rsq; 
end

disp(['saving participant' num2str(ippant) "-" thisSpeed]);




ppantData_psychPlot.([speedNames{ispeed} '_xDim']) = plotOUT.xDim;
ppantData_psychPlot.([speedNames{ispeed} '_pFit']) = plotOUT.pfit;
ppantData_psychPlot.([speedNames{ispeed} '_xbins']) = plotOUT.xbins;
ppantData_psychPlot.([speedNames{ispeed} '_ybins']) = plotOUT.ybins;
end %End of Speed Loop

% %% sanity check:
% clf
% cols= {'k','b','r'};
% for ispeed=1:3
%    xDim = ppantData_psychPlot.([speedNames{ispeed} '_xDim']);
%    pFit= ppantData_psychPlot.([speedNames{ispeed} '_pFit']);
%    dX = ppantData_psychPlot.([speedNames{ispeed} '_xbins']);
%    dY = ppantData_psychPlot.([speedNames{ispeed} '_ybins']);
%    
%    hold on;
%    plot(xDim, pFit, '-','color',cols{ispeed});
%    plot(dX,dY,'o', 'color',cols{ispeed});
%    
%     
% end


%%
    
% save at ppant level:
cd(savedatadir)

save(allppants(ippant).name, 'ppantData_psychPlot', 'ppantData','-append');
%


% store ppant data as new row in all

allppantsData(ippant,:) = ppantData;
% bit clunky, but stash each field of ppant structure in group level
% structure:
tfields = fieldnames(ppantData_psychPlot);
for ifield= 1:length(tfields)
allppantsData_psychPlot(ippant).([tfields{ifield}])= ppantData_psychPlot.([tfields{ifield}]);
end

end
% End of Participant loop.
cd(savedatadir);
save('GFX_tables', 'allppantsData', 'allppantsData_psychPlot');


