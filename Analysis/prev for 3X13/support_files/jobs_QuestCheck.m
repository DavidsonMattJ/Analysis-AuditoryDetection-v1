
% cd("/MATLAB Drive/Auditory Detection Task 1/Test_Data_TrialSummary/")

cd("/MATLAB Drive/Auditory Detection Task 1/Participant_Data_TrialSummary")



%readfile = "gcTest_2023-07-14-11-26_trialsummary.csv";

%readfile = "gcTest_2023-07-18-04-18_trialsummary.csv"

% readfile = "gcTest_2023-07-19-06-23_trialsummary.csv"

% readfile = "gcTest_2023-07-20-12-00_trialsummary.csv";

% readfile = "gcTest_2023-08-01-12-05_trialsummary.csv";

%readfile = "gcTest_2023-08-01-12-05_trialsummary.csv";

% readfile = "gcTest_2023-08-02-09-58_trialsummary.csv";

%readfile = "gcTest_2023-08-02-04-19_trialsummary.csv";

%readfile = "gcTest_2023-08-03-12-32_trialsummary.csv";

%readfile = "gcTest_2023-08-03-04-01_trialsummary.csv";

%readfile = "gcTest_2023-08-03-04-13_trialsummary.csv"

%readfile = "gcTest_2023-08-03-05-00_trialsummary.csv";

%readfile = "gcTest_2023-08-03-05-37_framebyframe.csv";

%readfile = "gcTest_2023-08-03-07-30_trialsummary.csv";

%readfile = "gcTest_2023-08-03-07-30_trialsummary.csv";

%readfile = "gcTestCP_2023-08-03-09-21_trialsummary.csv";

%readfile = "JM_2023-08-08-10-00_trialsummary.csv";

%readfile = "CYL_2023-08-09-02-10_trialsummary.csv";

readfile = "SF_2023-08-10-09-22_trialsummary.csv";

%readfile = "SK_2023-08-10-03-29_trialsummary.csv";

%readfile = "RH_2023-08-17-09-14_trialsummary.csv";


open(readfile);

        opts=detectImportOptions(readfile);
        %Defines the row location of channel variable name
        opts.VariableNamesLine = 1;
        %Specifies that the data is comma seperated
        opts.Delimiter =','; %Specifies that the data is comma seperated
        %Read the table
        mytab = readtable(readfile,opts, 'ReadVariableNames', true);
       
%%  CLEAN ROWS UNINTENTIONALLY COLLECTED AT START OF EACH TRIAL

% Remove rows where targRT (Target Response time) is less than 500ms (No
% stimuli are presented in this time)

cleanedData = mytab(mytab.targRT > 0.5, :);



%% Filter to only trials with Tone Present


data = cleanedData(cleanedData.tonePresent == 1, :);
%tonePresIdx = find(mytab.tonePresent == 1);
%data = mytab(tonePresIdx, :);

% Remove non responses

% All indexes with a response
responded = find(data.targResponse~=(-10));
respDat = data(responded, :);

data = respDat;


slowtrials = find(data.walkSpeed==1);
naturaltrials = find(data.walkSpeed==2);

slowDat = data(slowtrials, :);
normDat = data(naturaltrials, :);

%
usetrials ={slowtrials,naturaltrials};
speedDat = {slowDat,normDat};
heading = {'Slow Walk', 'Normal Walk'};
for ispeed=1:2
toneAmps = data.toneAmp(usetrials{ispeed});
subplot(3,2,ispeed)
plot(1:length(toneAmps), toneAmps, 'o-r');
titleText = sprintf('%s %s', 'Threshold Intensity - ', heading{ispeed});
title(titleText);
xlabel('Stimulus Presentation Number');
ylabel('Tone Intensity (dB relative)');


% Accuracy across all trials
overallAcc = height(data.correctResponse(usetrials{ispeed}))/ length((usetrials{ispeed}));



useDat = speedDat{ispeed};
nPres = height(useDat);


if ispeed==1
    overallAccuracySlow = overallAcc;
   
elseif ispeed ==2
    overallAccuracyFast = overallAcc;
end

%runningScore = 0;

runningAccuracy= zeros(1,nPres);
rollingAccuracy= zeros(1,nPres);


colToSum = useDat.correctResponse; % This data has already been restricted to cases where Tone is present so correctResponse == 1 means Correct Detect

for itone=1:nPres 

    %runningScore = runningScore + speedDat.correctResponse(itone);
    accumScore = sum(colToSum(1:itone));

    runningAccuracy(itone) = accumScore/itone;
    
    if (itone>= 11)
        lower = itone-10;
        %upper = itone + 20;
       
        % take average of the last 10 data points
    rollScore = sum(colToSum(lower+1:itone))
    rollingAccuracy(itone) = rollScore/(10);  
    else
       rollingAccuracy(itone)  = accumScore/itone; 
% For first 10 presentations there is not enough data to average over 10
% traisl so use running average for those trials

    end
end
%subplot(2,2,ispeed+2)
subplot(3,2,ispeed+2)
 plot(1:length(runningAccuracy), runningAccuracy, 'o-b');
 titleText2 = sprintf('%s %s', 'Running Accuracy - ', heading{ispeed})

 title(titleText2);
ylim([.5 1])
ylabel('Proportion of Tones Detected')
xlabel('Presentation Number')


subplot(3,2,ispeed+4)
plot(1:length(rollingAccuracy), rollingAccuracy, 'o-b');
titleText3 = sprintf('%s %s', 'Rolling Accuracy - ', heading{ispeed})

title(titleText3);
ylim([.5 1])
ylabel('Proportion of Tones Detected')
xlabel('Presentation Number')
end

shg


%% Check Accuracy

% Slow Data

% Tone Present, correct responses


% Tone Absent, correct responses



%Indexes where tonePresent matches targResponse
congruentIdx = find(respDat.targResponse == respDat.tonePresent);
numberCongruent = height(congruentIdx);
congruentDat = respDat(congruentIdx, :);

%incongruent = find(respDat.targResponse ~= respDat.tonePresent);

%Indexes with correctResponse = 1

correctIdx = find(respDat.correctResponse==1);
numberCorrect = height(correctIdx)
correctDat = respDat(correctIdx, :);

intersect(congruentDat,correctDat);

inconsistency = setdiff(congruentDat,correctDat);

%%  

% qstep0dat =  [];
% qstep1dat = [];
% qstep2dat =  [];
% qstep3dat =  [];
% qstep4dat =  [];
% qstep5dat =  [];
% qstep6dat =  [];



%Slow Walk Speed
%slowDat = respDat(respDat.walkSpeed == 1, :);
%Average Simulus intensity


% % Normal Walk Speed
% normDat = respDat(respDat.walkSpeed == 2, :);
% 
% useDat = {slowDat,normDat}
% speedStep = {slowStep(), normStep()};
% 
% for iSpeed = 1:2
% 
%     for iStep = 0:6
%     speedStep[i] = normDat({iSpeed}.qStep == iStep, :);
%     end
% 
% 
% end
% 
% 
% 

% Find rows where tonPresent = 1 but toneAmp = -100 (shouldn't be possible)
anomaly = mytab(mytab.tonePresent == 1, :);
%update to only toneAmp = -100
anomaly = anomaly(anomaly.toneAmp == -100, :);

nAnom = height(anomaly)
