

% CALLED FROM JOBS_import_preprocess.m


%% THIS SECTION UPDATES FILE PATH FOR DESKTOP OR ONLINE VERSION
% 
% MatlabOnline = 0
% 
% if MatlabOnline == 0
%     rootPath = 'C:/Users/mobil/Gabriel/MATLAB_Drive';
% else
%     rootPath = 'MATLAB DRIVE'
% end

%% THIS SECTION BEGINS THE HEAD TRACK ANALYSIS

cd(savedatadir);
allppants = dir([pwd filesep 'p_*.mat']);

%%
% load data, do x y z

for ippant = 1:length(allppants)

cd(savedatadir);


load(allppants(ippant).name, 'framebyframe_table');

filenameThisPpant = allppants(ippant).name; 

disp("Processing Head Tracking Data for Participant " +ippant);
 
mytab = framebyframe_table;

nTotalRows = height(mytab);

ppantName = string(mytab.participant(1));
subjectID = mytab.participant{1};


%%

% Set up an array to hold tables of quatile times.  Added to each Step
quantileTimes = {};

% quantileTimes = array2table(quantileTimes, 'VariableNames', {'Participant', 'Trial', 'Step', 'qstartTime', 'Walk Phase'});

%FOR EACH TRIAL

%for i= 10:15
for i= 1:max(mytab.trial)

    idxsThisTrial = find(mytab.trial == i);

    dataThisTrial = mytab(idxsThisTrial, :);

    tabHeadPosY = dataThisTrial(:,'head_Y');
    rawHeadPosY = table2array(tabHeadPosY);
    tabTime = dataThisTrial(:,'time');
    Time = table2array(tabTime);

    % Find the local maxima
    [pks_raw,locs_raw] = findpeaks(rawHeadPosY);
    pksTimes_raw = Time(locs_raw);

    %Find Gradient of floor (or floor calibration error)
    y2 = pks_raw(end);
    x2 = pksTimes_raw(end);

    y1 = pks_raw(10);    % Use 10th maxima as data at start is messy
    x1 = pksTimes_raw(10);

    floorGradient = (y2-y1)/(x2-x1);

    % Adjust data to correct for floor slope (data plus a line with inverse
    % gradient to the floor slope to balance).
    % (Facilitates filtering nonsensical maxima and minima due to small
    % fluctuations)

    Offset = -floorGradient * Time;
    HeadPosY = rawHeadPosY + Offset;

    % Find the local maxima of adjusted data
    [pks,locs] = findpeaks(HeadPosY, "MinPeakProminence",0.01);
    pksTimes = Time(locs);

    % Find the local minima of adjusted data
    [mins, locs_mins] = findpeaks(-HeadPosY, "MinPeakProminence",0.01);
    minsTimes = Time(locs_mins);
    mins = -mins;

%     % Plot Data - Head Y Position throughout trial
%     figure;
%     title('Trial', i);   % TO DO - ADD PARTICIPANT LABEL ONCE PPANT LOOP ESTABLISHED
% 
%     hold on
%     plot(Time, HeadPosY)
% 
%     % Add markers at the maxima locations
%     scatter(pksTimes,pks,'r','filled');
%     scatter(minsTimes,mins,'g','filled');



    % FOR EACH STEP
    for j=1:(length(locs_mins)-1)

        % Find Time Stamps for start and end of gait cycle (One step,
        % Minima to Minima
        tmin(j) = Time(locs_mins(j));  % Start of step
        tmin(j+1) = Time(locs_mins(j+1));  % End of Step

        % Divide step into 4 equal time divisions
        q1start = tmin(j);

        % q1end = tmin(j) + (tmin(j+1)-tmin(j))/4;

        q2start = tmin(j) + (tmin(j+1)-tmin(j))/4;

        %q2end = tmin(j) + (tmin(j+1)-tmin(j))/2;

        q3start = tmin(j) + (tmin(j+1)-tmin(j))/2;

        %q3end = tmin(j) + (tmin(j+1)-tmin(j))*(3/4);

        q4start = tmin(j) + (tmin(j+1)-tmin(j))*(3/4);

        q4end = tmin(j+1);

        %Add the qStart times to an array so they can be added as rows
        %sequentially

        qStart = nan(4,1);

        qStart(1) = q1start;
        qStart(2) = q2start;
        qStart(3) = q3start;
        qStart(4) = q4start;

        % Assign Walk Phase Values to table based on divisions above

        % For all rows in this step where Trial Time is greater than or
        % equal to q1Start and less that q2Start set Walk Phase = 1

        % For all rows in this step where Trial Time is greater than or
        % equal to q2Start and less that q3Start set Walk Phase = 2

        % For all rows in this step where Trial Time is greater than or
        % equal to q3Start and less that q4Start set Walk Phase = 3

        % For all rows in this step where Trial Time is greater than or
        % equal to q1Start and less that q4end (same as q1Start of next loop) set Walk Phase = 4
        %%

        % Create table of Quantile start times of this Step
        %qTimesThisStep = nan(4,5);
        %qTimesThisStep = array2table(qTimesThisStep, 'VariableNames', {'Participant', 'Trial', 'Step', 'qstartTime', 'Walk Phase'});

        qTimesThisStep = table('Size',[4 5],'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, 'VariableNames', {'Participant', 'Trial', 'Step', 'qstartTime', 'Walk Phase'});

        % qTimesThisStep.Properties.VariableTypes{'Participant'} = 'string';
        qTimesThisStep.Participant(1:end) = ppantName;
        qTimesThisStep.Trial(1:end) = i;
        qTimesThisStep.Step(1:end) = j;


        % FOR EACH QUANTILE OF EACH STEP
        
        for k=1:4      % 4 slices of each step

            qTimesThisStep.qstartTime(k) = qStart(k);
            qTimesThisStep.("Walk Phase")(k) = k;

        end
        
        % End of Quantile For loop

   % Add Table of Quantile Time Data for this step to overall list of
   % Quantile times.

        quantileTimes{end+1} = qTimesThisStep;

%         %% SHADE QUARTILES
% 
%         % Set y boundaries for shading of the quartiles on plot
%         yHigh = max(pks) + (max(pks)-min(mins))*0.1;
%         yLow = min(mins) - (max(pks)-min(mins))*0.1;
% 
%         % Shade all times in Quartile 1 of Walk Phase Red
%         xfill1 = [q1start q2start q2start q1start];
%         yfill1 = [yLow yLow yHigh yHigh];
%         fill(xfill1, yfill1, 'r', 'FaceAlpha', 0.1, "EdgeColor","none") % creates a filled region with red color and 20% transparency
% 
%         % Shade all times in Quartile 2 of Walk Phase Green
%         xfill2 = [q2start q3start q3start q2start];
%         yfill2 = [yLow yLow yHigh yHigh];
%         fill(xfill2, yfill2, 'g', 'FaceAlpha', 0.1, "EdgeColor","none")
% 
%         % Shade all times in Quartile 3 of Walk Phase Blue
%         xfill3 = [q3start q4start q4start q3start];
%         yfill3 = [yLow yLow yHigh yHigh];
%         fill(xfill3, yfill3, 'b', 'FaceAlpha', 0.1, "EdgeColor","none")
% 
%         % Shade all times in Quartile 4 of Walk Phase Yellow
%         xfill4 = [q4start q4end q4end q4start];
%         yfill4 = [yLow yLow yHigh yHigh];
%         fill(xfill4, yfill4, 'y', 'FaceAlpha', 0.1, "EdgeColor","none")


    end
% End of Step For Loop


    % Concatenate the table of this Step onto the overall Quantile Times Table
    %quantileTimes = vertcat(quantileTimes : qTimesThisStep)
    %idxTarg = find()

    
end
% End of Trial For Loop

allQuantileTimes = vertcat(quantileTimes{:});



%SAVING
quantiledata_table = allQuantileTimes;

disp(['saving data for participant ' num2str(ippant)]);

cd(savedatadir);
save(filenameThisPpant, 'quantiledata_table', '-append');


% save(savefilename, ...
%     'Participant', 'Trial', 'Step', 'qstartTime', 'Walk Phase');


%End of Per Participant For Loop
end





