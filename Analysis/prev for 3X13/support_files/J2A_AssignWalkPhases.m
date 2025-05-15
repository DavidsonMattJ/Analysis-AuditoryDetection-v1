
% RUN HEAD TRACKING ANALYSIS SCRIPT BEFORE THIS SCRIPT

% THIS SCRIPT IS CALLED FROM JOBS_import_preprocess




%% ASSIGN WALK PHASES SCRIPT BEGINS HERE

% This script will take the time stamps of events in the summary files
% and cross reference them with quantile times data to determine which walk
% phase the participant was in when each event occurred.



% updated March 2023 (MD) to assign gait quartiles to visual and auditory
%and September 2023 (GC) to assign to Auditory only
% events

% Open participant Summary
% Find and open the corresesponding participant Quantile data
% 

cd(savedatadir);
allppants = dir([pwd filesep 'p_*.mat']);

%% Load Each Participants Head Tracking Quantile Data

% per participant, load the data, perform some analyses...

%for ippant = 1:length(allppants)  
for ippant = 1 :length(allppants)  % Reinstate ppant loop after testing
    

cd(savedatadir);

%Load Quantile data for participant

% load(allppants(ippant).name, 'quantiledata_table');
% disp(["Loading Quantile Data for Participant " num2str(ippant)]);
% quantiletab = quantiledata_table;

subjectID = allppants(ippant).name;

%Load Summary data for participant
load(allppants(ippant).name, 'summary_table', 'subjID', 'HeadPos');


disp(['Loading Summary Data for Participant ' num2str(ippant)]);
 
nRowsSummary = size(summary_table,1);


%STEPS
%For each line in the summary file (event) take the trial number and event time.
%Find the data in the quantile times for that trial
%Find quantile times in this trial that came before the event qStart < tEvent
%Find maximum value of the remaing times
%Take the Walk phase of this line and assign it to the walk phase for this
%event in the the summary table.

%For each line in the summary file (event) take the trial number and event time.

for iEvent = 1:nRowsSummary
    
%     disp('Starting Loop'+iEvent);
    
    itrial = summary_table.trial(iEvent); % no longer +1 (already performed in import job)
    
    
    %eventTime = summary_table.stimSeqStartTime(iEvent); % first stim, regardless of modality.
    %audTime = summary_table.audCodedTime(iEvent);
    %visTime = summary_table.vistargOnset(iEvent);   Not used in this
    %version

    
   %Used event time as the logic that follows makes most sense that way
    eventTime = summary_table.targOnset(iEvent);
    respTime = summary_table.targRT(iEvent);



    % skip this trial if previously IDd as requiring exclusion:
    skip=0;
    rejTrials_AUD_detect_v1;
    if skip 
        continue % next event if this trial should be ignored.
    end
    
      
%        saveFlds_pcnt= {'trgO_gPcnt','respO_gPcnt'};
       saveFlds_type= {'trgO', 'respO'};

       % some misalignment to be skipped:
       if itrial > length(HeadPos)
           disp(['WARNING! more trials in Summary Data than Head tracking, skipping',....
               'for this participant: ' subjID])
       continue
       end
           
    
    if HeadPos(itrial).blockType==0
        % place nan instead.
        summary_table.([saveFlds_type{1} '_gPcnt'])(iEvent) = nan;
        summary_table.([saveFlds_type{2} '_gPcnt'])(iEvent) = nan;
        summary_table.(['trgO_gStart_samp'])(iEvent) = nan;
        summary_table.(['trgO_gFin_samp'])(iEvent) = nan;                          
        summary_table.(['respO_gStart_samp'])(iEvent) = nan;
        summary_table.(['respO_gFin_samp'])(iEvent) = nan;
        summary_table.(['trgO_gCount'])(iEvent) = nan; % indx of the gait it fell in.
        summary_table.(['respO_gCount'])(iEvent) = nan;
        continue % jump to next iter of for loop
    end
    
  %Find the data in the quantile times for that trial  
%     idxquantileDataThisTrial = find(quantiletab.Trial == trialThisEvent);
    
    trs_samp= HeadPos(itrial).Y_gait_troughs;
    trs_sec= HeadPos(itrial).Y_gait_troughs_sec;
    
    testTimes = trs_sec;
    
    %nearest gaitonset, repeat a process for different events.
    %testEvents= {eventTime, visTime, audTime, lastStim, respTime};
    
    testEvents= {eventTime, respTime};  %eventTime = targOnset



     % step through each of these gait events, (perform the same operation
     % each time), then save in a new column. 
     % We're converting the time(sec) onset of each event into a
     % percentage, relative to occurrence within the co-occuring
     % single-step cycle. (1-100%).
     
     % later scripts will bin into quantiles/ quintiles etc.
      

    
    for igEvent= 1:length(testEvents)
        thisEvent= testEvents{igEvent};
        
    gO = find(testTimes<thisEvent, 1, 'last');
    
        if isempty(gO) || gO ==length(testTimes)      % in this case, we are in the last step of the trial, so skip.
            disp(['Warning! no gait info for iEvent: ' num2str(iEvent) '(miss)' ]);
            
            % place nan for the gait pcnt, and start/end samps
            summary_table.([saveFlds_type{igEvent} '_gPcnt'])(iEvent) = nan;
            summary_table.([saveFlds_type{igEvent} '_gStart_samp'])(iEvent) = nan;
            summary_table.([saveFlds_type{igEvent} '_gFin_samp'])(iEvent) = nan;
            summary_table.([saveFlds_type{igEvent} '_gCount'])(iEvent) = nan;
            
        else
            
            % event as percentage:
            % extract the event as % step.[0 100]
            gaitsamps =trs_samp(gO):trs_samp(gO+1); % avoid counting edge bins twice.
            gaitTimes = HeadPos(itrial).time(gaitsamps);
            resizeT= imresize(gaitTimes', [1,100]);
            %             gPcnt = dsearchn(resizeT', thisEvent);
            
            % or take as proportion of total.
            tDur = gaitTimes(end)- gaitTimes(1);
            tE= thisEvent-gaitTimes(1);
            gPcnt= round((tE/tDur)*100);
            
            % > add to table:
            summary_table.([saveFlds_type{igEvent}  '_gPcnt'])(iEvent) = gPcnt;
            
           %> also add the start and end samples for this gait 
           % (used in next analyses).
           gStart = gaitsamps(1);
           gEnd = gaitsamps(end);
           if gStart ==0 || gEnd==0
               error('check code');
           end
%store
                summary_table.([saveFlds_type{igEvent} '_gStart_samp'])(iEvent) = gStart;
                summary_table.([saveFlds_type{igEvent} '_gFin_samp'])(iEvent) = gEnd;
                summary_table.([saveFlds_type{igEvent} '_gCount'])(iEvent) = gO;
               

        end

    end % each gevent (columns to add)
       
end % each row of summarytable. 

 disp(['Loop for ppant' num2str(ippant ) ', event '  num2str(iEvent)  ' COMPLETE']); 
% End of per event loop


% correct for no respsonses :


%% critical for later stages, remove the data for 'bad' trials.
   allbadtrials = badtrials; % pulled from rejTrials_detectv3.
   allts = summary_table.trial;
   remtrials = ismember(allts, allbadtrials);
   summary_table(remtrials,:)=[];
   
  save(allppants(ippant).name, 'summary_table', '-append');
    
end
%End of Participant loop.



% 