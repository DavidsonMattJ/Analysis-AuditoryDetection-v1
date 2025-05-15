% jobs_AssignStridephases

% j2C_append_doubGCData 
% building off the back of j2- > adding an extra column for the
% classification of an event in either a LRL or RLR step sequence. 
% necessary for showing double step proportions of DVs etc.

% append event percentages, relative to doubGait.

%laptop:
cd(savedatadir);
allppants = dir([pwd filesep 'p_*.mat']);
%%


for ippant= 1:length(allppants)
    cd(savedatadir)    %%load data from import job.
    load(allppants(ippant).name, ...
        'HeadPos', 'summary_table', 'subjID');
    savename = allppants(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % per event, find the relevant gait. then add the gPcnt relative
    % to LRL anr RLR sequence (if possible).

    allevents = size(summary_table,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary.
        if strcmp(summary_table.isStationary{ievent}, 'True')
            %>>>>> add this new info to our table,
            %trgOnset:
%           
            summary_table.trgO_gPcnt_LRL(ievent)= nan;
            summary_table.trgO_gPcnt_RLR(ievent)= nan;
            
            
            %response
            summary_table.respO_gPcnt_LRL(ievent)= nan;
                summary_table.respO_gPcnt_RLR(ievent)= nan;


            continue
        else
            % check if this trial is flagged for rejection (visual inspect
            % output of j1).
            skip=0;
            itrial = summary_table.trial(ievent);
            rejTrials_AUD_detect_v1; %toggles skip based on bad trial ID
            if skip ==1
                continue;
            end

            %for each event, conver the target onset time, and reaction onset
            %time, to a % of a double gait cycle.

            tTrial = summary_table.trial(ievent);

            trs = HeadPos(tTrial).Y_gait_troughs;
            trs_sec = HeadPos(tTrial).Y_gait_troughs_sec;
            pks = HeadPos(tTrial).Y_gait_peaks;
            trialTime = HeadPos(tTrial).time;
            eventcolumns = {'targOnset', 'targRT'};
            savecols = {'trgO', 'respO'};

            for ieventtype=1:2
                %%Target onsets first:
                % find the appropriate gait. Final trough that the event is later
                % than:
                thisEvent = summary_table.([ eventcolumns{ieventtype} ])(ievent);
                
               gindx = find(thisEvent>trs_sec,1,'last');


                %don't fill for misses, or before/after first/last step in a trial.
                if thisEvent ==0 || isempty(gindx) || gindx== length(trs)  %
                    summary_table.([savecols{ieventtype}  '_gPcnt_LRL'])(ievent)= nan;
                    summary_table.([savecols{ieventtype}  '_gPcnt_RLR'])(ievent)= nan;
                   

                else % extract info we need:

                    % we have the gIndx, so look at prev and current step.
                    %
                    currentFt = summary_table.([savecols{ieventtype} '_gFoot'])(ievent); % note this was provided in j2B
                    if strcmp(currentFt, 'LR')
                        currDft = 'RLR';
                        prevDft = 'LRL';
                    else
                        currDft = 'LRL';
                        prevDft= 'RLR';
                    end

                    skipCurr=0;
                    skipPrev=0;
                    % extract the event as % double gait.[0 100]
                    
                    % might not have enough steps. 
                    if gindx>=2
                    prevSamps= trs(gindx-1):trs(gindx+1);
                    else 
                        prevSamps=nan;
                        skipPrev=1;
                    end
                    % same for current step
                    if gindx+2 <= length(trs)
                    currSamps= trs(gindx):trs(gindx+2);
                    else
                        currSamps=nan;
                        skipCurr=1;
                    end
%                     
%                     old way:
%                     resizeT= imresize(gaitTimes', [1,100]);
%                     gPcnt = dsearchn(resizeT', thisEvent);
                        useSamps = {currSamps, prevSamps};
                        doubgPcnts=nan(1,2);
                        for iwin=1:2
                            if skipCurr==1 && iwin==1
                                continue % pass over.
                            end
                            if skipPrev==1 && iwin==2
                                continue %pass aswell
                            end
                            gaitSamps= useSamps{iwin};
                            gaitTimes = HeadPos(tTrial).time(gaitSamps);
                        
                            % or take as proportion of total.
                            tDur = gaitTimes(end)- gaitTimes(1);
                            tE= thisEvent-gaitTimes(1);
                            gPcnt= round((tE/tDur)*100);
                            %
                            if gPcnt== 100
                                gPcnt=1;
                            end
                            doubgPcnts(iwin)= gPcnt;

                        end

                    %>>>>> add this new info to our table,
                    summary_table.([savecols{ieventtype}  '_gPcnt_' currDft])(ievent)= doubgPcnts(1);
                    summary_table.([savecols{ieventtype}  '_gPcnt_' prevDft])(ievent)= doubgPcnts(2);
                   


                end % if not misses
            end % targ and response onset

        end % if not practice

    end % each row in table (event)

% debug. Show trgOs
% ptrgOs = trial_summaryTable.trgO_gPcnt;
% ptrgOs = ptrgOs(~isnan(ptrgOs));
% hist(ptrgOs, 100);

%% include zscored version of RTs:
% allRTs = summary_table.clickRT;
%  userow = find(~isnan(allRTs));
% zRTs = zscore(allRTs(userow));
%    %place in table as new columnL
%    trial_summaryTable.z_clickRT = nan(size(trial_summaryTable,1),1);
%    trial_summaryTable.z_clickRT(userow) = zRTs;
% 


   %% critical for later stages, remove the data for 'bad' trials.
   allbadtrials = badtrials; % pulled from rejTrials_detectv3.
   allts = summary_table.trial;
   remtrials = ismember(allts, allbadtrials);
   summary_table(remtrials,:)=[];
%%
%%
   disp(['Finished appending DOUBLE gait percentage data for ... ' subjID]);
   save(savename, 'summary_table','-append');
end % participant

