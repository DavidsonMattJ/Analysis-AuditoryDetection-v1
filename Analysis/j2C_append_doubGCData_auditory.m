% j2C_append_doubGCData_gabor 
% building off the back of j2- > adding an extra column for the
% classification of an event in either a LRL or RLR step sequence. 
% necessary for showing double step proportions of DVs etc.

% append event percentages, relative to doubGait.

%laptop:
set_myOnlineDirectories_AUD;
cd(procdatadir)
% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%


for ippant= 1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % per event, find the relevant gait. then add the gPcnt relative
    % to LRL anr RLR sequence (if possible).

    allevents = size(trial_summaryTable,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary.
        if trial_summaryTable.isPrac(ievent) || trial_summaryTable.isStationary(ievent)
            %>>>>> add this new info to our table,
            %trgOnset:
%           
            trial_summaryTable.trgO_gPcnt_LRL(ievent)= nan;
            trial_summaryTable.trgO_gPcnt_RLR(ievent)= nan;
            
            
            %response
            trial_summaryTable.respO_gPcnt_LRL(ievent)= nan;
                trial_summaryTable.respO_gPcnt_RLR(ievent)= nan;


            continue
        else
            % check if this trial is flagged for rejection (visual inspect
            % output of j1).
            skip=0;
            itrial = trial_summaryTable.trial(ievent);
            rejTrials_detectvAud    ;
            if skip ==1
                continue;
            end

            %for each event, conver the target onset time, and reaction onset
            %time, to a % of a double gait cycle.

            tTrial = trial_summaryTable.trial(ievent);

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
                thisEvent = trial_summaryTable.([ eventcolumns{ieventtype} ])(ievent);
                
               gindx = find(thisEvent>trs_sec,1,'last');


                %don't fill for misses, or before/after first/last step in a trial.
                if thisEvent ==0 || isempty(gindx) || gindx== length(trs)  %
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt_LRL'])(ievent)= nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt_RLR'])(ievent)= nan;
                   

                else % extract info we need:

                    % we have the gIndx, so look at prev and current step.
                    %
                    currentFt = trial_summaryTable.([savecols{ieventtype} '_gFoot'])(ievent); % note this was provided in j2B
                    if strcmp(currentFt, 'LR');
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

                            % note that the result can now be 0 ! which 
                            % messes with later calculations. 
                            if ismember(gPcnt, [0,1,100])
                                % flip coin and give to either 1 or 100.
                                if randi([0 1])
                                    gPcnt=1;
                                else
                                    gPcnt=100;
                                end
                            end

                            
                            doubgPcnts(iwin)= gPcnt;

                        end

                    %>>>>> add this new info to our table,
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt_' currDft])(ievent)= doubgPcnts(1);
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt_' prevDft])(ievent)= doubgPcnts(2);
                   


                end % if not misses
            end % targ and response onset

        end % if not practice


    %safety check - should only have a LRL classification if the there is a
    %corresponding LR or RL to enable comparison later.
    if isnan(trial_summaryTable.trgO_gPcnt(ievent)) &&  ~isnan(trial_summaryTable.trgO_gPcnt_LRL(ievent));

        disp('check for errros');
    end
    end % each row in table (event)

% debug. Show trgOs
% ptrgOs = trial_summaryTable.trgO_gPcnt;
% ptrgOs = ptrgOs(~isnan(ptrgOs));
% hist(ptrgOs, 100);

%% include zscored version of RTs:
allRTs = trial_summaryTable.clickRT;
 userow = find(~isnan(allRTs));
zRTs = zscore(allRTs(userow));
   %place in table as new columnL
   trial_summaryTable.z_clickRT = nan(size(trial_summaryTable,1),1);
   trial_summaryTable.z_clickRT(userow) = zRTs;



   %% critical for later stages, remove the data for 'bad' trials.
   allbadtrials = badtrials; % pulled from rejTrials_detectv3.
   allts = trial_summaryTable.trial;
   remtrials = ismember(allts, allbadtrials);
   trial_summaryTable(remtrials,:)=[];

   %% quick sanity check:
% 
% allgO = trial_summaryTable.trgO_gPcnt;
% allgO_doubLRL = trial_summaryTable.trgO_gPcnt_LRL;
% allgO_doubRLR = trial_summaryTable.trgO_gPcnt_RLR;
% 
% LRtrials = find(strcmp(trial_summaryTable.trgO_gFoot, 'LR'));
% RLtrials = find(strcmp(trial_summaryTable.trgO_gFoot, 'RL'));
% bothtrials = [LRtrials;RLtrials];
% 
% clf;
% subplot(241);
% histogram(allgO(LRtrials),200); title(['LR n=' num2str(length(LRtrials))]);
% subplot(242);
% histogram(allgO(RLtrials),200); title(['LR n=' num2str(length(RLtrials))]);
% subplot(243); % combined.
% histogram(allgO(bothtrials),200); title(['LR+RL n=' num2str(length(bothtrials))]);
% subplot(2,4,4)
% %all:
% histogram(allgO,200); title(['column n=' num2str(length(allgO))]);
% % % now see if the stride makes sense:
% subplot(2,4,5);
% histogram(allgO_doubLRL,200);title(['LRL n=' num2str(length(find(~isnan(allgO_doubLRL))))]);
% subplot(2,4,6);
% histogram(allgO_doubRLR,200);title(['RRL n=' num2str(length(find(~isnan(allgO_doubRLR))))]);
% 
% shg
%% display counts atedges (for debugging).
disp(['counts for LRL g1 = ' num2str(length(find(trial_summaryTable.trgO_gPcnt_LRL==1)))]);
disp(['counts for LRL g100 = ' num2str(length(find(trial_summaryTable.trgO_gPcnt_LRL==100)))]);
disp(['counts for RLR g1 = ' num2str(length(find(trial_summaryTable.trgO_gPcnt_RLR==1)))]);
disp(['counts for RLR g100 = ' num2str(length(find(trial_summaryTable.trgO_gPcnt_RLR==100)))]);
%%
   disp(['Finished appending DOUBLE gait percentage data for ... ' subjID]);
   save(savename, 'trial_summaryTable','-append');
end % participant

