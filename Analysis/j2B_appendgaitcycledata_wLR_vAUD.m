% j2B_appendgaitcycleData_LRft_vAUD

% Here we will load the summary data table, and add event markers  new columns, for when events
% happened relative to % gait cycle completion.
% For now, events of interest are target onset, and response onset. 

% jobs:
% -  appends whether L/R ft (now in j2B).


% % % % Auditory task VERSION % % % 

%laptop:
set_myOnlineDirectories_AUD;     
cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%


for ippant =1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos',  'trial_summaryTable', 'trialInfo','subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.


    allevents = size(trial_summaryTable,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary.
        if trial_summaryTable.isPrac(ievent) || trial_summaryTable.isStationary(ievent)
            
            %>>>>> add this new info to our table,            
            trial_summaryTable.trgO_gFoot(ievent) = {'nan'};
             trial_summaryTable.respO_gFoot(ievent) = {'nan'};
            
            continue
        else
            % check if this trial is flagged for rejection (visual inspect
            % output of j1).

            itrial = trial_summaryTable.trial(ievent);

            skip=0;
            rejTrials_detectvAud; %toggles skip based on bad trial ID
            if skip ==1
                continue;
            end

            %for each event, determine whether the gait cycle is Left-Right
            %ft, or Right-Left ft. This is done by observing the sway (left
            % to right), in head position.
            
             

            % Head position data:
            tmpPos=  squeeze(HeadPos(itrial).Y);
            tmpSway = squeeze(HeadPos(itrial).Z);
            tmpwalkDir = squeeze(HeadPos(itrial).X);
            
            trialTime = HeadPos(itrial).time;
        
            % quick classification:
            if mod(itrial,2)~=0 % odd numbers (walking toward computer)
                % Then descending on X, (i.e. first trajectory), more positive z values
                % are left side of the body.
                walkDir= 'IN';
                ylab='Right - > Left';
            elseif mod(itrial,2) ==0
                % ascending (the return trajectory), pos z values are RHS
                walkDir= 'OUT';

                ylab='Left - > Right';
            end
            

            % note that for each event (row int the table), we could have
            % separate gaits for target onset and response. 
            % So calculate both
            savecols = {'trgO', 'respO'};
            for ieventtype=1:2

            gStart = trial_summaryTable.([savecols{ieventtype} '_gStart_samp'])(ievent); 
            gEnd = trial_summaryTable.([savecols{ieventtype} '_gFin_samp'])(ievent);
           
            % if no response though, continue
            if any(isnan([gStart gEnd]))
                continue
            else
                % is the z value increasing or decreasing, relative to feet placement?

                %
                midlineS = mean(tmpSway([gStart, gEnd]));
                meanS = mean(tmpSway(gStart:gEnd));

                gaitSway = tmpSway(gStart:gEnd);

                if strcmp(walkDir, 'IN') && meanS>midlineS
                    % if we were approaching computer, and swinging positive
                    ft='LR'; % foot placement was right-left.
                elseif strcmp(walkDir, 'IN') && meanS<midlineS
                    %swinging opp direction.
                    ft= 'RL';
                elseif strcmp(walkDir, 'OUT') && meanS>midlineS
                    %if back to computer, positive swing is LR
                    ft= 'RL';
                elseif strcmp(walkDir, 'OUT') && meanS<midlineS
                    ft= 'LR';
                end


                %% sanity check:
%                 clf; subplot(211);
% 
%                 plot(trialTime, tmpSway);
%                 subplot(212);
%                 plot(trialTime, tmpSway, 'color', [.8 .8 .8]);
%                 % to ease interp, reorient to allocentric:
%                 if strcmp(walkDir, 'OUT')
%                     set(gca,'Ydir','reverse')
%                      ylabel('Right -Left')
%                 else
%                     ylabel(ylab)
%                 end
% 
%                 %overlay
%                 hold on;
%                 plot(trialTime(gStart:gEnd), gaitSway, 'k', 'linew', 2)
%                 hold on;
%                 plot(trialTime([gStart, gEnd]), tmpSway([gStart, gEnd]), 'ro')
%                 %add text to be sure
%                 title(['Walking ' walkDir ', gait: ' ft])
               

                %%
                %>>>>> add this new info to our table,


                trial_summaryTable.([savecols{ieventtype} '_gFoot'])(ievent) = {ft};
            end
            end % eventtype (targetonset or resp in gait).

        end % if not practice

    end % each row in table (event)

   disp(['Finished appending gait percentage data for ... ' subjID]);
   save(savename, 'trial_summaryTable','-append');
end % participant

