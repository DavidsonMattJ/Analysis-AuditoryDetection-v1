% jobs_Assign_stepSupport


% Here we will continue the previous job which allocated gait percentiles
% to target and response events.
% Here we add details, such as estimating which foot the person was
% standing on per step.

% -  appends whether L/R ft 


%%
cd(savedatadir)    %%load data from import job.
pfols = dir([pwd filesep 'p_*']);
%%
nsubs= length(pfols);
for ippant =1:nsubs
    cd(savedatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'summary_table', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.


    allevents = size(summary_table,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary.
        if strcmp(summary_table.isStationary{ievent}, 'True');
            
            %>>>>> add this new info to our table,            
            summary_table.trgO_gFoot(ievent) = {'nan'};
             summary_table.respO_gFoot(ievent) = {'nan'};
            
            continue
        else
            % check if this trial is flagged for rejection (visual inspect
            % output of j1).

            itrial = summary_table.trial(ievent);

            skip=0;
            rejTrials_AUD_detect_v1;
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

            gStart = summary_table.([savecols{ieventtype} '_gStart_samp'])(ievent); 
            gEnd = summary_table.([savecols{ieventtype} '_gFin_samp'])(ievent);
           
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


                summary_table.([savecols{ieventtype} '_gFoot'])(ievent) = {ft};
            end
            end % eventtype (targetonset or resp in gait).

        end % if not practice

    end % each row in table (event)

   disp(['Finished appending gait percentage data for ... ' subjID]);
   save(savename, 'summary_table','-append');
end % participant

