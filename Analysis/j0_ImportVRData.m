% Discrimination experiment Auditory -
%%  Import from csv. FramebyFrame, then summary data.

%%%%%% v4: AUDITORY data (Discrimination version) %%%%%%

set_myOnlineDirectories_AUD;

% note now updated for the split dataset (with and without eye).
sancheckplots=0; % toggle to show FA per trial etc.
cd(datadir)
    pfols = dir([pwd filesep '*framebyframe.csv']);

% show ppant numbers: in command window
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)


%revisit (mismatch) 

%% Per csv file, import and wrangle into Matlab Structures, and data matrices:
for ippant =1:length(pfols)
    cd(datadir)

    pfols = dir([pwd filesep '*framebyframe.csv']);

    %% load subject data as table.
    filename = pfols(ippant).name;
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(1)-1);
%%

    

    if ippant<10
        pnum = ['0' num2str(ippant)];
    else
        pnum = num2str(ippant);
    end
        savename= ['p_'  subjID  '_' pnum '_summary_data'];
    
        
    %read table
    cd(datadir);
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    opts = setvaropts(opts,"date", 'InputFormat', 'dd/Mm/uuuu');
    disp(['reading large frame x frame file now...']);
    T = readtable(filename,opts);
    ppant = T.participant{1};
    disp(['Preparing participant ' ppant ' ' pnum]);

%% first repair the column headers:
   %% Fix column labelling error (issue with recordData script in Unity)

            searchString = 'blockTypetime';

            for col = 1 : width(T) % Iterate through each column
                if strcmp(T.Properties.VariableNames{col}, searchString)
                    
                    % Reassign Variable Names
                    T.Properties.VariableNames{'blockTypetime'} = 'blockType';
                    T.Properties.VariableNames{'head_X'} = 'time';
                    T.Properties.VariableNames{'head_Y'} = 'head_X';
                    T.Properties.VariableNames{'head_Z'} = 'head_Y';
                    T.Properties.VariableNames{'Var9'} = 'head_Z';
                    
                    break; % Exit the loop once a match is found
                end
                

            end

    
  %% also create HeadPos structure used in later scripts:
            HeadPos=[];
            alltrials = unique(T.trial);
            
            
            for itrial= 1:length(alltrials)
                                
                rowIndx = find(T.trial==itrial);
                %store:
                HeadPos(itrial).Y =  T.head_Y(rowIndx);
                HeadPos(itrial).X =  T.head_X(rowIndx);
                HeadPos(itrial).Z =  T.head_Z(rowIndx);
                HeadPos(itrial).time =  T.time(rowIndx);
                HeadPos(itrial).blockType =  unique(T.blockType(rowIndx)); % just 1
                
            end
            
            
    %% use logical indexing to find all relevant info (in cells)
    
    cd(procdatadir)


    try save(savename, 'HeadPos',  'subjID', 'ppant', '-append');
    catch
        save(savename, 'HeadPos' ,'subjID', 'ppant');
    end

    %     end
    %
    %% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ now summary data



    cd(datadir)
    pfols= dir([pwd filesep '*trialsummary.csv']);
    nsubs= length(pfols);
    %
    filename = pfols(ippant).name;

    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);
    %read table
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    opts = setvaropts(opts,"date", 'InputFormat', 'dd/Mm/uuuu');
    T = readtable(filename,opts);
    rawSummary_table = T;
    
    disp(['Preparing participant ' T.participant{1}  ' ' pnum]);

%%
    %% FAlse alarms: repair FA data in table (happens if we have multiple FAs, extra columns are added, which need to be collapsed.
    ColIndex = find(strcmp(T.Properties.VariableNames, 'FA_rt'), 1);
    %repair string columns:
    for ix=ColIndex:size(T,2)
        tmp = table2array(T(:, ix));
        if iscell(tmp)
            % where is there non-zero data, to convert?
            replaceRows= find(~cellfun(@isempty,tmp));
            %convert each row
            newD=zeros(size(tmp,1),1);
            for ir=1:length(replaceRows)
                user = replaceRows(ir);
                newD(user) = str2num(cell2mat(tmp(user )));
            end

            %% change table column type:
            colnam = T.Properties.VariableNames(ix);
            T.(colnam{:}) = newD;
            %             T(:,ix) =table(newD');

        end


    end
    %% extract Target onsets per trial (as table).
    %% and Targ RTs, contrast values if they exist.
    frametrials = length(alltrials); 
    alltrials = unique(T.trial); 
    disp([subjID  'has ' num2str(length(alltrials)) '(summary)'])


    if frametrials~=length(alltrials)
        error([ subjID ': trial mismatch!'])
    end

    %we want to repair the table as we go.
    trial_summaryTable= T;
    %remove some cols we don't need:
    trial_summaryTable.date=[];

    try trial_summaryTable.qStep=[]; % not all have it for some reason.
    catch
    end
    %rename correctResponse to targCor, and tonePresent to signalPresent  (match this/later scripts).
    colid1= find(contains(trial_summaryTable.Properties.VariableNames, 'correctResponse'));
    trial_summaryTable.Properties.VariableNames{colid1}= 'targCor';
    colid2= find(contains(trial_summaryTable.Properties.VariableNames, 'tonePresent'));
    trial_summaryTable.Properties.VariableNames{colid2}= 'signalPresent';
    %%

    %%There's a bug.
    % Some trials start with a RT on frame 1 (first 11ms). Remove theses
    % data rows.
    a= find(trial_summaryTable.targRT < .5);
    b= find(trial_summaryTable.targRT>0); % no responses are -10. (so this finds those RTs betwen 0 and .02
    remrows = intersect(a,b);
    trial_summaryTable(remrows,:)=[];
    


    for itrial= 1:length(alltrials)
        thistrial= alltrials(itrial);
        relvrows = find(trial_summaryTable.trial ==thistrial); % unity idx at zero.

        %create vectors for storage:
        tOnsets = trial_summaryTable.targOnset(relvrows);
        tRTs = trial_summaryTable.targRT(relvrows);
        tCor = trial_summaryTable.targCor(relvrows);
        tFAs= table2array(trial_summaryTable(relvrows(end), ColIndex:end));
        tFAs= tFAs(~isnan(tFAs));



        %% Reaction time region (tidy / reclassify some values).

        RTs = (tRTs - tOnsets);
        %         %note that negative RTs, indicate that no response was recorded:
        tOmit = find(RTs<=0);
        if ~isempty(tOmit)
            %             tCor(tOmit) = NaN; % don't count those incorrects, as a mis identification.
            RTs(tOmit)=NaN; % remove no respnse
        end

        % now we can add this RT information to our revised table
        trial_summaryTable.clickRT(relvrows) = RTs;


        %% we also want to reject very short RTs (reclassify as a FA).

        if any(find(RTs<0.15))
            disp(['Trial: ' num2str(thistrial) ' Suspicious RT']);
            checkRT= find(RTs<.15);
            %debug to check:

            % reclassify data (mark as incorrect, and record as FA).
            for ifalseRT = 1:length(checkRT)
                indx = relvrows(checkRT(ifalseRT));
                %repair table to avoid counting twice! 

                trial_summaryTable.FA_rt(indx) =  trial_summaryTable.targRT(indx);
                trial_summaryTable.clickRT(indx) = NaN;
                trial_summaryTable.targOnset(indx)=NaN;
                trial_summaryTable.signalPresent(indx)= NaN;
                trial_summaryTable(indx, 11:16)= table(NaN);
                

                %debug add to fig:
                if sancheckplots
                    clf; title(['Trial: ' num2str(thistrial) 'early RT?']); hold on;
                    plot(trialInfo(itrial).times, trialInfo(itrial).targstate); hold on;
                    %
                    plot([tRTs(checkRT(ifalseRT)), tRTs(checkRT(ifalseRT))], [0 1], 'om')


                    % add each 
                    plotRTs= tRTs(tRTs>0);
                    for ip=1:length(plotRTs)
                        plot([plotRTs(ip) plotRTs(ip)], [0 .5], 'k', 'linew',2);
                    end
                    shg;
%             pause; % wait for click
                end
            end
%             
        % update for next catch:
        tOnsets = trial_summaryTable.targOnset(relvrows);
        tRTs = trial_summaryTable.targRT(relvrows);

        end
        %%
        % seems some FA are missing, do a quick check to reclassify double
        % responses and determine FAs.
        % frequently multiple triggers are pulled (could be self
        % corrections?).
        % find the *second* trigger within a window, and convert to a FA.
        % these are identfiable as duplicate entries to the same target onset.
        if any(diff(tOnsets)==0)

           
            if sancheckplots
                %plot
                clf;
                %head data
                plot(trialInfo(itrial).times, HeadPos(itrial).Y); hold on;;

                yyaxis right
                %trg and rts:
                for itrgo = 1:length(tOnsets)
                    plot([tOnsets(itrgo) tOnsets(itrgo)], [0, 1], 'k-o')
                    text(tOnsets(itrgo), 1, num2str(itrgo))                  
                end
                for irt = 1:length(tRTs)
                    plot([tRTs(irt) tRTs(irt)], [0, 1], 'r-o')
                    text(tRTs(irt), 1, num2str(irt))
                end

                ylim([0 2]);
                shg
            end




            % reclassify data (mark as incorrect, and record as FA).
            checkRT= find(diff(tOnsets)==0);
            %debug to check:

            % reclassify data (mark as incorrect, and record as FA).
            for ifalseRT = 1:length(checkRT)
                indx = relvrows(checkRT(ifalseRT));
                %repair table to avoid counting twice!
                % note we are assuming the second response is the error:
                repindx= indx+1;

                trial_summaryTable.FA_rt(repindx) =  trial_summaryTable.targRT(repindx);
                trial_summaryTable.clickRT(repindx) = NaN;
                trial_summaryTable.targOnset(repindx)=NaN;
                trial_summaryTable.signalPresent(repindx)= NaN;
                trial_summaryTable(repindx, 11:16)= table(NaN);
            end
  % update for next catch:
        tOnsets = trial_summaryTable.targOnset(relvrows);
        tRTs = trial_summaryTable.targRT(relvrows);
        end  %any duplicates

        % end RT region


    end
    %repair trial indexing:
    trial_summaryTable.trial =  trial_summaryTable.trial +1;
    trial_summaryTable.block =  trial_summaryTable.block +1;
    trial_summaryTable.trialID =  trial_summaryTable.trialID +1;

    % remove '0' time stamp for any empty FAs cells.
    tFAs =find(trial_summaryTable.FA_rt==0);
    trial_summaryTable.FA_rt(tFAs)=nan;

    %also convert cells to easier format (only some ppants).
    isstatall= trial_summaryTable.isStationary;

    if iscell(isstatall)
        isstatint = contains(isstatall, 'TRUE', 'IgnoreCase',true);
        trial_summaryTable.isStationary = isstatint;
    end

% change all -10 to nan (these were no response recorded.
chRT= find(trial_summaryTable.targRT==-10);
trial_summaryTable.targRT(chRT)=nan;

    % add a column for easy categorisation (H, M, FA, CR).= 1,2,3,4
    %
    presData = trial_summaryTable.signalPresent;
    sigPresent = find(trial_summaryTable.signalPresent);
    sigAbsent = find(trial_summaryTable.signalPresent==0);
    respPres = find(trial_summaryTable.targResponse ==1);
    respAbs = find(trial_summaryTable.targResponse ==0);

    Hits = intersect(sigPresent, respPres);
    Misses = intersect(sigPresent, respAbs);
    FalseAlarms = intersect(sigAbsent, respPres);
    CorrectRejections = intersect(sigAbsent, respAbs);

    trial_summaryTable.SDTcat(Hits) = 1;
    trial_summaryTable.SDTcat(Misses) = 2;
    trial_summaryTable.SDTcat(FalseAlarms) = 3;
    trial_summaryTable.SDTcat(CorrectRejections) = 4;

    %save this structure for later analysis per gait-cycle:
    disp(['Saving trial summary data ... ' subjID]);
    cd(procdatadir)

    save(savename, 'trial_summaryTable','-append');


end % participant
% end % with out without.
