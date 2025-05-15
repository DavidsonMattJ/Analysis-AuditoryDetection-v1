
% This script loads the raw csvs and saves as a matlab table. 
% -- called from JOBS_import_preprocess;

%% THIS SECTION UPDATES FILE PATH FOR DESKTOP OR ONLINE VERSION

set_myOnlineDirectories_AUD;

cd(rawdatadir)
allppants_frame = dir([pwd filesep '*framebyframe.csv']);
allppants_summary = dir([pwd filesep '*trialsummary.csv']);
    pfols = dir([pwd filesep '*framebyframe.csv']);

% show ppant numbers: in command window
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)

%% load and save frame by frame as table.
%% summary
usefiles = {allppants_frame, allppants_summary};

for ifiletype=2%1:2

    
        readfiles = usefiles{ifiletype};

    for ippant = 1:length(readfiles)
        cd(rawdatadir);
        readfile = readfiles(ippant).name;

        opts=detectImportOptions(readfile);
        %Defines the row location of channel variable name
        opts.VariableNamesLine = 1;
        %Specifies that the data is comma seperated
        opts.Delimiter =','; %Specifies that the data is comma seperated
        %Read the table
        mytab = readtable(readfile,opts, 'ReadVariableNames', true);
        %initial and index.


    if ippant<10
        pnum = ['0' num2str(ippant)];
    else
        pnum = num2str(ippant);
    end
        saveID = ['p_' mytab.participant{1} '_' pnum];

        cd(savedatadir);
    subjID = readfile(1:2);

        disp(['Saving file ' num2str(ifiletype) ' for ppant ' num2str(ippant)]);
        if ifiletype==1
            

            %% Fix column labelling error (issue with recordData script in Unity)

            searchString = 'blockTypetime';

            for col = 1 : width(mytab) % Iterate through each column
                if strcmp(mytab.Properties.VariableNames{col}, searchString)
                    
                    % Reassign Variable Names
                    mytab.Properties.VariableNames{'blockTypetime'} = 'blockType';
                    mytab.Properties.VariableNames{'head_X'} = 'time';
                    mytab.Properties.VariableNames{'head_Y'} = 'head_X';
                    mytab.Properties.VariableNames{'head_Z'} = 'head_Y';
                    mytab.Properties.VariableNames{'Var9'} = 'head_Z';
                    
                    break; % Exit the loop once a match is found
                end
                

            end
%% also repair trial index now (makes everything easier later).
            mytab.trial = mytab.trial+1;

            
            
            
            
            
            
            framebyframe_table= mytab;
            
            disp(['saving framebyframe data for participant ' num2str(ippant)]);
            save([saveID '_data'], 'framebyframe_table', 'subjID');
            
            
            
            
            %% also create HeadPos structure used in later scripts:
            HeadPos=[];
            alltrials = unique(framebyframe_table.trial);
            if min(alltrials)~=1
                error('check code!');
            end
            
            for itrial= 1:length(alltrials)
                                
                rowIndx = find(framebyframe_table.trial==itrial);
                %store:
                HeadPos(itrial).Y =  framebyframe_table.head_Y(rowIndx);
                HeadPos(itrial).X =  framebyframe_table.head_X(rowIndx);
                HeadPos(itrial).Z =  framebyframe_table.head_Z(rowIndx);
                HeadPos(itrial).time =  framebyframe_table.time(rowIndx);
                HeadPos(itrial).blockType =  unique(framebyframe_table.blockType(rowIndx)); % just 1
                
            end
            
            % save:
            save([saveID '_data'], 'HeadPos', '-append');
            
        else
            
            %%  CLEAN ROWS UNINTENTIONALLY COLLECTED AT START OF EACH TRIAL

            % Remove rows where targRT (Target Response time) is less than 500ms (No
            % stimuli are presented in this time)
 
            
            mytab= mytab(mytab.targRT > 0.5, :);

           %% also repair trial index:
            mytab.trial = mytab.trial+1;
            
            
            % what to do with the omitted case i.e. RTs that exceeded are
            % response window (900 ms) are coded as -10 (neither correct or
            % error). 
            % can change to error (0), or omit?
            
            omitRespRows = find(mytab.targResponse==-10);
            
            mytab.correctResponse(omitRespRows)=nan;
            
             % add a column for easy categorisation (H, M, FA, CR).= 1,2,3,4
    %
    presData = mytab.tonePresent;
    sigPresent = find(mytab.tonePresent);
    sigAbsent = find(mytab.tonePresent==0);
    respPres = find(mytab.targResponse ==1);
    respAbs = find(mytab.targResponse ==0);

    Hits = intersect(sigPresent, respPres);
    Misses = intersect(sigPresent, respAbs);
    FalseAlarms = intersect(sigAbsent, respPres);
    CorrectRejections = intersect(sigAbsent, respAbs);

    mytab.SDTcat(Hits) = 1;
    mytab.SDTcat(Misses) = 2;
    mytab.SDTcat(FalseAlarms) = 3;
    mytab.SDTcat(CorrectRejections) = 4;
            % change all RTs to nan 
            summary_table= mytab;
            disp(['saving summary table data for participant ' num2str(ippant)]);
            save([saveID '_data'], 'summary_table', '-append');

        end

        
    end % participants


   

end % filetype

%% 


