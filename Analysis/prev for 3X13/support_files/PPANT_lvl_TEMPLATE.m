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
%     rootPath = 'MATLAB DRIVE'
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
BFPpath = "/MATLAB Drive/Auditory Detection Task 1/Support_Scripts/";
addpath(BFPpath);


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

%for ippant= 1:length(allppants)
for ippant= 1
    
    cd(savedatadir);
    load(allppants(ippant).name, 'summary_table');
    
    disp("Processing Participant Level Data for Participant "+ippant);
    
    
    mytab = summary_table;
    
   
    %% change no responses to nans (to avoid stuffing up averages).
    missedRTs = mytab.targRT==-10;
    missedrts_indx = find(missedRTs);
    mytab.targRT(missedrts_indx)=nan;

    %Filter Data to only walking trials
    walkingIdx = find(mytab.blockType~=0);
    walkingData = mytab(walkingIdx,:);
    
    % create slow data and fast data.
    slowtrials = walkingData.blockType==1;
    fasttrials = walkingData.blockType==2;
    
    slowData = walkingData(slowtrials,:);
    fastData = walkingData(fasttrials,:);
 













    





    %%
    % SOAs = unique(walkingData.SOA);
    % nSOAs = length(SOAs);
    % 
    % %preallocate for all data (walk speeds combined)
    % propSameSOA_all = nan(1,nSOAs);
    % RTfromstim1_perSOA_all= nan(1,nSOAs);
    % RTfromstim2_perSOA = nan(1,nSOAs);
    % nThisSOA = nan(1,nSOAs);
    

%%





    % for our 3 data types (all data, slow ,and fast).
    speedData ={walkingData, slowData, fastData};
    
    save_speed ={'_all', '_slow', '_fast'};
    
    for iCond= 1:length(speedData)
        
        useData = speedData{iCond};
        
        % For this speed condition, create a matrix to hold SOA data for
        % each gait phase
        
        %Preallocating tables for holding Propsame By Gate and SOA
        nanSetup = nan(nGaitPhases,nSOAs);
        
        %         propSameByGaitandSOA = array2table(nanSetup, 'VariableNames', string(SOAs)','RowNames', {'GaitPhase 1', 'GaitPhase 2', 'GaitPhase 3', 'GaitPhase 4'});
        %         propSameByGaitandSOA_visAligned=  propSameByGaitandSOA;
        %         propSameByGaitandSOA_audAligned=  propSameByGaitandSOA;
        %         propSameByGaitandSOA_respAligned= propSameByGaitandSOA;
        
        propSame=[]; % proportion same per SOA
        rt1=[]; % RT from stim 1, per SOA.
        rt2=[]; % RT from stim 2, per SOA.
        nThisSOA=[];
        
        % >> store:
        ppantData_SOAs.SOAs = SOAs;
        
        for iSOA=1:nSOAs
            
            % Number of targets with SOA(i)
            thisSOA = SOAs(iSOA);
            % trial indx.
            trialsWithSOA= find(useData.SOA== thisSOA);
            dataThisSOA = useData(trialsWithSOA, :);
            
            nThisSOA(iSOA) = length(trialsWithSOA);
            
            % Number of SOA(i) trials with 'SAME' response
            SameSOA = find(useData.responseType(trialsWithSOA)==1); % 1 is Same, 2 is Different, 0 is absent
            nSameSOA = length(SameSOA);
            
            % Proportion of Same Responses on those trials
            %             propSame(iSOA) = nSameSOA/nThisSOA(iSOA);
            
            %rt for this SOA.
            %             rt1(iSOA)  = mean(useData.RTfromStim1(trialsWithSOA), 'omitnan');
            %             rt2(iSOA)  = mean(useData.RTfromStim2(trialsWithSOA), 'omitnan');
            
            %>> store, note we include the speed suffix based on condition:
            ppantData_SOAs.(['nTrialsSOA' save_speed{iCond}])(iSOA)          = length(trialsWithSOA);
            ppantData_SOAs.(['rt1' save_speed{iCond}])(iSOA)                 = mean(useData.RTfromStim1(trialsWithSOA), 'omitnan');
            ppantData_SOAs.(['rt2' save_speed{iCond}])(iSOA)                 = mean(useData.RTfromStim2(trialsWithSOA), 'omitnan');
            ppantData_SOAs.(['propSame' save_speed{iCond}])(iSOA)            = nSameSOA/nThisSOA(iSOA);
            
            % adapted to account for alignment with different stim:
            % first case is aligning to the first stim in the AV
            % sequence
            % second case is aligning to visual onset
            % third is aligning to auditory onset
            % fourth is alligning to the response.
            
            % consider aligning to second stim also (will be a mix).
            
            alignAt={'_First', '_Vis', '_Aud','_Last', '_Resp'};
            for ialign= 1:length(alignAt) % each event type above:
                
                usefield = ['StepPcnt' alignAt{ialign}];
                
                for iGaitPhase = 1:nGaitPhases
                    
                    % find trials that have a gpcnt for this quantile
                    % member
                    findpcnts = gQindx{iGaitPhase};
                    
                    %intersection of this SOA, and whichever gait quantile:
                    idxThisSOAthisGaitPhase = find(ismember(dataThisSOA.(usefield), findpcnts));
                    
                    %                     idxThisSOAthisGaitPhase = find(dataThisSOA.(usefield)== iGaitPhase);
                    % using the specific field below:
                    %                     idxThisSOAthisGaitPhase = find(dataThisSOA.(usefield)== iGaitPhase);
                    
                    dataThisSOAthisGaitPhase = dataThisSOA(idxThisSOAthisGaitPhase,:);
                    
                    nThisSOAthisGaitPhase = height(idxThisSOAthisGaitPhase);
                    
                    idxSameThisSOAthisGaitPhase = find(dataThisSOAthisGaitPhase.responseType == 1);
                    
                    nSameThisSOAthisGaitPhase = height(idxSameThisSOAthisGaitPhase);
                    
                    propSameThisSOAthisGaitPhase = nSameThisSOAthisGaitPhase / nThisSOAthisGaitPhase;
                    
                    %propSameByGaitandSOA(iGaitPhase,iSOA) = propSameThisSOAthisGaitPhase;
                    
                    %>> store, note we include the speed suffix based on condition:
                    ppantData_SOAs.(['nThisSOAthisGaitPhase' save_speed{iCond}])(iGaitPhase,iSOA)= nThisSOAthisGaitPhase;
                    ppantData_SOAs.(['propSameByGait' save_speed{iCond} alignAt{ialign} 'Aligned'])(iGaitPhase, iSOA)= propSameThisSOAthisGaitPhase;
                    
                    
                    if any(isnan(propSameThisSOAthisGaitPhase(:)))
                        
%                         pause
                    end
                    
                end %End gate Phase Loop
                
                %while within the SOA and respAlign loop, calc percentwise
                %propSame:
                alignedPropSame=nan(1,100);
                N_alignedPropSame=nan(1,100);
                usefield = ['StepPcnt' alignAt{ialign}];

            % for each pcntge, calculate prop same:
            for ipcnt=1:100
                pcnt_idx= find(ismember(dataThisSOA.(usefield), ipcnt));
                % what prop same this idx?
                resps= dataThisSOA.responseType(pcnt_idx);
                nSame = find(resps==1);
                propSame_idx = length(nSame)/length(resps);
                alignedPropSame(ipcnt)= propSame_idx;
                N_alignedPropSame(ipcnt)= length(resps);

            end


            ppantData_SOAs.(['propSameByGaitPcnt_perSOA_' save_speed{iCond} alignAt{ialign} 'Aligned'])(iSOA,:) = alignedPropSame;
            ppantData_SOAs.(['nStimByGaitPcnt_perSOA_' save_speed{iCond} alignAt{ialign} 'Aligned'])(iSOA,:) = N_alignedPropSame;
       
              
            end % each align type: rename as we go.
        end  % End of SOA loop
        
        %
        % outside of SOA loop, calculate the proportion same over the
        % gait-percentages. per response alignment.
        for ialign= 1:length(alignAt) % each event type above:
            alignedPropSame=nan(1,100);
            N_alignedPropSame=nan(1,100);
            usefield = ['StepPcnt' alignAt{ialign}];

            % for each pcntge, calculate prop same:
            for ipcnt=1:100
                pcnt_idx= find(ismember(useData.(usefield), ipcnt));
                % what prop same this idx?
                resps= useData.responseType(pcnt_idx);
                nSame = find(resps==1);
                propSame_idx = length(nSame)/length(resps);
                alignedPropSame(ipcnt)= propSame_idx;
                N_alignedPropSame(ipcnt)= length(resps);

            end


            ppantData_SOAs.(['propSameByGaitPcnt' save_speed{iCond} alignAt{ialign} 'Aligned']) = alignedPropSame;
            ppantData_SOAs.(['nStimByGaitPcnt' save_speed{iCond} alignAt{ialign} 'Aligned']) = N_alignedPropSame;
        end
    end % End of Speed condition loop - for all data, slow data only, fast data only.
    
    disp([' saving participant  ' num2str(ippant) '(synch functions)'])
    % append data to save file:
    subjID= summary_table.participant{1};
    save(allppants(ippant).name, 'ppantData_SOAs','-append');
    
    
end
% End of Participant loop.
