% jobs_concatenate_groupSOA

% -- called from JOBS_import_preprocess;


% Gather output from jobs_calculateParticipantlevelSynchfn
%


%specify save directory:
% savedatadir='/MATLAB Drive/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

allppants = dir([pwd filesep 'p_*data.mat']);
%%

% useppants= [1:16,18:19]; % 17 and 20 have poor gaussian fits.

useppants=1:length(allppants);
nsubs= length(useppants);
GFX_Data_SOAs=[];
for ippant= 1:nsubs
    
    loadp = useppants(ippant);
    cd(savedatadir);
    load(allppants(loadp).name, 'ppantData_SOAs');
    %'propSame_all', 'propSame_slow', 'propSame_fast',...
    %             'RTfromStim1_all', 'RTfromStim1_slow', 'RTfromStim1_fast',...
    %             'RTfromStim2_all', 'RTfromStim2_slow', 'RTfromStim2_fast',...
    %             'nperSOA_all', 'nperSOA_slow', 'nperSOA_fast', 'subjID', 'SOAs', 'propSamebyGaitandSOA_all',...
    %             );
    
    if ippant==1
        allflds = fields(ppantData_SOAs);
        
    end
    
    % GFX_table.propSame_all(ippant,:) = propSame_all;
    disp(['concatenating participant ' num2str(ippant)])
    for ippantfield= 1:length(allflds)
        tmpdata = ppantData_SOAs.([allflds{ippantfield}]);
        if size(tmpdata,1)==1 % array
            GFX_Data_SOAs.([allflds{ippantfield}])(ippant,:) = tmpdata;
        else % matrix:
            GFX_Data_SOAs.([allflds{ippantfield}])(ippant,:,:) = tmpdata;
            
        end
    end % each field
end % ppant
%
save('GFX_SOAdata' , 'GFX_Data_SOAs');
%
%%


%

