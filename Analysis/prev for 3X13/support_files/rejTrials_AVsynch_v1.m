% rejTrials_AVsynch_v1


% certain bad trials (identified in figure folder).
% look at 'TrialHeadYPeaks', to see if the tracking dropped out on any
% trials. Or if participants were not walking smoothly.

% abbreviations for rejected trials:
% s = poor signal quality (drop-outs/discontinuities)
% g= poor gait extraction (head tracking unclear).

% step in to reject particular ppant+trial combos.
badtrials=[];
switch subjectID
    case 'p_AL_01_data.mat'
       badtrials=[];
    case 'p_AN_02_data.mat'
        badtrials=[];
    case 'p_CK_03_data.mat'
        badtrials=[];
        
    case 'p_CW_04_data.mat'
        badtrials=[];
    case 'p_DL_05_data.mat' % query REJ. Large slope artefact.
        badtrials=9; 
    case 'p_EO_06_data.mat'
        badtrials=[];
    case 'p_FS_07_data.mat'
        badtrials=[6,45,189,196];
        
    case 'p_GHN_08_data.mat'
        badtrials=174;
    case 'p_JKE_09_data.mat'       
        badtrials=[50,53,128, 169];
    case 'p_JRL_10_data.mat'
        badtrials=[];
    case 'p_JT_11_data.mat'
        badtrials=[5,6,23:30, 43,136,139,152, 173];
    case 'p_LDM_12_data.mat'
        badtrials=[22,36,53113,135];
    case 'p_LS_13_data.mat'
        badtrials=11;
    case 'p_MMM_14_data.mat'
        badtrials=[6,111, 177];
    case 'p_MTT_15_data.mat'
        badtrials=[];
    case 'p_RJ_16_data.mat'
        
       badtrials=[];
    case 'p_SB_17_data.mat' % query reject. Very large slope. Less data (n=148) 
        badtrials=[];
    case 'p_SK_18_data.mat'
        badtrials=[9, 84];
    case 'p_TK_19_data.mat' % query REJ, very large slope (throughout).
        badtrials=[6];
        
        % (HAVE REJECTED): low trials, nans in gait SOA.    
    case 'p_WW_20_data.mat' % query REJ, less data (n=145).   
        badtrials=[5,97, 141, 143:145];
    
    case 'p_XH_21_data.mat'
        badrials=[48];

    case 'p_AC_01_data.mat'
        badtrials=[54,66];

    case 'p_AK_02_data.mat'
        badtrials=[]; % no bad trials.
    case 'p_DG_03_data.mat'
        badtrials=[]; % none
        
        
end

    
%%
if ismember(itrial,badtrials)
    disp(['Skipping bad trial ' num2str(itrial) ' for ' subjID]);
    skip=1;
end
%%