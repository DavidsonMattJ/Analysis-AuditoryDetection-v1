% rejTrials_AUD_detect_v1


% certain bad trials (identified in figure folder).
% look at 'TrialHeadYPeaks', to see if the tracking dropped out on any
% trials. Or if participants were not walking smoothly.

% abbreviations for rejected trials:
% s = poor signal quality (drop-outs/discontinuities)
% g= poor gait extraction (head tracking unclear).

% step in to reject particular ppant+trial combos.
badtrials=[];
switch subjID
    case 'AC'
        badtrials=[200]; %g
    case 'AH'
        badtrials=[]; % none
        
    case 'AW'
        badtrials=[]; % none
    
    case 'CH'
        badtrials=[73];% g
    case 'CY'
        badtrials=200; % g
    
    case 'CZ'       
        badtrials=[5];% g
    case 'EC'
        badtrials=[]; % none (199 trials total)
    case 'GW'
        badtrials=[]; % none (199 trials total)
    case 'HD'
        badtrials=[]; % none
    case 'JK'
    badtrials=[170,200]; % g g (strange double recording in 200)
    
    case 'JL'
        badtrials=[]; % none
        
    case 'JM'
        badtrials=[200]; % g (same strange extra trials in 200)
   
    case 'KH'
        badtrials=[]; % none
    
    case 'KN' % Strong slope (198 trials total )
        badtrials=[]; % none
    
    case 'LN'
        badtrials=[]; % none
        
    case 'LZ'
        badtrials=[55]; % (199 trials total)
        
    case 'MC'
        badtrials=35; %g/s
    case 'MF'
        badtrials=[];  % none
    case 'NH'
        badtrials=[]; % none
    case 'RH'
        badtrials=[126]; %g/s (199 trials total)
        
    case 'RW'
        badtrials=[]; % none;
    case 'SF'
        badtrials=[138]; % s
    
    case 'SK'
    badtrials=[]; %none (198 trials total)
    case 'SM'
        badtrials=[]; % none
    case 'TP'
        badtrials=[]; % none
    case 'YB'
        badtrials=[]; % none
        
        
end

    
%%
if ismember(itrial,badtrials)
    disp(['Skipping bad trial ' num2str(itrial) ' for ' subjID]);
    skip=1;
end
%%