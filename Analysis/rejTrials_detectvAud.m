% rejTrials_detectvAud% rejTrials_detectvGabor%

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
        badtrials = []; % none.
    case 'AD'
        badtrials = []; % none.
    case 'AGH'
        badtrials = []; % none.
    case 'AH'
        badtrials = []; % none.
    case 'AW'
        badtrials = []; % none.
    case 'CHL'
        badtrials= [72];
    case 'CYL'
        badtrials=[43];
    case 'CZ'
        badtrials=[];
    case 'EC'
        badtrials=[];

    case 'EM'
        badtrials=[];
    case 'EWC'
        badtrials=[];
    case 'GW'
        badtrials=[];

    case 'HB'
        badtrials=[15];

    case 'HD'
        badtrials=[];

    case 'JH'
        badtrials=[60, 158,160,163,164,166];

    case 'JK'
        badtrials=[169];
    case 'JL'
        badtrials=[];
    case 'JM'
        badtrials=[];

    case 'JO'
        badtrials=[];
    case 'JT'
        badtrials=[60];

    case 'KHS'
        badtrials=[];
    case 'KN'
        badtrials=[108,143];

        case 'LN'
        badtrials=[];
        
    case 'LZ'
        badtrials=[54];

        case 'MC'
        badtrials=[34];

        case 'MF'
        badtrials=[125];
        
    case 'MS'
        badtrials=[];
    case 'NHGT'
        badtrials=[];
        
    case 'PC'
        badtrials=[191];
    
    case 'RH'
        badtrials=[125,126];
    case 'RW'
        badtrials=[];

    case 'SC'
        badtrials=[]; % good examples of asymmetric steps/head Y
    case 'SF'
        badtrials=[];
    
    case 'SK'
        badtrials=[98,100];

    case 'SM'
        badtrials=[];
    case 'SP'
        badtrials=[];
    case 'TF'
        badtrials=[184,194];

    case 'TPS'
        badtrials=[];
    case 'YB'
        badtrials=[150];
        

end



%%
if ismember(itrial,badtrials)
    disp(['Skipping bad trial ' num2str(itrial) ' for ' subjID]);
    skip=1;
end
%%