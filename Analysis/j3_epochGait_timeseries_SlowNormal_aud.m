% j3_epochGait_timeseries_SlowNormal

% using the trough information in the summary table (or HeadPos), we will
% now epoch the raw, resampled, and normalized versions of the gait
% time-series (for later plots).

%to do: get response gaits as well (at the moment only gaits when a target
%was presented).


% Auditory version % Two walking speeds
%laptop:

set_myOnlineDirectories_AUD;

cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%%
 

resampSize=100; 


for ippant =1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j3 gait time series data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.



    %create a large empty matrix to staw raw gait, and resampled sizes
    %also:

    gait_ts_raw= zeros(1,500); % in samples
    gait_ts_resamp= zeros(1,100); % resampled to 100%
    [trialallocation, gaitDuration,gaitIdx, gaitSamps, gaitStart, gaitFin, walkSpeed, SDTcat]= deal(nan);
    gaitFeet= {'L'};
    % including walk speed and result (H,M,FA,CR).
    gait_ts_gData= table(trialallocation, walkSpeed, SDTcat, gaitIdx, gaitFeet,  gaitStart,gaitSamps,gaitDuration);
    

    % lets also grab the double GC:
    doubgait_ts_raw= zeros(1,500);
    doubgait_ts_resamp= zeros(1,200);

    gcounter=1;

%find all gait onsets (in samples)
gOnsets = find(~isnan(trial_summaryTable.trgO_gStart_samp));

% % This works for single Gaits only:
    for igait = 1:length(gOnsets)

        rowIndx = gOnsets(igait);
        itrial = trial_summaryTable.trial(rowIndx);
        
        %check if this trial should be skipped:
        skip=0;
        rejTrials_detectvAud;
        if skip ==1
            continue;
        end
        % continue:
        gaitStFin = [trial_summaryTable.trgO_gStart_samp(rowIndx),trial_summaryTable.trgO_gFin_samp(rowIndx)];
        gaitDur = trial_summaryTable.trgO_gFin_sec(rowIndx)- trial_summaryTable.trgO_gStart_sec(rowIndx);
        gFt = trial_summaryTable.trgO_gFoot(rowIndx);

        % % % now using a detrended version (smoother plots).
        % for each gait we have, extract the raw head time series,
       

% detrend the gait to avoid slant effects?
headY = detrend(HeadPos(itrial).Y);
% headY = (HeadPos(itrial).Y);
rawHead = headY(gaitStFin(1):gaitStFin(2));
% Eye_Y = detrend(EyePos(itrial).Y_clean);
% rawEye_Y= Eye_Y(gaitStFin(1):gaitStFin(2));

% while here, also extract eye origin and gaze direction.
% %  rawHead = HeadPos(itrial).Y(gaitStFin(1):gaitStFin(2));        
% rawGaze_Y = EyeDir(itrial).Y_clean(gaitStFin(1):gaitStFin(2));
% rawGaze_Z = EyeDir(itrial).Z_clean(gaitStFin(1):gaitStFin(2));
% % rawEye_Y = EyePos(itrial).Y_clean(gaitStFin(1):gaitStFin(2));
% rawEye_Z = EyePos(itrial).Z_clean(gaitStFin(1):gaitStFin(2));

        %also resample along x dimension:
        resampHead = imresize(rawHead', [1,resampSize]);
%         resampGY= imresize(rawGaze_Y', [1,resampSize]);
%         resampGZ= imresize(rawGaze_Z', [1,resampSize]);
%         resampOY= imresize(rawEye_Y', [1,resampSize]);
%         resampOZ= imresize(rawEye_Z', [1,resampSize]);

        
        
        %gaitcount(within trial)
        gIdx = find(HeadPos(itrial).Y_gait_troughs==(gaitStFin(1)));
        
        %>>>>> store:
        gait_ts_raw(igait,1:length(rawHead)) = rawHead;
        gait_ts_resamp(igait,:) = resampHead;
%         gait_ts_eyedirY_resamp(igait,:)= resampGY;
%         gait_ts_eyedirZ_resamp(igait,:)= resampGZ;
%         gait_ts_eyeposY_resamp(igait,:)= resampOY;
%         gait_ts_eyeposZ_resamp(igait,:)= resampOZ;
% 
%         gait_ts_gData.trialallocation(igait) = itrial;
        gait_ts_gData.gaitDuration(igait)= gaitDur;
        gait_ts_gData.gaitFeet(igait) = gFt;
        gait_ts_gData.gaitIdx(igait) = gIdx;
        gait_ts_gData.gaitStart(igait) = gaitStFin(1);
        gait_ts_gData.gaitSamps(igait) = length(rawHead);

%         table(trialallocation, walkSpeed, SDTcat, gaitIdx, gaitFeet,  gaitStart,gaitSamps,gaitDuration);
    gait_ts_gData.trialallocation(igait)= itrial;
    gait_ts_gData.SDTcat(igait)= trial_summaryTable.SDTcat(rowIndx);
    gait_ts_gData.walkSpeed(igait)= trial_summaryTable.walkSpeed(rowIndx);
        %prefill next row to suppress output warnings in command window:
        gait_ts_gData(igait+1,:) = gait_ts_gData(igait,:) ;
    end
    %remove last row (autofilled)
    gait_ts_gData(igait+1,:)=[];
    %convert zeros to nan to tidy plot:
    gait_ts_raw(gait_ts_raw==0)=nan;
    gait_ts_resamp(gait_ts_resamp==0)=nan;

%% % now also epoch the doubgc Head data, by iterating over all gaits and epoching n+1.
% 
% sltrials = find(gait_ts_gData.walkSpeed==1);
% 
% ntrials = find(gait_ts_gData.walkSpeed==2);
% 
% clf;
% t1 = squeeze(mean(gait_ts_raw(sltrials,:),1));
% 
% t2 = squeeze(mean(gait_ts_raw(ntrials,:),1));
% 
% plot(t1); hold on; plot(t2, 'r');
% shg
%%

  for igait = 1:length(gOnsets)

        rowIndx = gOnsets(igait);
        itrial = trial_summaryTable.trial(rowIndx);
        
        %check if this trial should be skipped:
        skip=0;
        rejTrials_detectvAud; %toggles skip based on bad trial ID (itrial)
        if skip ==1
            continue;
        end
        % continue by looking for gait +1.

        gaitIdx = trial_summaryTable.trgO_gCount(rowIndx);
        allGs = HeadPos(itrial).Y_gait_troughs;
        %select either the current or prev started step for this doubGC.
        try
            gaitStFin= allGs(gaitIdx):allGs(gaitIdx+2);
        catch
            gaitStFin= allGs(gaitIdx-1):allGs(gaitIdx+1);
        end
            
%      

% detrend the gait to avoid slant effects?
headY = detrend(HeadPos(itrial).Y);
% headY = (HeadPos(itrial).Y);
rawHead = headY(gaitStFin);

        %also resample along x dimension:
        resampHead = imresize(rawHead', [1,200]); % double the length.
 %>>>>> store:
        doubgait_ts_raw(igait,1:length(rawHead)) = rawHead;
        doubgait_ts_resamp(igait,:) = resampHead;





  end
    %convert zeros to nan to tidy plot:
    doubgait_ts_raw(gait_ts_raw==0)=nan;
    doubgait_ts_resamp(gait_ts_resamp==0)=nan;





disp(['saving gait and double gait time series data for ' subjID])
cd(procdatadir);
save(savename, 'gait_ts_raw', 'gait_ts_resamp', 'gait_ts_gData', ...
    'doubgait_ts_raw', 'doubgait_ts_resamp','-append');



end % per ppant

