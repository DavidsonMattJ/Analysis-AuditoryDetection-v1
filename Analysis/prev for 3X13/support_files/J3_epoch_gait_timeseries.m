% j3_epochGait_timeseries

% using the trough information in the summary table (or HeadPos), we will
% now epoch the raw, resampled, and normalized versions of the gait
% time-series (for later plots).

%to do: get response gaits as well (at the moment only gaits when a target
%was presented).



resampSize=100; 

cd(savedatadir);
pfols = dir([pwd filesep 'p_*']);
%%
for ippant =1:length(pfols)
    cd(savedatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'summary_table', 'subjID')        
    savename = pfols(ippant).name;
    disp(['Preparing j3 gait time series data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.


    allevents = size(summary_table,1);

    %create a large empty matrix to staw raw gait, and resampled sizes
    %also:

    gait_ts_raw= zeros(1,500); % in samples
    [gait_ts_resamp, gait_ts_eyedirY_resamp,...
        gait_ts_eyedirZ_resamp,gait_ts_eyeposY_resamp,...
        gait_ts_eyeposZ_resamp] = deal(zeros(1,100)); % resampled to 100%


    % lets also grab the double GC:
    doubgait_ts_raw= zeros(1,500);
    doubgait_ts_resamp= zeros(1,200);

    [trialallocation, gaitDuration,gaitIdx, gaitSamps, gaitStart, gaitFin, walkSpeed]= deal(nan);
    gaitFeet= {'L'};
    
    % create an empty table with useful information, so that we can access
    % appropriately later
    gait_ts_gData= table(trialallocation,gaitIdx, walkSpeed, gaitFeet,  gaitStart,gaitSamps,gaitDuration);
    
%find all gait onsets (in samples)
gOnsets = find(~isnan(summary_table.trgO_gStart_samp));

% % This works for single Gaits only:
    for igait = 1:length(gOnsets)

        rowIndx = gOnsets(igait);
        itrial = summary_table.trial(rowIndx);
        
        %check if this trial should be skipped:
        skip=0;
        rejTrials_AUD_detect_v1; %toggles skip based on bad trial ID (itrial)
        if skip ==1
            continue;
        end
        % continue:
        gaitStFin = [summary_table.trgO_gStart_samp(rowIndx),summary_table.trgO_gFin_samp(rowIndx)];
        gaitDur = HeadPos(itrial).time(gaitStFin(2) - gaitStFin(1)); 
        gFt = summary_table.trgO_gFoot(rowIndx);

        % % % now using a detrended version (smoother plots).
        % for each gait we have, extract the raw head time series,
       

% detrend the gait to avoid slant effects?
headY = detrend(HeadPos(itrial).Y);
% headY = (HeadPos(itrial).Y);
rawHead = headY(gaitStFin(1):gaitStFin(2));

        %also resample along x dimension:
        resampHead = imresize(rawHead', [1,resampSize]);

        
        
        %gaitcount(within trial)
        gIdx = find(HeadPos(itrial).Y_gait_troughs==(gaitStFin(1)));
        
        %>>>>> store:
        gait_ts_raw(igait,1:length(rawHead)) = rawHead;
        gait_ts_resamp(igait,:) = resampHead;
    
        %store data for this gait event, in our table
        gait_ts_gData.trialallocation(igait) = itrial;
        gait_ts_gData.walkSpeed(igait) = HeadPos(itrial).blockType;  
        gait_ts_gData.gaitDuration(igait)= gaitDur;
        gait_ts_gData.gaitFeet(igait) = gFt;
        gait_ts_gData.gaitIdx(igait) = gIdx;
        gait_ts_gData.gaitStart(igait) = gaitStFin(1);
        gait_ts_gData.gaitSamps(igait) = length(rawHead);
        %prefill next row to suppress output warnings in command window:
        gait_ts_gData(igait+1,:) = gait_ts_gData(igait,:) ;
    end
    %remove last row (autofilled)
    gait_ts_gData(igait+1,:)=[];
    %convert zeros to nan to tidy plot:
    gait_ts_raw(gait_ts_raw==0)=nan;
    gait_ts_resamp(gait_ts_resamp==0)=nan;

%% % now also epoch the doubgc Head data, by iterating over all gaits and epoching n+1.



  for igait = 1:length(gOnsets)

        rowIndx = gOnsets(igait);
        itrial = summary_table.trial(rowIndx);
        
        %check if this trial should be skipped:
        skip=0;
        rejTrials_AUD_detect_v1; %toggles skip based on bad trial ID (itrial)
        if skip ==1
            continue;
        end
        % continue by looking for gait +1.

        gaitIdx = summary_table.trgO_gCount(rowIndx);
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
cd(savedatadir);
save(savename, 'gait_ts_raw', 'gait_ts_resamp', 'gait_ts_gData', ...
    'doubgait_ts_raw', 'doubgait_ts_resamp','-append');

    
%%  sanity check  - plot results
 slowSteps = find(gait_ts_gData.walkSpeed==1);
 natSteps = find(gait_ts_gData.walkSpeed==2);
 
 %plot both to check for bugs:
 useData = {gait_ts_raw, doubgait_ts_raw};
 clf;
 for isteplength=1:2
 subplot(1,2,isteplength);
 %slow,natural, single step:
 slowS= nanmean(useData{isteplength}(slowSteps,:),1);
  natS= nanmean(useData{isteplength}(natSteps,:),1);
  plot(slowS, 'b'); hold on
  plot(natS, 'r'); hold on
  xlim([0 60*isteplength]);
  ylabel('detrended head height')
  legend('slow', 'natural')
 end
     
     %%
 

%     %% sanity check to find long gaits in trials:
%     clf;
%     subplot(211)
%     plot(gait_ts_raw');
%     subplot(212);
%       plot(gait_ts_resamp');
% %%
%     %
%     hold on;
%     for igait = 1:size(gait_ts_raw,1);
%         %plot trial number
%         %plot at end point:
%         tmp = gait_ts_raw(igait,:);
%         lastsamp= max(find(~isnan(tmp)));
%         %if large add text.
%         if lastsamp >75
%             text(lastsamp, tmp(lastsamp), num2str([gait_ts_trialallocation(igait)]))
%         end
%     end % add gait text if long:

end % per ppant

