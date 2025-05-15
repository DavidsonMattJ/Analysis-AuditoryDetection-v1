% j5_createNull_testGFX_aud
%% similar to the concat job (j4), except we disregard gait position information,
% and instead sample uniformly from all positions to create a null distribution.
%  [nbins x nPerm] per subject and datatype.

% Auditory
clear all; close all;

% this creates null data, per participant, requires about 10 s per participant.

set_myOnlineDirectories_AUD;

% show ppant numbers:
cd(procdatadir);
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%


%% concat data:
%preallocate storage:

gaittypes = {'single gait' , 'double gait'};
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preallocate storage:

%use same spacing as observed data, for the bins:
%add extra fields for the shuffled data.
cd([procdatadir filesep 'GFX'])
load('GFX_Data_inGaits', 'GFX_TargPos_nullData', 'GFX_RespPos_nullData','pidx1', 'pidx2');

% GFX_TargPos_nullData=[]; % store the results per subj, for stats.
% GFX_RespPos_nullData=[];

nPerm = 1000; % how many times to perform the resampling?

for ippant =1:length(pfols)
    tic

    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'trial_summaryTable','subjID');
    Ppant_gaitData= trial_summaryTable;

    Ltrials= strcmp(trial_summaryTable.trgO_gFoot, 'LR');

    Rtrials= strcmp(trial_summaryTable.trgO_gFoot, 'RL');

    disp(['constructing null distributions for ' subjID]);

    %%


    useindexes= {pidx1, pidx2};

    for nGait=2%1:2
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
        ppantData= Ppant_gaitData;

        if nGait==1
            useL=find(Ltrials);
            useR=find(Rtrials);

            % split by foot (LR, RL, all combined)
        else

            % note that here we have data in separate columns (so find non-nan
            % rows.)
            useL =find(~isnan(trial_summaryTable.trgO_gPcnt_LRL));
            useR =find(~isnan(trial_summaryTable.trgO_gPcnt_RLR));

        end
        % also the combined case:
        useAll= unique([useL; useR]);

        trialstoIndex ={useL, useR, useAll};

        %index for slow /nat speed trials.
        SlowWalk = find(trial_summaryTable.walkSpeed==1);
        NormWalk = find(trial_summaryTable.walkSpeed==2);
        speedIndex= {SlowWalk, NormWalk, sort([SlowWalk;NormWalk])};



        % split by foot (LR, RL, all combined)
        % !the combined case counts everything twice.

        for iSpeed = 1:3
            for iLR=1:3

                useft= find(trialstoIndex{iLR});
                useSpd = speedIndex{iSpeed};
                uset = intersect(useft, useSpd);

                % note that if we have the double gc and combined case.
                % also need the following intersection:
                if iLR==3 && nGait==2
                    Lcomb = intersect(useL, useSpd);
                    Rcomb = intersect(useR, useSpd);

                end

                %% restruct to relevant data (using index for trgO and respO)
                % which restricted dataset are we using?
                if nGait==1

                    %Onsets:
                    tmp_targOnsets = ppantData.trgO_gPcnt(uset);
                    tmp_respOnsets = ppantData.respO_gPcnt(uset);

                    %DVs:
                    % Order of DVs, for each search IV above:
                    tmp_Correct = ppantData.targCor(uset);
                    tmp_clickRT = ppantData.clickRT(uset);
                    % tmp_clickRT = ppantData.z_clickRT(uset);

                    tmp_SDTcat = ppantData.SDTcat(uset);

                elseif nGait==2  % stride based (LRL, or RLR)

                    %Onsets
                    if iLR==1
                        tmp_targOnsets = ppantData.trgO_gPcnt_LRL(uset); %LRL onset
                        tmp_respOnsets = ppantData.respO_gPcnt_LRL(uset);
                    elseif iLR==2
                        tmp_targOnsets = ppantData.trgO_gPcnt_RLR(uset); %RLR onset
                        tmp_respOnsets = ppantData.respO_gPcnt_RLR(uset);
                    elseif iLR==3
                        tmp_targOnsets= [ppantData.trgO_gPcnt_LRL(Lcomb); ppantData.trgO_gPcnt_RLR(Rcomb)];
                        tmp_respOnsets= [ppantData.respO_gPcnt_LRL(Lcomb); ppantData.respO_gPcnt_RLR(Rcomb)];
                    end


                    % DVs unchanged:
                    % Order of DVs, for each search IV above:
                    tmp_Correct = ppantData.targCor(uset);
                    tmp_clickRT = ppantData.clickRT(uset);
                    % tmp_clickRT = ppantData.z_clickRT(uset);
                    tmp_SDTcat = ppantData.SDTcat(uset);

                    %xcept in combined case, where we need to extract
                    %from two columns:
                    %overwrite the DVs:
                    if iLR==3
                        % special case, use all data (both columns stacked).
                        tmp_Correct = [ppantData.targCor(Lcomb);ppantData.targCor(Rcomb)];
                        tmp_clickRT = [ppantData.clickRT(Lcomb);ppantData.clickRT(Rcomb)];
                        % tmp_clickRT = [ppantData.z_clickRT(Lcomb);ppantData.z_clickRT(Rcomb)];
                        tmp_SDTcat = [ppantData.SDTcat(Lcomb);ppantData.SDTcat(Rcomb)];
                    end
                end
                %% %%%%%%%%%%%%%%%%
                %% for gait percent 1-100, calculate DVs
                % - acc, rt, dprime, crit etc.

                % then calculate binnedversions

                %preallocate (nPerm x pidx):
                [targOns_Acc_shuff, respOns_Acc_shuff,...
                    targOns_RTs_shuff, respOns_RTs_shuff,...
                    targOns_dprime_shuff, respOns_dprime_shuff,...
                    targOns_crit_shuff, respOns_crit_shuff,...
                    targOns_Counts_shuff, respOns_Counts_shuff]  =deal(nan(nPerm, pidx(end)));

                % for the binned dprime and criterion, we need to store
                % the counts of H,M,DA,CR, at each indx
                [ targOns_SDTcat_shuff,respOns_SDTcat_shuff] = deal(nan(1,pidx(end),4)); % 1,2,3,4 for H,M,FA,CR


                % note that now, we will sample randomly from all
                % positions, based on the number in real data.

                % the easiest way is to simply find all the data in our
                % restricted set of trials, shuffle, and place back in the
                % same index.

                % then repeat the shuffle and main analysis nPerm times for
                %tmp_targOnsets
                %tmp_respOnsets
                %tmp_Correct
                %tmp_clickRT
                %tmp_SDTcat
                %%%%%%%%%%%%%%%%%%

                %
                realOnsets= tmp_targOnsets;
                for iperm=1:nPerm

                    %shuffle the order:
                    permTargOnsets = shuffle(tmp_targOnsets);
                    permRespOnsets = shuffle(tmp_respOnsets);

                    % now repeat the analysis from j4, storing result per
                    % shuffle:


                    for ip=1:pidx(end)
                        %use the location in the shuffled index:
                        randtargIDX = find(permTargOnsets==ip);
                        randrespIDX = find(permRespOnsets==ip);

                        % Accuracy:
                        targOns_Acc_shuff(iperm,ip) = nanmean(tmp_Correct(randtargIDX));
                        respOns_Acc_shuff(iperm,ip) = nanmean(tmp_Correct(randrespIDX));

                        %RT:
                        targOns_RTs_shuff(iperm,ip) = nanmean(tmp_clickRT(randtargIDX));
                        respOns_RTs_shuff(iperm,ip) = nanmean(tmp_clickRT(randrespIDX));

                        %dprime and criterion

                        useOnsets= {randtargIDX, randrespIDX};

                        for itargResp=1:2
                            onsetidx = useOnsets{itargResp};
                            %store counts for calculations:
                            Hitn = length(find(tmp_SDTcat(onsetidx)==1));
                            Missn = length(find(tmp_SDTcat(onsetidx)==2));

                            FAn = length(find(tmp_SDTcat(onsetidx)==3));
                            CRn = length(find(tmp_SDTcat(onsetidx)==4));

                            HitRate = Hitn / (Hitn+Missn);
                            FArate = FAn / (FAn + CRn);
                            % calc:
                            dp = norminv(HitRate)- norminv(FArate);
                            crit = -0.5*(norminv(HitRate)+ norminv(FArate));

                            if itargResp==1
                                targOns_dprime_shuff(iperm,ip) = dp;
                                targOns_crit_shuff(iperm,ip) = crit;
                                targOns_SDTcat_shuff(iperm,ip,1) = Hitn;
                                targOns_SDTcat_shuff(iperm,ip,2) = Missn;
                                targOns_SDTcat_shuff(iperm,ip,3) = FAn;
                                targOns_SDTcat_shuff(iperm,ip,4) = CRn;

                            else
                                respOns_dprime_shuff(iperm,ip) = dp;
                                respOns_crit_shuff(iperm,ip) = crit;
                                respOns_SDTcat_shuff(iperm,ip,1) = Hitn;
                                respOns_SDTcat_shuff(iperm,ip,2) = Missn;
                                respOns_SDTcat_shuff(iperm,ip,3) = FAn;
                                respOns_SDTcat_shuff(iperm,ip,4) = CRn;
                            end

                        end

                        % !!  to store null counts counts, see below (slightly different case)


                    end % ip (1-100)

                end % perm (nPerm)
                %%
                %% %%%%%%%%%%%%%%%%
                % Step through slightly different analysis for shuffling
                % counts per timepoint.


                % 2 main analysis (for counts).

                %-targOnset
                %-clickOnset

                %work out what the count per actual IDX was, then shuffle
                %these counts.

                %%%%%%%%%%%%%%%%%%


                %% store the histcounts, for target contrast level per pos.

                % start with the real onsets (not shuffled) for each type
                % (target and response)

                usePos= {tmp_targOnsets, tmp_respOnsets};
                for itype=1:2
                    useposdata = usePos{itype};

                    % raw hit counts, per gaitsample:
                    real_counts= histcounts(useposdata, pidx(end));

                    %critical step -convert 0 to nans. this avoids averaging
                    %problems (as above, in data anaysis)
                    % real_counts(real_counts==0)=NaN;

                    %now create a null, by shuffling these counts per perm.

                    % preallocate:
                    outgoingnull = nan(nPerm, pidx(end));
                    for iperm = 1:nPerm


                        a = rand(1,length(real_counts));
                        [~, i] = sort(a);

                        outgoingnull(iperm,:) =real_counts(i);
                    end

                    switch itype
                        case 1
                            targOns_Counts_shuff = outgoingnull;
                        case 2
                            respOns_Counts_shuff = outgoingnull;
                    end
                end % both types.


                %% %%%%%%%%%%%%%%%%%%%
                % BINNING
                % Step through all shuffled data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%
                % This first version simply averages the data within a bin.
                %
                % this will be compared to the specific SDT case of
                % recalculating SDT metrics based on raw data in the bin
                % (to avoid floor/ceiling effects for example).


                BinMeData = {targOns_Acc_shuff, respOns_Acc_shuff, ...
                    targOns_RTs_shuff, respOns_RTs_shuff,...
                    targOns_dprime_shuff,respOns_dprime_shuff,... % expect lots of nans for data paucity.
                    targOns_crit_shuff,respOns_crit_shuff,... % expect lots of nans for data paucity.
                    targOns_Counts_shuff, respOns_Counts_shuff};

                % note dprime and criterion were included for comparison
                % with the alternative below (calculating based on data
                % within bin, not averaging over bin indexes).

                for itype = 1:length(BinMeData)

                    usedata = BinMeData{itype};
                    outgoingbinned=nan(nPerm, length(pidx)-1);

                    for ibin=1:length(pidx)-1
                        idx = pidx(ibin):pidx(ibin+1);

                        outgoingbinned(:,ibin) = nanmean(usedata(:,idx),2);

                    end

                    %rename this time. to avoid headaches
                    switch itype
                        case 1
                            targOns_Acc_bin= outgoingbinned;
                        case 2
                            respOns_Acc_bin= outgoingbinned;
                        case 3
                            targOns_RT_bin= outgoingbinned;
                        case 4
                            respOns_RT_bin= outgoingbinned;
                        case 5
                            targOns_dP_bin= outgoingbinned;
                        case 6
                            respOns_dP_bin= outgoingbinned;
                        case 7
                            targOns_Crit_bin= outgoingbinned;
                        case 8
                            respOns_Crit_bin= outgoingbinned;
                        case 9
                            targOns_Counts_bin= outgoingbinned;
                        case 10
                            respOns_Counts_bin= outgoingbinned;
                    end
                end % per type


                %% %%%%%%%%%%%%%%%%%%%
                % Step through all shuffled data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%

                %Acc, Counts, resp Counts,
                BinMeData = {targOns_Acc_shuff, respOns_Acc_shuff, ...
                    targOns_RTs_shuff, respOns_RTs_shuff,...
                    targOns_dprime_shuff,respOns_dprime_shuff,... % expect lots of nans for data paucity.
                    targOns_crit_shuff,respOns_crit_shuff,... % expect lots of nans for data paucity.
                    targOns_Counts_shuff, respOns_Counts_shuff};

                for itype = 1:length(BinMeData)

                    usedata = BinMeData{itype};
                    outgoingbinned=nan(nPerm, length(pidx)-1);

                    for ibin=1:length(pidx)-1
                        idx = pidx(ibin):pidx(ibin+1);

                        outgoingbinned(:,ibin) = nanmean(usedata(:,idx),2);

                    end


                    %rename this time. to avoid headaches
                    switch itype
                        case 1
                            targOns_Acc_bin_shuff= outgoingbinned;
                        case 2
                            respOns_Acc_bin_shuff= outgoingbinned;
                        case 3
                            targOns_RT_bin_shuff= outgoingbinned;
                        case 4
                            respOns_RT_bin_shuff= outgoingbinned;
                        case 5
                            targOns_dP_bin_shuff= outgoingbinned;
                        case 6
                            respOns_dP_bin_shuff= outgoingbinned;
                        case 7
                            targOns_Crit_bin_shuff= outgoingbinned;
                        case 8
                            respOns_Crit_bin_shuff= outgoingbinned;
                        case 9
                            targOns_Counts_bin_shuff= outgoingbinned;
                        case 10
                            respOns_Counts_bin_shuff= outgoingbinned;
                    end
                end % per type


                %% to shuffle the dprime and criterion results, select 10 pos at random (non consecutive).
                % special cases for sensitivity and criterion. We will compute
                % the sensitivity and criterion based on all data within the
                % bin.
                BinMeData = {targOns_SDTcat_shuff, respOns_SDTcat_shuff};

                for itype = 1:2% targ ons, resp ons

                    usedata = BinMeData{itype};


                    [dprime_bin,criterion_bin, HitRate_bin, FARate_bin]=...
                        deal(nan(nPerm,length(pidx)-1));

                    for ibin=1:length(pidx)-1
                        idx = pidx(ibin):pidx(ibin+1);

                        HITn_bin = sum(usedata(:,idx,1),2); % sum counts over bin range
                        MISSn_bin = sum(usedata(:,idx,2),2); % sum counts over bin range
                        FAn_bin = sum(usedata(:,idx,3),2); % sum counts over bin range
                        CRn_bin = sum(usedata(:,idx,4),2); % sum counts over bin range

                        Hratebin = HITn_bin./ (HITn_bin + MISSn_bin);
                        FAratebin = FAn_bin./ (FAn_bin + CRn_bin);


                        %remove 0 and 1 floor and ceiling effects
                        %(stuffs calcs).

                        FAratebin(FAratebin==0)=.001;

                        Hratebin(Hratebin==0)=.001;

                        Hratebin(Hratebin==1)=.999;

                        FAratebin(FAratebin==1)=.999;


                        % alternatively:
                        % Correct for Hit Rate = 1 and FA rate = 0 in bins
                        % by adjusting only the extreme values by replacing
                        % rates of 0 with 0.5/n and rates of 1 with (n−0.5)/n
                        % where n is the total number of targets for condition
                        %  in the bin (HR) or total number of responses in the
                        % bin (FA) (Macmillan & Kaplan, 1985).

                        % FAratebin(FAratebin==0) = 0.5/(FAn_bin + CRn_bin);
                        % Hratebin(Hratebin==0) = 0.5/(HITn_bin + MISSn_bin);
                        %
                        %
                        % FAratebin(FAratebin==1) = (FAn_bin + CRn_bin)-0.5/(FAn_bin + CRn_bin);
                        % Hratebin(Hratebin==1) = (HITn_bin + MISSn_bin)-0.5/(HITn_bin + MISSn_bin);


                        %% store

                        HitRate_bin(:,ibin) = Hratebin;
                        FARate_bin(:,ibin) = FAratebin;

                        %calculate dprime and criterion, per bin.
                        dprime_bin(:,ibin) = norminv(Hratebin)- norminv(FAratebin);
                        criterion_bin(:,ibin) =  -0.5*(norminv(Hratebin)+ norminv(FAratebin));




                        if isinf(dprime_bin(:,ibin))
                            disp('check')
                        end

                    end % for ibin

                    if itype==1

                        targOns_dprime_bin_shuff = dprime_bin;
                        targOns_criterion_bin_shuff =criterion_bin;

                        targOns_HitRate_bin_shuff = HitRate_bin;
                        targOns_FARate_bin_shuff= FARate_bin;

                    elseif itype==2

                        respOns_dprime_bin_shuff = dprime_bin;
                        respOns_criterion_bin_shuff =criterion_bin;

                        respOns_HitRate_bin_shuff = HitRate_bin;
                        respOns_FARate_bin_shuff= FARate_bin;

                    end


                end % itype (targ, resp)

                %% store for all ppants, preparing to save:
                gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:

               
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'counts']) = targOns_Counts_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'Acc']) = targOns_Acc_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'rts']) = targOns_RTs_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'dprime']) =targOns_dprime_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'crit']) =targOns_crit_shuff;
                
                %binned version:
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts']) = targOns_Counts_bin_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc']) = targOns_Acc_bin_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts']) = targOns_RT_bin_shuff;
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime2']) = targOns_dP_bin_shuff; % averaged within bin, '2' since the second preference
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit2']) = targOns_Crit_bin_shuff; % averaged within bin
                                
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime']) = targOns_dprime_bin_shuff; % calculated for bin
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit']) = targOns_criterion_bin_shuff; % calculated for bin               
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR']) = targOns_HitRate_bin_shuff; % calculated for bin
                GFX_TargPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA']) = targOns_FARate_bin_shuff; % calculated for bin
                
                %also using response position as index:
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'counts']) = respOns_Counts_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'Acc']) = respOns_Acc_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'rts']) = respOns_RTs_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'dprime']) =respOns_dprime_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'crit']) =respOns_crit_shuff;
                
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts']) = respOns_Counts_bin_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc']) = respOns_Acc_bin_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts']) = respOns_RT_bin_shuff;                
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime2']) = respOns_dP_bin_shuff; % averaged within bin, '2' since the second preference
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit2']) = respOns_Crit_bin_shuff; % averaged within bin
                
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime']) = respOns_dprime_bin_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit']) = respOns_criterion_bin_shuff;               
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR']) = respOns_HitRate_bin_shuff;
                GFX_RespPos_nullData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA']) = respOns_FARate_bin_shuff;
                %                 GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = respRT_bin_shuff;


            end % iLR

            disp(['finished speed ' num2str(iSpeed) ' of 3']);
        end % iSpeed
    end %ngait
    toc
end % ppant

%%
cd([procdatadir filesep 'GFX']);

disp('resaving GFX data with shuffle attached...');

save('GFX_Data_inGaits', ...
    'GFX_TargPos_nullData','GFX_RespPos_nullData',...
    '-append');