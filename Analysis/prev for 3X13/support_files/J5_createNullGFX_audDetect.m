% j5_createNullGFX_audDetect
%% similar to the concat job (j4), except we disregard gait position information, 
 % and instead sample uniformly from all positions to create a null distribution.
%  [nbins x nPerm] per subject and datatype.
%% show ppant numbers:
cd(savedatadir);
pfols = dir([pwd filesep '*_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
procdatadir= savedatadir; % compatability below.

%% concat data:
%preallocate storage:

gaittypes = {'single gait' , 'double gait'};
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %preallocate storage:
   
     %use same spacing as observed data, for the bins:
     %add extra fields for the shuffled data.
     cd([procdatadir filesep 'GFX'])
    load('GFX_Data_inGaits', 'GFX_TargPosData', 'GFX_RespPosData','pidx1', 'pidx2');
    
%     GFX_TargPos_nullData=[]; % store the results per subj, for stats.
%     GFX_RespPos_nullData=[];
    
    nPerm = 1000; % how many times to perform the resampling?
    
    for ippant =1:length(pfols)
        
        cd(procdatadir)    %%load data from import job.
        load(pfols(ippant).name, ...
            'summary_table','subjID');
    % should have been done earlier, but make a column of clickRT for
    % compatability below:
    summary_table.clickRT = summary_table.targRT - summary_table.targOnset;
    
        Ppant_gaitData= summary_table;
    
        Ltrials= strcmp(summary_table.trgO_gFoot, 'LR');
    
        Rtrials= strcmp(summary_table.trgO_gFoot, 'RL');
       
        disp(['constructing null distributions for ' subjID]);
 
        %%

        useindexes= {pidx1, pidx2};
        
        for nGait=2% 
    
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
        ppantData= Ppant_gaitData;
        if nGait==1
        useL=Ltrials;
        useR=Rtrials;
        useAll = 1:length(Ltrials);
        else
            
            useL =~isnan(summary_table.trgO_gPcnt_LRL);
            useR =~isnan(summary_table.trgO_gPcnt_RLR);
            useAll= 1:length(Ltrials);
            
        end
        
       trialstoIndex ={find(useL), find(useR), useAll};
        
        
        % also include a sub index for different walk speeds.
        %in summary table, index of different speeds:
        slowspeedTrials = find(summary_table.walkSpeed==1);
        naturalspeedTrials = find(summary_table.walkSpeed==2);
        % all trials will be their combination.
        allspd = [slowspeedTrials; naturalspeedTrials];
        speedstoIndex={slowspeedTrials, naturalspeedTrials,sort(allspd) };
        
        
        
        for iSpeed=1:3        
            for iLR=3%1:3
            
                useft= trialstoIndex{iLR};
                usespd = speedstoIndex{iSpeed};
                uset = intersect(useft,usespd); % trials that satisfy both conditions.
        
                % note that if we have the double gc and combined case.
                % also need the following intersection:
                if iLR==3 && nGait==2
                    Lcomb = intersect(find(useL), usespd);
                    
                    Rcomb = intersect(find(useR), usespd);
                    
                end
                
                %% %%%%%%%%%%%%%%%% 
                % Step through different data types :
                
                % note that now, we will sample randomly from all
                % positions, based on the number in real data.
                
               % Main analysis focus on performance at a single gait percent
                %- Relative to target onset in gait:
                % -- Accuracy
                % -- RT
                % -- dprime
                % -- criterion
                %- Relative to response onset in gait:
                % -- Accuracy
                % -- RT
                % -- dprime
                % -- criterion
        
                 if nGait==1 % single step
                    
                    searchPosIDX = {ppantData.trgO_gPcnt(uset), ... % trg onset as %
                        ppantData.respO_gPcnt(uset)}; % resp O as pcnt}
                
                else %(Much messier, but due to overlap in events occurring across 2 gaits.)
                    %use strides:
                    strideSuffix = {'LRL', 'RLR'};
                    if iLR==1
                        searchPosIDX = {ppantData.trgO_gPcnt_LRL(uset), ... % trg onset as %
                            ppantData.respO_gPcnt_LRL(uset)}; % resp onset as %
                    elseif iLR==2
                        searchPosIDX = {ppantData.trgO_gPcnt_RLR(uset), ... % trg onset as %
                            ppantData.respO_gPcnt_RLR(uset)}; % resp onset as %
                        
                    elseif iLR==3
                        % both. Problem is they can have identical row info.
                        
                        searchPosIDX = {[ppantData.trgO_gPcnt_LRL(Lcomb); ppantData.trgO_gPcnt_RLR(Rcomb)], ... % trg onset as %
                            [ppantData.respO_gPcnt_LRL(Lcomb); ppantData.respO_gPcnt_RLR(Rcomb)]};
                    end
                end
                
                % we have different operations to perform:
                
                % Order of DVs, for each search IV above:
                searchPosDVs = {ppantData.correctResponse(uset),... % acc
                    ppantData.clickRT(uset),...  %RT
                    ppantData.SDTcat(uset),... % for dprine
                    ppantData.SDTcat(uset)}; % for criterion
                
                % %except in this particular case:
                if iLR==3 && nGait==2
                    
                    %and define the DV of interest:
                    searchPosDVs = {[ppantData.correctResponse(Lcomb);ppantData.correctResponse(Rcomb)], ... %Acc, RT,
                        [ppantData.clickRT(Lcomb);ppantData.clickRT(Rcomb)],...
                        [ppantData.SDTcat(Lcomb);ppantData.SDTcat(Rcomb)],...
                        [ppantData.SDTcat(Lcomb);ppantData.SDTcat(Rcomb)]};
                    
                end
                
                
                %also store the HR and FAR per pos, for later binning.
                 outgoing_SDT= nan(2, pidx(end));
                 
                 for iAlignto= 1:2 % targOnset or respOnset.
                     
                     for iAnalysis=1:4
                         %preallocate:
                         outgoingnull = nan(nPerm,pidx(end)); %
                        allpos = searchPosIDX{iAlignto};
                        allDVs = searchPosDVs{iAnalysis};
                        null_CVs = nan(3,nPerm); % 5, median, and 95% intervals.
                    
                        allpos_real = allpos(~isnan(allpos));
                    allDVs_real = allDVs(~isnan(allpos));
                
                % for actual gait pos, select another at random:
                for ip=1:pidx(end)     
                     
%                     % we may want to match the N datapoints:
%                     % how many observed data points this index?
%                      actualIDX=find(allpos_real==ip);
%                      actualN = length(actualIDX);
%                      actualObs = nanmean(allRTs_real(actualIDX));

                        tmpstor=nan(1, nPerm);                     
                    for iperm = 1:nPerm
                        % instead of the actual index (ip), replace with an indices at
                        % random from all possible:
                        useIDX = randi(pidx(end),1);
                        
%                         thisrand = randi(length(allpos_real),[1,actualN]);
                        thisrand= find(allpos_real==useIDX);
                        
                        useDVs = allDVs_real(thisrand);
                        
                         %% Compute accuracy.
                            if iAnalysis==1
%                                 
                                
                                % remove nans from length calcs.
                                tmpResp = useDVs;
                                useResp = tmpResp(~isnan(tmpResp));
                                tmpstor(iperm)= sum(useResp, 'omitnan') / length(useResp);
                                
                                
                                %% take mean for these  trials (if RT):
                            elseif iAnalysis==2
                                
                                tmpstor(iperm) = mean(useDVs,1, 'omitnan');
                                
                            elseif iAnalysis ==3
                                %compute sensitivity
                                %sensitivity from HR and FAR.
                                %HR is HR/ H+M;
                                % FA is FA/ FA+CR;
                                tmpResp= useDVs;
                                
                                Hitn = length(find(tmpResp==1));
                                Missn = length(find(tmpResp==2));
                                
                                FAn = length(find(tmpResp==3));
                                CRn = length(find(tmpResp==4));
                                
                                HitRate = Hitn / (Hitn+Missn);
                                FArate = FAn / (FAn + CRn);
                                dp = norminv(HitRate)- norminv(FArate);
                                tmpstor(iperm) = dp;
%                                 outgoing_SDT(1,ip)= HitRate;
%                                 outgoing_SDT(2,ip)= FArate;
                                
                            elseif iAnalysis==5 % compute criterion
                                tmpResp= useDVs;
                                
                                Hitn = length(find(tmpResp==1));
                                Missn = length(find(tmpResp==2));
                                
                                FAn = length(find(tmpResp==3));
                                CRn = length(find(tmpResp==4));
                                
                                HitRate = Hitn / (Hitn+Missn);
                                FArate = FAn / (FAn + CRn);
                                
                                crit = -0.5*(norminv(HitRate)+ norminv(FArate));
                                tmpstore(iperm) = crit;
                                
                                
                                
                            end % ianalysis
                            
                    end % for nPerm
%                    
                    %store null
                    outgoingnull(:, ip) = tmpstor;
                end
                   % now rename
                        switch iAnalysis
                            case 1
                                outgoing_Acc= outgoingnull;
                            case 2
                                outgoing_RT= outgoingnull;
                            case 3
                                outgoing_dprime= outgoingnull;
                                
                            case 4
                                outgoing_crit = outgoingnull;
                        end
                        
                    end  % ianalysis
                    
                    if iAlignto==1 % rename once more:
                        
                        targOns_Acc_null= outgoing_Acc;
                        targOns_RTs_null= outgoing_RT;
                        targOns_dprime_null = outgoing_dprime;
                        targOns_crit_null = outgoing_crit;
%                         targOns_HR = outgoing_SDT(1,:);
%                         targOns_FAR = outgoing_SDT(2,:);
                        
                    elseif iAlignto==2
                        
                        respOns_Acc_null= outgoing_Acc;
                        respOns_RTs_null= outgoing_RT;
                        respOns_dprime_null = outgoing_dprime;
                        respOns_crit_null = outgoing_crit;
%                         respOns_HR = outgoing_SDT(1,:);
%                         respOns_FAR = outgoing_SDT(2,:);
                        
                        
                    end
                    
                end %iAlign to.
                
                
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
                usePos= {ppantData.trgO_gPcnt(uset),...% gait pcnt for target onset.                
                ppantData.respO_gPcnt(uset)}; % gaitpcnt for resp onset.
            
            for itype=1:2
                useposdata = usePos{itype};
                
                % raw hit counts, per gaitsample:
                real_counts= histcounts(useposdata, pidx(end));
                
                %critical step -convert 0 to nans. this avoids averaging
                %problems (as above, in data anaysis)
                
                real_counts(real_counts==0)=NaN;
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
                        targOns_Counts_null = outgoingnull;
                    case 2
                        respOns_Counts_null = outgoingnull;
                end
            end % both types.
            
                %% %%%%%%%%%%%%%%%%%%%
                % Step through all shuffled data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%
               
                binmedata = {targOns_RTs_null, respOns_RTs_null, targOns_Acc_null,...
                    targOns_Counts_null, respOns_Counts_null};
                
                for itype = 1:length(binmedata)
                   
                   usedata = binmedata{itype};                   
                   outgoingbinned=nan(nPerm, length(pidx)-1);
                   
                   for ibin=1:length(pidx)-1
                       idx = pidx(ibin):pidx(ibin+1);
                       
                       outgoingbinned(:,ibin) = nanmean(usedata(:,idx),2);
                       
                   end
                    
                   %rename this time. to avoid headaches
                   switch itype
                       case 1
                           trgRT_bin_shuff= outgoingbinned;
                       case 2
                           respRT_bin_shuff= outgoingbinned;
                       case 3
                           trgAcc_bin_shuff= outgoingbinned;
                       case 4
                           trgCounts_bin_shuff = outgoingbinned;
                       case 5
                           respCounts_bin_shuff = outgoingbinned;
                   end
                end % per type
                %% extra handling for binned SDT output:
                   %% % special cases for sensitivity and criterion. We will compute
            % the sensitivity and criterion based on all data within the
            % bin.

            if nGait==1
            searchPosIDX = {ppantData.trgO_gPcnt(uset),...
                ppantData.respO_gPcnt(uset)}; % targ onset as %% targ onset as %
                

            else %(Much messier, but due to overlap in events occurring across 2 gaits.)
                if iLR==1
                    searchPosIDX = {ppantData.trgO_gPcnt_LRL(uset),...% trg onset as %
                        ppantData.respO_gPcnt_LRL(uset)};  % resp O
                elseif iLR==2
                    searchPosIDX = {ppantData.trgO_gPcnt_RLR(uset),...% trg onset as %
                        ppantData.respO_gPcnt_RLR(uset)};  % resp O
                    
                elseif iLR==3
                    % both. Problem is they can have identical row info.
                    searchPosIDX = {[ppantData.trgO_gPcnt_LRL(Lcomb); ppantData.trgO_gPcnt_RLR(Rcomb)],...
                        [ppantData.respO_gPcnt_LRL(Lcomb); ppantData.respO_gPcnt_RLR(Rcomb)]};

                end
            end


            searchPosDVs= ppantData.SDTcat(uset); % for both dprime and criterion


            % %except in this particular case:
            if iLR==3 && nGait==2

                %and define the DV of interest:
                searchPosDVs = [ppantData.SDTcat(Lcomb); ppantData.SDTcat(Rcomb)];
                    
            end

            % accumulate all data within a relevant bin.

            
            datatobin = searchPosDVs;

            for iAlignto=1:2
                
                [dprime_bin, criterion_bin,...
                    HitRate_bin, FARate_bin]=deal(nan(1,length(pidx)-1));               
                
                allpos = searchPosIDX{iAlignto};
                
            for ibin=1:length(pidx)-1
                idx = pidx(ibin):pidx(ibin+1);
                % for all index points, collect the data we need.
                HITn_bin=[];
                MISSn_bin=[];
                FAn_bin=[];
                CRn_bin=[];
                
                for ipos = 1:length(idx)
                    searchidx= idx(ipos);
                    
                    %find all rows with this onset:
                    useb = find(allpos==searchidx);

                    HITn_bin(ipos) = length(find(datatobin(useb)==1));

                    MISSn_bin(ipos) = length(find(datatobin(useb)==2));

                    FAn_bin(ipos) = length(find(datatobin(useb)==3));

                    CRn_bin(ipos) = length(find(datatobin(useb)==4));

                end
                %total for this bin:
                Hratebin = sum(HITn_bin)/ (sum(HITn_bin) + sum(MISSn_bin));
                FAratebin = sum(FAn_bin)/ (sum(FAn_bin) + sum(CRn_bin));

                if FAratebin==0 % HACK! 
                    FAratebin=.001;
                end
                if Hratebin==0
                    Hratebin=.001;
                end

                if Hratebin==1
                    Hratebin=.999;
                end
                if FAratebin==1
                    FAratebin=.999;
                end

                %calculate dprime and criterion, per bin.
                dprime_bin(ibin) = norminv(Hratebin)- norminv(FAratebin);
                criterion_bin(ibin) =  -0.5*(norminv(Hratebin)+ norminv(FAratebin));

                HitRate_bin(ibin) = Hratebin;
                FARate_bin(ibin) = FAratebin;
                
                if isinf(dprime_bin(ibin))
                    disp('check')
                end
            end % bin
            
            if iAlignto==1
                targOns_dprime_bin = dprime_bin;
                targOns_crit_bin = criterion_bin;
                targOns_HR_bin = HitRate_bin;
                targOns_FA_bin= FARate_bin;
            else
                respOns_dprime_bin = dprime_bin;
                respOns_crit_bin = criterion_bin;
                respOns_HR_bin = HitRate_bin;
                respOns_FA_bin= FARate_bin;
            end
            
            end % iAlign
               %% store for all ppants, preparing to save:
                gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
                
                %add binned, shuffled vers:
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts_shuff']) = trgCounts_bin_shuff;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc_shuff']) = trgAcc_bin_shuff;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = trgRT_bin_shuff;
                
                %also using click pos as index:                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts_shuff']) = respCounts_bin_shuff;           
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = respRT_bin_shuff;
                
                
            end % iLR
        end % iSpeed
        end
 

    
    end % ppant
    
    %%
    cd([procdatadir filesep 'GFX']);
    
    disp('resaving GFX data with shuffle attached...');
    save('GFX_Data_inGaits', ...
      'GFX_TargPosData','GFX_RespPosData',...
        '-append');