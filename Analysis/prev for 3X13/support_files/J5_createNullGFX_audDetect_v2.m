% J5_createNullGFX_audDetect_v2
% here cleaning up the order of operations, to perform a series of
% computations per shuffle (rather than shuffling for every computation, as
% was done previously).in j5_createNullGFX_audDetect
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
       
        disp(['constructing null distributions for ' subjID ', ' num2str(ippant) '/' num2str(length(pfols))]);
 
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
                    Lcomb = intersect(find(useL), usespd); %LRL at a particular speed.
                    
                    Rcomb = intersect(find(useR), usespd); % RLR at a particular speed.
                    
                end
                %% restruct to relevant data (using index for trgO and respO)
                % which restricted dataset are we using?
                    if nGait==1
                        %IVs
                        tmp_targOnsets = ppantData.trgO_gPcnt(uset);
                        tmp_respOnsets = ppantData.respO_gPcnt(uset);
                         % Order of DVs, for each search IV above:
                         tmp_Correct = ppantData.correctResponse(uset);
                         tmp_clickRT = ppantData.clickRT(uset);
                         tmp_SDTcat = ppantData.SDTcat(uset);
                         
                    elseif nGait==2 && iLR<3
                        
                        if iLR==1
                            tmp_targOnsets = ppantData.trgO_gPcnt_LRL(uset); %LRL onset
                            tmp_respOnsets = ppantData.respO_gPcnt_LRL(uset);
                        elseif iLR==2
                            tmp_targOnsets = ppantData.trgO_gPcnt_RLR(uset); %RLR onset
                            tmp_respOnsets = ppantData.respO_gPcnt_RLR(uset);
                        end
                        % DVs unchanged:
                        % Order of DVs, for each search IV above:
                        tmp_Correct = ppantData.correctResponse(uset);
                        tmp_clickRT = ppantData.clickRT(uset);
                        tmp_SDTcat = ppantData.SDTcat(uset);
                        
                        
                    elseif nGait==2 && iLR==3 % special case, use all data (both columns stacked).
                        
                        tmp_targOnsets= [ppantData.trgO_gPcnt_LRL(Lcomb); ppantData.trgO_gPcnt_RLR(Rcomb)];                        
                        tmp_respOnsets= [ppantData.respO_gPcnt_LRL(Lcomb); ppantData.respO_gPcnt_RLR(Rcomb)];
                        
                        tmp_Correct = [ppantData.correctResponse(Lcomb);ppantData.correctResponse(Rcomb)];
                        tmp_clickRT = [ppantData.clickRT(Lcomb);ppantData.clickRT(Rcomb)];
                        tmp_SDTcat = [ppantData.SDTcat(Lcomb);ppantData.SDTcat(Rcomb)];
                        
                    end
                %% Here is the trick, for each gait x speed x LR combination, create the null dataset of 
                %new position indices.
                % one for target onset, one for response onset, everything
                % else flows from that.
                
                % pseudo-code:
                % for gait percent 1-100, calculate DVs
                % - acc, rt, dprime, crit etc.
                 
                % then calculate binnedversions
                
                 % for actual gait pos, select another at random:
                 
                 %preallocate:
                 [targOns_Acc_null, respOns_Acc_null,...
                     targOns_RTs_null, respOns_RTs_null,...
                     targOns_dprime_null, respOns_dprime_null,...
                     targOns_crit_null, respOns_crit_null,...
                     targOns_Counts_null, respOns_Counts_null]  =deal(nan(nPerm, pidx(end)));
                 
                 % for the binned dprime and criterion, we need to store
                 % the counts of H,M,DA,CR, at each indx, and shuffle:
                 [targOns_SDTcat_null,respOns_SDTcat_null] = deal(nan(nPerm,pidx(end),4)); % 1,2,3,4 for H,M,FA,CR
                 
                 
                 for iperm = 1:nPerm
                    
                     for ip=1:pidx(end)  
                    
                    
                    %instead of using ip, select at random.
                    useIDX = randi(pidx(end),1); 
                    
                    %find row numbers:
                    thisrand_targO= find(tmp_targOnsets==useIDX);
                    thisrand_respO= find(tmp_respOnsets==useIDX);
                        
                    %now find allocate data at useIDX at ip (1-100%).
                    % Accuracy:
                    targOns_Acc_null(iperm,ip) = nanmean(tmp_Correct(thisrand_targO));
                    respOns_Acc_null(iperm,ip) = nanmean(tmp_Correct(thisrand_respO));
                    
                    %RT:
                    targOns_RTs_null(iperm,ip) = nanmean(tmp_clickRT(thisrand_targO));
                    respOns_RTs_null(iperm,ip) = nanmean(tmp_clickRT(thisrand_respO));
                    
                    %dprime and criterion
                     
                    useOnsets= {thisrand_targO, thisrand_respO};
                                
                        for itargResp=1:2
                            onsetidx = useOnsets{itargResp};
                            %counts for calcuulations:
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
                                    targOns_dprime_null(iperm,ip) = dp;
                                    targOns_crit_null(iperm,ip) = crit;
                                    targOns_SDTcat_null(iperm,ip,1) = Hitn;
                                    targOns_SDTcat_null(iperm,ip,2) = Missn;
                                    targOns_SDTcat_null(iperm,ip,3) = FAn;
                                    targOns_SDTcat_null(iperm,ip,4) = CRn;
                                    
                                else
                                    respOns_dprime_null(iperm,ip) = dp;
                                    respOns_crit_null(iperm,ip) = crit;
                                    respOns_SDTcat_null(iperm,ip,1) = Hitn;
                                    respOns_SDTcat_null(iperm,ip,2) = Missn;
                                    respOns_SDTcat_null(iperm,ip,3) = FAn;
                                    respOns_SDTcat_null(iperm,ip,4) = CRn;
                                end
                                
                        end
                        
                        % store counts as well:
                        targOns_Counts_null(iperm, ip) = length(thisrand_targO);
                        respOns_Counts_null(iperm, ip) = length(thisrand_respO);
                        
                     end % ip (1-100)
                     %
                     if mod(iperm,100)==0
                         disp(['fin perm ' num2str(iperm) ' of ' num2str(nPerm)]);
                     end
                 end % perm (1-nPerm)
                
               
            
                %% %%%%%%%%%%%%%%%%%%%
                % Step through all shuffled data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%
               % we can bin using averaging if not the SDT summary
               
                
               binmedata = {targOns_Acc_null, respOns_Acc_null, ...
                   targOns_RTs_null, respOns_RTs_null,...
                   targOns_dprime_null,respOns_dprime_null,... % expect lots of nans for data paucity.
                   targOns_crit_null,respOns_crit_null,... % expect lots of nans for data paucity.
                    targOns_Counts_null, respOns_Counts_null};
                % note dprime and criterion were included for comparison
                % with the alternative below (calculating based on data
                % within bin, not averaging over bin indexes).
                
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
                           trgAcc_bin_shuff= outgoingbinned;
                       case 2
                           respAcc_bin_shuff= outgoingbinned;
                       case 3
                           trgRT_bin_shuff= outgoingbinned;
                       case 4
                           respRT_bin_shuff= outgoingbinned;
                       case 5
                           trgdP_bin_shuff = outgoingbinned;
                       case 6
                           respdP_bin_shuff = outgoingbinned;
                       case 7
                           trgCrit_bin_shuff = outgoingbinned;
                       case 8
                           respCrit_bin_shuff = outgoingbinned;
                       case 9
                           trgCounts_bin_shuff = outgoingbinned;
                       case 10
                           respCounts_bin_shuff = outgoingbinned;
                   end
                end % per type
                
                % 
                
                
                %% extra handling for binned SDT output:
                   binmedata = {targOns_SDTcat_null, respOns_SDTcat_null};
                   
                   for itype = 1:length(binmedata)
                       
                       usedata = binmedata{itype};
                       
                       
                       [dprime_bin,criterion_bin, HitRate_bin, FARate_bin]=...
                           deal(nan(nPerm,length(pidx)-1));
                       
                       
                       for ibin=1:length(pidx)-1
                           idx = pidx(ibin):pidx(ibin+1);
                           
                           HITn_bin = sum(usedata(:,idx,1),2); % sum counts over bin range
                           MISSn_bin = sum(usedata(:,idx,2),2); % sum counts over bin range
                           FAn_bin = sum(usedata(:,idx,3),2); % sum counts over bin range
                           CRn_bin = sum(usedata(:,idx,4),2); % sum counts over bin range
                           
                           %                     HITn_bin(ipos) = length(find(datatobin(useb)==1));
                           %
                           %                     MISSn_bin(ipos) = length(find(datatobin(useb)==2));
                           %
                           %                     FAn_bin(ipos) = length(find(datatobin(useb)==3));
                           %
                           %                     CRn_bin(ipos) = length(find(datatobin(useb)==4));
                           Hratebin = HITn_bin./ (HITn_bin + MISSn_bin);
                           FAratebin = FAn_bin./ (FAn_bin + CRn_bin);
                           
                           
                           %remove 0 and 1 floor and ceiling effects
                           %(stuffs calcs).
                               FAratebin(FAratebin==0)=.001;
                           
                               Hratebin(Hratebin==0)=.001;
                           
                               Hratebin(Hratebin==1)=.999;
                           
                               FAratebin(FAratebin==1)=.999;
                           
                           
                           %calculate dprime and criterion, per bin.
                           dprime_bin(:,ibin) = norminv(Hratebin)- norminv(FAratebin);
                           criterion_bin(:,ibin) =  -0.5*(norminv(Hratebin)+ norminv(FAratebin));
                           
                           HitRate_bin(:,ibin) = Hratebin;
                           FARate_bin(:,ibin) = FAratebin;
                           
                           
                           if isinf(dprime_bin(ibin))
                               disp('check')
                           end
                           
                       end % for ibin
                       if itype==1
                           
                           targOns_dprime_bin = dprime_bin;
                           targOns_criterion_bin =criterion_bin;
                           
                           targOns_HitRate_bin = HitRate_bin;
                           targOns_FARate_bin= FARate_bin;
                           
                       elseif itype==2
                           
                           respOns_dprime_bin = dprime_bin;
                           respOns_criterion_bin =criterion_bin;
                           
                           respOns_HitRate_bin = HitRate_bin;
                           respOns_FARate_bin= FARate_bin;
                           
                       end
                       
                       
                   end % both onset types
                   
                
           
               %% store for all ppants, preparing to save:
                gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
                
                %add binned, shuffled vers to existing data structure
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts_shuff']) = trgCounts_bin_shuff;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc_shuff']) = trgAcc_bin_shuff;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = trgRT_bin_shuff;
               %we have 2 versions of the SDT output, saving ver 2 for now:
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime_shuff']) = targOns_dprime_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit_shuff']) = targOns_criterion_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR_shuff']) = targOns_HitRate_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA_shuff']) = targOns_FARate_bin;
                
                %also using click pos as index:                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts_shuff']) = respCounts_bin_shuff;           
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc_shuff']) = respAcc_bin_shuff;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = respRT_bin_shuff;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime_shuff']) = respOns_dprime_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit_shuff']) = respOns_criterion_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR_shuff']) = respOns_HitRate_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA_shuff']) = respOns_FARate_bin;
                
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