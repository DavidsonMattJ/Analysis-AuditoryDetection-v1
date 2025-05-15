% j4_crunchPFX_concatGFX


% Here we load all data from previous analyses, preparing to save before
% plotting in next round of jobs.

%%%%%% Auditory detection v1.



cd(savedatadir);
pfols= dir([pwd filesep 'p_*']);
procdatadir= savedatadir; % for compatability below.

%% concat data:
%preallocate storage:

pidx1=ceil(linspace(1,100,21)); % length n-1
pidx2=ceil(linspace(1,100,41));%
useindexes = {pidx1, pidx2};
gaittypes = {'single gait' , 'double gait'};
%

GFX_headY=[];

GFX_TargPosData=[];
GFX_RespPosData=[];
subjIDs=cell([1, length(pfols)]);

avWalkParams=[]; % added post hoc to calculate average stride length (duration sec).

for ippant =1:length(pfols)
    cd(procdatadir)
    load(pfols(ippant).name, 'subjID', 'summary_table',...
        'gait_ts_gData', 'gait_ts_resamp',...
        'doubgait_ts_resamp');
    
    
    
    % should have been done earlier, but make a column of clickRT for
    % compatability below:
    summary_table.clickRT = summary_table.targRT - summary_table.targOnset;
    
    subjIDs{ippant} = subjID;
    
    disp(['concatenating subject ' subjID]);
    
    % first retrieve index for separate gaits (L/Right feet).
    Ltrials= strcmp(summary_table.trgO_gFoot, 'LR');
    Rtrials= strcmp(summary_table.trgO_gFoot, 'RL');
    
    % mean head pos (per speed):
    slowsteps = find(gait_ts_gData.walkSpeed==1);
    naturalsteps = find(gait_ts_gData.walkSpeed==2);
    GFX_headY(ippant,1).gc = mean(gait_ts_resamp(slowsteps,:),1, 'omitnan');
    GFX_headY(ippant,2).gc = mean(gait_ts_resamp(naturalsteps,:),1, 'omitnan');
    GFX_headY(ippant,3).gc = mean(gait_ts_resamp,1, 'omitnan'); % all steps
    
    GFX_headY(ippant,1).doubgc = mean(doubgait_ts_resamp(slowsteps,:),1, 'omitnan');
    GFX_headY(ippant,2).doubgc = mean(doubgait_ts_resamp(naturalsteps,:),1, 'omitnan');
    GFX_headY(ippant,3).doubgc = mean(doubgait_ts_resamp,1, 'omitnan'); % all steps
    
    Ltrials_ts= strcmp(gait_ts_gData.gaitFeet, 'LR');
    Rtrials_ts= strcmp(gait_ts_gData.gaitFeet, 'RL');
    avWalkParams(ippant)= mean(gait_ts_gData.gaitDuration(Ltrials_ts)) + mean(gait_ts_gData.gaitDuration(Rtrials_ts));
    
    
    %% also calculate binned versions
    % These take the mean over an index range, per gait cycle position
    
    for nGait=1:2  % single step or stride based analysis.
        
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
        ppantData= summary_table;
        
        if nGait==1
            useL=Ltrials;
            useR=Rtrials;
            useAll = 1:length(Ltrials);
            
            % split by foot (LR, RL, all combined)
        else
            
            
            %%% % % % ! ! ! ! ! ! ! UNFINISHED
            % note that here we have data in separate columns (so find non-nan
            % rows.)
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
        
        % different combinations to use for indexing:
        for iSpeed=1:3        % third is combined.
            for iLR=1:3 %L, R, combined
                
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
                % Step through different data types (DVs and event positions),
                % and analyses:
                
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
                
                
                % to add:  FA, omissions etc.?
                
                %%%%%%%%%%%%%%%%%%
                % for each of these analysis  define the onsets we will aggregate:
                
                %i.e.  target onsets or response onsets within a gait?
                
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
                        
                %% loop through the index analysis:
                for iAlignto= 1:2 % targOnset or respOnset.
                    
                    for iAnalysis=1:4
                        %preallocate:
                        outgoing = nan(1,length(pidx)); %
                        allpos = searchPosIDX{iAlignto};
                        allDVs = searchPosDVs{iAnalysis};
                        
                        
                        
                        
                        % for actual gait pos, select reldata
                        for ip=1:pidx(end)
                            useb = find(allpos==ip);
                            
                            %% Compute accuracy.
                            if iAnalysis==1
                                
                                
                                % remove nans from length calcs.
                                tmpResp = allDVs(useb);
                                useResp = tmpResp(~isnan(tmpResp));
                                outgoing(ip)= sum(useResp, 'omitnan') / length(useResp);
                                
                                
                                %% take mean for these  trials (if RT):
                            elseif iAnalysis==2
                                
                                outgoing(ip) = mean(allDVs(useb),1, 'omitnan');
                                
                            elseif iAnalysis ==3
                                %compute sensitivity
                                %sensitivity from HR and FAR.
                                %HR is HR/ H+M;
                                % FA is FA/ FA+CR;
                                tmpResp= allDVs(useb);
                                
                                Hitn = length(find(tmpResp==1));
                                Missn = length(find(tmpResp==2));
                                
                                FAn = length(find(tmpResp==3));
                                CRn = length(find(tmpResp==4));
                                
                                HitRate = Hitn / (Hitn+Missn);
                                FArate = FAn / (FAn + CRn);
                                dp = norminv(HitRate)- norminv(FArate);
                                outgoing(ip) = dp;
                                outgoing_SDT(1,ip)= HitRate;
                                outgoing_SDT(2,ip)= FArate;
                                
                            elseif iAnalysis==5 % compute criterion
                                tmpResp= allDVs(useb);
                                
                                Hitn = length(find(tmpResp==1));
                                Missn = length(find(tmpResp==2));
                                
                                FAn = length(find(tmpResp==3));
                                CRn = length(find(tmpResp==4));
                                
                                HitRate = Hitn / (Hitn+Missn);
                                FArate = FAn / (FAn + CRn);
                                
                                crit = -0.5*(norminv(HitRate)+ norminv(FArate));
                                outgoing(ip) = crit;
                                
                                
                                
                            end % ianalysis
                            
                        end % pidx loop
                        
                        % now rename
                        switch iAnalysis
                            case 1
                                outgoing_Acc= outgoing;
                            case 2
                                outgoing_RT= outgoing;
                            case 3
                                outgoing_dprime= outgoing;
                                
                            case 4
                                outgoing_crit = outgoing;
                        end
                        
                    end  % ianalysis
                    
                    if iAlignto==1 % rename once more:
                        
                        targOns_Acc= outgoing_Acc;
                        targOns_RTs= outgoing_RT;
                        targOns_dprime = outgoing_dprime;
                        targOns_crit = outgoing_crit;
%                         targOns_HR = outgoing_SDT(1,:);
%                         targOns_FAR = outgoing_SDT(2,:);
                        
                    elseif iAlignto==2
                        
                        respOns_Acc= outgoing_Acc;
                        respOns_RTs= outgoing_RT;
                        respOns_dprime = outgoing_dprime;
                        respOns_crit = outgoing_crit;
%                         respOns_HR = outgoing_SDT(1,:);
%                         respOns_FAR = outgoing_SDT(2,:);
                        
                        
                    end
                    
                end %iAlign to.
                
                
                
                
                
                %% store histcounts, per target onset pos, and resp onset pos.
                tPos = ppantData.trgO_gPcnt(uset);
                clkPos= ppantData.respO_gPcnt(uset); % gait pcnt for response onset.
                
                % % correct for special case:
                if nGait==2
                    if iLR==1
                        tPos = ppantData.trgO_gPcnt_LRL(uset);
                        clkPos= ppantData.respO_gPcnt_LRL(uset);
                    elseif iLR==2
                        tPos = ppantData.trgO_gPcnt_RLR(uset);
                        clkPos= ppantData.respO_gPcnt_RLR(uset);
                    elseif iLR==3
                        tPos = [ppantData.trgO_gPcnt_LRL(uset),ppantData.trgO_gPcnt_RLR(uset)];
                        clkPos= [ppantData.respO_gPcnt_LRL(useL);ppantData.respO_gPcnt_RLR(useR)];
                    end
                    
                end
                % % % !!
                
                usetypes = {tPos, clkPos};
                for iAnalysis=1:2
                    datain = usetypes{iAnalysis};
                    % remove zero pos points before counting:
                    datain(datain==0)=[];
                    datain(isnan(datain))=[];
                    % raw hit counts, per gaitsample:
                    tmp_counts = histcounts(datain, pidx(end));
                    
                    %critical step -convert 0 to nans. this avoids averaging
                    %problems.
                    tmp= tmp_counts; % trgO_counts, / clkOcounts saved below.
                    tmp(tmp==0)=NaN;
                    tmptoAvg=tmp;
                    
                    % rename
                    if iAnalysis==1
                        trgO_counts = tmp_counts;
                        trgOtoAvg = tmptoAvg; % with nans removed.
                    elseif iAnalysis==2
                        respO_counts = tmp_counts;
                        respOtoAvg = tmptoAvg;
                    end
                end
                %% %%%%%%%%%%%%%%%%%%%
                % Step through all data types, performing binning based on
                % length of pidx.
                %%%%%%%%%%%%%%%%%%%
                %trg onset counts
                %resp onset counts
                %acc per trg onset
                %rt per trg onset
                %acc per resp onset
                %rt per resp onset
                % add:
                
                usedata= {trgOtoAvg, respOtoAvg, ... % using ver with NaNs removed
                    targOns_Acc,targOns_RTs, ...                    
                    respOns_Acc,respOns_RTs};
                
                for iAnalysis= 1:length(usedata)
                    
                    %% for all types, take mean over bin indices.
                    out_bin=nan(1,length(pidx)-1);
                    datatobin = usedata{iAnalysis};
                    for ibin=1:length(pidx)-1
                        idx = pidx(ibin):pidx(ibin+1);
                        
                        out_bin(ibin) = mean(datatobin(idx), 'omitnan');
                        
                    end
                    
                    switch iAnalysis
                        case 1
                            trgCounts_bin = out_bin;
                        case 2
                            respCounts_bin = out_bin;
                        case 3
                            targOns_Acc_bin = out_bin;
                            
                        case 4
                            targOns_RT_bin= out_bin;
                        case 5
                            respOns_Acc_bin = out_bin;
                        case 6
                            respOns_RT_bin = out_bin;
                    end
                    
                end
                
                
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
                %% %%%%%%%%%%%%%%%%%%%%%%%%% now save:
                gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
                
                %add targ data to one matrix, resp data to another:
                %raw (no binning)
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'counts']) = trgO_counts;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'Acc']) = targOns_Acc;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'rts']) = targOns_RTs;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'dprime']) =targOns_dprime;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'crit']) =targOns_crit;
                
                %binned version:
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts']) = trgCounts_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc']) = targOns_Acc_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts']) = targOns_RT_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime']) = targOns_dprime_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit']) = targOns_crit_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR']) = targOns_HR_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA']) = targOns_FA_bin;
                
                %also using response position as index:
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'counts']) = respO_counts;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'Acc']) = respOns_Acc;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'rts']) = respOns_RTs;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'dprime']) =respOns_dprime;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'crit']) =respOns_crit;
                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts']) = respCounts_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc']) = respOns_Acc_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts']) = respOns_RT_bin;                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime']) = respOns_dprime_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit']) = respOns_crit_bin;               
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR']) = respOns_HR_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA']) = respOns_FA_bin;
                
                
            end % iLR
        end % ispeed
    end % nGaits.
    
    % % after all said and done, concatenate the two single gaits to make a
    % % stride (for plots).
    % % store as a new double gc (i.e. a stride).
    %
    % useD = {GFX_TargPosData, GFX_RespPosData};
    % %store field names, but only first time through, otherwise we will append
    % %too many times.
    % if ippant==1
    % allfieldsTarg = fieldnames(GFX_TargPosData);
    % allfieldsResp= fieldnames(GFX_RespPosData);
    % end
    %
    % for idatatype = 1:2
    %
    %     if idatatype==1
    %         allfields= allfieldsTarg;
    %     else
    %         allfields= allfieldsResp;
    %     end
    %
    %
    %     for ifield = 1:length(allfields)
    %         % per field, store with new field name prefix (doubgc_) and
    %         % horizonrally concatenate the data
    %         if ifield==4
    %             pause(.1)
    %         end
    %         useD{idatatype}(ippant,4).(['doub' allfields{ifield}]) = ...
    %             [useD{idatatype}(ippant,1).([allfields{ifield}]), ... %LR
    %             useD{idatatype}(ippant,2).([allfields{ifield}])]; %RL
    %     end
    %
    %
    %
    % end
    %  %rename to save  with new info
    % GFX_TargPosData = useD{1};
    % GFX_RespPosData= useD{2};
    %
end % ppant
%% save GFX
try cd([procdatadir filesep 'GFX']);
catch
    mkdir([procdatadir filesep 'GFX']);
cd([procdatadir filesep 'GFX']);
end
    disp('Saving GFX');
save('GFX_Data_inGaits', ...
    'GFX_headY', 'GFX_TargPosData','GFX_RespPosData',...
    'subjIDs', 'pidx1', 'pidx2', 'gaittypes');%, '-append');

