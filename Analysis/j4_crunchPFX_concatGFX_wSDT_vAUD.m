% j4_crunchPFX_concatGFX_wSDT_vAUD

% Here we load each participants data table, and perform a series of for each index
% in either the step, stride, or other breakdowns (speed, left vs right, etc). 

%%%%%% Auditory discrimination project %%%%%%



set_myOnlineDirectories_AUD
% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%


%% concat data:
%preallocate storage:
% how many bins?
% pidx1=ceil(linspace(1,100,21)); % length n-1
% pidx2= ceil(linspace(1,100,41));%

pidx1=ceil(linspace(1,100,11)); % length n-1 (so 10 for a single step)
pidx2= ceil(linspace(1,100,21));% (now 20 for a stride)

useindexes = {pidx1, pidx2};

gaittypes = {'single gait' , 'double gait'};
%

GFX_headY=[];
GFX_TargPosData=[];
GFX_RespPosData=[];
GFX_TargParam=[];
subjIDs=cell([1, length(pfols)]);
%%
for ippant =1:nsubs
    cd(procdatadir)
    load(pfols(ippant).name, 'subjID', 'trial_summaryTable',...
        'gait_ts_gData', 'gait_ts_resamp', 'doubgait_ts_resamp');
%%

    subjIDs{ippant} = subjID;

    disp(['concatenating subject ' subjID]);

    % first retrieve index for separate gaits (L/Right feet).
    %logical indexing to be used below:
    %only applies to single steps:
    Ltrials= strcmp(trial_summaryTable.trgO_gFoot, 'LR');
    Rtrials= strcmp(trial_summaryTable.trgO_gFoot, 'RL');

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
    
    ppantData= trial_summaryTable;
    
        
    % Include an index for different walk speeds.
    % 
    slowspeedTrials = find(ppantData.walkSpeed==1);
    naturalspeedTrials = find(ppantData.walkSpeed==2);
    % all trials will be their combination.
    allspd = [slowspeedTrials; naturalspeedTrials];

    speedstoIndex={slowspeedTrials, naturalspeedTrials,sort(unique(allspd)) };



    %% Now cycle through the table:
    %pseudo code:
    % 1) restrict eligible events (rows in table), to some combination of:
    % -  gait(step vs stride)
    % - speed (slow vs normal vs combined)
    % - stance (left support, right support, combined)

    % 2) with this restricted, perfrom an order of operations for all
    % target/response onsets at each position in the 1-100 index over
    % step/stride.
    % - counts (target onset/response onset)
    % - accuracy
    % - RT
        
for    nGait=1:2 % single or dual steps.

    pidx=useindexes{nGait}; % indices for binning (single or double gc).
  
    %restrict to which rows based on single or dual steps?
    if nGait==1
        useL=find(Ltrials);
        useR=find(Rtrials);
    else
        useL =find(~isnan(ppantData.trgO_gPcnt_LRL));
        useR =find(~isnan(ppantData.trgO_gPcnt_RLR));       
    end
        % also the combined case:
        useAll= unique([useL; useR]);
        
        %store for later use:
       trialstoIndex ={useL, useR, useAll};
        
        
        % debugging: 
        pc=1; clf;

%% now cycle through other IV combinations and restrict data range as appropriate:
        for iSpeed=1:3     %slow, normal, combined   
            for iLR=1:3  % left, right, combined
           
                useft= trialstoIndex{iLR};
                usespd = speedstoIndex{iSpeed};
                uset = intersect(useft,usespd); % trials that satisfy all conditions.
        
                % note that if we have the double gc and combined case.
                % also need the following intersection:
                if iLR==3 && nGait==2
                    Lcomb = intersect(useL, usespd); %LRL at a particular speed.                    
                    Rcomb = intersect(useR, usespd); % RLR at a particular speed.                    
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

                %% for gait percent 1-100, calculate DVs
                % - acc, rt, dprime, crit etc.
                 
                % then calculate binnedversions
                
                 %preallocate:
                 [targOns_Acc, respOns_Acc,...
                     targOns_RTs, respOns_RTs,...
                     targOns_dprime, respOns_dprime,...
                     targOns_crit, respOns_crit,...
                     targOns_Counts, respOns_Counts]  =deal(nan(1, pidx(end)));
                 
                 % for the binned dprime and criterion, we need to store
                 % the counts of H,M,DA,CR, at each indx
                 [ targOns_SDTcat,respOns_SDTcat] = deal(nan(1,pidx(end),4)); % 1,2,3,4 for H,M,FA,CR
                 
                 % %debug plotting, show the histograms for targ counts as we
                 % %go.
                 % 
                 % subplot(3,3, pc);
                 % histogram(tmp_targOnsets,100);
                 % title(['iLR(' num2str(iLR) '), nG(' ...
                 %     num2str(nGait) '), n=' num2str(length(tmp_targOnsets)) ]);
                 % pc=pc+1;
                 % shg

                     for ip=1:pidx(end)  
                    
                        useIDX = ip; %1-100?
                    
                    %find row numbers:
                    thisIDX_targO= find(tmp_targOnsets==useIDX);
                    thisIDX_respO= find(tmp_respOnsets==useIDX);
                   
                    
                    %now find allocate data at useIDX at ip (1-100%).
                    % Accuracy:
                    targOns_Acc(1,ip) = nanmean(tmp_Correct(thisIDX_targO));
                    respOns_Acc(1,ip) = nanmean(tmp_Correct(thisIDX_respO));
                    
                    %RT:
                    targOns_RTs(1,ip) = nanmean(tmp_clickRT(thisIDX_targO));
                    respOns_RTs(1,ip) = nanmean(tmp_clickRT(thisIDX_respO));
                    
                    %dprime and criterion
                     
                    useOnsets= {thisIDX_targO, thisIDX_respO};
                                
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
                                    targOns_dprime(1,ip) = dp;
                                    targOns_crit(1,ip) = crit;
                                    targOns_SDTcat(1,ip,1) = Hitn;
                                    targOns_SDTcat(1,ip,2) = Missn;
                                    targOns_SDTcat(1,ip,3) = FAn;
                                    targOns_SDTcat(1,ip,4) = CRn;
                                    
                                else
                                    respOns_dprime(1,ip) = dp;
                                    respOns_crit(1,ip) = crit;
                                    respOns_SDTcat(1,ip,1) = Hitn;
                                    respOns_SDTcat(1,ip,2) = Missn;
                                    respOns_SDTcat(1,ip,3) = FAn;
                                    respOns_SDTcat(1,ip,4) = CRn;
                                end
                                
                        end
                        
                        % store counts as well:
                        targOns_Counts(1, ip) = length(thisIDX_targO);
                        respOns_Counts(1, ip) = length(thisIDX_respO);
                        
                     end % ip (1-100)
                
            
                %% %%%%%%%%%%%%%%%%%%%
                % Step through all  data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%
               % This first version simply averages the data within a bin.
%               
                % this will be compared to the specific SDT case of
                % recalculating SDT metrics based on raw data in the bin
                % (to avoid floor/ceiling effects for example).
               
                
               BinMeData = {targOns_Acc, respOns_Acc, ...
                   targOns_RTs, respOns_RTs,...
                   targOns_dprime,respOns_dprime,... % expect lots of nans for data paucity.
                   targOns_crit,respOns_crit,... % expect lots of nans for data paucity.
                    targOns_Counts, respOns_Counts};

                % note dprime and criterion were included for comparison
                % with the alternative below (calculating based on data
                % within bin, not averaging over bin indexes).
                
                for itype = 1:length(BinMeData)
                   
                   usedata = BinMeData{itype};                   
                   outgoingbinned=nan(1, length(pidx)-1);
                   
                   for ibin=1:length(pidx)-1
                       idx = pidx(ibin):pidx(ibin+1);
                       
                       outgoingbinned(1,ibin) = nanmean(usedata(1,idx),2);
                       
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
                
                % 
                
                
                %% extra handling for binned SDT output:
                   BinMeData = {targOns_SDTcat, respOns_SDTcat};
                   
                   for itype = 1:2% targ ons, resp ons
                       
                       usedata = BinMeData{itype};
                       
                       
                       [dprime_bin,criterion_bin, HitRate_bin, FARate_bin]=...
                           deal(nan(1,length(pidx)-1));
                       
                       
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
                
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'counts']) = targOns_Counts;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'Acc']) = targOns_Acc;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'rts']) = targOns_RTs;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'dprime']) =targOns_dprime;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'crit']) =targOns_crit;
                
                %binned version:
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts']) = targOns_Counts_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc']) = targOns_Acc_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts']) = targOns_RT_bin;
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime2']) = targOns_dP_bin; % averaged within bin, '2' since the second preference
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit2']) = targOns_Crit_bin; % averaged within bin
                                
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime']) = targOns_dprime_bin; % calculated for bin
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit']) = targOns_criterion_bin; % calculated for bin               
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR']) = targOns_HitRate_bin; % calculated for bin
                GFX_TargPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA']) = targOns_FARate_bin; % calculated for bin
                
                %also using response position as index:
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'counts']) = respOns_Counts;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'Acc']) = respOns_Acc;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'rts']) = respOns_RTs;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'dprime']) =respOns_dprime;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'crit']) =respOns_crit;
                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_counts']) = respOns_Counts_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_Acc']) = respOns_Acc_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_rts']) = respOns_RT_bin;                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime2']) = respOns_dP_bin; % averaged within bin, '2' since the second preference
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit2']) = respOns_Crit_bin; % averaged within bin
                
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_dprime']) = respOns_dprime_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_crit']) = respOns_criterion_bin;               
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_HR']) = respOns_HitRate_bin;
                GFX_RespPosData(ippant,iSpeed,iLR).([gaitnames{nGait} 'binned_FA']) = respOns_FARate_bin;
                
                
            end % iLR
        end % iSpeed
        end %nGait
 

    
    end % ppant
    
%% save GFX
cd([procdatadir filesep 'GFX']);
disp('Saving GFX');
save('GFX_Data_inGaits', ...
    'GFX_headY', 'GFX_TargPosData','GFX_RespPosData',...
    'GFX_TargParam',...
    'subjIDs', 'pidx1', 'pidx2', 'gaittypes','-append');
%%
% %% plot Hit rate per orientation:
% gfx_counts=[];
% gfx_HR=[];
% for ippant= 1:size(GFX_TargParam,2)
% 
% gfx_counts(ippant,:) = GFX_TargParam(ippant).orientation_Histcounts_binned;
% gfx_HR(ippant,:) = GFX_TargParam(ippant).orientation_HR_binned;
% end
% 
% %%
% clf;
% 
% mB = nanmean(gfx_counts,1);
% stE= CousineauSEM(gfx_counts);
% 
% clf;
% subplot(211)
% bar(1:length(mB), mB); hold on;
% errorbar(1:length(mB), mB, stE,'linestyle','none','color','k');
% ylabel('target counts');
% xlabel('Orientation degrees (1-180 binned)');
% ylim([0 7])
% % set(gca,'XTickLabel',{'1-10','11-20','21-30','31-40','41-50','51-60','61-70','71'})
% %
% hold on;
% yyaxis right
% % subplot(212)
% mHR = nanmean(gfx_HR);
% errHR = CousineauSEM(gfx_HR);
% 
% sh=shadedErrorBar(1:length(mHR), mHR, errHR,{'color','r'},1)
% set(gca,'YColor','r');
% ylabel('Hit Rate');
% ylim([0 .6])

