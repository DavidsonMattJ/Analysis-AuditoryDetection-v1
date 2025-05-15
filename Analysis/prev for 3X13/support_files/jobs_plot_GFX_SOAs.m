% jobs_plot_GFX_SOAs

% -- called from JOBS_import_preprocess;




%////  This script is old, now deprecated. Use the gauss instead.

%% PLOT GFX SOAs BEGINS HERE

% this script loads the presaved participant level data (SOA variables).
% makes a plot per participant, comparing slow and fast walk conditions.


% - load data
% - extract trials per SOA type
% - calculate proportion same
% save
cd(savedatadir);

% groupdata = dir([pwd filesep 'p_*.mat']);


%%
% load group level data:

load('GFX_SOAdata.mat');



[nppants, nGatePhases, nSOAs] = size(GFX_Data_SOAs.nThisSOAthisGaitPhase_all,2);

% set up a figure

%%
clf;
% 
GFX_propSame_slow= GFX_Data_SOAs.propSame_slow;
ppantMax= max(GFX_propSame_slow,[],2);
ppantMax_forDiv = repmat(ppantMax, [1,7])
normData_slow = GFX_propSame_slow ./ppantMax_forDiv;



GFX_propSame_fast= GFX_Data_SOAs.propSame_fast;
ppantMax= max(GFX_propSame_fast,[],2);
ppantMax_forDiv = repmat(ppantMax, [1,7])
normData_fast = GFX_propSame_fast ./ppantMax_forDiv;


%%  group averages:
figure(1); clf;
pM= mean(GFX_propSame_slow,1);
pS= plot(SOAs, pM, 'bo-', 'MarkerFaceColor', 'b'); hold on;
% add individual points?
hold on;
% plot(SOAs, GFX_propSame_slow, '.b');
% 
stE= CousineauSEM(GFX_propSame_slow);
% stE= std(GFX_propSame_slow,0,1) ./ sqrt(size(GFX_propSame_slow,1));
hold on
errorbar(SOAs, pM, stE, 'LineStyle','none','color', 'b');
%
 pM=mean(GFX_propSame_fast,1);

stE= CousineauSEM(GFX_propSame_fast);
% stE= std(GFX_propSame_fast,0,1) ./ sqrt(size(GFX_propSame_fast,1));

pF=plot(SOAs, pM, 'r-o'); hold on;
hold on;
errorbar(SOAs, pM, stE, 'LineStyle','none', 'color', 'r');
title('Group Level Mean Proportion of "SAME" Responses')
set(gca, 'XTicklabel', SOAs);
ylabel('Proportion of "SAME" Responses');
xlabel('Stimulus Onset Asynchrony(Seconds)')
legend([pS, pF,], {'slow walk', 'normal walk'});


% can probably just fit a gussian at this stage:

%% normalized averages:
figure(2)

pM= mean(normData_slow,1);
plot(pM, 'b-o'); hold on;

stE= std(normData_slow,0,1) ./ sqrt(size(normData_slow,1));
hold on
errorbar(1:7, pM, stE, 'LineStyle','none','color', 'b');
%
pM=mean(normData_fast,1);
stE= std(normData_fast,0,1) ./ sqrt(size(normData_fast,1));

plot(pM, 'r-o'); hold on;
hold on;
errorbar(1:7, pM, stE, 'LineStyle','none', 'color', 'r');
title('Group Level Mean Proportion of "SAME" Responses (Normalized)')

set(gca, 'XTicklabel', SOAs);
ylabel('Proportion of "SAME" Responses');
xlabel('Stimulus Onset Asynchrony(Seconds)');

%% Take Mean of Prop Same across group, all walk Speeds, per SOA

nSOAs = length(SOAs);
GroupMeanSame_all = nan(2,7);

GFX_propSame_all= GFX_Data_SOAs.propSame_all;

for iSOA = 1:nSOAs

    GroupMeanSame_allThisSOA = mean(GFX_propSame_all(:,iSOA));
    
    GroupMeanSame_all(2,iSOA) = GroupMeanSame_allThisSOA;
    GroupMeanSame_all(1,iSOA) = SOAs(iSOA);

end




%% Take Mean of Prop Same across group, Slow Walk, per SOA

nSOAs = length(SOAs);
GroupMeanSame_slow  = nan(2,7);

for iSOA = 1:nSOAs

    GroupMeanSame_slowThisSOA = mean(GFX_propSame_slow(:,iSOA));
    
    GroupMeanSame_slow(2,iSOA) = GroupMeanSame_slowThisSOA;
    GroupMeanSame_slow(1,iSOA) = SOAs(iSOA);

end

%% Take Mean of Prop Same across group, Fast Walk, per SOA

nSOAs = length(SOAs);
GroupMeanSame_fast  = nan(2,nSOAs); % mean, std dev

for iSOA = 1:nSOAs

    GroupMeanSame_fastThisSOA = mean(GFX_propSame_fast(:,iSOA));
    
    GroupMeanSame_fast(2,iSOA) = GroupMeanSame_fastThisSOA;
    GroupMeanSame_fast(1,iSOA) = SOAs(iSOA);

end


%% Take Group Mean RTfromStim1, across all speeds


Speeddiff = GroupMeanSame_fast - GroupMeanSame_slow;

%% Mean of in-Participant change in PropSame by Walk Speed

GroupMeanChangebyWalkSpeed = nan(1,nSOAs);

for iSOA = 1:nSOAs
    GroupMeanChangebyWalkSpeedThisSOA = mean(GFX_ChangeAcrossSpeeds(:,iSOA));

    GroupMeanChangebyWalkSpeed(1,iSOA) = GroupMeanChangebyWalkSpeedThisSOA

end


%% Mean of Abs of in-Participant change in PropSame by Walk Speed

GroupMeanAbsChangebyWalkSpeed = nan(1,nSOAs);

for iSOA = 1:nSOAs
    GroupMeanAbsChangebyWalkSpeedThisSOA = mean(GFX_AbsChangeAcrossSpeeds(:,iSOA));

    GroupMeanAbsChangebyWalkSpeed(1,iSOA) = GroupMeanAbsChangebyWalkSpeedThisSOA

end


%% Means of Response times
FocusVariables = {'GroupMeanRT1_all' , 'GroupMeanRT1_slow', 'GroupMeanRT1_fast', 'GroupMeanRT2_all' , 'GroupMeanRT2_slow', 'GroupMeanRT2_fast' }

nrowsGroupRt = length(FocusVariables);

%GroupRTData = 

GroupRTData = nan(nrowsGroupRt,nSOAs)
nVariables = length(FocusVariables);
VariableLabels = [nVariables,1];
