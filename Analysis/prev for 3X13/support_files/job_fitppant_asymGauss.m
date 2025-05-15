% job_fitppant_asymGauss
% this script loads individual subject data and fits a gaussian to the
% synchrony function for various conditions. 

% called from the JOBS_import_preprocess.m index script.
% test the gaussian fit parameter recovery.




% savedatadir= string(rootPath) + '/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

allppants = dir([pwd filesep 'p_*.mat']);
%%
%%
% Number of slices of Gait Cycle - Quantiles set at 4 in head
% tracking script

job=[];
job.fitper_participant=1;

job.plotPPantFitts=0;

% >>>>>
% now fitting an asymmetric gaussian to obtain width of LHS and RHS:


NLIN= 1; % 1 for NLINFIT; 0 or FMINSEARCH; 
              % Note: Only NLINFIT gives CIs and takes Weights

% e.g.:
% 
% startingVals = [ 0  50  100   1 ]; % .1];
% if NLIN
%     [ estimates,  resid, jacob, covarEst mse ] = nlinfit(SOAs, pSync,  @trySkewedGaussLR_nlin,startingVals );
%     ci = nlparci(estimates, resid, 'jacobian' , jacob); % calculate 95% confidence limits around estimates.
% else
%     options = optimset; 
%     estimates = fminsearch(@trySkewedGaussLR,startingVals,options,SOAs ,pSync);    
% end
% estimates = round(estimates*1000)/1000;
% 
% a = estimates(1); % Mean
% b = estimates(2); % SD of left half
% c = estimates(3); % SD of right half



figure(1); clf;
set(gcf, 'units', 'normalized', 'position', [.05 .05 .8 .8])



if job.fitper_participant==1
for ippant= 1:length(allppants)

    
    
% default fit params:    

cd(savedatadir);
load(allppants(ippant).name, 'ppantData_SOAs');

%store all the data types we will fit in a matrix , to easily cycle
%through
% note that all synch functions have the same number of data points: (7)

myGaussMatrix_bySpeed=[];
myGaussMatrix_bySpeed(1,:)= ppantData_SOAs.propSame_all;
myGaussMatrix_bySpeed(2,:)= ppantData_SOAs.propSame_slow;
myGaussMatrix_bySpeed(3,:)= ppantData_SOAs.propSame_fast;


% keep track of data in each location:
fitsTitle={'pSame_all', 'pSame_slow', 'pSame_fast'};

% next bit is confusing, so will make a matrix with different dims per relevant
% condition.
nsubs= length(allppants);

[nGaits, nSOAs]= size(ppantData_SOAs.nThisSOAthisGaitPhase_all);



speedFlds= {'_all', '_slow', '_fast'};
alignFlds= {'_FirstAligned', '_VisAligned', '_AudAligned', '_LastAligned', '_RespAligned'};
nspeeds= 3;
nalignments= length(alignFlds);
SOAs = ppantData_SOAs.SOAs;
SOAs_ms= SOAs*1000';
% populate with the relevant data:
myGaussMatrix_bySpeedbyResp= zeros(nspeeds, nalignments, nGaits, nSOAs);

for ispeed= 1:nspeeds
   for ialign=1:nalignments

       datain = ppantData_SOAs.(['propSameByGait' speedFlds{ispeed} alignFlds{ialign}]);
       
%speed x align x cond 
    myGaussMatrix_bySpeedbyResp(ispeed, ialign, :,:)= datain;
   end
end


% prepare output structs
participant_asymGauss_Fits_bySpeed=[];
participant_asymGauss_Fits_bySpeedbyAlignbyGait=[];

for igspeed= 1:size(myGaussMatrix_bySpeed,1) % 

%pseudo code: interpolate the observed data to 1000 points, to improve the
%fit. 
%fit a gaussian using the above model, retain the parameters within this
%loop and save.

thisData = myGaussMatrix_bySpeed(igspeed,:);

% as per DA discussion, normalize all to 1 to increase fit strength:
%norm between 0 and 1:

% thisData_scld=  (thisData - min(thisData)) / (max(thisData) - min(thisData));
%actually just norm height to 1, so that we dont shrink the width when
%scaling to zero:
thisData_upper=  thisData./ max(thisData);

thisData= thisData_upper;
% %%
% fnan = find(isnan(thisData));
% if ~isempty(fnan)
% thisData(fnan)= (thisData(fnan-1) + thisData(fnan+1))/2;
% end
% here we perform the fit:, note the bounds may need to be adjusted for
% accurate LHS/RHS fit. 

    
% startingVals = [ 0  50  100   1 ]; % .1];
% if NLIN
%     [ estimates,  resid, jacob, covarEst mse ] = nlinfit(SOAs, pSync,  @trySkewedGaussLR_nlin,startingVals );
%     ci = nlparci(estimates, resid, 'jacobian' , jacob); % calculate 95% confidence limits around estimates.
% else
%     options = optimset; 
%     estimates = fminsearch(@trySkewedGaussLR,startingVals,options,SOAs ,pSync);    
% end
% estimates = round(estimates*1000)/1000;
% 
% a = estimates(1); % Mean
% b = estimates(2); % SD of left half
% c = estimates(3); % SD of right half
%%
startingVals = [ 0  50  100   1 ]; % .1];
[ estimates,  resid, jacob, covarEst mse ] = nlinfit(SOAs_ms', thisData,  @trySkewedGaussLR_nlin,startingVals );
estimates = round(estimates*1000)/1000;
% sanity check (plot fit on top of observed data).
% clf;
% plot(SOAs_ms, thisData);
% hold on;
fittedC =  trySkewedGaussLR_nlin(estimates, SOAs_ms');
% hold on;
% plot(SOAs_ms, fittedC, 'r');
%%
participant_asymGauss_Fits_bySpeed(igspeed).datais= fitsTitle{igspeed};
participant_asymGauss_Fits_bySpeed(igspeed).gaussfit = [estimates];
participant_asymGauss_Fits_bySpeed(igspeed).gaussData_Y = fittedC;
participant_asymGauss_Fits_bySpeed(igspeed).gaussData_X = SOAs_ms;
participant_asymGauss_Fits_bySpeed(igspeed).observedData = thisData;


end % iGspeed
%% now perform fits on the extra gait split:


for igspeed= 1:size(myGaussMatrix_bySpeedbyResp,1) % 
for igalign = 1:nalignments
    for igait= 1:size(myGaussMatrix_bySpeedbyResp,3)
%pseudo code: interpolate the observed data to 1000 points, to improve the
%fit. 
%fit a gaussian using the above model, retain the parameters within this
%loop and save.

thisData = squeeze(myGaussMatrix_bySpeedbyResp(igspeed, igalign, igait,:));

% as per DA discussion, normalize all to 1 to increase fit strength:

% thisData=  (thisData - min(thisData)) / (max(thisData) - min(thisData));

thisData= thisData./max(thisData);


% hack: interp any nans.
% fnan = find(isnan(thisData));
% if ~isempty(fnan) && fnan~=1
% thisData(fnan)= (thisData(fnan-1) + thisData(fnan+1))/2;
% elseif  
    
% thisData(fnan)= (thisData(fnan-1) + thisData(fnan+1))/2;

startingVals = [ 0  50  100   1 ]; % .1];
[ estimates,  resid, jacob, covarEst mse ] = nlinfit(SOAs_ms', thisData',  @trySkewedGaussLR_nlin,startingVals );
estimates = round(estimates*1000)/1000;
%% sanity check (plot fit on top of observed data).
% clf;
% plot(SOAs_ms, thisData);
% hold on;
% hold on;
% plot(SOAs_ms, fittedC, 'r');


fittedC =  trySkewedGaussLR_nlin(estimates, SOAs_ms');

if any(abs(estimates)>1000)
    disp('check code, poor fit')
    estimates= nan(1,length(estimates));
    fittedC = nan(1, length(SOAs));

end

%%
disp(['estimates for ppant ' num2str(ippant) ' ' num2str(estimates)])
participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).datais= ['pSame' speedFlds{igspeed} alignFlds{igalign}];
participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).gaussfit = estimates;
participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).gaussData_Y = fittedC;
participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).gaussData_X = SOAs_ms;
participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).observedData = thisData;


    end % gait 
end %align
end % speed





save(allppants(ippant).name,'participant_asymGauss_Fits_bySpeed', 'participant_asymGauss_Fits_bySpeedbyAlignbyGait', '-append');
            disp(['finished all fits for ppant ' num2str(ippant)])
end % ppant
end % fit job

%%
if job.plotPPantFitts==1
%%
speedFlds= {'_all', '_slow', '_fast'};
alignFlds= {'_FirstAligned', '_VisAligned', '_AudAligned', '_LastAligned', '_RespAligned'};
%%
    for ippant = 1:length(allppants)
    clf;
        set(gcf, 'units', 'normalized', 'position', [0.05 0.05  .95 .95], 'color', 'w');

        cd(savedatadir);
load(allppants(ippant).name, 'participant_asymGauss_Fits_bySpeed','participant_asymGauss_Fits_bySpeedbyAlignbyGait', 'subjID')

cd(figdir)
cd('Participant asym Gaussian fits (all)');
 
%%
 SOAs= participant_asymGauss_Fits_bySpeed(1).gaussData_X; % in ms
 clf
 speedCols= {'k', 'b', 'r'};
for igspeed= 1:length(participant_asymGauss_Fits_bySpeed)
%%
    subplot(1,2,1); hold on;
   thisData= participant_asymGauss_Fits_bySpeed(igspeed).observedData;
    plot(SOAs, thisData, 'k-o', 'MarkerFaceColor', speedCols{igspeed});
    title(['Each speed, observed data, normalized']);
    hold on;
    ylim([0 1])
    %fix ylim
    yl= get(gca, 'ylim');
   xlim([SOAs(1) SOAs(end)])
    subplot(1,2,2); hold on;
    yf = participant_asymGauss_Fits_bySpeed(igspeed).gaussData_Y;
    g=plot(SOAs, yf, '-', 'LineWidth',2, 'color', speedCols{igspeed}); % fit.
%     legend(g, 'fit')
    
    ylim([0 1])
    
    % hold on, plot Mu. 
    mu = participant_asymGauss_Fits_bySpeed(igspeed).gaussfit(1);
        sigm_LHS = participant_asymGauss_Fits_bySpeed(igspeed).gaussfit(2);
        sigm_RHS = participant_asymGauss_Fits_bySpeed(igspeed).gaussfit(3);
        ytop = max(yf);
        hold on; 
        %
        % plot gauss mean:
        plot([mu mu], ylim, 'color', speedCols{igspeed})
        % plot sig (at Half max).
        plot([mu-sigm_LHS mu], [ytop/2+.05*igspeed ,ytop/2 + .05*igspeed ], 'linew', 1,'color', speedCols{igspeed})
        plot([mu mu+sigm_RHS], [ytop/2+.05*igspeed ,ytop/2 + .05*igspeed ], 'linew',4,'color', speedCols{igspeed})
        
        xlim([SOAs(1) SOAs(end)])
    
end % igfit.
shg
%%
print('-dpng', ['Participant ' num2str(ippant) ' ' subjID '_asymGauss by speed']);
%%

    
    
for igspeed= 1:size(participantGauss_Fits_bySpeedbyAlignbyGait,1) % 3 speeds.
        clf;
    scounter=1;
    for igalign= 1:nalignments % alignments

    
        for igait= 1:size(participant_asymGauss_Fits_bySpeedbyAlignbyGait,3) % 5 quintiles.
        
subplot(nalignments,nGaits,scounter); cla
   thisData= participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed, igalign, igait).observedData;
    plot(SOAs, thisData, 'k-o', 'MarkerFaceColor', 'k');
    
    ylim([0 1]);
    
    title(['Speed: ' speedFlds{igspeed} 'aligned to ' alignFlds{igalign}], 'interpreter', 'none')

    if igait==5
    text(0, .35, ['Sp:' num2str(igspeed)], 'HorizontalAlignment', 'center');
    text(0, .2, ['Al:' num2str(igalign)], 'HorizontalAlignment', 'center');
    text(0, .05, ['G:' num2str(igait)], 'HorizontalAlignment', 'center');
    hold on;
    end
    %fix ylim
    hold on;
    yl= get(gca, 'ylim');
    
    yf= participant_asymGauss_Fits_bySpeedbyAlignbyGait(igspeed, igalign, igait).gaussData_Y;
    if ~isempty(yf)
    g=plot(SOAs, yf, ':', 'LineWidth',2,'color', speedCols{igspeed}); % fit.
%     legend(g, 'fit')
    ylim(yl);
    end
    ylim([0 1]);
    scounter=scounter+1;
        end %igait. 
    end % aling
    
print('-dpng', ['Participant ' num2str(ippant) ' ' subjID '_asymGauss by DV alignment, speed' speedFlds{igspeed} ]);
end % igspeed





%
% title(['Participant '  num2str(ippant) ', ' allppants(ippant).name])
% text()
    end
    %%

end
% % Define a function to find the FWHM of an approximately Gaussian curve
% function fwhm = calcFWHM(x, y)
%     
% % for fzero to work , the function must cross zero line.
% y= y-mean(y);
% % Find the maximum value and corresponding x-coordinate
%     ymax = max(y);
%     xmax = x(y == ymax);
%     
%     % Define a function to find the roots of
%     halfmax = @(xx) abs(ymax/2 - interp1(x, y, xx, 'spline'));
%     
%     % Find the left and right roots of the halfmax function
%     left_root = fzero(halfmax, [min(x), xmax]);
%     right_root = fzero(halfmax, [xmax, max(x)]);
%     
%     % Calculate the FWHM as the difference between the left and right roots
%     fwhm = abs(right_root - left_root);
%     
% end