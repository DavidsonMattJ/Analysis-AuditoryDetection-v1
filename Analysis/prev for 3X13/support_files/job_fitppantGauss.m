% job_fitppantGauss
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

job.plotPPantFitts=1;

% >>>>>


% including option to use FWHM (since some synch functions are asymmetric).


%create a handle to a function called gaus that takes input parameters:
% x, array of x values
% mu, mean
% sig, SD
% amp, pos or negative curve
% vo, vertical offset.

% gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
% gaussEqn = 'a*exp(-((x-b)^2/(2*c^2)))+d';


%simpler ver (2 free params):

gaus = @(x,mu,sig)exp(-(((x-mu).^2)/(2*sig.^2)));
gaussEqn = 'exp(-((x-b)^2/(2*c^2)))';


figure(1); clf;
set(gcf, 'units', 'normalized', 'position', [.05 .05 .8 .8])


% note because we know the shape of the gaussians (approx), we can help the
% fit by specifying some lower bounds

foptions =fitoptions(gaussEqn);
% for each coefficient, apply some sensible bounds.
% order is a, b, c, d (amp, mean, SD, offset)


if job.fitper_participant==1
for ippant= 1:length(allppants)

    
    
% default fit params:    
% foptions.Lower= [0, -.3, 0,0];
% foptions.Upper= [10, .3, 20,1];

foptions.Lower= [-.3, .1];
foptions.Upper= [.3, 1];
foptions.StartPoint= [.2, 10];
foptions.MaxIter= 2000;

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
participantGauss_Fits_bySpeed=[];
participantGauss_Fits_bySpeedbyAlignbyGait=[];

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
%%
fnan = find(isnan(thisData));
if ~isempty(fnan)
thisData(fnan)= (thisData(fnan-1) + thisData(fnan+1))/2;
end

% as a bit of a hack, linearly interpolate the data points

% interpAt = linspace(SOAs(1), SOAs(end), 1000);
% propInterp= interp1(SOAs, thisData, interpAt);

% fix to floor:
interpAt = linspace(-1, 1, 1000);
propInterp= interp1([-1,SOAs',1], [0, thisData, 0], interpAt);

% [x1,x2,fw]= halfmaxpoints(interpAt, propInterp);
% fw= calcFWHM(interpAt, propInterp);

% perform fit :
    fit1 = fit(interpAt', propInterp', gaussEqn, foptions);

    %plot the results of the fit on top of the data:
    hold on;
%     fit_amp= fit1.a;

    fit_mu= fit1.b;

    fit_sig= fit1.c;

%     fit_vo= fit1.d;
%     yf = gaus(interpAt, fit_mu, fit_sig, fit_amp, fit_vo);

    yf = gaus(interpAt, fit_mu, fit_sig);
    
participantGauss_Fits_bySpeed(igspeed).datais= fitsTitle{igspeed};
participantGauss_Fits_bySpeed(igspeed).gaussfit = [fit_mu, fit_sig];
participantGauss_Fits_bySpeed(igspeed).gaussData_Y = yf;
participantGauss_Fits_bySpeed(igspeed).gaussData_X = interpAt;
participantGauss_Fits_bySpeed(igspeed).observedData = thisData;


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
fnan = find(isnan(thisData));
if ~isempty(fnan)
thisData(fnan)= (thisData(fnan-1) + thisData(fnan+1))/2;
end
% as a bit of a hack, linearly interpolate the data points
% interpAt = linspace(SOAs(1), SOAs(end), 1000);
% propInterp= interp1(SOAs, thisData, interpAt);

%fix to flor (improves fits)
 interpAt = linspace(-1, 1, 1000);
propInterp= interp1([-1,SOAs',1], [0, thisData', 0], interpAt);


% perform fit :
    fit1 = fit(interpAt', propInterp', gaussEqn, foptions);

    %plot the results of the fit on top of the data:
    hold on;
%     fit_amp= fit1.a;

    fit_mu= fit1.b;

    fit_sig= fit1.c;

    if fit_sig<.1
        error(' check fit for ppant ') % will be the case for delta functions.
    end
%     fit_vo= fit1.d;
%     yf = gaus(interpAt, fit_mu, fit_sig, fit_amp, fit_vo);

    yf = gaus(interpAt, fit_mu, fit_sig);
    
participantGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).datais= ['pSame' speedFlds{igspeed} alignFlds{igalign}];
participantGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).gaussfit = [fit_mu, fit_sig];
participantGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).gaussData_Y = yf;
participantGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).gaussData_X = interpAt;
participantGauss_Fits_bySpeedbyAlignbyGait(igspeed,igalign,igait).observedData = thisData;


    end % gait 
end %align
end % speed





save(allppants(ippant).name,'participantGauss_Fits_bySpeed', 'participantGauss_Fits_bySpeedbyAlignbyGait', '-append');
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
load(allppants(ippant).name, 'participantGauss_Fits_bySpeed','participantGauss_Fits_bySpeedbyAlignbyGait', 'subjID')

cd(figdir)
cd('Participant Gaussian fits (all)');
 
%%
 SOAs= [-0.32, -.14, -.05, .04, .13, .22, .4];
 clf
 speedCols= {'k', 'b', 'r'};
for igspeed= 1:length(participantGauss_Fits_bySpeed)

    subplot(1,2,1); hold on;
   thisData= participantGauss_Fits_bySpeed(igspeed).observedData;
    plot(SOAs, thisData, 'k-o', 'MarkerFaceColor', speedCols{igspeed});
    title(['Each speed, observed data, normalized']);
    hold on;
    ylim([0 1])
    %fix ylim
    yl= get(gca, 'ylim');
   xlim([SOAs(1) SOAs(end)])
    subplot(1,2,2); hold on;
    interpAt = participantGauss_Fits_bySpeed(igspeed).gaussData_X;
    % shrink size of interp to make the plots faster.
   
%     propInterp= interp1(SOAs, thisData, interpAt);

%     plot(interpAt, propInterp, 'r:'); % interpolated
   
    yf= participantGauss_Fits_bySpeed(igspeed).gaussData_Y;
    g=plot(interpAt, yf, '-', 'LineWidth',2, 'color', speedCols{igspeed}); % fit.
%     legend(g, 'fit')
    
    ylim([0 1])
    
    % hold on, plot Mu. 
    mu = participantGauss_Fits_bySpeed(igspeed).gaussfit(1);
        sigm = participantGauss_Fits_bySpeed(igspeed).gaussfit(2);
        ytop = dsearchn(interpAt', mu');
        hold on; 
        % plot gauss mean:
        plot([mu mu], [0 yf(ytop)], 'color', speedCols{igspeed})
        % plot sig (at Half max).
        plot([mu-sigm/2 mu+sigm/2], [yf(ytop)/2+.05*igspeed ,yf(ytop)/2 + .05*igspeed ], 'linew', 1,'color', speedCols{igspeed})
        
        xlim([SOAs(1) SOAs(end)])
    
end % igfit.
shg
%%
print('-dpng', ['Participant ' num2str(ippant) ' ' subjID '_Gauss by speed']);
%%

    
    
for igspeed= 1:size(participantGauss_Fits_bySpeedbyAlignbyGait,1) % 3 speeds.
        clf;
    scounter=1;
    for igalign= 1:nalignments % alignments

    
        for igait= 1:size(participantGauss_Fits_bySpeedbyAlignbyGait,3) % 5 quintiles.
        
subplot(nalignments,nGaits,scounter); cla
   thisData= participantGauss_Fits_bySpeedbyAlignbyGait(igspeed, igalign, igait).observedData;
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
    yl= get(gca, 'ylim');
    yyaxis right
    interpAt = participantGauss_Fits_bySpeedbyAlignbyGait(igspeed, igalign, igait).gaussData_X;
%     propInterp= interp1(SOAs, thisData, interpAt);

%     plot(interpAt, propInterp, 'r:'); % interpolated
   
    yf= participantGauss_Fits_bySpeedbyAlignbyGait(igspeed, igalign, igait).gaussData_Y;
    g=plot(interpAt, yf, '-', 'LineWidth',2,'color', speedCols{igspeed}); % fit.
%     legend(g, 'fit')
    ylim(yl);
    ylim([0 1]);
    scounter=scounter+1;
        end %igait. 
    end % aling
    
print('-dpng', ['Participant ' num2str(ippant) ' ' subjID '_Gauss by DV alignment, speed' speedFlds{igspeed} ]);
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