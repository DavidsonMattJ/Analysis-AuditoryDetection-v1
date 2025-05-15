
% j9_testfourier_Obs_NullPFX_vAUD (slow)
%
%This job precalculates the fourier fits for our data, per participant, and compares to the shuffled dataset
% Note that the previous jobs (j4_concatGFX and j5_createNull) need to be
% completed.

%%%%%% Auditory. %%%%%%%%
% now with optimisations to use parrallel processing and reduce computation
% time.
clear all; close all;
set_myOnlineDirectories_AUD;


% this script can be used to just quickly test the observed result across frequencies,
% or to also test the full permutation/null distribution, to create CVs for significance testing
% at each frequency:

testFull_perm = 1;

usebin=1; %to do: adapt script to also fit fourier to raw data.

% Currently running the full perm on Target counts (double GC), to aid
% participant exclusion:
% have completed for ispeed ==1 (target counts, gait=2). elapsed time with
% parfor for N=40 ppants = 50 minutes. (.:. 1 speed per train trip).

%attempting speed 2 overnight..

%% determine numcores available, leave 1 for normal processing
% Check if parallel pool exists
pool = gcp('nocreate');  % Get current parallel pool without creating one

if isempty(pool)
    % No parallel pool exists, create one
    disp('Creating new parallel pool...');
    parpool('local', min(12, feature('numcores')-1));
else
    % Pool already exists, display info about it
    disp(['Parallel pool already exists with ' num2str(pool.NumWorkers) ' workers']);
end

%%
cd([procdatadir filesep 'GFX'])
% There are various data types we will want to test
load('GFX_Data_inGaits.mat');
load('GFX_Data_inGaits_FourierFits.mat') ; % we will append to the presaved, 
%in case completing at separate times.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   There are various datatypes, gaits, etc to step through (and save).
% PFX_FourierNull=[]; % on first run,


%how many shuffs?
nPerm = size(GFX_RespPos_nullData(1,1,1).doubgc_Acc,1);
% will need to run overnight (very slow) on full_perm

testDVs = {'Accuracy', 'RT','dprime','crit','HR','FA','counts','counts'};

gaitfield = {'gc', 'doubgc'};
binfield = {'','_binned'};
tic
for testtype= 7%1:8

    dataIN = GFX_TargPosData;
    shuffIN= GFX_TargPos_nullData;
    typeOnset = 'Target';
    typeDV = testDVs{testtype};

    if testtype==8
        dataIN = GFX_RespPosData;
        shuffIN= GFX_RespPos_nullData;
        typeOnset = 'Response';
        typeDV = testDVs{testtype};

    end


    usebin=1; % USE BINNED VERSIONS OF DATA (for now).
    % binfield = {'','_binned'};

    
    %%
    cfg=[];
    cfg.subjIDs = subjIDs;
    cfg.type = typeOnset;
    cfg.DV = typeDV;
%     cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel = 'GFX'; % plot separate figures per participant
    cfg.norm=0; % already z scored, so don't tweak.
    cfg.ylims = [-.15 .15]; % if norm =0;
    cfg.normtype= 'relative';
    %%
    
    %just one gait at a time (its a slow process).
   for  nGaits_toPlot=2%1:2
    
       
    if nGaits_toPlot==1
        
        pidx= cfg.pidx1;
        ftnames= {'LR', 'RL', 'combined'};
    else
        pidx= cfg.pidx2;
        ftnames= {'LRL', 'RLR', 'combined'};
    end
    
    %note that pidx is adjusted to all datapoints, if not using the  bin.
    if usebin==0
        pidx=1:size(ppantData,2);
    end
    
    %%
   
    iLR=3; % use combined data (not separating feet)
   
    ppantData=[];
    plotHead=[];
    shuffData=[];
    
    %which field of datastructure to extract/plot?
    if strcmp(cfg.DV, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
        
    elseif strcmp(cfg.DV, 'Accuracy')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_Acc'];
        
    else % for counts HR,FA,crit,drime, the names match
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_' typeDV ];
       
    
    end
    


    for iSpeed=2%1:3 %slow, natural, combined.
    %collate data:

    for isub= 1:size(dataIN,1)
        
        ppantData(isub,:)= dataIN(isub,iSpeed,iLR).(usefield);
        if testFull_perm==1
        shuffData(isub,:,:) = shuffIN(isub,iSpeed,iLR).(usefield);
        end

    end
    %% if normON , normalize as appropriate
%     
%     if cfg.norm==1
%         pM = nanmean(ppantData,2);
%         meanVals= repmat(pM, 1, size(ppantData,2));
%         
%         
%         if strcmp(cfg.normtype, 'absolute')
%             data = ppantData - meanVals;
%         elseif strcmp(cfg.normtype, 'relative')
%             data = ppantData  ./ meanVals;
%             data=data-1;
%         elseif strcmp(cfg.normtype, 'relchange')
%             data = (ppantData  - meanVals) ./ meanVals;
%         elseif strcmp(cfg.normtype, 'normchange')
%             data = (ppantData  - meanVals) ./ (ppantData + meanVals);
%         elseif strcmp(cfg.normtype, 'db')
%             data = 10*log10(ppantData  ./ meanVals);
%         end
%         
%         ppantData= data;
%     end
    
    
    %x axis:          %approx centre point of the binns.
    mdiff = round(mean(diff(pidx)./2));
    xvec = pidx(1:end-1) + mdiff;
    
    %% extract fourier fits per shuffled series:
    % Declaring the type of fit.
    FitType = 'fourier1';
    % Creating and showing a table array to specify bounds.
    CoeffNames = coeffnames(fittype(FitType));
    %%
    %set bounds for w
    CoeffBounds = array2table([-Inf(1,length(CoeffNames));...
        Inf(1,length(CoeffNames))],'RowNames',...
        ["lower bound", "upper bound"],'VariableNames',CoeffNames);
    %%
    % Specifying bounds according to the position shown by the table.
    % e.g. to force fit with w ~ 1.545, we ChangeBound of w parameter
    % CoeffBounds.w(1) = 1.54;
    % CoeffBounds.w(2) = 1.55;
    %
    Hzspace = [0.01:.2:10];
    
    
    %preallocate:

    fits_Rsquared_obsrvd = nan(1, length(Hzspace));
    if testFull_perm==1
    fits_Rsquared_shuff = nan(size(shuffData,2), length(Hzspace));
    fits_Rsquared_shuffCV = nan(3, length(Hzspace));
    end
    
%     meanShuff = squeeze(mean(shuffData,1));
%     gM=squeeze(nanmean(ppantData));
    
    %best fit (use model params below):
    gM= squeeze(mean(ppantData,1));
    f = fit(xvec', gM', 'fourier1'); %unbounded
 
    
    %% loop over participants (slow)!, pull out the computations which 
    % don't need to be performed each loop:
    
    %e.g. predetermine the fit options per frequency to test.
    %%
    fitoptsperHz=[];
    for ifreq= 1:length(Hzspace)
        % include period and Rsquared
        %treat max xvec as our full 'period'
        %             Hzapp = xvec(end)/ (2*pi/(f.w));
        testw = 2*pi*Hzspace(ifreq)/xvec(end);
        
        CoeffBounds.w(1) = testw;
        CoeffBounds.w(2) = testw;
        
        %update fit opts settings
        
        %Update Fit Options setting.
        FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
            'Upper',table2array(CoeffBounds(2,:)));
        
        fitoptsperHz{ifreq}= FitOpts;
    end
    %%


    for ippant = 1:size(ppantData,1)
    % step through w, forcing fit at particular periods, by updating the
    % bounds in fit options
    % 
    % 
    
            for ifreq= 1:length(Hzspace)
        
        usefit = fitoptsperHz{ifreq};
        %first test this period on observed data
        %set last coefficient value in the fourier model (w) to this value:
        
        [tmpf,gof] = fit(xvec', ppantData(ippant,:)', 'fourier1', usefit );
        
        % how good/bad was the fit on observed data?
        fits_Rsquared_obsrvd(1,ifreq) = gof.rsquare;
        %%
        if testFull_perm==1
            
        % for each shuffle as well, calc preallocate and calc the fits.
        permarray=[];
        parfor iperm = 1:nPerm
            % try
            [tmpf,gof] = fit(xvec', squeeze(shuffData(ippant,iperm,:)), 'fourier1', usefit);
            fits_Rsquared_shuff(iperm, ifreq) = gof.rsquare;
            % catch
                % disp([' skipping perm ' num2str(iperm) ' ' cfg.type ' ' cfg.DV ...
            % ', gaits(' num2str(nGaits_toPlot) ') likely a NaN'])
            % end
        end % per perm

        
        % per freq, store the 95%CI. of Rsq values(plotted below)
        fits_Rsquared_shuffCV(:,ifreq) = quantile(fits_Rsquared_shuff(:,ifreq), [.05, .5, .95]);
        end
        toc
    disp(['Fin freq ' num2str(ifreq)])
    end
    
    %store:
    PFX_FourierNull(ippant, iSpeed).([cfg.type 'Ons_' usefield '_fitsRsq_Obs']) = fits_Rsquared_obsrvd;
    PFX_FourierNull(ippant,iSpeed).([cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV']) = fits_Rsquared_shuffCV;
    
    disp(['finished gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV 'ippant: ' num2str(ippant) ])
    end % participant

    disp(['finished gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV ', speed ' num2str(iSpeed)])
    end
end % nGaits

end %itype

%%
save('GFX_Data_inGaits_FourierFits', 'PFX_FourierNull', 'Hzspace');%,'-append')
