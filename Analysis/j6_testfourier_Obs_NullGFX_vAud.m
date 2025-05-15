
% j8_testfourier&NullGFX_vAUD
%
%This job precalculates the fourier fits for our data, based on the shuffled dataset
% Note that the previous jobs (j4_concatGFX and j7_createNull) need to be
% completed.

set_myOnlineDirectories_AUD;


% this script can be used to just quickly test the GFX across frequencies,
% or to also test the full perm to create null distributions, calculating 
% the fit per shuffle as well.:

testFull_perm = 0;


%% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%


%% There are various data types we will want to test
cd([procdatadir filesep 'GFX'])
load('GFX_Data_inGaits.mat');
load('GFX_Data_inGaits_FourierFits.mat') ; % we will append to the presaved, 
%in case completing at separate times.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   There are various datatypes, gaits, etc to step through (and save).
% GFX_FourierNull=[];
%%

testDVs = {'Accuracy', 'RT','dprime','crit','HR','FA','counts','counts'};

for testtype= 1:length(testDVs)

    usebin=1; % USE BINNED VERSIONS OF DATA (for now).
    
    %test DV relative to target/response onset:
    dataIN = GFX_TargPosData;
    typeOnset = 'Target';
    typeDV = testDVs{testtype};
    
    if testtype==8
        dataIN = GFX_RespPosData;
        typeOnset = 'Response';
        typeDV = testDVs{testtype};

    end

   
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
    
    %just one gait at a time (its a slow process). restrict to 2 gaits for
    %better Fourier fits.

   for  nGaits_toPlot=2%:2
    
    % plot_FourierFit(cfg,dataIN);
    %%
   
   
        iLR=3;
    gaitfield = {'gc', 'doubgc'};
    binfield = {'','_binned'};
    
    ppantData=[];
    
    shuffData=[];
    
    %which field of datastructure to plot?
    if strcmp(cfg.DV, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
     
    elseif strcmp(cfg.DV, 'Accuracy')
        usefield = [gaitfield{nGaits_toPlot} '_binned_Acc'];
    
    elseif  strcmp(cfg.DV, 'crit')
        usefield = [gaitfield{nGaits_toPlot} '_binned_crit'];
    else
          usefield = [gaitfield{nGaits_toPlot} '_binned_' cfg.DV];
    end
    

    % we also have 2 speeds:
    for iSpeed=1:3%1:3


    %collate data:
    for isub= 1:size(dataIN,1)
        
        ppantData(isub,:)= dataIN(isub,iSpeed, iLR).(usefield);
        if testFull_perm==1
        shuffData(isub,:,:) = dataIN(isub,iSpeed, iLR).([usefield '_shuff']);
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
    %% other specs:
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
    % perW=per
    fits_Rsquared_obsrvd = nan(1, length(Hzspace));
    
    fits_Rsquared_shuffCV = nan(3, length(Hzspace));
    if testFull_perm==1
    fits_Rsquared_shuff = nan(size(shuffData,2), length(Hzspace));
    meanShuff = squeeze(nanmean(shuffData,1));
    end

    gM=squeeze(nanmean(ppantData));
    
    %best fit (use model params below):
    f = fit(xvec', gM', 'fourier1'); %unbounded
    
    % step through w, forcing fit at particular periods, by updating the
    % bounds in fit options.
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
        
        %first test this period on observed data
        %set last coefficient value in the fourier model (w) to this value:
        
        [tmpf,gof] = fit(xvec', gM', 'fourier1', FitOpts);
        % how good/bad was the fit?
        fits_Rsquared_obsrvd(1,ifreq) = gof.rsquare;
%%         for each shuffle as well, calc the fits.        
        if testFull_perm==1 

        for iperm = 1:size(meanShuff)
            try
            [tmpf,gof] = fit(xvec', squeeze(meanShuff(iperm,:))', 'fourier1', FitOpts);
            fits_Rsquared_shuff(iperm,ifreq) = gof.rsquare;
            catch
                disp([' skipping perm ' num2str(iperm) ' ' cfg.type ' ' cfg.DV ...
            ', gaits(' num2str(nGaits_toPlot) ') likely a NaN'])
            end
        end % per freq
        
%         per freq, store the 95%CI. of Rsq values(plotted below)
        fits_Rsquared_shuffCV(:,ifreq) = quantile(fits_Rsquared_shuff(:,ifreq), [.05, .5, .95]);
        disp(['fin perm for freq ' num2str(ifreq) ' / ' num2str(length(Hzspace)) ' gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV])
        end
    end
    
    %store:
    GFX_FourierNull(iSpeed).([cfg.type 'Ons_' usefield '_fitsRsq_Obs']) = fits_Rsquared_obsrvd;
    GFX_FourierNull(iSpeed).([cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV']) = fits_Rsquared_shuffCV;
%     
    disp(['finished gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV  ', speed ' num2str(iSpeed)])
    end % iSpeed
   end % gaits
end %itype


% GFX_FourierNull(iSpeed).nShuff = size(meanShuff,1);
%%
save('GFX_Data_inGaits_FourierFits', 'GFX_FourierNull', 'Hzspace','-append');