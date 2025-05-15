
% J7_testfourier_Obs_Null_GFX_AudDetect
%
%This job precalculates the fourier fits for our data, based on the shuffled dataset
% Note that the previous jobs (j4_concatGFX and j7_createNull) need to be
% completed.
%% show ppant numbers:
cd(savedatadir);
pfols = dir([pwd filesep '*_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
procdatadir= savedatadir; % compatability below.


%% There are various data types we will want to test
cd([procdatadir filesep 'GFX'])
load('GFX_Data_inGaits.mat');
load('GFX_Data_inGaits_FourierFits.mat') ; % we will append to the presaved, 
%in case completing at separate times.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   There are various datatypes, gaits, etc to step through (and save).
% GFX_FourierNull=[];

for testtype= 3%[1,3,6,7]

    usebin=1; % USE BINNED VERSIONS OF DATA (for now).
    switch testtype
        case 1
            %test RT relative to target onset:
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV = 'RT';
        case 2
            %%
            %test RT relative to response onset:
            dataIN = GFX_RespPosData;
            typeOnset = 'Response';
            typeDV = 'RT';
            %%
        case 3
            %test ACC relative to target onset.
            
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV = 'Accuracy';
            
        case 4
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_RespPosData;
            typeOnset = 'Response';
            typeDV='Counts';
        case 5
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV='Counts';
            
        case 6
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV='dprime';
        case 7
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV='criterion';
            
            case 8
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV='HR';
            
            case 9
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV='FA';
            
            
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
    cfg.usebin=1;
    %%
    
    %just one gait at a time (its a slow process).
    for  nGaits_toPlot=2%1:2
        
        % plot_FourierFit(cfg,dataIN);
        %%
        
        
        iLR=3;
        gaitfield = {'gc', 'doubgc'};
        binfield = {'','_binned'};
        
        tmpfield = [gaitfield{nGaits_toPlot} binfield{cfg.usebin+1}];
        
        ppantData=[];
        
        shuffData=[];
    
        
        if strcmp(cfg.DV, 'RT')
            usefield = [tmpfield '_rts'];
        elseif strcmp(cfg.DV, 'Accuracy')
            usefield = [tmpfield '_Acc'];
        elseif strcmp(cfg.DV, 'dprime')
            usefield = [tmpfield '_dprime'];
        elseif strcmp(cfg.DV, 'criterion')
            usefield = [tmpfield '_crit'];
        elseif strcmp(cfg.DV, 'HR')
            usefield = [tmpfield '_HR'];
        elseif strcmp(cfg.DV, 'FA')
            usefield = [tmpfield '_FA'];
        elseif strcmp(cfg.DV, 'Counts')
            usefield = [tmpfield '_counts'];
        end
        
        
        for ispeed=3%1:3
            %collate data:
            for isub= 1:size(dataIN,1)
                
                ppantData(isub,:)= dataIN(isub,ispeed,iLR).([usefield]);
                shuffData(isub,:,:) = dataIN(isub,ispeed,iLR).([usefield '_shuff']);
                
            end
            %% if normON , normalize as appropriate
            
            if cfg.norm==1
                pM = nanmean(ppantData,2);
                meanVals= repmat(pM, 1, size(ppantData,2));
                
                
                if strcmp(cfg.normtype, 'absolute')
                    data = ppantData - meanVals;
                elseif strcmp(cfg.normtype, 'relative')
                    data = ppantData  ./ meanVals;
                    data=data-1;
                elseif strcmp(cfg.normtype, 'relchange')
                    data = (ppantData  - meanVals) ./ meanVals;
                elseif strcmp(cfg.normtype, 'normchange')
                    data = (ppantData  - meanVals) ./ (ppantData + meanVals);
                elseif strcmp(cfg.normtype, 'db')
                    data = 10*log10(ppantData  ./ meanVals);
                end
                
                ppantData= data;
            end
            
            
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
    fits_Rsquared_shuff = nan(size(shuffData,2), length(Hzspace));
    fits_Rsquared_shuffCV = nan(3, length(Hzspace));
    meanShuff = squeeze(mean(shuffData,1));
   
   %remove any infs
   tmpinf = isinf(ppantData);
   ppantData(tmpinf)=nan;

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
%          testw = 2*pi*Hzspace(ifreq)/100;
        
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
        for iperm = 1:500%size(meanShuff,1)
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
        disp(['fin perm for freq ' num2str(ifreq) ' / ' num2str(length(Hzspace)) ' gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV ',speed: ' num2str(ispeed)])
    end
    
    % sanity check.
%     clf;
%     plot(Hzspace, fits_Rsquared_obsrvd)


    %store:
    GFX_FourierNull(ispeed).([cfg.type 'Ons_' usefield '_fitsRsq_Obs']) = fits_Rsquared_obsrvd;
    GFX_FourierNull(ispeed).([cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV']) = fits_Rsquared_shuffCV;
    end
    disp(['finished gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV])
end % nGaits

end %itype
GFX_FourierNull(ispeed).nShuff = size(meanShuff,1);
%%
save('GFX_Data_inGaits_FourierFits', 'GFX_FourierNull', 'Hzspace','-append')