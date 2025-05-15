function plot_GFX_groupFits_vAUD_MSver2(cfg,testData, GFX_FourierNull)
% % % helper function to plot a MS ready summary of the group fits to
% detection performance.

% will include an example trial, and target onset coutns over the gait.

% !! only the gait (linked steps),in this version:

% called from j9__Plot_Data_FourierFits_GFX

% % % % set some parameters:

GFX_headY = cfg.HeadData;

gaitfield = {'gc', 'doubgc'};
gaitprint = {'step', 'stride'};
binfield = {'','_binned'};
nGaits_toPlot= cfg.nGaits_toPlot;
usebin=cfg.usebin;
ppantData=[];
plotHead=[];
% 3 colours 1 per DV
barCols= {'b', 'r', 'm', 'c', 'y'};
speedCols= [.2 .2 1; 1 .2 .2];
% fntsize=12;

nDVs= length(testData);
iLR=3;

xindexes = {cfg.pidx1, cfg.pidx2};

figure(1);
clf;
set(gcf, 'units', 'normalized', 'position', [.01 .01 .95 .95], 'color', 'w');

% first plot the example trial (function below).
subplot(2,4,1:2);

plotextrial(cfg);
% % now plot the target onsets for the next panels:

dataIN= testData{1}; % target onset info:
ylabis = ['Target counts'];



%%
for igait= 2%1:2 % both gait versions:
    
    usefield = [gaitfield{igait} '_binned_counts'];
    pidx = xindexes{igait};
    
    for iSpeed=1:2
        %collate data across subjs:
        ppantData=[];
        plotHead=[]
        for isub= 1:size(dataIN,1)
            
            ppantData(isub,:)= dataIN(isub,iSpeed,iLR).(usefield);
            plotHead(isub,:) = GFX_headY(isub,iSpeed).([gaitfield{igait}]);
            
        end
        
        % set up plot:
        subplot(2,4,3)
        hold on
        % yyaxis right
        stEH= CousineauSEM(plotHead);
        pH= nanmean(plotHead,1);
        
        %resize for plot?
        pHs= imresize(pH, [1,100]);
        stEHs= imresize(stEH, [1,100]);
        sh= shadedErrorBar(1:100, pHs, stEHs,{'color',speedCols(iSpeed,:),'linew',1}) ;
        sh.mainLine.LineStyle='-';
        sh.mainLine.LineWidth = 4;
        legp = sh.mainLine;
        set(gca,'fontsize',cfg.fntsize, 'YColor', 'k');
        
        % ylim([-.2 .05]) % if detrended
        ylabel('detrended Head height')
        xlabel( '% stride-cycle');
        
        %% add swing/stance phases?
        if iSpeed==2
            markersAt = ceil([15, 75 115 175]/2); % note that 40-60 is the usual stance-swing split.
            yat = get(gca, 'ylim');
            %make room for data:
            ylim(yat)
            ylim([yat(1)-.005 yat(2)])
            Yverts = [yat(1) - .005, yat(1)];
            Xverts = [ 0, markersAt, 100];
            pos1 = [0 Yverts(1) markersAt(1) .005]; % x y width height
            pos2 = [markersAt(1) Yverts(1) markersAt(2)-markersAt(1) .005];
            pos3 = [markersAt(2) Yverts(1) markersAt(3)-markersAt(2), .005];
            pos4 = [markersAt(3) Yverts(1) markersAt(4)-markersAt(3), .005];
            pos5 = [markersAt(4) Yverts(1) 100-markersAt(4), .005];
            
            
            h=rectangle('Position', pos1, 'FaceColor',[.8 .8 .8]);
            h=rectangle('Position', pos2, 'FaceColor','w');
            h=rectangle('Position', pos3, 'FaceColor',[.8 .8 .8]);
            h=rectangle('Position', pos4, 'FaceColor','w');
            h=rectangle('Position', pos5, 'FaceColor',[.8 .8 .8]);
            
            text(mean(markersAt(1:2)), Yverts(1), 'Swing', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize',12);
            text(mean(markersAt(2:3)), Yverts(1), 'Stance', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize',12);
            
            text(mean(markersAt(3:4)), Yverts(1), 'Swing', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize',12);
            %
        end
        %%
        
        % set xticks at approx centre points of the bins.
        mdiff = round(mean(diff(pidx)./2));
        xvec = pidx(1:end-1) + mdiff; % xaxis vector
        midp=xvec(ceil(length(xvec)/2)); % mid point:
        
        
        %preset fit range:
        
        
        %restrict to 0 - 10 cps: for fits
        Hzspace = [0.01:.2:10];
        
        testw_low = 2*pi*Hzspace(1)/xvec(end);
        testw_high = 2*pi*6/xvec(end);
        
        FitType = 'fourier1';
        % Creating and showing a table array to specify bounds.
        CoeffNames = coeffnames(fittype(FitType));
        
        %set bounds for w (when fitting)
        CoeffBounds = array2table([-Inf(1,length(CoeffNames));...
            Inf(1,length(CoeffNames))],'RowNames',...
            ["lower bound", "upper bound"],'VariableNames',CoeffNames);
        CoeffBounds.w(1) = testw_low;
        CoeffBounds.w(2) = testw_high;
        
        %update fit opts settings
        
        %Update Fit Options setting.
        FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
            'Upper',table2array(CoeffBounds(2,:)));
        %%
        
        
        
        subplot(4,4,4 +4*(iSpeed-1));
        [bh, Fh]=bar_wFourier(xvec, ppantData, FitOpts);
        
        % tidy with some specifics:
        Fh(2).Color= 'k';
        Fh(2).Visible= 'off'
        ylabel(ylabis)
        set(gca,'fontsize', cfg.fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})
        xlabel(['Target onset (% ' gaitprint{igait} '-cycle )']);%
        lg= get(gca, 'legend');
        lg.Visible= 'off';
        % ylim([0 10*igait])
        
        gM= nanmean(ppantData,1);
        sdrange = max(gM) - min(gM);
        ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
        
        %% apply best fit (overlay)
        [f,gof]= fit(xvec',gM',  'fourier1', FitOpts);
        
        % plot:
        hold on;
        %             yyaxis right
        h=plot(f, xvec, gM);%,
        h(2).LineWidth = 2;
        h(2).Color = speedCols(iSpeed,:);
        %
        %treat max xvec as our full 'period'
        fitperiod = f.w;
        %convert to period per samples.
        % include period and Rsquared
        %treat max xvec as our full 'period'
        Hzapp = xvec(end)/ (2*pi/(f.w));
        legdetails = [sprintf('%.2f', Hzapp) '(cpg), R^2 = ' sprintf('%.2f', gof.rsquare) ];
        
        % legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];
        
        
        ylabel(['Counts']);%
        legend off
        if iSpeed==2
            xlabel(['Auditory onset (% ' gaitprint{igait} '-cycle )']);%
        end
%         legend(h(2), legdetails, 'fontsize', cfg.fntsize-2, 'autoupdate', 'off', 'Location', 'NorthEast')
    end
end


shg









%% now cycle through data types and plot each time:

%
% nDVs
% subspots1= [nDVs+1:(2*nDVs)+1];
% subspots2 = subspots1+nDVs+1;
% subspots= [subspots1;subspots2];

%sub(3,6,ip); starting second row.

subspots=[4,5,7,8];
ic=1;
for igait=2%1:2
    
    for iSpeed=3%1:2
        pidx = xindexes{igait};
        iLR=3;
        nGaits_toPlot= igait;
        
        for id= 1:length(testData)
            
            
            dataIN= testData{id};
            
            
            %% wrangle data to plot:
            
            %which field of the datastructure to plot?
            if strcmp(cfg.DV(id), 'RT')
                usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
                
            elseif strcmp(cfg.DV(id), 'Accuracy')
                usefield = [gaitfield{nGaits_toPlot} '_binned_Acc'];
                
                
            elseif strcmp(cfg.DV(id), 'Counts')
                usefield = [gaitfield{nGaits_toPlot} '_binned_counts'];
                
            elseif strcmp(cfg.DV(id), 'criterion')
                usefield = [gaitfield{nGaits_toPlot} '_binned_crit'];
                
                
                
            elseif strcmp(cfg.DV(id), 'dprime')
                usefield = [gaitfield{nGaits_toPlot} '_binned_dprime'];
                
                
            elseif strcmp(cfg.DV(id), 'FA')
                usefield = [gaitfield{nGaits_toPlot} '_binned_FA'];
                
            elseif strcmp(cfg.DV(id), 'HR')
                usefield = [gaitfield{nGaits_toPlot} '_binned_HR'];
                
            end
            ylabis= [cfg.DV(id)];
            
            ppantData=[];
            plotHead=[];
            %collate data across subjs:
            for isub= 1:size(dataIN,1)
                
                ppantData(isub,:)= dataIN(isub,iSpeed,iLR).(usefield);
                
                plotHead(isub,:) = GFX_headY(isub,iSpeed).([gaitfield{nGaits_toPlot}]);
            end
            
            
            %% other specs:
            if nGaits_toPlot==1
                
                pidx= cfg.pidx1;
                ftnames= {'LR', 'RL', 'combined'};
            else
                pidx= cfg.pidx2;
                ftnames= {'LRL', 'RLR', 'combined'};
            end
            
            %note that pidx is adjusted to all datapoints, if not using the  binned version.
            if usebin==0
                pidx=1:size(ppantData,2);
            end
            
            % set xticks at approx centre points of the bins.
            mdiff = round(mean(diff(pidx)./2));
            xvec = pidx(1:end-1) + mdiff;
            
            
            
            
            
            
            
            subplot(2,3, subspots(ic))
            
            [bh, Fh]=bar_wFourier(xvec, ppantData,FitOpts);
            
            gM = squeeze(mean(ppantData));
            stE = CousineauSEM(ppantData);
            hold on;
            % finely sampled bar, each gait "%" point.
            bh=bar(xvec, gM);
            hold on;
            errorbar(xvec, gM, stE, ...
                'color', 'k',...
                'linestyle', 'none',...
                'linew', 2);
            bh.FaceColor = [.9 .9 .9];
            bh.EdgeColor = barCols{id};
            
            sdrange = max(gM) - min(gM);
            ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
            ylsat = get(gca, 'ylim');
            
            
            midp=xvec(ceil(length(xvec)/2));
            
            %% apply best fit (overlay)
            
            if any(isnan(gM))
                tmpN= find(isnan(gM));
                gM(tmpN) = nanmean(gM);
            end
            
            
            %first test this period on observed data
            %set last coefficient value in the fourier model (w) to this value:
            
            [f,gof] = fit(xvec', gM', 'fourier1', FitOpts);
            % [f,gof]= fit(xvec',gM',  'fourier1');
            
            % plot:
            hold on;
            %             yyaxis right
            h=plot(f, xvec, gM);%,
            h(2).LineWidth = 2;
            h(2).Color = barCols{id};
            %%
            %treat max xvec as our full 'period'
            fitperiod = f.w;
            %convert to period per samples.
            % include period and Rsquared
            %treat max xvec as our full 'period'
            Hzapp = xvec(end)/ (2*pi/(f.w));
            legdetails = [sprintf('%.2f', Hzapp) '(cps), R^2 = ' sprintf('%.2f', gof.rsquare) ];
            
            % legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];
            
            
            legend(h(2), legdetails, 'fontsize', cfg.fntsize-2, 'autoupdate', 'off', 'Location', 'NorthEast')
            
            ylabel(ylabis)
            set(gca,'fontsize', cfg.fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})
            
            xlabel([ 'Auditory onset (% ' gaitprint{nGaits_toPlot} '-cycle )']);%
            
            
            
            
            
            %Hzspace% loaded above
            fits_Rsquared_obsrvd=GFX_FourierNull(iSpeed).([cfg.type{id} 'Ons_' usefield '_fitsRsq_Obs']);
            if cfg.plotShuff
                fits_Rsquared_shuffCV=GFX_FourierNull(iSpeed).([cfg.type{id} 'Ons_' usefield '_fitsRsq_ShuffCV']);
            end
            
            if iSpeed<3
            subplot(2,3,6+3*(iSpeed-1));
            else
                subplot(2,3,6);
            end
                
            hold on;
            
            % if the first, plot a grey version and withold handle for the legend
            if id==1 && cfg.plotShuff==1
                pG= plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:),':','linew',2,'color', [.3 .3 .3]);
            end
            pO=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_obsrvd,'color', barCols{id}, 'linew', 3);
            ylabel('R^2');
            xlabel(['Frequency (cycles-per-' gaitprint{igait} ')'])
            
            
            legHZ(id) = pO;
            if cfg.plotShuff
                hold on;
                ph=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:), ':', 'linew', 2, 'color', barCols{id});
            end
            set(gca,'fontsize', cfg.fntsize)
            % title('(forced) Fits per frequency', 'fontsize', 18)
            xlim([.5 10])
            xlim([.5 5])
            ylim([0 .9]);
            %% plot below patch?
            % create the patch
            if cfg.plotShuff
                x=cfg.Hzspace(1:length(fits_Rsquared_obsrvd));
                y=fits_Rsquared_shuffCV(3,:);
                patch([x fliplr(x)], [zeros(size(x)) fliplr(y)], [.8 .8 .8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
            
            
            
            
            if (id==3 || id==6) && cfg.plotShuff==1
                legend(pG, '95% CI')
            end
            
            
            
            
            
            ic=ic+1;
            
        end % id
    end % iSpeed
end % igait

%%
% overlay head position?
% hold on
% yyaxis right
% stEH= CousineauSEM(plotHead);
% pH= nanmean(plotHead,1);
% sh= shadedErrorBar(1:size(plotHead,2), pH, stEH,'k',1) ;
% sh.mainLine.LineStyle='-';
% sh.mainLine.LineWidth = 4;
% legp = sh.mainLine;
% set(gca,'ytick', [], 'fontsize',fntsize, 'YColor', 'k');
%
% ylim([-.2 .05])
%>>>>> END DV with fits

%%  also show fit strength across frequencies:
%% extracted fourier fits per shuffled series (calculated in j8_' loaded above)

ylim([0 .9])











end
function plotextrial(cfg)

setmydirs_detectv3;
cd(procdatadir);
subjID= 'ABB';
itrial= 124;
cd(procdatadir)
lfile  = dir([pwd filesep subjID '*']);
load(lfile.name)
fntsize= cfg.fntsize;

%
plot(trialInfo(124).times, HeadPos(124).Y, 'k', 'LineWidth',2);
% axis tight;
xlabel('Trial time (sec)');
ylabel('Head height (m)');
set(gca,'fontsize', fntsize);
ylim([1.69 1.78])
%
% tOnsets = HeadPos(124).tr
relT= find(trial_summaryTable.trial == itrial);
tOnsets = trial_summaryTable.targOnset(relT);
tResult = trial_summaryTable.targCor(relT);
tResult= [1,4,3,1,2,1];
hold on;
trialTimes = trialInfo(itrial).times;
HeadData= HeadPos(itrial).Y;
txtHeight = 1.7;
for itarg = 1:length(tOnsets)
    
    pAt = dsearchn(trialTimes, tOnsets(itarg)');
    
    if tResult(itarg)==1% Hit.
        %     ht=text(tOnsets(itarg), txtHeight, 'H', 'HorizontalAlignment', 'center')
        hl=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'r');
        
        plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'k')
        
    elseif tResult(itarg)==2 % Miss
        
        %     mt=text(tOnsets(itarg), txtHeight, 'M', 'HorizontalAlignment','center')
        ml=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'w');
        plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'w')
    
    elseif tResult(itarg)==3 % FA
        
        %     mt=text(tOnsets(itarg), txtHeight, 'M', 'HorizontalAlignment','center')
       fal=plot(tOnsets(itarg), HeadData(pAt), 'bo', 'MarkerFaceColor', 'w');
        plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'w')
    elseif tResult(itarg)==4 % CR
        
        %     mt=text(tOnsets(itarg), txtHeight, 'M', 'HorizontalAlignment','center')
        crl=plot(tOnsets(itarg), HeadData(pAt), 'bo', 'MarkerFaceColor', 'b');
        plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'b')
    end
    
    % connect to base
    plot([tOnsets(itarg) tOnsets(itarg)], [1.695 1.72], 'k-')
    % plot(tOnsets(itarg), 1.7, 'dk', 'MarkerFaceColor', 'k')
end
legend([hl,ml,fal,crl], {'Hit', 'Miss', 'FA','CR'})

text(0.1, 1.71, {['Auditory'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', 15);
xlim([0 9])
box off
set(gca,'fontsize',15)
end



% %
function [bh, Fh]= bar_wFourier(xvec, ppantData,FitOpts)

% called to plot the nice bar charts with fourier overlayed. returns handle
% for bar (bh) and Fourier fit (Fh)

gM = squeeze(mean(ppantData));
stE = CousineauSEM(ppantData);
hold on;
% finely sampled bar, each gait "%" point.
bh=bar(xvec, gM);
hold on;
errorbar(xvec, gM, stE, ...
    'color', 'k',...
    'linestyle', 'none',...
    'linew', 2);
bh.FaceColor = [.9 .9 .9];
% bh.EdgeColor = barCols{id};

sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylsat = get(gca, 'ylim');


midp=xvec(ceil(length(xvec)/2));

if any(isnan(gM))
    tmpN= find(isnan(gM));
    gM(tmpN) = nanmean(gM);
end

%% apply best fit (overlay)
[f,gof]= fit(xvec',gM',  'fourier1', FitOpts);

% plot:
hold on;
%             yyaxis right
Fh=plot(f, xvec, gM);%,
Fh(2).LineWidth = 2;
%treat max xvec as our full 'period'
% fitperiod = f.w;
%convert to period per samples.
% include period and Rsquared
%treat max xvec as our full 'period'
% Hzapp = xvec(end)/ (2*pi/(f.w));
% legdetails = [sprintf('%.2f', Hzapp) ' cycles per, R^2 = ' sprintf('%.2f', gof.rsquare) ];
% legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];

% legend(Fh(2), legdetails, 'fontsize', 8, 'autoupdate', 'off', 'Location', 'NorthEast')


end