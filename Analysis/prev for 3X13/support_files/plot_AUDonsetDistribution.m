function  plot_AUDonsetDistribution(dataIN, cfg)
% helper function to plot either the distribution of all target onsets, or
% all response onsets (clicks), relative to gait "%"

% called from the script
% plot_ReactionTime_

GFX_headY = cfg.HeadData;
usecols = {[0 .7 0], [.7 0 0], [.7 0 .7]}; % R Gr Prp

%Note that this plots both target and response onset relative to gait.
% For the response ver, we also overlay all False Alarms recorded
useFA =strcmp(cfg.type, 'Response');
if cfg.usebin
    useFA=0; % dont plot FA on binned x axis
end

% set colours for diff speeds:
speedCols = {'b', 'r'};
figure(1); clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0 0 .9 .9]);
nsubs = length(cfg.subjIDs);
%% set up shortcuts for accessing data:
binfields = {'', '_binned'};
usefield = [binfields{cfg.usebin+1} '_counts'];
usegaitnames= {'gc','doubgc'};

if strcmp(cfg.plotlevel, 'PFX')

    for ippant = 1:nsubs
        clf;
        pc=1; % plot counter
        psubj= cfg.subjIDs{ippant}; % print ppid.
        % both this and the next use the same figure function:
        pspots=1:6;

        for nGaits_toPlot=1:2

            if nGaits_toPlot==1
                pidx= cfg.pidx1;
                ftnames= {'LR', 'RL', 'combined'};
            else
                pidx= cfg.pidx2;
              
            end

            legp=[]; % for legend
            for iLR=1:3

                % % note that for doubgc, we skip to index 4.
                

                    ppantData= dataIN(ippant,iLR).([usegaitnames{nGaits_toPlot} usefield]);
                    plotHead = GFX_headY(ippant).(usegaitnames{nGaits_toPlot});

              

                %x axis:
                if cfg.usebin

                    %x axis:          %approx centre point of the binns.
                    mdiff = round(mean(diff(pidx)./2));
                    xvec = pidx(1:end-1) + mdiff;
                else  % note that if we aren't using the binned versions, use full xaxis.
                    %         (overwrite the above)
                    xvec = 1:pidx(end);
                end

                %%

              
                subplot(2,3,pspots(pc))
                hold on;
                yyaxis left

                % finely sampled bar, each gait "%" point.
                bh=bar(xvec, ppantData);
                bh.FaceColor = speedCols{iSpeed};
                legp(iLR)= bh;

                ylabel([cfg.type ' onset [counts]']);
                ylim([0 10*nGaits_toPlot])
                %% add FA if response type:

%                 if useFA
%                     bh=bar(xvec, faData);
%                     bh.FaceColor = 'm';
%                 end
                %%


                yyaxis right
                pHs= imresize(plotHead, [1,100]);
                ph=plot(pHs, ['k-o'], 'linew', 3); hold on
                set(gca,'ytick', []);
                if cfg.usebin==0
                    title([psubj ' (N ' num2str(sum(ppantData, 'omitnan')) ') ' ftnames{iLR}], 'interpreter', 'none');
                else
                    title([psubj ' ' ftnames{iLR}], 'interpreter', 'none');
                end

                midp=xvec(ceil(length(xvec)/2));
                set(gca,'fontsize', 15, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})

                xlabel([ '% of gait-cycle ']);%

                %                 ylim([0 max(plotHead)]);
                pc=pc+1; %plotcounter


                legend(bh, ['Total from average ~ ' num2str(sum(ppantData)) ])
                % pspots=pspots+1;
            end % i LR

        end % nGaits.
        %%
        cd([cfg.figdir filesep  cfg.type ' onset distribution'])
        shg

        print([psubj ' ' cfg.type ' onset distribution binned(' num2str(cfg.usebin) ')'],'-dpng');
    end % ppant
else % plot level GFX


    clf;
    pc=1; % plot counter
    psubj= 'GFX';
    % both this and the next use the same figure function:
    pspots=1:3;

    for nGaits_toPlot=1:2

        if nGaits_toPlot==1
            pidx= cfg.pidx1;
            ftnames= {'LR', 'RL', 'combined'};
        else
            pidx= cfg.pidx2;
            ftnames= {'LRL', 'RLR', 'combined'};
        end

        legp=[]; % for legend
        for iSpeed=1:2
        iLR=3 % focusing on all steps for now.

            % % note that for doubgc, we skip to index 4.
            GFXppantData=[];
            GFXheadData=[];
          
                for ippant= 1:nsubs
                    GFXppantData(ippant,:)= dataIN(ippant,iSpeed,iLR).([usegaitnames{nGaits_toPlot} usefield]);
                    GFXheadData(ippant,:) = GFX_headY(ippant,iSpeed).([usegaitnames{nGaits_toPlot}]);
                end
          



            %x axis:
            if cfg.usebin

                %x axis:          %approx centre point of the binns.
                mdiff = round(mean(diff(pidx)./2));
                xvec = pidx(1:end-1) + mdiff;
            else  % note that if we aren't using the binned versions, use full xaxis.
                %         (overwrite the above)
                xvec = 1:pidx(end);
            end

            %%
                subplot(2,2,pc)
           
            hold on;
            yyaxis left

            % finely sampled bar, each gait "%" point.
            GFXdata = squeeze(mean(GFXppantData,1));
            stE= CousineauSEM(GFXppantData);
            bh=bar(xvec, GFXdata);
            bh.FaceColor = speedCols{iSpeed};
            bh.FaceAlpha = .2;
            legp(iLR)= bh;
            hold on;
            errorbar(xvec, GFXdata, stE, 'LineStyle', 'none', 'Color', 'k');

            ylabel([cfg.type ' onset [counts]']);
            %% add FA if response type:

%             if useFA
%                 bh=bar(xvec, faData);
%                 bh.FaceColor = 'm';
%             end
            %%


            yyaxis right
            mH = squeeze(mean(GFXheadData,1));
            stE = CousineauSEM(GFXheadData);


           mHs= imresize(mH, [1,100]);

           stEs= imresize(stE, [1,100]);
            ph= shadedErrorBar(1:length(mHs), mHs, stEs, [],1);
%             ph=plot(plotHead, ['k-o'], 'linew', 3); hold on
            
            set(gca,'ytick', []);
            if cfg.usebin==0
                title([psubj ' (N ' num2str(sum(GFXppantData(:), 'omitnan')) ') ' ftnames{iLR}], 'interpreter', 'none');
            else
                title([psubj ' ' ftnames{iLR}], 'interpreter', 'none');
            end

            midp=xvec(ceil(length(xvec)/2));
            set(gca,'fontsize', 15, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})

            xlabel([ '% of gait-cycle ']);%
            %                 ylim([0 max(plotHead)]);
            pc=pc+1; %plotcounter

            % pspots=pspots+1;
        end % i LR

    end % nGaits.
    %%
    try cd([cfg.figdir filesep  cfg.type ' onset distribution'])
    catch 
        mkdir([cfg.figdir filesep  cfg.type ' onset distribution'])
    cd([cfg.figdir filesep  cfg.type ' onset distribution']);
    end
    shg

    print([psubj ' ' cfg.type ' onset distribution binned(' num2str(cfg.usebin) ')'],'-dpng');
end % plotlevel



end