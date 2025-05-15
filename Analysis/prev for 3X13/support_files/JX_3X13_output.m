%JX_3X13_output

% this loads the group data saved in GFX_tables, created within 
%jobs_calculatePpantPsychoemtric. 


% it plots some basic output for the assingment

%Fig 1: 
% overall Acc, RT, HitR, FAR, dp, criterion (slow and natural speeds).


%Fig 2: psychometric fits (on average)

%Fig 3: *another script -> DVs over the stride-cycle.


%%
cd('/Users/matthewdavidson/Documents/MATLAB/Auditory Detection while Walking/Processed_Data');
load('GFX_tables.mat');

%%

% Figure 1.
speedstoPlot= {'All', 'Slow', 'Natural'};
DVstoPlot ={'Accuracy', 'RT', 'HitRate', 'FAR', 'dprime', 'criterion'};

ipanel=1;
clf


for iDV = 1:length(DVstoPlot)
       
        data1 = allppantsData.([DVstoPlot{iDV} '_Slow']);
        data2 = allppantsData.([DVstoPlot{iDV} '_Natural']);
        
        barD = [data1, data2];
        
        mB= nanmean(barD,1);
%         stE = CousineauSEM(barD);
        stE = std(barD)./sqrt(size(barD,1));
        
        subplot(2,3,ipanel)
        bh=bar(1:2, mB);
        bh.FaceColor='flat';
        bh.CData(1,:)  =[0,0,1];
        bh.CData(2,:)  =[1,0,0];
        bh.FaceAlpha=.4;
        hold on;
        % add ppantmarkers
        %jitter xpos?
        
        
        scatter(repmat(1,nppants), data1,40,'k');
        scatter(repmat(2,nppants), data2,40,'k');
        
        
        errorbar(1:2, mB,stE,'linestyle', 'none', 'color', 'k','linew',2);
        ylabel(DVstoPlot{iDV});
        
        nppants= length(data1);
        
        %adjust range,
        minY = min(barD(:));
        maxY = max(barD(:));
        mnY = mean(barD(:));
        yR = maxY-minY;
        ylim([mnY - yR , mnY+yR]);
        xlabel('Walk speed');
        set(gca,'XTickLabel', {'Slow', 'Natural'})
        
        
        % ttest
        [h,p,ci,stat] = ttest(data1,data2);
        plot([1,2],[mnY+yR*.65 mnY+yR*.65], 'k-') 
        if p<=.05
            
            text(1.5, mnY+yR*.75, ['p = ' sprintf('%.3f',p)], 'HorizontalAlignment', 'center','fontsize',15);
            
            % disp output to write up.
%             disp(['t(' stat.df
        else
            text(1.5, mnY+yR*.75, ['ns'], 'HorizontalAlignment', 'center','fontsize',15);
        end
            
        set(gca,'fontsize', 15)
        
        ipanel=ipanel+1;
end

set(gcf,'color','w');
    
        