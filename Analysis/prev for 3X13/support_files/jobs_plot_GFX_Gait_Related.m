% jobs_plot_GFX_Gait_Related

% -- called from JOBS_import_preprocess;


%% THIS SECTION UPDATES FILE PATH FOR DESKTOP OR ONLINE VERSION

MatlabOnline = 0

if MatlabOnline == 0
    rootPath = 'C:/Users/mobil/Gabriel/MATLAB_Drive';
else
    rootPath = 'MATLAB DRIVE'
end

%% PLOT GFX GAIT RELATED BEGINS HERE

% this script loads the presaved participant level data (SOA variables).
% makes a plot per participant, comparing slow and fast walk conditions.

% - load data
% - extract trials per SOA type
% - calculate proportion same
% - save

%specify save directory:
savedatadir = string(rootPath)+'/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

% groupdata = dir([pwd filesep 'p_*.mat']);

figdir = string(rootPath)+'/AV Synchrony Exp/Figures';

%%
% load group level data:

load('GFX_SOAdata.mat');



nGaitPhases = 4;
nSOAs = length(SOAs);



%%  Calculate Group Means of PropSame by Gait Phase per SOA


speed = {'all','slow','fast'};

for icond = 1:length(speed);
    
    % ALL SPEEDS
    
    % Create an array to hold Group Mean Data in
    GFX_Means_Propsame_GaitPhase_SOA = nan(4,7);
    
    
    
    %Load the GFX Data for This Speed condition
    if icond == 1
        GFX_propSamebyGaitandSOA = GFX_propSamebyGaitandSOA_all;
    elseif icond == 2
        GFX_propSamebyGaitandSOA = GFX_propSamebyGaitandSOA_slow;
    elseif icond == 3
        GFX_propSamebyGaitandSOA = GFX_propSamebyGaitandSOA_fast;
    end
    
    
    for iGate = 1:nGaitPhases
        
        disp('Processing Gait ' + string(iGate));
        
        for iSOA = 1:nSOAs
            
            GFX_Means_Propsame_GaitPhase_SOA(iGate, iSOA) = mean(GFX_propSamebyGaitandSOA(:,iGate,iSOA), 'omitnan');
            %meanThisGaitandSOA = mean(GFX_propSamebyGaitandSOA(:,iGate,iSOA));      
            
%             plot(1:4, GFX_Means_Propsame_GaitPhase_SOA, styleThisLine);
            
        end
        
        % Fit Gaussian
        
         
        % perform fit :
%     fit1 = fit(X', yData', gaussEqn);
        
    end % End of per Gait Phase Loop
    
    if icond == 1
        GFX_Means_Propsame_GaitPhase_SOA_all = GFX_Means_Propsame_GaitPhase_SOA
        
    elseif icond == 2
        GFX_Means_Propsame_GaitPhase_SOA_slow = GFX_Means_Propsame_GaitPhase_SOA
        
    elseif icond == 3
        GFX_Means_Propsame_GaitPhase_SOA_fast = GFX_Means_Propsame_GaitPhase_SOA
    end
    
    
    %% Group Means of Synchrony Functions per Gait Cycle
    %%
    % for iVariable = 1:nVariables
    %
    %     ThisVariable = FocusVariables(iVariable);
    %     disp(ThisVariable);
    %     GroupRTData(iVariable,1)
    %     VariableLabels(iVariable) = ThisVariable;
    %
    %
    %     for iSOA = 1:nSOAs
    %
    %         GroupMeanThisVariableThisSOA = mean(ThisVariable(:,iSOA));
    %
    %         GroupMeanSame_all(2,iSOA) = GroupMeanSame_allThisSOA;
    %         GroupMeanSame_all(1,iSOA) = SOAs(iSOA);
    %
    %     end
    
    %%  Plot Group Level Synchrony Curves by Gait Phase and Fit Gaussians
    

        %% For Regular Gaussians
           %%  Plot Group Level Synchrony Curves by Gait Phase and Fit Gaussians
    
    % Regular Gaussian
    gaussEqn = 'a*exp(-((x-b)^2/(2*c^2)))+d';   % Where a = amplitude, b= mean, c= Standard Deviation, and d is v0 (the vertical offset)
    
   
    
    nvar = 4 % Variables a, b, c and d
    
    gaussians = nan(4,4)
    gaussians = array2table(gaussians, 'RowNames', {'GaitPhase 1', 'GaitPhase 2',...
        'GaitPhase 3', 'GaitPhase 4'}, 'VariableNames', {'Amplitude', 'Mean', 'StdDev', 'v0'});
    
    clf
    for ig=1:nGaitPhases
        
        
        %Plot original (unnormalised) Synch Curves
        
        plot(SOAs, GFX_Means_Propsame_GaitPhase_SOA(ig,:));
        
        hold on
        
        
        %FIT GAUSSIAN
        
        %create a handle to a function called gaus that takes input parameters:
        % x, array of x values
        % mu, mean
        % sig, SD
        % amp, pos or negative curve
        % vo, vertical offset.
        
        % Regular Gaussian
        gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
                
       % Vector for this SOA
        y = GFX_Means_Propsame_GaitPhase_SOA(ig,:) 
       
        
        %Plot Normalised Synch Functions
        plot(SOAs, y);
        
        %Take transpose of vector to
        y = y'            % Fit function requires a column vector argument.
        
        % Fit Gaussian
        startPoint = [0.1 0.85 0.9 0.2] %startPoint = [mean, amplitude, standard deviation, v0]
        [fitThisGaitPhase, gof] = fit(SOAs, y, gaussEqn, 'StartPoint', startPoint); %gof is goodness of fit
                
        %Copy details of Gaussians to the gaussian table
        gaussians.Amplitude(ig) = fitThisGaitPhase.a
        gaussians.Mean(ig) = fitThisGaitPhase.b
        gaussians.StdDev(ig) = fitThisGaitPhase.c
        gaussians.v0(ig) = fitThisGaitPhase.d

        
        
        %plot the results of the fit on top of the data:
        hold on;
        
        fit_x = SOAs
        
        fit_amp= fitThisGaitPhase.a;

        fit_mu= fitThisGaitPhase.b;

        fit_sig= fitThisGaitPhase.c;

        fit_vo= fitThisGaitPhase.d;

 %plot this gaussian:
        yf = gaus(fit_x, fit_mu, fit_sig, fit_amp, fit_vo); % This envokes the function handle above.
        plot(SOAs, yf, 'k-', 'LineWidth',1 )
        
end

        
        %% For Normalised Gaussians
        
%         %%  Plot Group Level Synchrony Curves by Gait Phase and Fit Gaussians
%     
%     % Regular Gaussian
%     %gaussEqn = 'a*exp(-((x-b)^2/(2*c^2)))+d';   % Where a = amplitude, b= mean, c= Standard Deviation, and d is v0 (the vertical offset)
%     
%     % Normalised Gaussian
%     gaussEqn = 'exp(-((x-b)^2/(2*c^2)))';   % Where amplitude is 1, b= mean, c= Standard Deviation, and d is v0 (the vertical offset) = 0)
%     
%     nvar = 2 % Variables a, b, c and d
%     
%     gaussians = nan(4,2)
%     gaussians = array2table(gaussians, 'RowNames', {'GaitPhase 1', 'GaitPhase 2',...
%         'GaitPhase 3', 'GaitPhase 4'}, 'VariableNames', {'Mean', 'StdDev'});
%     
%     clf
%     for ig=1:nGaitPhases
%         
%         %%Plot original (unnormalised) Synch Curves
%         %plot(SOAs, GFX_Means_Propsame_GaitPhase_SOA(ig,:));
%         
%         hold on
%         
%         
%         %FIT GAUSSIAN
%         
%         %create a handle to a function called gaus that takes input parameters:
%         % x, array of x values
%         % mu, mean
%         % sig, SD
%         % amp, pos or negative curve
%         % vo, vertical offset.
%         
%         % Regular Gaussian
%         %gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
%         
%         
%         %Gaussian for normalised data (amp =1, vo = 0, simplifies things to
%         %2 variables)
%         gaus = @(x,mu,sig)exp(-(((x-mu).^2)/(2*sig.^2)));
%          
%         
%        % Vector for this SOA
%         y = GFX_Means_Propsame_GaitPhase_SOA(ig,:) 
%         
%         yMax = max(y);
%         
%         % Normalise 
%         yNorm = y / yMax
%         
%         %Plot Normalised Synch Functions
%         plot(SOAs, yNorm);
%         
%         %Take transpose of vector to
%         yNorm = yNorm'            % Fit function requires a column vector argument.
%         
%         % Fit Gaussian
%         fitThisGaitPhase = fit(SOAs, yNorm, gaussEqn);
%                 
%         %Copy details of Gaussians to the gaussian table
%         
%         gaussians.Mean(ig) = fitThisGaitPhase.b
%         gaussians.StdDev(ig) = fitThisGaitPhase.c
%         
% 
%         
%         
%         hold on;
%         
%         fit_x = SOAs
%         
% 
%         fit_mu= fitThisGaitPhase.b;
% 
%         fit_sig= fitThisGaitPhase.c;
% 
%         
%         %plot the results of the fit on top of the data:
%         hold on;
%         
%         fit_x = SOAs
%         
%         %fit_amp= fitThisGaitPhase.a;
% 
%         fit_mu= fitThisGaitPhase.b;
% 
%         fit_sig= fitThisGaitPhase.c;
% 
%         %fit_vo= fitThisGaitPhase.d;
% 
%         %plot this gaussian:
%         yf = gaus(fit_x, fit_mu, fit_sig); % This envokes the function handle above.
%         plot(SOAs, yf, 'k-', 'LineWidth',1 )
%         
%          end

 %%       
     
   
    
    legend('q1','q2','q3','q4');
    xlabel('SOA');
    ylabel('Prop. "Same" responses');
    
    thisspeedCond = speed(icond);
    titletextA = 'Group Level Synchrony Curves by Gait Phase. Walkspeed: ' + string(thisspeedCond);
    title(titletextA)
    shg
       
    cd(figdir);
    cd('SynchFunctionperGait')
    filenameSOAgPhase = 'GFX_SynchCurves_per_GaitPhase_' + string(speed(icond));
    print(filenameSOAgPhase, '-dpng');
 
  hold off    
   
  %Save Gaussians info
  
      if icond == 1
        gaussians_all = gaussians;
        %save('GFX_SOAdata.mat', 'gaussians_all',);
        
    elseif icond == 2
        gaussians_slow = gaussians;
        %save('GFX_SOAdata.mat', 'gaussians_slow');
        
    elseif icond == 3
        gaussians_fast = gaussians;
        %save('GFX_SOAdata.mat', 'gaussians_fast');    
      end   %End of If statements selecting which walk speed to save gaussians under.
 
    %% Plot Group Level PropSame by Gait Cycle per SOA
    
    clf
    for isoa=1:nSOAs
    
       % For Regular Gaussians
       plot(1:4, GFX_Means_Propsame_GaitPhase_SOA(:,isoa));
       
       % For Normalised
       %plot(1:4, GFX_Means_Propsame_GaitPhase_SOA(:,isoa));
       
       hold on
    
    end
    
    legend(num2str(SOAs));

    titletextB = 'Group Level PropSame by Gait Cycle per SOA. Walkspeed: ' + string(thisspeedCond);
    title(titletextB);
    
    xlabel('Gait Phase');
    ylabel('Prop. "Same" responses');

    shg
    
    cd(figdir);
    cd('PropSamebyGaitperSOA');
    
    filenamePhaseSOAg = 'GFX_Mean_by_GaitPhase_perSOA' + string(speed(icond));
    print(filenamePhaseSOAg, '-dpng');
    
    
end % End of Speed condition Loop