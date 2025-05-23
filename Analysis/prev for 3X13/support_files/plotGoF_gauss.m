function statsOUT =plotGoF_gauss(plotData, plot_xloc)
%compare linear and quadratic GoF statistics:

%
%% test regression statistics.

%convert to table:
nquintiles=size(plotData,2);
dataR = repmat([1:nquintiles], [size(plotData,1),1]); % quintile labels.
subjectID = repmat([1:size(plotData,1)]', [1,size(plotData,2)]); 
%%
testd = plotData; % alpha power.
newT = [subjectID(:), dataR(:), testd(:)];

tbA = array2table(newT);
tbA.Properties.VariableNames = {'subjectID', 'GaitQ', 'data'};

% perform stepwise regression.
% starting with a constant term, linear, interaction, then quadratic.
%%
% starting with constant term, , see if model fit improves with linear,
%
%%
UNterm1='data ~  1+(1|subjectID)'; % constant, with random intercepts per subj
UNterm2='data ~  GaitQ + (1|subjectID)'; % as above, + fixed effect of alpha

UNterm3='data ~  GaitQ^2 + (1|subjectID)';  % include quadratic term
% UNterm3='data ~  -Alpha + Alpha^2 + (1|subjectID) '; % as above, but pure quadratic effect (pure) 


baseMDL = fitlme(tbA, UNterm1);
linearMDL = fitlme(tbA, UNterm2);
quadMDL = fitlme(tbA, UNterm3); %'FitMethod', 'REML') - more conservative).

%%
%coeff beta estimates:
% e.g. 
%linearMDL.
c1=compare(baseMDL, linearMDL);
c2=compare(baseMDL, quadMDL);
c3=compare(linearMDL, quadMDL);

% for ease of writing up, summarise the stats for output
statsOUT=[];
statsOUT.compBase_v_Linear=c1;
statsOUT.compBase_v_Quad=c2;
statsOUT.compLinear_v_Quad=c3;

statsOUT.MDLs.baseMDL=baseMDL;
statsOUT.MDLs.linearMDL=linearMDL;
statsOUT.MDLs.quadMDL=quadMDL;

%% display result in comm window.

    LRstat = statsOUT.compBase_v_Linear.LRStat(2);
    pval = statsOUT.compBase_v_Linear.pValue(2);
    
    LRstat2 = statsOUT.compBase_v_Quad.LRStat(2);
    pval2 = statsOUT.compBase_v_Quad.pValue(2);
    
    LRstatq= statsOUT.compLinear_v_Quad.LRStat(2);  
    pvalq= statsOUT.compLinear_v_Quad.pValue(2);

disp('>>>>>>')
disp(['>>>>>> base vs linear LRtest = ' num2str(LRstat) ',p=' num2str(pval)]); 
disp(['>>>>>> base vs quad LRtest = ' num2str(LRstat2) ',p=' num2str(pval2)]); 
disp(['>>>>>> linear vs quadratic LRtest = ' num2str(LRstatq) ',p=' num2str(pvalq)]); 


%% add to plot the best fit (if sig)


mD= squeeze(nanmean(plotData,1));
[F1,G1] = fit([plot_xloc]', mD', 'poly1');
[F2,G2] = fit([plot_xloc]', mD', 'poly2');
%
%%
%Goodness of fits:
goodnessfitsare = [G1.adjrsquare, G2.adjrsquare];


%Output of Fits
Fitsare  = {F1, F2};%

%whichever was a better fit, now select that to plot on top of figure:
[~,polyu] = max(goodnessfitsare);
[~, polyu]= min([pval,pval2]);
% lgprint = ['best fit, adj. R^2 =' sprintf('%.3f',goodnessfitsare(polyu))];

%fit and evaluate best polynomial:
Fplot = Fitsare{polyu};

% return % if not plotting


canPlot=0;
% plot only if sig:
if polyu==1 && pval<.05 % sig linear fit compared to basic model.
    canPlot=1;
elseif polyu==2 && pval2 <.05;
    canPlot=1;
end

if pval<0.05 ||  pval2<0.05
    canPlot=1;
end
if canPlot
    
    %%
    hold on
    % hold on; plot best fit:
    pl=  plot(Fplot, plot_xloc, mD);
    pl(2).LineWidth = 3;
    pl(2).Color = 'k';
    % legend([pl(2)], { lgprint}, 'location', 'SouthEast');
    
    
    % %%% run sig tests?
    % [p, ~] = polyfit([1:size(mD,2)]', mD', 1); % linear
    % bestparam(1) = p(1);
    % [p, ~] = polyfit([1:size(mD,2)]', mD', 2); % quad
    % bestparam(2) = p(1);
    
mD= squeeze(nanmean(plotData,1));
[F1,G1] = fit(plot_xloc', mD', 'poly1');
[F2,G2] = fit(plot_xloc', mD', 'poly2');
%
%%
%Goodness of fits:
goodnessfitsare = [G1.adjrsquare, G2.adjrsquare];


%Output of Fits
Fitsare  = {F1, F2};%

%whichever was a better fit, now select that to plot on top of figure:
[~,polyu] = max(goodnessfitsare);

% lgprint = ['best fit, adj. R^2 =' sprintf('%.3f',goodnessfitsare(polyu))];

%fit and evaluate best polynomial:
Fplot = Fitsare{polyu};
% 
% hold on;
pl=  plot(Fplot, plot_xloc, mD);
pl(2).LineWidth = 3;
pl(2).Color = 'k';
% legend([pl(2)], { lgprint}, 'location', 'SouthEast');


% %%% run sig tests?
% [p, ~] = polyfit([1:size(mD,2)]', mD', 1); % linear
% bestparam(1) = p(1);
% [p, ~] = polyfit([1:size(mD,2)]', mD', 2); % quad
% bestparam(2) = p(1);

legend off
end


end