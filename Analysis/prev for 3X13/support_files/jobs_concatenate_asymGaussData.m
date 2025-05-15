% jobs_concatenate_asymGaussData




%specify save directory:
% savedatadir='/MATLAB Drive/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

allppants = dir([pwd filesep 'p_*data.mat']);


% useppants= [1:16,18:19]; % 17 and 20 have poor gaussian fits.
useppants=1:length(allppants);
nsubs= length(useppants);

%%

GFX_table= table();
% GFX_table

GFX_asymgauss_Cond_x_Yfit= zeros(nsubs, 3, 7); % dims: subjects, cond(all,slow,fast), sample.
GFX_asymgauss_Cond_xparam=zeros(nsubs, 3, 3); % last dim: mean and LHS , RHS width of gaussian fit.
% same per gait

GFX_asymgauss_Cond_x_Gait_x_Yfit= [];%zeros(nsubs, 3,4, 4,1000); % extra dims for cond, alignment, and gait quantile
GFX_asymgauss_Cond_x_Gait_param=[];%zeros(nsubs, 3,4,4,2); % last dim is params (mu sigma)  




for ippant= 1:length(useppants)

pINDX= useppants(ippant);
    cd(savedatadir);
    load(allppants(pINDX).name, 'participant_asymGauss_Fits_bySpeed','participant_asymGauss_Fits_bySpeedbyAlignbyGait');


    % GFX_table.propSame_all(ippant,:) = propSame_all;
    disp(['concatenating participant ' num2str(pINDX)])

    for icond= 1:3 % all, slow, fast:
        %first store the interpolated gaussian for plotting:
        GFX_asymgauss_Cond_x_Yfit(ippant,icond,:) = participant_asymGauss_Fits_bySpeed(icond).gaussData_Y; %

        % % % now store the parameters for specific between condition stats:
        GFX_asymgauss_Cond_xparam(ippant,icond,1)= participant_asymGauss_Fits_bySpeed(icond).gaussfit(1); % b is coefficient for mean.
        GFX_asymgauss_Cond_xparam(ippant,icond,2)= participant_asymGauss_Fits_bySpeed(icond).gaussfit(2); % c is coefficient for SD (LHS).
        GFX_asymgauss_Cond_xparam(ippant,icond,3)= participant_asymGauss_Fits_bySpeed(icond).gaussfit(3); % c is coefficient for SD (LHS).

    end


   for ispeed=1:size(participant_asymGauss_Fits_bySpeedbyAlignbyGait,1)
       for ialignment= 1:size(participant_asymGauss_Fits_bySpeedbyAlignbyGait,2)
           for igait= 1:size(participant_asymGauss_Fits_bySpeedbyAlignbyGait,3)
              
GFX_asymgauss_Cond_x_Gait_x_Yfit(ippant, ispeed,ialignment,igait,:) = participant_asymGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussData_Y;  

GFX_asymgauss_Cond_x_Gait_param(ippant, ispeed, ialignment, igait,1)= participant_asymGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussfit(1);  
GFX_asymgauss_Cond_x_Gait_param(ippant, ispeed, ialignment, igait,2)= participant_asymGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussfit(2);
GFX_asymgauss_Cond_x_Gait_param(ippant, ispeed, ialignment, igait,3)= participant_asymGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussfit(3);


           end
       end
   end

end % per participant


%% also store the X axis for plots.
gauss_xvec = participant_asymGauss_Fits_bySpeed(1).gaussData_X; % same for all conds.

% save
save('GFX_asymGaussData',...
    'GFX_asymgauss_Cond_x_Yfit', 'GFX_asymgauss_Cond_xparam',...
    'GFX_asymgauss_Cond_x_Gait_x_Yfit','GFX_asymgauss_Cond_x_Gait_param',...
    'gauss_xvec');