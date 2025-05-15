% jobs_concatenate_GaussData




%specify save directory:
% savedatadir='/MATLAB Drive/AV Synchrony Exp/Processed_Data';
cd(savedatadir);

allppants = dir([pwd filesep 'p_*data.mat']);


% useppants= [1:16,18:19]; % 17 and 20 have poor gaussian fits.
useppants=1:length(allppants)
nsubs= length(useppants);

%%

GFX_table= table();
% GFX_table

GFX_gauss_Cond_x_Yfit= zeros(nsubs, 3, 1000); % dims: subjects, cond(all,slow,fast), sample.
GFX_gauss_Cond_xparam=zeros(nsubs, 3, 2); % last dim: mean and width of gaussian fit.(and amp)
% same per gait

GFX_gauss_Cond_x_Gait_x_Yfit= [];%zeros(nsubs, 3,4, 4,1000); % extra dims for cond, alignment, and gait quantile
GFX_gauss_Cond_x_Gait_param=[];%zeros(nsubs, 3,4,4,2); % last dim is params (mu sigma)  




for ippant= 1:length(useppants)

pINDX= useppants(ippant);
    cd(savedatadir);
    load(allppants(pINDX).name, 'participantGauss_Fits_bySpeed','participantGauss_Fits_bySpeedbyAlignbyGait');


    % GFX_table.propSame_all(ippant,:) = propSame_all;
    disp(['concatenating participant ' num2str(pINDX)])

    for icond= 1:3 % all, slow, fast:
        %first store the interpolated gaussian for plotting:
        GFX_gauss_Cond_x_Yfit(ippant,icond,:) = participantGauss_Fits_bySpeed(icond).gaussData_Y; %

        % % % now store the parameters for specific between condition stats:
        GFX_gauss_Cond_xparam(ippant,icond,1)= participantGauss_Fits_bySpeed(icond).gaussfit(1); % b is coefficient for mean.
        GFX_gauss_Cond_xparam(ippant,icond,2)= participantGauss_Fits_bySpeed(icond).gaussfit(2); % c is coefficient for SD.
%         GFX_gauss_Cond_xparam(ippant,icond,3)= participantGauss_Fits_bySpeed(icond).gaussfit.a; % a is coefficient for gauss amp.
    end


   for ispeed=1:size(participantGauss_Fits_bySpeedbyAlignbyGait,1)
       for ialignment= 1:size(participantGauss_Fits_bySpeedbyAlignbyGait,2)
           for igait= 1:size(participantGauss_Fits_bySpeedbyAlignbyGait,3)
GFX_gauss_Cond_x_Gait_x_Yfit(ippant, ispeed,ialignment,igait,:) = participantGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussData_Y;  
GFX_gauss_Cond_x_Gait_param(ippant, ispeed, ialignment, igait,1)= participantGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussfit(1);  
GFX_gauss_Cond_x_Gait_param(ippant, ispeed, ialignment, igait,2)= participantGauss_Fits_bySpeedbyAlignbyGait(ispeed, ialignment, igait,:).gaussfit(2);


           end
       end
   end

end % per participant


%% also store the X axis for plots.
gauss_xvec = participantGauss_Fits_bySpeed(1).gaussData_X; % same for all conds.

% save
save('GFX_GaussData',...
    'GFX_gauss_Cond_x_Yfit', 'GFX_gauss_Cond_xparam',...
    'GFX_gauss_Cond_x_Gait_x_Yfit','GFX_gauss_Cond_x_Gait_param',...
    "gauss_xvec");