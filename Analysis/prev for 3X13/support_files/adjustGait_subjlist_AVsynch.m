%adjustGait_subjlist_AVsynch
%this script simply holds the hard-coded adjustments per individual

usetype=1; % default 
% usetype=1% details below. catches reduced step lengths for shorter

% participants.
% if strcmp(subjID, 'AL')
%     usetype=1;
% else % back to original:
%     
%     pkdist = ceil(0.17*Fs);
%     pkheight = .02;
%     
% end

%% change search params:

if usetype==1 %(default
    pkdist=15;
    pkheight = 0.005; % (m)
elseif usetype==2
    
    pkduration = 0.4; % min seconds between pks/troughs
    pkdist = ceil(pkduration*Fs); % 400ms. (36 samps).
    
    pkheight= .01;
end