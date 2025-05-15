%adjustGait_subjlist_vAUD
%this script simply holds the hard-coded adjustments per individual

usetype=1; % default
% usetype=1% details below. catches reduced step lengths for shorter

% participants.
if strcmp(subjID, 'xx')
    %do something
    usetype
end

%% change search params:
switch usetype
    case 1
        pkdist=15; %samples between
        pkheight=.005;
    case 2
        
        pkheight = .02;
    case 3
          pkdist = 15;
        pkheight = .015;

end

