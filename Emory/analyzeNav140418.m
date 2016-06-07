% analyzeNav140418.m
%
% uses data from parsed python files (.par)
% these files contain eye, position and velocity data all using the same
% timescale, but with samples of different data types falling at different
% timespoints; the timescale is resampled at 1000 Hz to match sampling rate 
% of neural recordings, and additional data points between samples are 
% filled in using interpolation.
%
% 140418 M. Jutras
%
% this code requires ft_nearest.m, inpaint_nans.m, ft_progress.m

% python samples position/direction at ~83 Hz; eye data @240

% timestamp values in panda log are already in ms

% specify the directory containing log files
sessionDir = 'R:\Buffalo Lab\Virtual Navigation\Free Foraging\Giz\session_1191';

% import panda log file info
[logtime_eye, id_eye, x_eye, y_eye] = textread(fullfile(sessionDir,'eye.par'), '%f%s%f%f'); % eye data
[logtime_ava, id_ava, x_ava, y_ava] = textread(fullfile(sessionDir,'avatar.par'), '%f%s%f%f'); % avatar (position) data


% eye data: timestamps, x and y coordinates in viewing field
time_eye = logtime_eye(strcmp(id_eye,'eye'));
xeye_real = x_eye(strcmp(id_eye,'eye'));
yeye_real = y_eye(strcmp(id_eye,'eye'));

% position data: timestamps, x and y coordinates in virtual environment
time_pos = logtime_ava(strcmp(id_ava,'pos'));
xpos_real = x_ava(strcmp(id_ava,'pos'));
ypos_real = y_ava(strcmp(id_ava,'pos'));

% direction data: timestamps and direction
time_dir = logtime_ava(strcmp(id_ava,'dir'));
dir_real = x_ava(strcmp(id_ava,'dir'));
% % if using interpolation, convert to degrees/radians after interpolating
% dir_degrees = mod(dir_real+90,360);
% dir_radians = dir_degrees * pi/180;

% velocity data
time_vel = logtime_ava(strcmp(id_ava,'speed'));
vel_real = x_ava(strcmp(id_ava,'speed'));

% start and duration of each trial (when bananas appear)
% the last instance of '100' will mark the end of the last trial
log_trl_starttime = logtime_eye(strcmp(id_eye,'100'));
log_trl_length = diff(log_trl_starttime);

%% get all the data into one structure

trlinf = [];
c=1;
ft_progress('init', 'etf',     'Please wait...');
for trllop = 2:length(log_trl_starttime)-1 % skip first trial, position data doesn't start right away
    
    ft_progress(trllop/(length(log_trl_starttime)-1), 'Processing event %d from %d', trllop, (length(log_trl_starttime)-1));
    
    trltim = [log_trl_starttime(trllop) log_trl_starttime(trllop+1)-1];
    
    % buffer the start and end times for each trial 100 samples on each
    % side (50 for position and direction)
    trlstart_eye = ft_nearest(time_eye,trltim(1))-100;
    trlend_eye = ft_nearest(time_eye,trltim(2))+100;
    
    trlstart_pos = ft_nearest(time_pos,trltim(1))-50;
    trlend_pos = ft_nearest(time_pos,trltim(2))+50;
    
    trlstart_dir = ft_nearest(time_dir,trltim(1))-50;
    trlend_dir = ft_nearest(time_dir,trltim(2))+50;
    
    trlstart_vel = ft_nearest(time_vel,trltim(1))-50;
    trlend_vel = ft_nearest(time_vel,trltim(2))+50;
    
    trlinf{c}.eyetime = time_eye(trlstart_eye:trlend_eye);
    trlinf{c}.postime = time_pos(trlstart_pos:trlend_pos);
    trlinf{c}.dirtime = time_dir(trlstart_dir:trlend_dir);
    trlinf{c}.veltime = time_vel(trlstart_vel:trlend_vel);
    trlinf{c}.eyedat = [xeye_real(trlstart_eye:trlend_eye)'; yeye_real(trlstart_eye:trlend_eye)'];
    trlinf{c}.posdat = [xpos_real(trlstart_pos:trlend_pos)'; ypos_real(trlstart_pos:trlend_pos)'];
    trlinf{c}.dirdat = dir_real(trlstart_dir:trlend_dir)';
    trlinf{c}.veldat = vel_real(trlstart_vel:trlend_vel)';
    trlinf{c}.trltim = trltim;

    c=c+1;
   
end
ft_progress('close')

clear logtime_* id_eye x_eye y_eye id_ava x_ava y_ava
clear time_* xeye_real yeye_real xpos_real ypos_real dir_real vel_real



%% resample to 1000 Hz and interpolate data points between samples

% the data structure created here is easily merged with the Fieldtrip
% (neural data analysis toolbox) data structure

clear data

data.time = [];
data.eyedat = [];
data.posdat = [];
data.dirdat = [];
data.veldat = [];

ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(trlinf)
    
    ft_progress(trllop/length(trlinf), 'Processing event %d from %d', trllop, length(trlinf));

    data.time{trllop} = trlinf{trllop}.trltim(1):trlinf{trllop}.trltim(2);
        
    % interpolate eye data
    eyetime_interp = trlinf{trllop}.eyetime(1):trlinf{trllop}.eyetime(end);
    eye_interp = nan(2,length(eyetime_interp));
    for k=1:length(eyetime_interp)
        if ismember(eyetime_interp(k),trlinf{trllop}.eyetime)
            eye_interp(:,k) = trlinf{trllop}.eyedat(:,find(trlinf{trllop}.eyetime==eyetime_interp(k),1,'last'));
        end
    end
    eye_interp=inpaint_nans(eye_interp,2);

    % interpolate position data
    postime_interp = trlinf{trllop}.postime(1):trlinf{trllop}.postime(end);
    pos_interp = nan(2,length(postime_interp));
    for k=1:length(postime_interp)
        if ismember(postime_interp(k),trlinf{trllop}.postime)
            pos_interp(:,k) = trlinf{trllop}.posdat(:,find(trlinf{trllop}.postime==postime_interp(k),1,'last'));
        end
    end
    pos_interp=inpaint_nans(pos_interp,2);

    % interpolate direction data
    dirtime_interp = trlinf{trllop}.dirtime(1):trlinf{trllop}.dirtime(end);
    dir_interp = nan(1,length(dirtime_interp));
    for k=1:length(dirtime_interp)
        if ismember(dirtime_interp(k),trlinf{trllop}.dirtime)
            dir_interp(1,k) = trlinf{trllop}.dirdat(1,find(trlinf{trllop}.dirtime==dirtime_interp(k),1,'last'));
        end
    end
    dir_interp=inpaint_nans(dir_interp,2);

    % interpolate velocity data
    veltime_interp = trlinf{trllop}.veltime(1):trlinf{trllop}.veltime(end);
    vel_interp = nan(1,length(veltime_interp));
    for k=1:length(veltime_interp)
        if ismember(veltime_interp(k),trlinf{trllop}.veltime)
            vel_interp(1,k) = trlinf{trllop}.veldat(1,find(trlinf{trllop}.veltime==veltime_interp(k),1,'last'));
        end
    end
    vel_interp=inpaint_nans(vel_interp,2);

    
    data.eyedat{trllop} = eye_interp(:,find(eyetime_interp==trlinf{trllop}.trltim(1)):find(eyetime_interp==trlinf{trllop}.trltim(2)));
    data.posdat{trllop} = pos_interp(:,find(postime_interp==trlinf{trllop}.trltim(1)):find(postime_interp==trlinf{trllop}.trltim(2)));
    data.dirdat{trllop} = dir_interp(find(dirtime_interp==trlinf{trllop}.trltim(1)):find(dirtime_interp==trlinf{trllop}.trltim(2)));
    data.veldat{trllop} = vel_interp(find(veltime_interp==trlinf{trllop}.trltim(1)):find(veltime_interp==trlinf{trllop}.trltim(2)));
    
    % additional step converts direction data to angle in radians
    data.dirdat{trllop} = mod(data.dirdat{trllop} + 90, 360) * pi/180;
    
end
ft_progress('close')


%% save the new data structure

save(fullfile(sessionDir,'navData.mat'),'data')

%% create quiver plots to visualize each trial

trlsel = 1;

[u,v] = pol2cart(data.dirdat{trlsel},data.veldat{trlsel});

figure
quiver(data.posdat{trlsel}(1,1:100:end),data.posdat{trlsel}(2,1:100:end),u(1:100:end),v(1:100:end))
xlim([-5 5])
ylim([-5 5])
title(['Trial ' num2str(trlsel)])


%%

trlsel = 1;

clear M1
figure
c=1;
for timlop = 1:50:length(data.time{trlsel})

    quiver(data.posdat{trlsel}(1,timlop),data.posdat{trlsel}(2,timlop),u(timlop),v(timlop))
    xlim([-5 5]);ylim([-5 5])
    
    M1(c)=getframe;
    
    c=c+1;
    
end

movie2avi(M1,fullfile(sessionDir,['trl' num2str(trlsel) '_nav.avi']),'FPS',20)

