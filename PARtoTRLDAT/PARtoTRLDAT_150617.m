function trldat = PARtoTRLDAT(sessionDir,runpar)

% imports data from parsed python files (.par)
% these files contain eye, position and velocity data all using the same
% timescale, but with samples of different data types falling at different
% timespoints; the timescale is resampled at 1 kHz to match sampling rate 
% of neural recordings, and additional data points between samples are 
% filled in using interpolation.
%
% also imports data related to target (fruit) placement and acquisition
% during each trial, including alpha (transparency) levels
%
% use the trldat structure to align the neural data using BRtoFT
%
% 150421 Mike Jutras
%
% NOTE:
% before running, must run Python script MakeParFiles in the directory containing 
% the Python log files
% (e.g. R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\JN_14_08_25_13_12)
%
% CHANGE LOG:
% 150617 - MJJ - Removed savDir, since I reorganized the directories where
% Matlab code/.MAT files were stored; trldat files now saved to the same
% directory where the Python log files are stored, and can be moved after
% they are created.
%

if nargin==1
    runpar = 0; % if set to 1, run interpolation using parallel processing
end

logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\';


% import panda log file info
disp('Loading avatar & eye PAR files')

[logtime_eye, id_eye, x_eye, y_eye] = textread(fullfile(logDir,sessionDir,'eyepos.par'), '%f%s%f%f'); % eye data
[logtime_ava, id_ava, x_ava, y_ava] = textread(fullfile(logDir,sessionDir,'avatar.par'), '%f%s%f%f'); % avatar (position) data

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

% velocity data
time_vel = logtime_ava(strcmp(id_ava,'speed'));
vel_real = x_ava(strcmp(id_ava,'speed'));

% start of each trial
log_trl_starttime = logtime_eye(strcmp(id_eye,'NewTrial'));


% get all the data into one structure
disp('Grouping data into trials')

ft_defaults

trlinf = [];
trltimarr = [];
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(log_trl_starttime)-1 % don't include the last trial, usually gets cut off
    
    ft_progress(trllop/(length(log_trl_starttime)-1), 'Processing trial %d from %d', trllop, (length(log_trl_starttime)-1));
    
    % trial time according to panda log file
    trltim = [log_trl_starttime(trllop) log_trl_starttime(trllop+1)-1];

    % eye data
    if ~isempty(find(time_eye<trltim(1),1,'last'))
        trlstart_eye = find(time_eye<trltim(1),1,'last');
    else
        trlstart_eye = ft_nearest(time_eye,trltim(1));
    end
    if ~isempty(find(time_eye>trltim(2),1,'first'))
        trlend_eye = find(time_eye>trltim(2),1,'first');
    else
        trlend_eye = ft_nearest(time_eye,trltim(2));
    end

    % position data
    if ~isempty(find(time_pos<trltim(1),1,'last'))
        trlstart_pos = find(time_pos<trltim(1),1,'last');
    else
        trlstart_pos = ft_nearest(time_pos,trltim(1));
    end
    if ~isempty(find(time_pos>trltim(2),1,'first'))
        trlend_pos = find(time_pos>trltim(2),1,'first');
    else
        trlend_pos = ft_nearest(time_pos,trltim(2));
    end

    % direction data
    if ~isempty(find(time_dir<trltim(1),1,'last'))
        trlstart_dir = find(time_dir<trltim(1),1,'last');
    else
        trlstart_dir = ft_nearest(time_dir,trltim(1));
    end
    if ~isempty(find(time_dir>trltim(2),1,'first'))
        trlend_dir = find(time_dir>trltim(2),1,'first');
    else
        trlend_dir = ft_nearest(time_dir,trltim(2));
    end

    % velocity data
    if ~isempty(find(time_vel<trltim(1),1,'last'))
        trlstart_vel = find(time_vel<trltim(1),1,'last');
    else
        trlstart_vel = ft_nearest(time_vel,trltim(1));
    end
    if ~isempty(find(time_vel>trltim(2),1,'first'))
        trlend_vel = find(time_vel>trltim(2),1,'first');
    else
        trlend_vel = ft_nearest(time_vel,trltim(2));
    end
    
    trltimarr(trllop,:,:) = [trltim; ...
        time_eye(trlstart_eye) time_eye(trlend_eye); ...
        time_pos(trlstart_pos) time_pos(trlend_pos); ...
        time_dir(trlstart_dir) time_dir(trlend_dir); ...
        time_vel(trlstart_vel) time_vel(trlend_vel)];
    
    trlinf{trllop}.eyetime = time_eye(trlstart_eye:trlend_eye);
    trlinf{trllop}.postime = time_pos(trlstart_pos:trlend_pos);
    trlinf{trllop}.dirtime = time_dir(trlstart_dir:trlend_dir);
    trlinf{trllop}.veltime = time_vel(trlstart_vel:trlend_vel);
    trlinf{trllop}.eyedat = [xeye_real(trlstart_eye:trlend_eye)'; yeye_real(trlstart_eye:trlend_eye)'];
    trlinf{trllop}.posdat = [xpos_real(trlstart_pos:trlend_pos)'; ypos_real(trlstart_pos:trlend_pos)'];
    trlinf{trllop}.dirdat = dir_real(trlstart_dir:trlend_dir)';
    trlinf{trllop}.veldat = vel_real(trlstart_vel:trlend_vel)';
    trlinf{trllop}.trltim = trltim;
   
end
ft_progress('close')

clear logtime_* id_eye x_eye y_eye id_ava x_ava y_ava
clear time_* xeye_real yeye_real xpos_real ypos_real dir_real vel_real


% resample to 1000 Hz and interpolate data points between samples
disp('Upsampling data points to 1000 Hz')

% the data structure created here is easily merged with the Fieldtrip
% (neural data analysis toolbox) data structure

clear dat*

dattime = cell(1,length(trlinf));
dateyedat = cell(1,length(trlinf));
datposdat = cell(1,length(trlinf));
datdirdat = cell(1,length(trlinf));
datveldat = cell(1,length(trlinf));

if runpar==1
    
    parfor trllop = 1:length(trlinf)

        % interpolate eye data
        eyetime_interp = trlinf{trllop}.eyetime(1):trlinf{trllop}.eyetime(end);
        eye_interp = nan(2,length(eyetime_interp));
        for k=1:length(trlinf{trllop}.eyetime)
            eye_interp(:,eyetime_interp==trlinf{trllop}.eyetime(k)) = trlinf{trllop}.eyedat(:,k);
        end
        eye_interp=fill_nans(eye_interp);
        
        % interpolate position data
        postime_interp = trlinf{trllop}.postime(1):trlinf{trllop}.postime(end);
        pos_interp = nan(2,length(postime_interp));
        for k=1:length(trlinf{trllop}.postime)
            pos_interp(:,postime_interp==trlinf{trllop}.postime(k)) = trlinf{trllop}.posdat(:,k);
        end
        pos_interp=fill_nans(pos_interp);
        
        % interpolate direction data
        dirtime_interp = trlinf{trllop}.dirtime(1):trlinf{trllop}.dirtime(end);
        dir_interp = nan(1,length(dirtime_interp));
        for k=1:length(trlinf{trllop}.dirtime)
            dir_interp(1,dirtime_interp==trlinf{trllop}.dirtime(k)) = trlinf{trllop}.dirdat(1,k);
        end
        dir_interp=fill_nans(dir_interp);
        
        % interpolate velocity data
        veltime_interp = trlinf{trllop}.veltime(1):trlinf{trllop}.veltime(end);
        vel_interp = nan(1,length(veltime_interp));
        for k=1:length(trlinf{trllop}.veltime)
            vel_interp(1,veltime_interp==trlinf{trllop}.veltime(k)) = trlinf{trllop}.veldat(1,k);
        end
        vel_interp=fill_nans(vel_interp);
        
        dateyedat{trllop} = eye_interp(:,find(eyetime_interp==trlinf{trllop}.trltim(1)):find(eyetime_interp==trlinf{trllop}.trltim(2)));
        datposdat{trllop} = pos_interp(:,find(postime_interp==trlinf{trllop}.trltim(1)):find(postime_interp==trlinf{trllop}.trltim(2)));
        datdirdat{trllop} = dir_interp(find(dirtime_interp==trlinf{trllop}.trltim(1)):find(dirtime_interp==trlinf{trllop}.trltim(2)));
        datveldat{trllop} = vel_interp(find(veltime_interp==trlinf{trllop}.trltim(1)):find(veltime_interp==trlinf{trllop}.trltim(2)));
        
        % additional step converts direction data to angle in radians
        datdirdat{trllop} = mod(datdirdat{trllop} + 90, 360) * pi/180;
        
        dattime{trllop} = trlinf{trllop}.trltim(1):trlinf{trllop}.trltim(2);
        
    end
        
else
    
    ft_progress('init', 'etf',     'Please wait...');
    for trllop = 1:length(trlinf)
        
        ft_progress(trllop/length(trlinf), 'Processing trial %d from %d', trllop, length(trlinf));
        
        % interpolate eye data
        eyetime_interp = trlinf{trllop}.eyetime(1):trlinf{trllop}.eyetime(end);
        eye_interp = nan(2,length(eyetime_interp));
        for k=1:length(trlinf{trllop}.eyetime)
            eye_interp(:,eyetime_interp==trlinf{trllop}.eyetime(k)) = trlinf{trllop}.eyedat(:,k);
        end
        eye_interp=fill_nans(eye_interp);
        
        % interpolate position data
        postime_interp = trlinf{trllop}.postime(1):trlinf{trllop}.postime(end);
        pos_interp = nan(2,length(postime_interp));
        for k=1:length(trlinf{trllop}.postime)
            pos_interp(:,postime_interp==trlinf{trllop}.postime(k)) = trlinf{trllop}.posdat(:,k);
        end
        pos_interp=fill_nans(pos_interp);
        
        % interpolate direction data
        dirtime_interp = trlinf{trllop}.dirtime(1):trlinf{trllop}.dirtime(end);
        dir_interp = nan(1,length(dirtime_interp));
        for k=1:length(trlinf{trllop}.dirtime)
            dir_interp(1,dirtime_interp==trlinf{trllop}.dirtime(k)) = trlinf{trllop}.dirdat(1,k);
        end
        dir_interp=fill_nans(dir_interp);
        
        % interpolate velocity data
        veltime_interp = trlinf{trllop}.veltime(1):trlinf{trllop}.veltime(end);
        vel_interp = nan(1,length(veltime_interp));
        for k=1:length(trlinf{trllop}.veltime)
            vel_interp(1,veltime_interp==trlinf{trllop}.veltime(k)) = trlinf{trllop}.veldat(1,k);
        end
        vel_interp=fill_nans(vel_interp);
        
        dateyedat{trllop} = eye_interp(:,find(eyetime_interp==trlinf{trllop}.trltim(1)):find(eyetime_interp==trlinf{trllop}.trltim(2)));
        datposdat{trllop} = pos_interp(:,find(postime_interp==trlinf{trllop}.trltim(1)):find(postime_interp==trlinf{trllop}.trltim(2)));
        datdirdat{trllop} = dir_interp(find(dirtime_interp==trlinf{trllop}.trltim(1)):find(dirtime_interp==trlinf{trllop}.trltim(2)));
        datveldat{trllop} = vel_interp(find(veltime_interp==trlinf{trllop}.trltim(1)):find(veltime_interp==trlinf{trllop}.trltim(2)));
        
        % additional step converts direction data to angle in radians
        datdirdat{trllop} = mod(datdirdat{trllop} + 90, 360) * pi/180;
        
        dattime{trllop} = trlinf{trllop}.trltim(1):trlinf{trllop}.trltim(2);
        
    end
    ft_progress('close')
    
end


clear trldat
trldat.time = dattime;
trldat.eyedat = dateyedat;
trldat.posdat = datposdat;
trldat.dirdat = datdirdat;
trldat.veldat = datveldat;

clear dattime dateyedat datposdat datdirdat datveldat

% % plot example trial(s) to double-check alignment
% trllop=2;
% figure;hold on;
% plot(trldat.time{trllop},trldat.posdat{trllop}');
% plot(trldat.time{trllop},trldat.dirdat{trllop}');
% plot(trldat.time{trllop},trldat.veldat{trllop}');
% plot(trldat.time{trllop},trldat.eyedat{trllop}')


% add fruit info
disp('Loading fruit PAR file')

% % appearance of new bananas: '100'
% log_trl_starttime = unique(logtime_eye(strcmp(id_eye,'100')));

[logtime_frt, id_frt, alph, x_frt, y_frt] = textread(fullfile(logDir,sessionDir,'fruit.par'), '%f%s%f%f%f'); % avatar (position) data

trldat.frtpos = cell(1);
trldat.frttim = cell(1);
trldat.alpha = cell(1);
trldat.newfrt = cell(1);
for trllop = 1:length(trldat.time)
    
    % logind: index of which data to use for the trial
    logind = logtime_frt>=trldat.time{trllop}(1) & logtime_frt<=trldat.time{trllop}(end);
    
    % fruit positions for each trial
    % include timestamps and 'alph' (0 for banana, 1 for cherry)
    frtposdum = [logtime_frt(logical(logind.*strcmp(id_frt,'pos'))) alph(logical(logind.*strcmp(id_frt,'pos'))) x_frt(logical(logind.*strcmp(id_frt,'pos'))) y_frt(logical(logind.*strcmp(id_frt,'pos')))];
    % some positions are all zeros (happens at beginning of session)
    trldat.frtpos{trllop} = frtposdum(frtposdum(:,3)~=0 | frtposdum(:,4)~=0,:);
    
    % timestamps for getting fruit for each trial
    trldat.frttim{trllop} = [logtime_frt(logical(logind.*strcmp(id_frt,'eaten'))) alph(logical(logind.*strcmp(id_frt,'eaten')))];
    
    % alpha levels with timestamps, for each trial
    trldat.alpha{trllop} = [logtime_frt(logical(logind.*strcmp(id_frt,'alpha'))) alph(logical(logind.*strcmp(id_frt,'alpha')))];
    % if no alpha specified at trial start, alpha is set to 0
    if isempty(find(logtime_frt==trldat.time{trllop}(1).*strcmp(id_frt,'alpha'),1))
        trldat.alpha{trllop} =  [trldat.time{trllop}(1) 0; trldat.alpha{trllop}];
    end
    
    trldat.newfrt{trllop} = logtime_frt(logical(logind.*strcmp(id_frt,'100')));

end


% save the behavioral data file

[~,sesnam]=fileparts(sessionDir);

save(fullfile(logDir,sessionDir,[sesnam '_trldat.mat']),'trldat')

disp(['Created ' fullfile(logDir,sessionDir,[sesnam '_trldat.mat'])])



function A = fill_nans(A)
% Replaces the nans in each column with previous non-nan values.

for ii = find(~isnan(A(1,:)))
    I = A(:,ii);
    subind = ii+1:find(~isnan(A(1,ii+1:end)),1,'first')+ii-1;
    A(:,subind) = repmat(I,1,length(subind));
end


