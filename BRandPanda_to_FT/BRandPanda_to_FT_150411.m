% BR&Panda_to_FT.m
%
% imports neural (LFP) data from Blackrock NS6 data (decimated to 1000 Hz)
%
% imports data from parsed python files (.par)
% these files contain eye, position and velocity data all using the same
% timescale, but with samples of different data types falling at different
% timespoints; the timescale is resampled at 1 kHz to match sampling rate 
% of neural recordings, and additional data points between samples are 
% filled in using interpolation.
%
% the data is integrated into Fieldtrip format for further analysis
%
% 141017 Mike Jutras (with parts adapted from quickLineup by Yoni Browning)
%
%
% this code requires ft_nearest.m, inpaint_nans.m, ft_progress.m
% before running, must run Python script MakeParFiles in the directory containing 
% the Python log files
% (e.g. R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\JN_14_08_25_13_12)

% Notes:
% revised in part from analyzeNav140418
% python samples position/direction at ~83 Hz; eye data @240
%
% timestamp values in panda log are already in ms (1 kHz)

% CHANGE LOG:
% 3/18/15 - MJJ - the data structure generated at the end now provides the
% correct info (and only the correct info, i.e. nothing extra from trldat)
% that Fieldtrip needs. "data.trial" includes eye and position data along
% with the neural data; these "analog" data are pulled from the Python log
% and interpolated to 1 kHz sampling. "data.time" is designated in units of
% seconds, with the first time sample, 0, corresponding with the first
% sample designated in data.sampleinfo.
% Also, uses neural data from NS6 file instead of NS2 file. NS2 data are
% delayed relative to NS6 data due to the filtering that we usually use.

% TO DO:
% - integrate NEV data, including trial event timestamps ("encodes")

%% specify filenames (need to do this manually for now)

% specify the directory containing log files
% BRnam = 'JN150209001'; sessionDir = 'JN_15_02_09\JN_15_02_09_15_34';
BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';


logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\2014';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
savDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\trial data';
decDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\NS6 - decimated';

% copy files to local drive if network transfer rate is slow
% BRDir = 'C:\Data\VR';
% savDir = 'C:\Data\VR';
% decDir = 'C:\Data\VR';

%% import panda log file info

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

% start and duration of each trial (when bananas appear)
% the last instance of '100' will mark the end of the last trial
log_trl_starttime = unique(logtime_eye(strcmp(id_eye,'100')));
log_trl_length = diff(log_trl_starttime);

%% get all the data into one structure

trlinf = [];
trltimarr = [];
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(log_trl_starttime)-1
    
    ft_progress(trllop/(length(log_trl_starttime)-1), 'Processing event %d from %d', trllop, (length(log_trl_starttime)-1));
    
    % trial time according to panda log file
    trltim = [log_trl_starttime(trllop) log_trl_starttime(trllop+1)-1];
    
%     % buffer the start and end times for each trial 100 samples on each side
%     trlstart_eye = ft_nearest(time_eye,trltim(1)-100);
%     trlend_eye = ft_nearest(time_eye,trltim(2)+100);
%     
%     trlstart_pos = ft_nearest(time_pos,trltim(1)-100);
%     trlend_pos = ft_nearest(time_pos,trltim(2)+100);
%     
%     trlstart_dir = ft_nearest(time_dir,trltim(1)-100);
%     trlend_dir = ft_nearest(time_dir,trltim(2)+100);
%     
%     trlstart_vel = ft_nearest(time_vel,trltim(1)-100);
%     trlend_vel = ft_nearest(time_vel,trltim(2)+100);

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


%% resample to 1000 Hz and interpolate data points between samples

% the data structure created here is easily merged with the Fieldtrip
% (neural data analysis toolbox) data structure

clear dat*

dattime = cell(1,length(trlinf));
dateyedat = cell(1,length(trlinf));
datposdat = cell(1,length(trlinf));
datdirdat = cell(1,length(trlinf));
datveldat = cell(1,length(trlinf));

t0 = clock;

parfor trllop = 1:length(trlinf)
        
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

    
    dateyedat{trllop} = eye_interp(:,find(eyetime_interp==trlinf{trllop}.trltim(1)):find(eyetime_interp==trlinf{trllop}.trltim(2)));
    datposdat{trllop} = pos_interp(:,find(postime_interp==trlinf{trllop}.trltim(1)):find(postime_interp==trlinf{trllop}.trltim(2)));
    datdirdat{trllop} = dir_interp(find(dirtime_interp==trlinf{trllop}.trltim(1)):find(dirtime_interp==trlinf{trllop}.trltim(2)));
    datveldat{trllop} = vel_interp(find(veltime_interp==trlinf{trllop}.trltim(1)):find(veltime_interp==trlinf{trllop}.trltim(2)));
    
    % additional step converts direction data to angle in radians
    datdirdat{trllop} = mod(datdirdat{trllop} + 90, 360) * pi/180;

    dattime{trllop} = trlinf{trllop}.trltim(1):trlinf{trllop}.trltim(2);
    
end

fprintf('%g\n',etime(clock,t0));

clear data
data.time = dattime;
data.eyedat = dateyedat;
data.posdat = datposdat;
data.dirdat = datdirdat;
data.veldat = datveldat;

clear dattime dateyedat datposdat datdirdat datveldat

% % plot example trial(s) to double-check alignment
% trllop=2;
% figure;hold on;
% plot(data.time{trllop},data.posdat{trllop}');
% plot(data.time{trllop},data.dirdat{trllop}');
% plot(data.time{trllop},data.veldat{trllop}');
% plot(data.time{trllop},data.eyedat{trllop}')

%% save the behavioral data file
% (if stopping here, otherwise proceed to create trial data file)

save(fullfile(savDir,[BRnam '_navdat.mat']),'data')


%% SKIP TO HERE: load the behavioral data file

load(fullfile(savDir,[BRnam '_navdat.mat']))

% rename data structure to avoid confusion with fieldtrip data
trldat = data;
clear data

% trldat now contains all of the relevant information from the python log
% file except for the banana/fruit locations and when the monkey got
% rewarded; we can add these now

%% add banana info

[logtime_ban, id_ban, ~, x_ban, y_ban] = textread(fullfile(logDir,sessionDir,'banana.par'), '%f%s%f%f%f'); % avatar (position) data
% [logtime_ban, id_ban, ~, x_ban, y_ban] = textread(fullfile(logDir,sessionDir,'fruit.par'), '%f%s%f%f%f'); % avatar (position) data

trldat.banpos = cell(1);
trldat.eattim = cell(1);
for trllop = 1:length(trldat.time)
    
    % logind: index of which data to use for the trial
    logind = logtime_ban>=trldat.time{trllop}(1) & logtime_ban<=trldat.time{trllop}(end);
    
    % banana positions for each trial
    trldat.banpos{trllop} = [x_ban(logical(logind.*strcmp(id_ban,'pos'))) y_ban(logical(logind.*strcmp(id_ban,'pos')))];
    
    % timestamps for getting bananas for each trial
    trldat.eattim{trllop} = logtime_ban(logical(logind.*strcmp(id_ban,'eaten')));
    
end


%% load the raw data (NS6, decimated version); align the NEV (timestamp) and NS (continuous) signals


NEV = openNEV(fullfile(BRDir,[BRnam '.nev']),'read'); % NEV file; contains event codes and timestamps
nevDTR = NEV.MetaTags.DateTimeRaw; % start time of NEV file
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000]; % convert start time to date vector

load(fullfile(decDir,[BRnam '_NS6_SF30.mat'])) % decimated NS6 file; 1000 Hz
ns6DTR = NS6.MetaTags.DateTimeRaw; % start time of NS6 file
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000]; % convert start time to date vector

NS2 = openNSx(fullfile(BRDir,[BRnam '.NS2']),'read','precision','double'); % NS2 file (contains eye & position data)
ns2DTR = NS2.MetaTags.DateTimeRaw; % start time of NS2 file
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000]; % convert start time to date vector

% yes, these three files all have different start times, by a millisecond or two.

NEVfs = 1/double(NEV.MetaTags.SampleRes); % sampling frequency of NEV file (should be 1/30000)
NS2fs = 1/NS2.MetaTags.SamplingFreq; % sampling frequency of NS2 file (should be 0.001)
DECfs = 0.001; % sampling frequency of decimated NS6 file
NS6fs = 1/NS6.MetaTags.SamplingFreq; % sampling frequency of NS6 file (should also be 1/30000)

% create time arrays for each file; timestamps match data samples
NS2ts = datenum(ns2datevec)*86400:NS2fs:datenum(ns2datevec)*86400+NS2.MetaTags.DataDurationSec-NS2fs;
NS6ts = datenum(ns6datevec)*86400:NS6fs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec-NS6fs;
NEVts = datenum(nevdatevec)*86400:NEVfs:datenum(nevdatevec)*86400+NEV.MetaTags.DataDurationSec-NEVfs;
DECts = datenum(ns6datevec)*86400:DECfs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec;

% check that length of time array matches # of datapoints:
length(NS6ts)==NS6.MetaTags.DataPoints
length(NEVts)==NEV.MetaTags.DataDuration
length(NS2ts)==NS2.MetaTags.DataPoints
length(DECts)==size(NS6.Data,2)

% align every time array to the first sample of the NS6 file;
% this brings everything onto the same time scale (time values line up
% across different data streams with different start points/sampling rates,
% everything based on the start of the NS6 file)
% [DECts and NS6ts will have the same start value, since the first sample is
% the same for both data streams]
NEVts = (NEVts-NS6ts(1))+NS6fs;
DECts = (DECts-NS6ts(1))+NS6fs;
NS2ts = (NS2ts-NS6ts(1))+NS6fs;
NS6ts = (NS6ts-NS6ts(1))+NS6fs;

% create label arrays for NS2 and NS6 channels
NS2lab = cell(1,length(NS2.ElectrodesInfo));
for lablop = 1:length(NS2.ElectrodesInfo)
    NS2lab{lablop} = NS2.ElectrodesInfo(lablop).Label;
end
NS6lab = cell(1,length(NS6.ElectrodesInfo));
for lablop = 1:length(NS6.ElectrodesInfo)
    NS6lab{lablop} = NS6.ElectrodesInfo(lablop).Label;
end


%% line up behavioral data with raw LFPs (adapted from quickLineup)

% find eye & position signals in NS2
eyeindx = find(strncmpi(NS2lab,'ainp1',5));
eyeindy = find(strncmpi(NS2lab,'ainp2',5));
posindx = find(strncmpi(NS2lab,'ainp3',5));
posindy = find(strncmpi(NS2lab,'ainp4',5));

eyedatNS2 = NS2.Data([eyeindx eyeindy],:);
posdatNS2 = NS2.Data([posindx posindy],:);

clear trial sampleinfo eyedatBR posdatBR
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(trldat.eyedat)
    ft_progress(trllop/length(trldat.eyedat), 'Processing event %d from %d', trllop, length(trldat.eyedat));
    
    [Xs, lag] = xcorr(eyedatNS2(1,:), trldat.eyedat{trllop}(1,:));
    [~,Is] = max(Xs);
    trl_start = lag(Is);
    trl_end = trl_start+length(trldat.eyedat{trllop})-1;
    
    trl_start_tim = NS2ts(trl_start);
    trl_end_tim = NS2ts(trl_end);
    [~,trl_start_ind] = min(abs(DECts-trl_start_tim));
    [~,trl_end_ind] = min(abs(DECts-trl_end_tim));
    
    % Save as a new file structure containing all the necessary info.
    trial{trllop} = double(NS6.Data(:,trl_start_ind:trl_end_ind));
    
    % append the eye & position data from the Python log file
    trial{trllop}(end+1:end+2,:) = trldat.eyedat{trllop};
    trial{trllop}(end+1:end+2,:) = trldat.posdat{trllop};
    
    % append the eye and position data from Blackrock file
    trial{trllop}(end+1:end+2,:) = eyedatNS2(:,trl_start:trl_end);
    trial{trllop}(end+1:end+2,:) = posdatNS2(:,trl_start:trl_end);
    
    % include both Blackrock and Python analog data for now, until we
    % decide which version has cleaner A2D conversion

    % corresponds to samples in decimated NS6 file
    sampleinfo(trllop,:) = [trl_start_ind trl_end_ind];
    
end
ft_progress('close')

% create the data structure for Fieldtrip
clear data
data.trial = trial;
data.sampleinfo = sampleinfo;
for trllop = 1:size(sampleinfo,1)
    data.time{trllop} = 0:0.001:(diff(sampleinfo(trllop,:),1,2)/1000);
end
data.fsample = 1000;
data.label = {'A01'; 'A02'; 'A03'; 'A04'; 'A05'; 'A06'; 'A07'; ...
    'A08'; 'A09'; 'A10'; 'A11'; 'A12'; 'B01'; 'B02'; 'B03'; 'B04'; ...
    'B05'; 'B06'; 'B07'; 'B08'; 'B09'; 'B10'; 'B11'; 'B12'; 'C01'; ...
    'C02'; 'C03'; 'C04'; 'C05'; 'C06'; 'C07'; 'C08'; 'C09'; 'C10'; ...
    'C11'; 'C12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'; 'eyeX_B'; ...
    'eyeY_B'; 'posX_B'; 'posY_B'};
clear trial sampleinfo eyedatNS2 posdatNS2


%% save the complete data file

save(fullfile(savDir,[BRnam '_NSdat.mat']),'data')


%% plot eye and position data trial by trial to ensure accuracy

close all
xpi = find(strcmp(data.label,'eyeX_P'));
xbi = find(strcmp(data.label,'eyeX_B'));
ppi = find(strcmp(data.label,'posX_P'));
pbi = find(strcmp(data.label,'posY_B')); % BR position data is rotated relative to Panda
for trllop = 1:length(data.trial)

    xp = data.trial{trllop}(xpi,:);
    xb = data.trial{trllop}(xbi,:);
    
    xpn = xp/(max(xp)-min(xp));
    xbn = xb/(max(xb)-min(xb));
    
    pp = data.trial{trllop}(ppi,:);
    pb = data.trial{trllop}(pbi,:);
    
    ppn = pp/(max(pp)-min(pp));
    pbn = pb/(max(pb)-min(pb));
    
    figure
    subplot(2,1,1);plot([xpn; xbn]')
    subplot(2,1,2);plot([ppn; pbn]')
    pause;close
    
end


%% compare timestamps when banana is eaten between NEV and Python logs

bandifmat = [];
for trllop = 1:length(data.trial)
    
    [~,NEVbnd1] = min(abs(NEVts-DECts(data.sampleinfo(trllop,1))));
    [~,NEVbnd2] = min(abs(NEVts-DECts(data.sampleinfo(trllop,2))));

    NEVtrlind = find(NEV.Data.SerialDigitalIO.TimeStamp>=NEVbnd1 & ...
        NEV.Data.SerialDigitalIO.TimeStamp<=NEVbnd2);
    
    % hopefully there are only 10 occurences of bananas being eaten!
    banind = find(NEV.Data.SerialDigitalIO.UnparsedData(NEVtrlind)==200);
    if length(banind)==10
        
        banTS = NEV.Data.SerialDigitalIO.TimeStamp(NEVtrlind(banind));
        bantim = NEVts(banTS);
        bandif = (bantim - DECts(data.sampleinfo(trllop,1)))*1000; % convert to ms
        
        bandifmat = [bandifmat; bandif' trldat.eattim{trllop}-trldat.time{trllop}(1)];
        
    end
    
end

median(diff(bandifmat,1,2)) % 66.0

figure;subplot(2,1,1);hold on
scatter(1:size(bandifmat,1),diff(bandifmat,1,2))
bandiftrlmat = [];
for k=1:length(data.trial)
    line([k*10+0.5 k*10+0.5],ylim,'Color','r')
    bandiftrlmat(:,k) = diff(bandifmat((k*10-9):k*10,:),1,2);
end
subplot(2,1,2)
scatter(1:10,mean(bandiftrlmat,2))

