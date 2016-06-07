% BR&Panda_to_FT.m
%
% imports neural (LFP) data from Blackrock NS6 data (decimated to 1000 Hz)
%
% imports data from parsed python files (.par)
% these files contain eye, position and velocity data all using the same
% timescale, but with samples of different data types falling at different
% timespoints; the timescale is resampled at 1000 Hz to match sampling rate 
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
% timestamp values in panda log are already in ms

%% specify filenames (need to do this manually for now)

% specify the directory containing log files
% BRnam = 'JN150209001'; sessionDir = 'JN_15_02_09\JN_15_02_09_15_34';
BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';


logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
savDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trial data';
decDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\NS6 - decimated';

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


%% load the raw data (NS6, decimated version); align the NEV (timestamp) and NS (continuous) signals

NS2 = openNSx(fullfile(BRDir,[BRnam '.NS2']),'read');
% NEV = openNEV(fullfile(BRDir,[BRnam '.nev'])); % integrate NEV later
% % NS6 = openNSx(fullfile(datdir,[BRnam '.ns6']));
% % load decimated MAT version of NS6 file instead of directly from NS6 file
load(fullfile(decDir,[BRnam '_NS6_SF30.mat']))

fs = 1/NS2.MetaTags.SamplingFreq;

% NS2_timestamp = datenum(ns2datevec)*86400:fs:datenum(ns2datevec)*86400+NS2.MetaTags.DataDurationSec-fs;
% NS6_timestamp = datenum(ns6datevec)*86400:fs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec-fs;


%% line up behavioral data with raw LFPs (adapted from quickLineup)


% find eye signals in NS2
eyeindx = nan;
eyeindy = nan;
posindx = nan;
posindy = nan;
lfpind = [];
for lablop = 1:length(NS2.ElectrodesInfo)
    if strncmpi(NS2.ElectrodesInfo(lablop).Label,'ainp1',5)
        eyeindx = lablop;
    elseif strncmpi(NS2.ElectrodesInfo(lablop).Label,'ainp2',5)
        eyeindy = lablop;
    elseif strncmpi(NS2.ElectrodesInfo(lablop).Label,'ainp3',5)
        posindx = lablop;
    elseif strncmpi(NS2.ElectrodesInfo(lablop).Label,'ainp4',5)
        posindy = lablop;
    elseif ismember(NS2.ElectrodesInfo(lablop).Label(1),['A' 'B' 'C'])
        lfpind = [lfpind; lablop];
    end
end

if isfield(data,'trial')
    data=rmfield(data,'trial');
end
if isfield(data,'eyedatBR')
    data=rmfield(data,'eyedatBR');
end
if isfield(data,'posdatBR')
    data=rmfield(data,'posdatBR');
end
if isfield(data,'trl')
    data=rmfield(data,'trl');
end

clear trial eyedatBR posdatBR trl
parfor trllop = 1:length(data.time)
    
    Xs = xcorr(data.eyedat{trllop}(1,:),double(NS2.Data(eyeindx,:)));
    [~,Is] = max(Xs);
    trl_start = round(length(Xs)/2-Is);
    trl_end = trl_start+length(data.time{trllop})-1;
    
    % Save as a new file structure containing all the necessary info.
    trial{trllop} = double(NS2.Data(lfpind,trl_start:trl_end));
    eyedatBR{trllop} = double(NS2.Data(eyeindx:eyeindy,trl_start:trl_end));
    posdatBR{trllop} = double(NS2.Data(posindx:posindy,trl_start:trl_end));
    
    trl(trllop,:) = [trl_start trl_end];
    
end

data.trial = trial;
data.eyedatBR = eyedatBR;
data.posdatBR = posdatBR;
data.sampleinfo = trl;
data.fsample = 1000;
clear trial eyedatBR posdatBR trl

% % plot eye and position data trial by trial to ensure accuracy
% close all
% for trllop = 1:length(data.eyedat)
%     
%     w=data.eyedat{trllop}(1,:);
%     ww=data.eyedatBR{trllop}(1,:);
%     x=w/(max(w)-min(w));
%     xx=ww/(max(ww)-min(ww));
%     
%     y=data.posdat{trllop}(1,:);
%     yy= data.posdatBR{trllop}(2,:); % BR position data is rotated relative to Panda
%     z=y/(max(y)-min(y));
%     zz=yy/(max(yy)-min(yy));
%     
%     figure
%     subplot(2,1,1);plot([x; xx]')
%     subplot(2,1,2);plot([z; zz]')
%     pause;close
%     
% end

% for k=1:size(data.trial{1},1)
%     data.label{k} = NS2.ElectrodesInfo(k).Label(1:find(double(NS2.ElectrodesInfo(k).Label)==0,1,'first')-1);
% end
data.label = {'A01'; 'A02'; 'A03'; 'A04'; 'A05'; 'A06'; 'A07'; ...
    'A08'; 'A09'; 'A10'; 'A11'; 'A12'; 'B01'; 'B02'; 'B03'; 'B04'; ...
    'B05'; 'B06'; 'B07'; 'B08'; 'B09'; 'B10'; 'B11'; 'B12'; 'C01'; ...
    'C02'; 'C03'; 'C04'; 'C05'; 'C06'; 'C07'; 'C08'; 'C09'; 'C10'; ...
    'C11'; 'C12'};


%% save the complete data file

save(fullfile(savDir,[BRnam '_trldat.mat']),'data')


%% SKIP TO HERE: load the complete data file

load(fullfile(savDir,[BRnam '_trldat.mat']))


%% create event and time arrays for the trial events

encode = double(NEV.Data.SerialDigitalIO.UnparsedData);
ts = double(NEV.Data.SerialDigitalIO.TimeStamp)/30;

event_arr = [];
time_arr = [];

c=0;

for rptlop = 1:length(encode)
    if ismember(encode(rptlop),1000:(1000+length(data.trial)))
        if c==0
            c=c+1;
            event_dum = encode(rptlop);
            time_dum = ts(rptlop);
        else
            if (encode(rptlop)-event_dum(1))==1
                if ~isempty(event_arr)
                    j = length(event_dum);
                    k = size(event_arr,1);
                    event_arr = [event_arr; nan(max([j k])-k,size(event_arr,2))];
                    time_arr = [time_arr; nan(max([j k])-k,size(time_arr,2))];
                    event_arr(:,c) = [event_dum; nan(max([j k])-j,1)];
                    time_arr(:,c) = [time_dum; nan(max([j k])-j,1)];
                else
                    event_arr(:,c) = event_dum;
                    time_arr(:,c) = time_dum;
                end
                c=c+1;
                event_dum = encode(rptlop);
                time_dum = ts(rptlop);
            else
                event_dum = [event_dum; encode(rptlop)];
                time_dum = [time_dum; ts(rptlop)];
            end
        end
    else
        event_dum = [event_dum; encode(rptlop)];
        time_dum = [time_dum; ts(rptlop)];
    end
end

