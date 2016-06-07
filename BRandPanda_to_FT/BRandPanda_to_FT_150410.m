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
% BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
% savDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\trial data';
% decDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\NS6 - decimated';

% copy files to local drive if network transfer rate is slow
BRDir = 'C:\Data\VR';
savDir = 'C:\Data\VR';
decDir = 'C:\Data\VR';

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

%% create trials based on the NEV file

% % for JN150209001, fix weird event where 1000 appears twice instead of 1001
% NEV.Data.SerialDigitalIO.UnparsedData(2202)=1001;

c=1;
m=1000; % trial onset codes start at 1000 and progress incrementally
clear NEVper
NEVper.val={};
NEVper.tim={};
for vallop = 1:length(NEV.Data.SerialDigitalIO.UnparsedData)
% for vallop = 2142:length(NEV.Data.SerialDigitalIO.UnparsedData) % for JN150209001
    if NEV.Data.SerialDigitalIO.UnparsedData(vallop)==m
        
        m=m+1;
        for k=1:find(NEV.Data.SerialDigitalIO.UnparsedData(vallop:end)==m,1,'first')-1
            
            % build NEVper.val, each cell contains one trial, each row is an event
            NEVper.val{c}(k,1) = NEV.Data.SerialDigitalIO.UnparsedData(vallop+k-1);
            
            % NEVevtTSind: NEV timestamp corresponding with the event code
            NEVevtTSind = NEV.Data.SerialDigitalIO.TimeStamp(vallop+k-1);
            
            % actual timestamp is given in "NS6 clock time"
            NEVper.tim{c}(k,1) = NEVts(NEVevtTSind);
            
        end
        
        c=c+1;
        vallop = vallop + find(NEV.Data.SerialDigitalIO.UnparsedData(vallop:end)==m,1,'first') - 1;
        
    end
end

% the number of trials (cells in NEVper.val and NEVper.tim) will hopefully
% match the number of trials in trldat

%% set up trial arrays with neural data, eye and position data from blackrock files

% find eye & position signals in NS2
eyeindx = find(strncmpi(NS2lab,'ainp1',5));
eyeindy = find(strncmpi(NS2lab,'ainp2',5));
posindx = find(strncmpi(NS2lab,'ainp3',5));
posindy = find(strncmpi(NS2lab,'ainp4',5));

% build a trial array using the timestamp markers in NEVper
clear trial sampleinfo eyedatBR posdatBR
for trllop = 1:length(NEVper.val)
    
    onset = NEVper.tim{trllop}(1);
    
    % find indices of nearest timepoint to trial onset in DEC and NS2
    [~,DEConsind] = min(abs(DECts-onset));
    [~,NS2onsind] = min(abs(NS2ts-onset));
    
    % trial: neural data from decimated NS6 file
    trial{trllop} = double(NS6.Data(:,DEConsind:DEConsind+length(trldat.time{trllop})-1));

    % build arrays for eye and position data from NS2 file
    eyedatBR{trllop} = NS2.Data(eyeindx:eyeindy,NS2onsind:NS2onsind+length(trldat.time{trllop})-1);
    posdatBR{trllop} = NS2.Data(posindx:posindy,NS2onsind:NS2onsind+length(trldat.time{trllop})-1);
    
    % sampleinfo comes into use later, in Fieldtrip format
    sampleinfo(trllop,:) = [DEConsind DEConsind+length(trldat.time{trllop})-1];

end

%% plot some data from an example trial to compare python log data to blackrock data

% look at avatar position when collecting bananas during trial #1, compared
% to position of bananas

% start with python log data
figure; hold on
% plot banana positions
scatter(trldat.banpos{1}(:,1),trldat.banpos{1}(:,2))
% plot avatar position when acquiring bananas
eatpos = [];
for banlop = 1:length(trldat.eattim{1})
    eatpos(banlop,:) = trldat.posdat{1}(:,trldat.time{1}==trldat.eattim{1}(banlop))';
end
scatter(eatpos(:,1),eatpos(:,2),'r')

% avatar position when acquiring bananas is slightly offset from actual
% banana position; this is ok, because acquisition occurs when collision
% spheres around avatar and banana touch

% now look at blackrock data
trltim = DECts(sampleinfo(1,1)):0.001:DECts(sampleinfo(1,2));
eattim = NEVper.tim{1}(NEVper.val{1}==200);
eatposBR = [];
for banlop = 1:length(eattim)
    [~,banind] = min(abs(trltim-eattim(banlop)));
    eatposBR(banlop,:) = posdatBR{1}(:,banind)';
end
figure;scatter(eatposBR(:,2),eatposBR(:,1),'r')

% try to line them up?
figure;hold on
scatter(eatposBR(:,2)/2600,eatposBR(:,1)/2600,'b')
scatter(eatpos(:,1),eatpos(:,2),'r')

% % do the same "actual" position vs. "eat" position, for trial 20
% % (doesn't work because different scale)
% figure;scatter(NEVper.val{20}(2:2:20),NEVper.val{20}(3:2:21))
% trltim = DECts(sampleinfo(20,1)):0.001:DECts(sampleinfo(20,2));
% eattim = NEVper.tim{20}(NEVper.val{20}==200);
% eatposBR = [];
% for banlop = 1:length(eattim)
%     [~,banind] = min(abs(trltim-eattim(banlop)));
%     eatposBR(banlop,:) = posdatBR{20}(:,banind)';
% end
% hold on;scatter(eatposBR(:,2),eatposBR(:,1),'r')

% plot "eat times" for each data stream
figure;hold on
scatter(1:length(eattim),(eattim-NEVper.tim{1}(1))*1000,'b')
scatter(1:length(trldat.eattim{1}),trldat.eattim{1}-trldat.time{1}(1),'r')

% blackrock "eat times" occur consistently later than python "eat times"
figure;hist(((eattim-NEVper.tim{1}(1))*1000)-(trldat.eattim{1}-trldat.time{1}(1)))
% ...suggesting that there is a lag from the event occurring from the
% monkey's perspective and when the event is marked in blackrock?


% plot the normalized position information from both data streams on the same plot
BRposnorm = [posdatBR{1}(1,:)/(max(posdatBR{1}(1,:))-min(posdatBR{1}(1,:)));
    posdatBR{1}(2,:)/(max(posdatBR{1}(2,:))-min(posdatBR{1}(2,:)))];
posnorm = [trldat.posdat{1}(2,:)/(max(trldat.posdat{1}(2,:))-min(trldat.posdat{1}(2,:)));
    trldat.posdat{1}(1,:)/(max(trldat.posdat{1}(1,:))-min(trldat.posdat{1}(1,:)))];
figure;hold on
plot(BRposnorm(1,:),'b')
plot(posnorm(1,:),'r')
% blackrock position data lag behind python by about 74 ms?
[Xs, lag] = xcorr(posnorm(1,:),BRposnorm(1,:));
[~,Is] = max(Xs);
datoff = lag(Is);
% cross correlation isn't showing the lag but it's definitely there!

% plot the eye data for the two data streams
eyedatBR_norm = [eyedatBR{1}(1,:)/(max(eyedatBR{1}(1,:))-min(eyedatBR{1}(1,:)));
    eyedatBR{1}(2,:)/(max(eyedatBR{1}(2,:))-min(eyedatBR{1}(2,:)))];
eyedat_norm = [trldat.eyedat{1}(1,:)/(max(trldat.eyedat{1}(1,:))-min(trldat.eyedat{1}(1,:)));
    trldat.eyedat{1}(2,:)/(max(trldat.eyedat{1}(2,:))-min(trldat.eyedat{1}(2,:)))];
figure;plot(eyedatBR_norm','b');hold on;plot(eyedat_norm','r')
% blackrock eye data lag behind python data by about 120 ms?
[Xs, lag] = xcorr(eyedatBR_norm(1,:),eyedat_norm(1,:));
[~,Is] = max(Xs);
datoff = lag(Is); % this shows lag to be 100 ms
% so, if the trial onset markers are simultaneous, then the blackrock
% eye data are delayed with respect to the python eye data
% or... if the eye data streams are in sync, then the trial onset markers
% are off... with the blackrock trial onset marker occuring earlier than
% the python marker... which doesn't make sense

% try another trial, #30
eyedatBR_norm = [eyedatBR{30}(1,:)/(max(eyedatBR{30}(1,:))-min(eyedatBR{30}(1,:)));
    eyedatBR{30}(2,:)/(max(eyedatBR{30}(2,:))-min(eyedatBR{30}(2,:)))];
eyedat_norm = [trldat.eyedat{30}(1,:)/(max(trldat.eyedat{30}(1,:))-min(trldat.eyedat{30}(1,:)));
    trldat.eyedat{30}(2,:)/(max(trldat.eyedat{30}(2,:))-min(trldat.eyedat{30}(2,:)))];
figure;plot(eyedatBR_norm','b');hold on;plot(eyedat_norm','r')
[Xs, lag] = xcorr(eyedatBR_norm(1,:),eyedat_norm(1,:));
[~,Is] = max(Xs);
datoff = lag(Is); % this shows lag to be only 26 ms

% do every trial!
datoff_pos = []; datoff_eye = [];
for trllop = 1:length(eyedatBR)
    
    BRposnorm = posdatBR{trllop}(1,:)/(max(posdatBR{trllop}(1,:))-min(posdatBR{trllop}(1,:)));
    posnorm = trldat.posdat{trllop}(2,:)/(max(trldat.posdat{trllop}(2,:))-min(trldat.posdat{trllop}(2,:)));
    [Xs, lag] = xcorr(BRposnorm,posnorm);
    [~,Is] = max(Xs);
    datoff_pos(trllop) = lag(Is);
    
    eyedatBR_norm = eyedatBR{trllop}(1,:)/(max(eyedatBR{trllop}(1,:))-min(eyedatBR{trllop}(1,:)));
    eyedat_norm = trldat.eyedat{trllop}(1,:)/(max(trldat.eyedat{trllop}(1,:))-min(trldat.eyedat{trllop}(1,:)));
    [Xs, lag] = xcorr(eyedatBR_norm,eyedat_norm);
    [~,Is] = max(Xs);
    datoff_eye(trllop) = lag(Is);
    
end
figure;hist(datoff_pos,20) % median position data lag is 0, with one outlier
figure;hist(datoff_eye,20) % median eye data lag is 25 ms

% so if the eye data were in sync, that suggests that the position data
% streams are off, since the lag isn't as bad... it's more likely that the
% eye data are off, and that the eye data in the blackrock NS2 file are
% inaccurate, which is bad!

%% compare time between trial start and getting bananas, between data streams

PYdel = [];
BRdel = [];
for trllop = 2:length(trldat.time)
    PYdel = [PYdel; trldat.eattim{trllop}-trldat.time{trllop}(1)];
    BRdel = [BRdel; (NEVper.tim{trllop}(NEVper.val{trllop}==200)-NEVper.tim{trllop}(1))*1000];
end
figure;hist(diff([PYdel BRdel],1,2),100) % median difference: -42 ms, python has longer delays


%% calculate saccade-aligned neural signal, using the two different eye signals

eyesigBR = NS2.Data(eyeindx:eyeindy,:);
eyesigPY = []; eyesigPY_tim = []; eyesigPY_trlind = [];
for trllop = 1:length(trldat.eyedat)
    eyesigPY = [eyesigPY trldat.eyedat{trllop}];
    eyesigPY_tim = [eyesigPY_tim trldat.time{trllop}];
    eyesigPY_trlind = [eyesigPY_trlind ones(size(trldat.time{trllop}))*trllop];
end

fltord = 40;
lowpasfrq = 80;
limBR = 100000;
limPY = 2000;
artpadding = 0.01;
Fs = 1000; % sampling rate in Hz

%low pass filter the eye position data_eye
nyqfrq = Fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
e_x_lowpassfilterBR=filtfilt(flt,1, eyesigBR(1,:));
e_y_lowpassfilterBR=filtfilt(flt,1, eyesigBR(2,:));
e_x_lowpassfilterPY=filtfilt(flt,1, eyesigPY(1,:));
e_y_lowpassfilterPY=filtfilt(flt,1, eyesigPY(2,:));

%differentiate and multiply with sampling rate to get velocity as deg/sec
x_vBR = diff(e_x_lowpassfilterBR) .* Fs;
y_vBR = diff(e_y_lowpassfilterBR) .* Fs;
x_vPY = diff(e_x_lowpassfilterPY) .* Fs;
y_vPY = diff(e_y_lowpassfilterPY) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
velBR = abs(complex(x_vBR,y_vBR));
velPY = abs(complex(x_vPY,y_vPY));


%detect saccade begins and saccade ends - blackrock
sacbeg = find(diff(velBR > limBR) > 0);
sacend = find(diff(velBR > limBR) < 0);
if velBR(end)>limBR
    sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
end
if velBR(1)>limBR
    sacend = sacend(2:end);
end

if size(sacbeg,1)
    sacbeg=sacbeg';
    sacend=sacend';
end
    
% artifact = round([sacbeg(:) - artpadding*Fs sacend(:) + artpadding*Fs]);
artifactBR = round([sacbeg(:) sacend(:)]);

%detect saccade begins and saccade ends - python
sacbeg = find(diff(velPY > limPY) > 0);
sacend = find(diff(velPY > limPY) < 0);
if velPY(end)>limPY
    sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
end
if velPY(1)>limPY
    sacend = sacend(2:end);
end

if size(sacbeg,1)
    sacbeg=sacbeg';
    sacend=sacend';
end
    
% artifact = round([sacbeg(:) - artpadding*Fs sacend(:) + artpadding*Fs]);
artifactPY = round([sacbeg(:) sacend(:)]);

% figure;hold on
% plot(1:40000,eyesigBR(:,1:40000)')
% for k=1:find(artifactBR(:,2)>40000,1,'first')
%     line([artifactBR(k,1) artifactBR(k,1)],ylim,'Color','g')
%     line([artifactBR(k,2) artifactBR(k,2)],ylim,'Color','r')
% end
% 
% figure;hold on
% plot(1:40000,eyesigPY(:,1:40000)')
% for k=1:find(artifactPY(:,2)>40000,1,'first')
%     line([artifactPY(k,1) artifactPY(k,1)],ylim,'Color','g')
%     line([artifactPY(k,2) artifactPY(k,2)],ylim,'Color','r')
% end

% merge artifacts when inter-artifact interval is less than 10 ms

sacarrBR = [];
iai=(artifactBR(2:end,1)-artifactBR(1:end-1,2))>10;
sacdumcol1 = artifactBR(2:end,1);
sacdumcol2 = artifactBR(1:end-1,2);
sacarrBR(:,1) = [artifactBR(1,1); sacdumcol1(iai,1)];
sacarrBR(:,2) = [sacdumcol2(iai,1); artifactBR(end,2)];

sacarrPY = [];
iai=(artifactPY(2:end,1)-artifactPY(1:end-1,2))>10;
sacdumcol1 = artifactPY(2:end,1);
sacdumcol2 = artifactPY(1:end-1,2);
sacarrPY(:,1) = [artifactPY(1,1); sacdumcol1(iai,1)];
sacarrPY(:,2) = [sacdumcol2(iai,1); artifactPY(end,2)];

% % find saccade start/end in task time
% sacdum = [NS2_timestamp(artifact(:,1))' NS2_timestamp(artifact(:,2))'];

%%%%%%%%%%%%%%%%%%%%%%%
% find saccade start/end in task time
NS6sacarr = [NS2ts(sacarrBR(:,1))' NS2ts(sacarrBR(:,2))'];

sacdatA = nan(size(sacarr,1),12,501);
sacdatB = nan(size(sacarr,1),12,501);
sacdatC = nan(size(sacarr,1),12,501);

c=1;

t0 = clock;
ft_progress('init', 'etf',     'Please wait...');
for saclop=1:size(NS2sacarr,1)-1
    
    ft_progress(saclop/(size(NS2sacarr,1)), 'Processing event %d from %d', saclop, size(NS2sacarr,1));

    if NS2sacarr(saclop,1)>NS6_timestamp(1)+0.2
        ts1 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)-0.2);
        ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop+1,1));
        if length(ts1:ts2)<=501
            sacdatA(c,:,1:length(ts1:ts2)) = NS6a.Data(:,ts1:ts2);
            sacdatB(c,:,1:length(ts1:ts2)) = NS6b.Data(:,ts1:ts2);
            sacdatC(c,:,1:length(ts1:ts2)) = NS6c.Data(:,ts1:ts2);
        else
            ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)+0.3);
            sacdatA(c,:,:) = NS6a.Data(:,ts1:ts2);
            sacdatB(c,:,:) = NS6b.Data(:,ts1:ts2);
            sacdatC(c,:,:) = NS6c.Data(:,ts1:ts2);
        end
        c=c+1;
    end
    
end
ft_progress('close')
fprintf('total analysis time: %g\n',etime(clock,t0));

trlcnt = nan(1,size(sacdatA,3));
for timlop=1:size(squeeze(sacdatA(1,:,:)),2)
    trlcnt(timlop) = length(find(~isnan(squeeze(sacdatA(1,:,timlop)))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% create the data structure for Fieldtrip

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
    'C11'; 'C12'; 'eyeX';'eyeY';'posX';'posY'};
clear trial trl

% % plot eye and position data trial by trial to ensure accuracy
% close all
% for trllop = 1:length(eyedat)
%     
%     w=eyedat{trllop}(1,:);
%     ww=eyedatBR{trllop}(1,:);
%     x=w/(max(w)-min(w));
%     xx=ww/(max(ww)-min(ww));
%     
%     y=posdat{trllop}(1,:);
%     yy= posdatBR{trllop}(2,:); % BR position data is rotated relative to Panda
%     z=y/(max(y)-min(y));
%     zz=yy/(max(yy)-min(yy));
%     
%     figure
%     subplot(2,1,1);plot([x; xx]')
%     subplot(2,1,2);plot([z; zz]')
%     pause;close
%     
% end


%% save the complete data file

save(fullfile(savDir,[BRnam '_NSdat.mat']),'data')


