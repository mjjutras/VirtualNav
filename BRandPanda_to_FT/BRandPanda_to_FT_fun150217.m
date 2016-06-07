% BR&Panda_to_FT.m
%
% create a Fieldtrip-format data structure using the following method:
% 
% - import neural (LFP) data from Blackrock NS6 file (decimated to 1000 Hz)
% - import eye data from Blackrock NS2 file
% - import data from parsed python files (.par)
%   these files contain eye, position and velocity data all using the same
%   timescale, but with samples of different data types falling at different
%   timespoints; the timescale is resampled at 1000 Hz to match sampling rate 
%   of neural recordings, and additional data points between samples are 
%   filled in using interpolation.
%
%
% 141017 Mike Jutras (with parts adapted from quickLineup by Yoni Browning)
%
%
% before running, must run Python script MakeParFiles in the directory containing 
% the Python log files
% (e.g. R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\JN_14_08_25_13_12)

% Notes:
% revised in part from analyzeNav140418 (MJ)
% python samples position/direction at ~83 Hz; eye data @240
%
% timestamp values in panda log are already in ms

%% specify filenames (need to do this manually for now)

% specify the directory containing log files
BRnam = 'JN150209001'; sessionDir = 'JN_15_02_09\JN_15_02_09_15_34';
% BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';


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
for trllop = 1:length(log_trl_starttime)-1
    
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
NEV = openNEV(fullfile(BRDir,[BRnam '.nev']));
% % NS6 = openNSx(fullfile(datdir,[BRnam '.ns6']));
% % load decimated MAT version of NS6 file instead of directly from NS6 file
load(fullfile(decDir,[BRnam '_NS6_SF30.mat']))

ns2DTR = NS2.MetaTags.DateTimeRaw;
% ns6DTR = NS6.MetaTags.DateTimeRaw; % (NS6 same as NS2, but Timestamp indicates slightly delayed start time)
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
% ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

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
data.trl = trl;
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

for k=1:size(data.trial{1},1)
    data.label{k} = NS2.ElectrodesInfo(k).Label(1:find(double(NS2.ElectrodesInfo(k).Label)==0,1,'first')-1);
end

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

% function i = nearest(array, val);
% 
% if val>max(array)
%   % return the last occurence of the nearest number
%   [dum, i] = max(flipud(array));
%   i = length(array) + 1 - i;
% else
%   % return the first occurence of the nearest number
%   [mindist, i] = min(abs(array(:) - val));
% end

% function NS6 = decimateNS6_30kto1k(BRnam,fildir,savdir,chans);
% 
% % BRnam = 'JN140825010';
% % fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
% % savdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\NS6 - decimated';
% % chans = 1:36;
% 
% if exist('openNSx','file')
%     NS6 = openNSx(fullfile(fildir,[BRnam '.ns6']),'read','channels',chans,'skipfactor',30);
%     save(fullfile(savdir,[BRnam '_NS6_SF30.mat']),'NS6')
%     disp([fullfile(savdir,[BRnam '_NS6_SF30.mat']) ' successfully saved'])
% else
%     disp('Cannot find openNSx in MATLAB path; cannot create decimated data file')
% end

% function B=inpaint_nans(A,method)
% % INPAINT_NANS: in-paints over nans in an array
% % usage: B=INPAINT_NANS(A)          % default method
% % usage: B=INPAINT_NANS(A,method)   % specify method used
% %
% % Solves approximation to one of several pdes to
% % interpolate and extrapolate holes in an array
% %
% % arguments (input):
% %   A - nxm array with some NaNs to be filled in
% %
% %   method - (OPTIONAL) scalar numeric flag - specifies
% %       which approach (or physical metaphor to use
% %       for the interpolation.) All methods are capable
% %       of extrapolation, some are better than others.
% %       There are also speed differences, as well as
% %       accuracy differences for smooth surfaces.
% %
% %       methods {0,1,2} use a simple plate metaphor.
% %       method  3 uses a better plate equation,
% %                 but may be much slower and uses
% %                 more memory.
% %       method  4 uses a spring metaphor.
% %       method  5 is an 8 neighbor average, with no
% %                 rationale behind it compared to the
% %                 other methods. I do not recommend
% %                 its use.
% %
% %       method == 0 --> (DEFAULT) see method 1, but
% %         this method does not build as large of a
% %         linear system in the case of only a few
% %         NaNs in a large array.
% %         Extrapolation behavior is linear.
% %
% %       method == 1 --> simple approach, applies del^2
% %         over the entire array, then drops those parts
% %         of the array which do not have any contact with
% %         NaNs. Uses a least squares approach, but it
% %         does not modify known values.
% %         In the case of small arrays, this method is
% %         quite fast as it does very little extra work.
% %         Extrapolation behavior is linear.
% %
% %       method == 2 --> uses del^2, but solving a direct
% %         linear system of equations for nan elements.
% %         This method will be the fastest possible for
% %         large systems since it uses the sparsest
% %         possible system of equations. Not a least
% %         squares approach, so it may be least robust
% %         to noise on the boundaries of any holes.
% %         This method will also be least able to
% %         interpolate accurately for smooth surfaces.
% %         Extrapolation behavior is linear.
% %
% %       method == 3 --+ See method 0, but uses del^4 for
% %         the interpolating operator. This may result
% %         in more accurate interpolations, at some cost
% %         in speed.
% %
% %       method == 4 --+ Uses a spring metaphor. Assumes
% %         springs (with a nominal length of zero)
% %         connect each node with every neighbor
% %         (horizontally, vertically and diagonally)
% %         Since each node tries to be like its neighbors,
% %         extrapolation is as a constant function where
% %         this is consistent with the neighboring nodes.
% %
% %       method == 5 --+ See method 2, but use an average
% %         of the 8 nearest neighbors to any element.
% %         This method is NOT recommended for use.
% %
% %
% % arguments (output):
% %   B - nxm array with NaNs replaced
% %
% %
% % Example:
% %  [x,y] = meshgrid(0:.01:1);
% %  z0 = exp(x+y);
% %  znan = z0;
% %  znan(20:50,40:70) = NaN;
% %  znan(30:90,5:10) = NaN;
% %  znan(70:75,40:90) = NaN;
% %
% %  z = inpaint_nans(znan);
% %
% %
% % See also: griddata, interp1
% %
% % Author: John D'Errico
% % e-mail address: woodchips@rochester.rr.com
% % Release: 2
% % Release date: 4/15/06
% 
% 
% % I always need to know which elements are NaN,
% % and what size the array is for any method
% [n,m]=size(A);
% A=A(:);
% nm=n*m;
% k=isnan(A(:));
% 
% % list the nodes which are known, and which will
% % be interpolated
% nan_list=find(k);
% known_list=find(~k);
% 
% % how many nans overall
% nan_count=length(nan_list);
% 
% % convert NaN indices to (r,c) form
% % nan_list==find(k) are the unrolled (linear) indices
% % (row,column) form
% [nr,nc]=ind2sub([n,m],nan_list);
% 
% % both forms of index in one array:
% % column 1 == unrolled index
% % column 2 == row index
% % column 3 == column index
% nan_list=[nan_list,nr,nc];
% 
% % supply default method
% if (nargin<2) || isempty(method)
%     method = 0;
% elseif ~ismember(method,0:5)
%     error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
% end
% 
% % for different methods
% switch method
%     case 0
%         % The same as method == 1, except only work on those
%         % elements which are NaN, or at least touch a NaN.
%         
%         % horizontal and vertical neighbors only
%         talks_to = [-1 0;0 -1;1 0;0 1];
%         neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
%         
%         % list of all nodes we have identified
%         all_list=[nan_list;neighbors_list];
%         
%         % generate sparse array with second partials on row
%         % variable for each element in either list, but only
%         % for those nodes which have a row index > 1 or < n
%         L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
%         nl=length(L);
%         if nl>0
%             fda=sparse(repmat(all_list(L,1),1,3), ...
%                 repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
%                 repmat([1 -2 1],nl,1),nm,nm);
%         else
%             fda=spalloc(n*m,n*m,size(all_list,1)*5);
%         end
%         
%         % 2nd partials on column index
%         L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(all_list(L,1),1,3), ...
%                 repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
%                 repmat([1 -2 1],nl,1),nm,nm);
%         end
%         
%         % eliminate knowns
%         rhs=-fda(:,known_list)*A(known_list);
%         k=find(any(fda(:,nan_list(:,1)),2));
%         
%         % and solve...
%         B=A;
%         B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
%         
%     case 1
%         % least squares approach with del^2. Build system
%         % for every array element as an unknown, and then
%         % eliminate those which are knowns.
%         
%         % Build sparse matrix approximating del^2 for
%         % every element in A.
%         % Compute finite difference for second partials
%         % on row variable first
%         [i,j]=ndgrid(2:(n-1),1:m);
%         ind=i(:)+(j(:)-1)*n;
%         np=(n-2)*m;
%         fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
%             repmat([1 -2 1],np,1),n*m,n*m);
%         
%         % now second partials on column variable
%         [i,j]=ndgrid(1:n,2:(m-1));
%         ind=i(:)+(j(:)-1)*n;
%         np=n*(m-2);
%         fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
%             repmat([1 -2 1],np,1),nm,nm);
%         
%         % eliminate knowns
%         rhs=-fda(:,known_list)*A(known_list);
%         k=find(any(fda(:,nan_list),2));
%         
%         % and solve...
%         B=A;
%         B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
%         
%     case 2
%         % Direct solve for del^2 BVP across holes
%         
%         % generate sparse array with second partials on row
%         % variable for each nan element, only for those nodes
%         % which have a row index > 1 or < n
%         L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
%         nl=length(L);
%         if nl>0
%             fda=sparse(repmat(nan_list(L,1),1,3), ...
%                 repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
%                 repmat([1 -2 1],nl,1),n*m,n*m);
%         else
%             fda=spalloc(n*m,n*m,size(nan_list,1)*5);
%         end
%         
%         % 2nd partials on column index
%         L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
%                 repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
%                 repmat([1 -2 1],nl,1),n*m,n*m);
%         end
%         
%         % fix boundary conditions at extreme corners
%         % of the array in case there were nans there
%         if ismember(1,nan_list(:,1))
%             fda(1,[1 2 n+1])=[-2 1 1];
%         end
%         if ismember(n,nan_list(:,1))
%             fda(n,[n, n-1,n+n])=[-2 1 1];
%         end
%         if ismember(nm-n+1,nan_list(:,1))
%             fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
%         end
%         if ismember(nm,nan_list(:,1))
%             fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
%         end
%         
%         % eliminate knowns
%         rhs=-fda(:,known_list)*A(known_list);
%         
%         % and solve...
%         B=A;
%         k=nan_list(:,1);
%         B(k)=fda(k,k)\rhs(k);
%         
%     case 3
%         % The same as method == 0, except uses del^4 as the
%         % interpolating operator.
%         
%         % del^4 template of neighbors
%         talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
%             0 1;0 2;1 -1;1 0;1 1;2 0];
%         neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
%         
%         % list of all nodes we have identified
%         all_list=[nan_list;neighbors_list];
%         
%         % generate sparse array with del^4, but only
%         % for those nodes which have a row & column index
%         % >= 3 or <= n-2
%         L = find( (all_list(:,2) >= 3) & ...
%             (all_list(:,2) <= (n-2)) & ...
%             (all_list(:,3) >= 3) & ...
%             (all_list(:,3) <= (m-2)));
%         nl=length(L);
%         if nl>0
%             % do the entire template at once
%             fda=sparse(repmat(all_list(L,1),1,13), ...
%                 repmat(all_list(L,1),1,13) + ...
%                 repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
%                 repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
%         else
%             fda=spalloc(n*m,n*m,size(all_list,1)*5);
%         end
%         
%         % on the boundaries, reduce the order around the edges
%         L = find((((all_list(:,2) == 2) | ...
%             (all_list(:,2) == (n-1))) & ...
%             (all_list(:,3) >= 2) & ...
%             (all_list(:,3) <= (m-1))) | ...
%             (((all_list(:,3) == 2) | ...
%             (all_list(:,3) == (m-1))) & ...
%             (all_list(:,2) >= 2) & ...
%             (all_list(:,2) <= (n-1))));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(all_list(L,1),1,5), ...
%                 repmat(all_list(L,1),1,5) + ...
%                 repmat([-n,-1,0,+1,n],nl,1), ...
%                 repmat([1 1 -4 1 1],nl,1),nm,nm);
%         end
%         
%         L = find( ((all_list(:,2) == 1) | ...
%             (all_list(:,2) == n)) & ...
%             (all_list(:,3) >= 2) & ...
%             (all_list(:,3) <= (m-1)));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(all_list(L,1),1,3), ...
%                 repmat(all_list(L,1),1,3) + ...
%                 repmat([-n,0,n],nl,1), ...
%                 repmat([1 -2 1],nl,1),nm,nm);
%         end
%         
%         L = find( ((all_list(:,3) == 1) | ...
%             (all_list(:,3) == m)) & ...
%             (all_list(:,2) >= 2) & ...
%             (all_list(:,2) <= (n-1)));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(all_list(L,1),1,3), ...
%                 repmat(all_list(L,1),1,3) + ...
%                 repmat([-1,0,1],nl,1), ...
%                 repmat([1 -2 1],nl,1),nm,nm);
%         end
%         
%         % eliminate knowns
%         rhs=-fda(:,known_list)*A(known_list);
%         k=find(any(fda(:,nan_list(:,1)),2));
%         
%         % and solve...
%         B=A;
%         B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
%         
%     case 4
%         % Spring analogy
%         % interpolating operator.
%         
%         % list of all springs between a node and a horizontal
%         % or vertical neighbor
%         hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
%         hv_springs=[];
%         for i=1:4
%             hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
%             k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
%             hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
%         end
%         
%         % delete replicate springs
%         hv_springs=unique(sort(hv_springs,2),'rows');
%         
%         % build sparse matrix of connections, springs
%         % connecting diagonal neighbors are weaker than
%         % the horizontal and vertical springs
%         nhv=size(hv_springs,1);
%         springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
%             repmat([1 -1],nhv,1),nhv,nm);
%         
%         % eliminate knowns
%         rhs=-springs(:,known_list)*A(known_list);
%         
%         % and solve...
%         B=A;
%         B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
%         
%     case 5
%         % Average of 8 nearest neighbors
%         
%         % generate sparse array to average 8 nearest neighbors
%         % for each nan element, be careful around edges
%         fda=spalloc(n*m,n*m,size(nan_list,1)*9);
%         
%         % -1,-1
%         L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % 0,-1
%         L = find(nan_list(:,3) > 1);
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % +1,-1
%         L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % -1,0
%         L = find(nan_list(:,2) > 1);
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % +1,0
%         L = find(nan_list(:,2) < n);
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % -1,+1
%         L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % 0,+1
%         L = find(nan_list(:,3) < m);
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % +1,+1
%         L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
%         nl=length(L);
%         if nl>0
%             fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
%                 repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
%                 repmat([1 -1],nl,1),n*m,n*m);
%         end
%         
%         % eliminate knowns
%         rhs=-fda(:,known_list)*A(known_list);
%         
%         % and solve...
%         B=A;
%         k=nan_list(:,1);
%         B(k)=fda(k,k)\rhs(k);
%         
% end
% 
% % all done, make sure that B is the same shape as
% % A was when we came in.
% B=reshape(B,n,m);
% 
% 
% % ====================================================
% %      end of main function
% % ====================================================
% % ====================================================
% %      begin subfunctions
% % ====================================================
% function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% % identify_neighbors: identifies all the neighbors of
% %   those nodes in nan_list, not including the nans
% %   themselves
% %
% % arguments (input):
% %  n,m - scalar - [n,m]=size(A), where A is the
% %      array to be interpolated
% %  nan_list - array - list of every nan element in A
% %      nan_list(i,1) == linear index of i'th nan element
% %      nan_list(i,2) == row index of i'th nan element
% %      nan_list(i,3) == column index of i'th nan element
% %  talks_to - px2 array - defines which nodes communicate
% %      with each other, i.e., which nodes are neighbors.
% %
% %      talks_to(i,1) - defines the offset in the row
% %                      dimension of a neighbor
% %      talks_to(i,2) - defines the offset in the column
% %                      dimension of a neighbor
% %
% %      For example, talks_to = [-1 0;0 -1;1 0;0 1]
% %      means that each node talks only to its immediate
% %      neighbors horizontally and vertically.
% %
% % arguments(output):
% %  neighbors_list - array - list of all neighbors of
% %      all the nodes in nan_list
% 
% if ~isempty(nan_list)
%     % use the definition of a neighbor in talks_to
%     nan_count=size(nan_list,1);
%     talk_count=size(talks_to,1);
%     
%     nn=zeros(nan_count*talk_count,2);
%     j=[1,nan_count];
%     for i=1:talk_count
%         nn(j(1):j(2),:)=nan_list(:,2:3) + ...
%             repmat(talks_to(i,:),nan_count,1);
%         j=j+nan_count;
%     end
%     
%     % drop those nodes which fall outside the bounds of the
%     % original array
%     L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m);
%     nn(L,:)=[];
%     
%     % form the same format 3 column array as nan_list
%     neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
%     
%     % delete replicates in the neighbors list
%     neighbors_list=unique(neighbors_list,'rows');
%     
%     % and delete those which are also in the list of NaNs.
%     neighbors_list=setdiff(neighbors_list,nan_list,'rows');
%     
% else
%     neighbors_list=[];
% end

