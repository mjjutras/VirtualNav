%% specify filenames (need to do this manually for now)

% specify the directory containing log files
% BRnam = 'JN150209001'; sessionDir = 'JN_15_02_09\JN_15_02_09_15_34';
BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';


% logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\2014';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
savDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trial data';
NS6Dir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\Data';

% copy files to local drive if network transfer rate is slow
logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\2014';
% logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data';
% savDir = 'C:\Data\VR';
decDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated';

%% get trial info from MAT files

load(fullfile(savDir,[BRnam '_navdat.mat']))
trldat = data;
load(fullfile(savDir,[BRnam '_NSdat.mat']))

%% open NS6 and NEV files, create timestamp arrays that align each
% data structure (from NSdat file) uses sampleinfo corresponding to the
% decimated NS6 file; time info from decimated NS6 file needs to be
% adjusted to match the raw version
% trldat (from navdat file) uses logtime from the python file

NS6 = openNSx(fullfile(NS6Dir,[BRnam '.ns6']),'noread');
NEV = openNEV(fullfile(BRDir,[BRnam '.nev']),'read');
% dec = load(fullfile(decDir,[BRnam '_NS6_SF30.mat']));

ns6DTR = NS6.MetaTags.DateTimeRaw;
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
% incorporate offset for NS6 (NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq)
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

NS6fs = 1/NS6.MetaTags.SamplingFreq;
NEVfs = 1/double(NEV.MetaTags.SampleRes);
DECfs = 0.001;

NS6ts = datenum(ns6datevec)*86400:NS6fs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec-NS6fs;
NEVts = datenum(nevdatevec)*86400:NEVfs:datenum(nevdatevec)*86400+NEV.MetaTags.DataDurationSec-NEVfs;
DECts = datenum(ns6datevec)*86400:DECfs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec;

% check that length of timestamp array matches # of datapoints:
length(NS6ts)==NS6.MetaTags.DataPoints
length(NEVts)==NEV.MetaTags.DataDuration

% align to first sample of NS6s; this brings everything onto the same time
% scale (time values line up across different data streams with different
% start points/sampling rates, everything based on the start of the NS6
% file)
NEVts = (NEVts-NS6ts(1))+NS6fs;
DECts = (DECts-NS6ts(1))+NS6fs;
NS6ts = (NS6ts-NS6ts(1))+NS6fs;

%% import python log file info: banana locations

[logtime_ban, id_ban, ~, x_ban, y_ban] = textread(fullfile(logDir,sessionDir,'banana.par'), '%f%s%f%f%f'); % avatar (position) data
% [logtime_ban, id_ban, ~, x_ban, y_ban] = textread(fullfile(logDir,sessionDir,'fruit.par'), '%f%s%f%f%f'); % avatar (position) data


%% create trials in NEV file

c=1;
m=1000;
clear NEVper
NEVper.val={};
NEVper.tim={};
for vallop = 1:length(NEV.Data.SerialDigitalIO.UnparsedData)
    if NEV.Data.SerialDigitalIO.UnparsedData(vallop)==m
        m=m+1;
        for k=1:find(NEV.Data.SerialDigitalIO.UnparsedData(vallop:end)==m,1,'first')-1
            NEVper.val{c}(k,1) = NEV.Data.SerialDigitalIO.UnparsedData(vallop+k-1);
            NEVevtTSind = NEV.Data.SerialDigitalIO.TimeStamp(vallop+k-1);
            NEVper.tim{c}(k,1) = NEVts(NEVevtTSind); % NS6 clock time
        end
        c=c+1;
        vallop = vallop + find(NEV.Data.SerialDigitalIO.UnparsedData(vallop:end)==m,1,'first') - 1;
    end
end


%%

for trllop = 1:length(data.trial)
    
    % difference between these two should be the same
    blkts = data.sampleinfo(trllop,1):data.sampleinfo(trllop,2);
    
    bantim_trl = logtime_ban(logtime_ban>=trlPyt(1)&logtime_ban<=trlPyt(2));
    id_trl = id_ban(logtime_ban>=trlPyt(1)&logtime_ban<=trlPyt(2));
    banx_trl = x_ban(logtime_ban>=trlPyt(1)&logtime_ban<=trlPyt(2));
    bany_trl = y_ban(logtime_ban>=trlPyt(1)&logtime_ban<=trlPyt(2));
    
    % find timestamps corresponding with banana eaten
    tsEat = bantim_trl(strcmp(id_trl,'eaten'));
    for k=1:length(tsEat)
        
        % index the timestamp in trldat.time for the trial (corresponds
        % with the event aligned to the Python eye data)
        [~,tsPytInd] = min(abs(trldat.time{trllop}-tsEat(k)));
        
        % find "NS6 time" for the event, using the data structure,
        % find index in the decimated NS6 file, use this to determine the
        % timestamp
        NS6tim = DECts(blkts(tsPytInd));
        
        % find the timestamp in NS6ts
        [~,NS6timInd] = min(abs(NS6ts-NS6tim));
        
        NEV.Data.SerialDigitalIO.TimeStampSec
    
    
    
    
%% determine avatar position when collecting bananas;
% compare between python log and blackrock data file

h=DECts(data.sampleinfo(1,1)):0.001:DECts(data.sampleinfo(1,2));
bantim = NEVper.tim{1}(NEVper.val{1}==200);
ban_x=[]; ban_y=[];
for banlop = 1:length(bantim)
    [~,banind] = min(abs(h-bantim(banlop)));
    ban_x(banlop) = data.trial{1}(39,banind);
    ban_y(banlop) = data.trial{1}(40,banind);
end
    