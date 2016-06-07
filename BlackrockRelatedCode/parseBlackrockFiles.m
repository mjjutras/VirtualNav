% MJJ 6/1/2016
% This code appears to have been created to "parse" the Matlab data in
% order to analyze periods during acute recordings where spikes were
% observed. These recordings were performed by lowering the electrodes
% through the hippocampus and pausing a specific electrode as soon as
% spiking activity was observed; the notes would specify the times that
% this would occur. The raw data could then be opened in KlustaKwik to
% perform sorting of the spike waveforms.

%% JN140620002

% Blackrock filename base
BRnam = 'JN140620002';
matdir = fullfile('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\parsed continuous data',BRnam); % directory to save the parsed data in MAT format
if ~exist(matdir)
    mkdir('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\parsed continuous data',BRnam)
end
datdir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';

% load the BR data
NS2 = openNSx(fullfile(datdir,[BRnam '.ns2']),'read');
NS6 = openNSx(fullfile(datdir,[BRnam '.ns6']),'read');
NEV = openNEV(fullfile(datdir,[BRnam '.nev']),'read');

ns2DTR = NS2.MetaTags.DateTimeRaw;
ns6DTR = NS6.MetaTags.DateTimeRaw;
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

%  date operations:
%  1 % day
%  1/24 % hour
%  1/1440 % minute (24*60)
%  1/86400 % second (24*60*60)

% find the last timepoint at which data is recorded in all three file (ID
% the file that cuts out the earliest)
[~,i]=min([NS2.MetaTags.DataDurationSec NS6.MetaTags.DataDurationSec NEV.MetaTags.DataDurationSec]);
if i==1
    lastsample_datevec = datevec(datenum(ns2datevec)+NS2.MetaTags.DataDurationSec/86400);
elseif i==2
    lastsample_datevec = datevec(datenum(ns6datevec)+NS6.MetaTags.DataDurationSec/86400);
elseif i==3
    lastsample_datevec = datevec(datenum(nevdatevec)+NEV.MetaTags.DataDurationSec/86400);
end

% parse each individual channel
for chnlop = 1:size(NS6.Data,1)

    clear spikestart_* spikeend_* data timestamp
    
    if chnlop == 1
        % Channel 1:
        % According to notes, spikes start at T15:37 (3:37 pm) (next line to verify, comment #1)
        % datestr(datenum(nevdatevec)+(NEV.Data.Comments.TimeStampSec(1)/86400))
        % Start sampling data at T15:35 pm
        spikestart_datevec = [nevdatevec(1:3) 15 35 0];
        spikestart_sec_ns2 = (datenum(spikestart_datevec)-datenum(ns2datevec))*86400;
        spikestart_sec_ns6 = (datenum(spikestart_datevec)-datenum(ns6datevec))*86400;
        spikestart_sec_nev = (datenum(spikestart_datevec)-datenum(nevdatevec))*86400;
        % End sampling on the last timestamp where sampling occurs in all
        % data formats
        spikeend_datevec = lastsample_datevec;
        spikeend_sec_ns2 = (datenum(lastsample_datevec)-datenum(ns2datevec))*86400;
        spikeend_sec_ns6 = (datenum(lastsample_datevec)-datenum(ns6datevec))*86400;
        spikeend_sec_nev = (datenum(lastsample_datevec)-datenum(nevdatevec))*86400;
    elseif chnlop == 2
        % Channel 2:
        % According to notes, spikes probably start ~ T15:25 (3:25 pm)
        spikestart_datevec = [nevdatevec(1:3) 15 25 0];
        spikestart_sec_ns2 = (datenum(spikestart_datevec)-datenum(ns2datevec))*86400;
        spikestart_sec_ns6 = (datenum(spikestart_datevec)-datenum(ns6datevec))*86400;
        spikestart_sec_nev = (datenum(spikestart_datevec)-datenum(nevdatevec))*86400;
        % End sampling on the last timestamp where sampling occurs in all
        % data formats
        spikeend_datevec = lastsample_datevec;
        spikeend_sec_ns2 = (datenum(lastsample_datevec)-datenum(ns2datevec))*86400;
        spikeend_sec_ns6 = (datenum(lastsample_datevec)-datenum(ns6datevec))*86400;
        spikeend_sec_nev = (datenum(lastsample_datevec)-datenum(nevdatevec))*86400;
    elseif chnlop == 3
        % Channel 3:
        % Start sampling on the first timestamp where sampling occurs in all
        % data formats [sampling rates])
        spikestart_datevec = datevec(max([datenum(ns2datevec) datenum(ns6datevec) datenum(nevdatevec)]));
        spikestart_sec_ns2 = (max([datenum(ns2datevec) datenum(ns6datevec) datenum(nevdatevec)])-datenum(ns2datevec))*86400;
        spikestart_sec_ns6 = (max([datenum(ns2datevec) datenum(ns6datevec) datenum(nevdatevec)])-datenum(ns6datevec))*86400;
        spikestart_sec_nev = (max([datenum(ns2datevec) datenum(ns6datevec) datenum(nevdatevec)])-datenum(nevdatevec))*86400;
        % According to notes, spikes end by T15:51 (next line to verify, comment #2)
        % datestr(datenum(nevdatevec)+(NEV.Data.Comments.TimeStampSec(2)/86400))
        % End sampling data at T15:51 pm
        spikeend_datevec = [nevdatevec(1:3) 15 51 0];
        spikeend_sec_ns2 = (datenum(spikeend_datevec)-datenum(ns2datevec))*86400;
        spikeend_sec_ns6 = (datenum(spikeend_datevec)-datenum(ns6datevec))*86400;
        spikeend_sec_nev = (datenum(spikeend_datevec)-datenum(nevdatevec))*86400;
    elseif chnlop == 4
        % Channel 4:
        % According to notes, spikes might have started around T15:50
        % Start sampling data at T15:50 pm
        spikestart_datevec = [nevdatevec(1:3) 15 50 0];
        spikestart_sec_ns2 = (datenum(spikestart_datevec)-datenum(ns2datevec))*86400;
        spikestart_sec_ns6 = (datenum(spikestart_datevec)-datenum(ns6datevec))*86400;
        spikestart_sec_nev = (datenum(spikestart_datevec)-datenum(nevdatevec))*86400;
        % End sampling on the last timestamp where sampling occurs in all
        % data formats
        spikeend_datevec = lastsample_datevec;
        spikeend_sec_ns2 = (datenum(lastsample_datevec)-datenum(ns2datevec))*86400;
        spikeend_sec_ns6 = (datenum(lastsample_datevec)-datenum(ns6datevec))*86400;
        spikeend_sec_nev = (datenum(lastsample_datevec)-datenum(nevdatevec))*86400;
    end
    
    % Create the raw data array and timestamp array
    % Timestamps: serial date number * 86400 (units: seconds)
    samp_start = round(spikestart_sec_ns6*NS6.MetaTags.SamplingFreq)+1;
    samp_end = round(spikeend_sec_ns6*NS6.MetaTags.SamplingFreq)+1;
    data = double(NS6.Data(chnlop,samp_start:samp_end));
    timestamp = datenum(spikestart_datevec)*86400:1/NS6.MetaTags.SamplingFreq:datenum(spikeend_datevec)*86400;

    % save the raw data and timestamps
    save(fullfile(matdir,['data_' NS6.ElectrodesInfo(chnlop).Label '.mat']), 'data')
    save(fullfile(matdir,['timestamp_' NS6.ElectrodesInfo(chnlop).Label '.mat']), 'timestamp')
    
end

%%
datdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\parsed continuous data\JN140620002';
load(fullfile([datdir '\chan4'],'data_chan4.mat'))
data=data';

fid=fopen(fullfile([datdir '\chan4'],'data_chan4.dat'),'wb');
count = fwrite(fid,data,'int16')
fclose('all')
