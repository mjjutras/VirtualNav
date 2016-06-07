% BRtoFT.m
%
% imports neural (LFP) data from Blackrock NS6 data (decimated to 1000 Hz)
% then uses trldat structure imported using PARtoTRLDAT to align the Python
% log info on a trial by trial basis
%
% the resulting data structure uses Fieldtrip format to facilitate analysis
%
% 150421 Mike Jutras (with parts adapted from quickLineup by Yoni Browning)
% (previous versions of code contained in BR&Panda_to_FT, analyzeNav140418)
%
%
% this code requires ft_progress.m (or you can comment out any line with a
% function starting with 'ft_')
%
% Notes:
% python samples position/direction at ~83 Hz; eye data @240
% timestamp values in panda log are already in ms (1 kHz)
%
% use of NEV file isn't required to align the data, but is a helpful check
% to verify the alignment of the behavioral events between the Python and
% Blackrock logs

%% specify filenames (need to do this manually for now)

% specify the directory containing log files
BRnam = 'JN150410003'; 
sessionDir = 'JN_BR_15_04_10\JN_BR_15_04_10_13_11';


logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
trlDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\trial data';
decDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated';

% copy files to local drive if network transfer rate is slow
% BRDir = 'C:\Data\VR';
% savDir = 'C:\Data\VR';
% decDir = 'C:\Data\VR';


%% load the behavioral data file

[~,sesnam]=fileparts(sessionDir);

load(fullfile(trlDir,[sesnam '_trldat.mat']))

% trldat now contains all of the relevant information from the python log


%% create array marking first "eat/reward" event in each trial

inieatarr = [];
timarr = [];
for trllop = 1:length(trldat.time)
    timarr = [timarr trldat.time{trllop}];
    dum = zeros(size(trldat.time{trllop}));
    dum(trldat.time{trllop}==trldat.frttim{trllop}(1)) = 1;
    inieatarr = [inieatarr dum];
end


%% load the raw data (NS6, decimated version); align the NEV (timestamp) and NS (continuous) signals

% open NEV file
NEV = openNEV(fullfile(BRDir,[BRnam '.nev']),'read'); % NEV file; contains event codes and timestamps

% open decimated NS6 file
load(fullfile(decDir,[BRnam '_NS6_SF30.mat'])) % decimated NS6 file; 1000 Hz

% open NS2 file
NS2 = openNSx(fullfile(BRDir,[BRnam '.NS2']),'read','precision','double'); % NS2 file (contains eye & position data)

NS2fs = 1/NS2.MetaTags.SamplingFreq; % sampling frequency of NS2 file (should be 0.001)
DECfs = 0.001; % sampling frequency of decimated NS6 file
NS6fs = 1/NS6.MetaTags.SamplingFreq; % sampling frequency of NS6 file (should also be 1/30000)

% find first non-zero value in NS6 file
nonzer = find(NS6.Data(1,:)~=0,1,'first');

% % this should be true
% length(nonzer:30:size(NS6b.Data,2))+(nonzer-1)==size(NS6.Data,2)

% the first nonzero value is the same for NS6b and decimated NS6
% I think this is a bug

% create time arrays for each file; timestamps match data samples
% NEV timestamps will be the same as those in NS6
NS6ts = NS6fs:NS6fs:(NEV.MetaTags.DataDuration*NS6fs);
DECts = [nan(1,nonzer-1) NS6ts(nonzer):DECfs:(size(NS2.Data,2)-nonzer)*NS2fs+NS6ts(nonzer)];
NS2ts = NS2fs:NS2fs:(size(NS2.Data,2)*NS2fs);

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

% for JN150209001, fix weird event where 1000 appears twice instead of 1001
if strcmp(BRnam,'JN150209001')
    NEV.Data.SerialDigitalIO.UnparsedData(2202)=1001;
end

c=1;
m=1000; % trial onset codes start at 1000 and progress incrementally
clear NEVper
NEVper.val={};
NEVper.tim={};
for vallop = 1:length(NEV.Data.SerialDigitalIO.UnparsedData)
% for vallop = 2142:length(NEV.Data.SerialDigitalIO.UnparsedData) % for JN150209001
    if NEV.Data.SerialDigitalIO.UnparsedData(vallop)==m
        
        m=m+1;
        for k=1:find(ismember(NEV.Data.SerialDigitalIO.UnparsedData(vallop+1:end),[1000 m]),1,'first')
            
            % build NEVper.val, each cell contains one trial, each row is an event
            NEVper.val{c}(k,1) = NEV.Data.SerialDigitalIO.UnparsedData(vallop+k-1);
            
            % NEVevtTSind: NEV timestamp corresponding with the event code
            NEVevtTSind = NEV.Data.SerialDigitalIO.TimeStamp(vallop+k-1);
            
            % actual timestamp is given in "NS6 clock time"
            NEVper.tim{c}(k,1) = NS6ts(NEVevtTSind);
            
        end
        
        if NEV.Data.SerialDigitalIO.UnparsedData(find(ismember(NEV.Data.SerialDigitalIO.UnparsedData(vallop+1:end),[1000 m]),1,'first') + vallop)==1000
            m=1000;
        end
        
        c=c+1;
        
        vallop = find(ismember(NEV.Data.SerialDigitalIO.UnparsedData(vallop+1:end),[1000 m]),1,'first') + vallop;
        
    end
end


%% get first "eat" event from each trial in NEV file

inieatarrNEV = zeros(size(DECts));
for trllop = 1:length(NEVper.val)
    
    if ~isempty(find(NEVper.val{trllop}==200,1,'first'))

        [~,inieatnevind] = min(abs(DECts-NEVper.tim{trllop}(find(NEVper.val{trllop}==200,1,'first'))));
        inieatarrNEV(inieatnevind) = 1;
        
    end
    
end

[Xs, lag] = xcorr(inieatarr, inieatarrNEV);
[~,Is] = max(Xs);
lag(Is)

% calculate lag between each Python event and NEV event, after accounting
% for total lag
tl1 = (1/1000:1/1000:length(inieatarr)/1000);
tl2 = (1/1000:1/1000:length(inieatarrNEV)/1000)+(lag(Is)/1000);
ts1 = tl1(logical(inieatarr));
ts2 = tl2(logical(inieatarrNEV));
lagarr = nan(length(ts1),2);
for trllop = 1:length(ts1)
    [~,dum] = min(abs(ts2-ts1(trllop)));
    lagarr(trllop,:) = [ts1(trllop) ts2(dum)];
end

%% line up behavioral data with raw LFPs (adapted from quickLineup)

% find eye & position signals in NS2
eyeindx = find(strncmpi(NS2lab,'ainp1',5));
eyeindy = find(strncmpi(NS2lab,'ainp2',5));
posindx = find(strncmpi(NS2lab,'ainp3',5));
posindy = find(strncmpi(NS2lab,'ainp4',5));

eyedatNS2 = NS2.Data([eyeindx eyeindy],:);
posdatNS2 = NS2.Data([posindx posindy],:);

ft_defaults

clear trial sampleinfo eyedatBR posdatBR
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(trldat.eyedat)
    ft_progress(trllop/length(trldat.eyedat), 'Processing event %d from %d', trllop, length(trldat.eyedat));
    
    [Xs1, lag] = xcorr(eyedatNS2(1,:), trldat.eyedat{trllop}(1,:));
    [Xs2, lag] = xcorr(eyedatNS2(2,:), trldat.eyedat{trllop}(2,:));
    [~,Is] = max(Xs1.*Xs2);
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

save(fullfile(trlDir,[sesnam '_NSdat.mat']),'data')


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

frtdifmat = [];
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(data.trial)
    ft_progress(trllop/length(data.trial), 'Processing event %d from %d', trllop, length(data.trial));
    
    [~,NEVbnd1] = min(abs(NS6ts-DECts(data.sampleinfo(trllop,1))));
    [~,NEVbnd2] = min(abs(NS6ts-DECts(data.sampleinfo(trllop,2))));

    NEVtrlind = find(NEV.Data.SerialDigitalIO.TimeStamp>=NEVbnd1 & ...
        NEV.Data.SerialDigitalIO.TimeStamp<=NEVbnd2);
    
    % hopefully there are only 10 occurences of bananas being eaten!
    frtind = find(NEV.Data.SerialDigitalIO.UnparsedData(NEVtrlind)==200);
    
    frtTS = NEV.Data.SerialDigitalIO.TimeStamp(NEVtrlind(frtind));
    frttim = NS6ts(frtTS);
    frtdif = (frttim - DECts(data.sampleinfo(trllop,1)))*1000; % convert to ms
    
    if length(trldat.frttim{trllop})==length(frtdif)
        frtdifmat = [frtdifmat; frtdif' trldat.frttim{trllop}-trldat.time{trllop}(1)];
    end
    
end
ft_progress('close')

median(diff(frtdifmat,1,2)) % 66.4

figure;subplot(2,1,1);hold on
scatter(1:size(frtdifmat,1),diff(frtdifmat,1,2))

