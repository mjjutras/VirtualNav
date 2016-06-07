function data = BRtoFT(sesDir,BRnam)

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
%
% CHANGE LOG:
% 4/25/15 MJJ - turned into a function; took out NEV and all related code
% (this was used to check alignment between reward events in the different
% signals; can check this after running)
% 4/26/15 MJJ - fixed eye data alignment, since there were quite a few
% instances where xcorr wasn't working on its own; the new way is slow and
% clunky but it seems to work.


logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
trlDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\trial data';
% decDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\NS6 - decimated';
% decDir = 'C:\Data\VR';
decDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated';

% load the behavioral data file; trldat contains info from python log
[~,sesnam]=fileparts(sesDir);
disp(['Loading ' fullfile(trlDir,[sesnam '_trldat.mat'])])
load(fullfile(trlDir,[sesnam '_trldat.mat']))

% open decimated NS6 file (1000 Hz)
disp(['Loading ' fullfile(decDir,[BRnam '_NS6_SF30.mat'])])
load(fullfile(decDir,[BRnam '_NS6_SF30.mat']))

% open NS2 file (contains eye & position data)
disp(['Loading ' fullfile(BRDir,[BRnam '.NS2'])])
NS2 = openNSx(fullfile(BRDir,[BRnam '.NS2']),'read','precision','double');

% create label arrays for NS2 and NS6 channels
NS2lab = cell(1,length(NS2.ElectrodesInfo));
for lablop = 1:length(NS2.ElectrodesInfo)
    NS2lab{lablop} = NS2.ElectrodesInfo(lablop).Label;
end
NS6lab = cell(1,length(NS6.ElectrodesInfo));
for lablop = 1:length(NS6.ElectrodesInfo)
    NS6lab{lablop} = NS6.ElectrodesInfo(lablop).Label;
end

% line up behavioral data with raw LFPs (adapted from quickLineup)

% find eye & position signals in NS2
eyeindx = find(strncmpi(NS2lab,'ainp1',5));
eyeindy = find(strncmpi(NS2lab,'ainp2',5));
posindx = find(strncmpi(NS2lab,'ainp3',5));
posindy = find(strncmpi(NS2lab,'ainp4',5));

eyedatNS2 = NS2.Data([eyeindx eyeindy],:);
posdatNS2 = NS2.Data([posindx posindy],:);

clear trial sampleinfo eyedatBR posdatBR
disp('Aligning eye data signals')
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(trldat.eyedat)
    ft_progress(trllop/length(trldat.eyedat), 'Processing event %d from %d', trllop, length(trldat.eyedat));
    
    if trllop==1
        [Xs1, lag] = xcorr(eyedatNS2(1,:), trldat.eyedat{trllop}(1,:));
        Xs2 = xcorr(eyedatNS2(1,:), trldat.eyedat{trllop}(1,300:end));
        Xs3 = xcorr(eyedatNS2(2,:), trldat.eyedat{trllop}(2,:));
        Xs4 = xcorr(eyedatNS2(2,:), trldat.eyedat{trllop}(2,300:end));
        [~,Is1] = max(Xs1.*(Xs1-Xs2).*Xs3.*(Xs3-Xs4));

        [~,locs] = findpeaks(Xs1.*Xs3,1000);
        [~,Is2] = min(abs(lag(round(locs*1000))-lag(Is1)));    

        trl_start = lag(round(locs(Is2)*1000))+1;
        trl_end = trl_start+length(trldat.eyedat{trllop})-1;
    else
        lagind1 = trl_start;
        lagind2 = trl_end+(size(trldat.eyedat{trllop}(1,:),2)*2);
        [Xs1, lag] = xcorr(eyedatNS2(1,lagind1:lagind2), trldat.eyedat{trllop}(1,:));
        Xs2 = xcorr(eyedatNS2(2,lagind1:lagind2), trldat.eyedat{trllop}(2,:));
        [~,Is1] = max(Xs1.*Xs2);

        trl_start = lag(Is1)+lagind1-1;
        trl_end = trl_start+length(trldat.eyedat{trllop})-1;
    end
    
    % Save as a new file structure containing all the necessary info.
    trial{trllop} = double(NS6.Data(:,trl_start:trl_end));
    
    % append the eye & position data from the Python log file
    trial{trllop}(end+1:end+2,:) = trldat.eyedat{trllop};
    trial{trllop}(end+1:end+2,:) = trldat.posdat{trllop};
    
    % append the eye and position data from Blackrock file
    trial{trllop}(end+1:end+2,:) = eyedatNS2(:,trl_start:trl_end);
    trial{trllop}(end+1:end+2,:) = posdatNS2(:,trl_start:trl_end);
    
    % include both Blackrock and Python analog data for now, until we
    % decide which version has cleaner A2D conversion

    % corresponds to samples in decimated NS6 file
    sampleinfo(trllop,:) = [trl_start trl_end];
    
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

