% getspikesfromNS6: uses the spike timestamps in the NEV file to get the
% raw data around each spike in the NS6 file, align across spikes to get an
% average raw signal
% 
% this verifies that the timing information is aligned correctly between
% the NEV and NS6 files for the purposes of designating event-aligned
% raw LFP activity
% 
% 150320 MJJ

%% specify filenames (need to do this manually for now)

% specify the directory containing log files
% BRnam = 'JN150209001'; sessionDir = 'JN_15_02_09\JN_15_02_09_15_34';
BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';


logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data';
% BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
% savDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\trial data';
% decDir = 'R:\Buffalo Lab\Virtual Navigation\NS6 - decimated';

% copy files to local drive if network transfer rate is slow
BRDir = 'C:\Data\VR';
savDir = 'C:\Data\VR';
decDir = 'C:\Data\VR';

%%

load(fullfile(savDir,[BRnam '_NSdat.mat']))
NEV = openNEV(fullfile(BRDir,[BRnam '.nev']),'read');
NS6 = openNSx(fullfile(BRDir,[BRnam '.ns6']),'read','c:1');

%%

ns6DTR = NS6.MetaTags.DateTimeRaw; % (NS6 same as NS2, but Timestamp indicates slightly delayed start time)
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

NS6fs = 1/NS6.MetaTags.SamplingFreq;
NEVfs = 1/double(NEV.MetaTags.SampleRes);

NS6ts = datenum(ns6datevec)*86400:NS6fs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec-NS6fs;
NEVts = datenum(nevdatevec)*86400:NEVfs:datenum(nevdatevec)*86400+NEV.MetaTags.DataDurationSec-NEVfs;

% check that these are 1:
length(NS6ts)==NS6.MetaTags.DataPoints
length(NEVts)==NEV.MetaTags.DataDuration

% align to first sample of NEV (should come before first sample of NS6)
NS6ts = NS6ts-NEVts(1);
NEVts = NEVts-NEVts(1);

% Beth saw units on A1,A3,B5,B8,C5 ([1 3 17 20 29])
% however, running
% "unique(NEV.Data.Spikes.Electrode(NEV.Data.Spikes.Unit==1))"
% gives spikes on [1 4 7 8 11 44]

% spikes on all electrodes
rawavgarr = nan(6,250);
spkwavarr = nan(6,48);
c=1;
for chnlop = [1 4 7 8 11 44]
    
    disp(['Channel ' num2str(chnlop)])
    
    spkind = NEV.Data.Spikes.Electrode==chnlop*NEV.Data.Spikes.Unit==1;
    spktim = NEV.Data.Spikes.TimeStamp(spkind);

    raw = int16(nan(length(spktim),250));
    ft_progress('init', 'etf', 'Please wait...');
    for k = 1:length(spktim)
    
        ft_progress(k/length(spktim), 'Processing channel %d from %d', k, length(spktim));
        
        [~,NS6spktim] = min(abs(NS6ts-NEVts(spktim(k))));

        dum = NS6.Data(1,(NS6spktim-125):(NS6spktim+124));

        raw(k,:) = dum-mean(dum);

    end
    ft_progress('close')
    
    rawavgarr(c,:) = mean(raw,1);
    
    spkwav = NEV.Data.Spikes.Waveform(:,spkind);
    spkwavarr(c,:) = mean(spkwav,2)';
    
    c=c+1;

end

save('C:\Data\VR\JN140825011_spikeraw.mat','rawavgarr','spkwavarr')

% figure;plot((-250:249)/30000,mean(raw,1))

