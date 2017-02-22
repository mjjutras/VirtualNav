% open the NEV and corresponding NS6 file
NEV = openNEV('C:\Data\Blackrock_VR\JN140520001.nev');
NS6 = openNSx('C:\Data\Blackrock_VR\JN140520001.ns6');

% apply >750Hz high-pass filter to the continuous (NS6) signal
Fbp = 750;
N = 4;
Fn = 30000/2;
[B, A] = butter(N, Fbp/Fn,'high');
sig = filtfilt(B, A, double(NS6.Data(3,:)));

% according to MetaTag info, the NEV file starts 2 ms before the NS6 file
% I initially thought this needed to be accounted for; however, after
% closer inspection it looks like the files are already fairly well aligned
% without accounting for this "offset"
filofs = NS6.MetaTags.DateTimeRaw(end)-NEV.MetaTags.DateTimeRaw(end);

% channel 3 had good spikes; these are unsorted, so just identify any
% waveform from channel 3 (background noise waveforms will be included)
% spktim will include the timestamps for these waveforms
spktim = NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==3);
% create datwin array using a 48-sample window in the filtered NS6 signal,
% with the start time taken from spktim (NEV timestamps)
datwin = nan(length(spktim),48);
for k=1:length(spktim)
    datwin(k,:) = sig(spktim(k):spktim(k)+47);
end

% plot a few sample spikes from NS6 and the same from NEV
figure;plot(1:48,datwin(2:8,:))
xlabel('sample #')
ylabel('Voltage (uV)')
title('Spike waveforms from NS6 file')
hh=find(NEV.Data.Spikes.Electrode==3);
figure;plot(1:48,NEV.Data.Spikes.Waveform(:,hh(2:8)))
xlabel('sample #')
ylabel('Voltage (uV)')
title('Spike waveforms from NEV file')

% plot the same sample spike from NEV and from NS6, on the same axes
figure;hold on
plot(1:48,datwin(2,:),'b')
plot(1:48,NEV.Data.Spikes.Waveform(:,hh(2)),'r')
scatter(1:48,datwin(2,:),'.b')
scatter(1:48,NEV.Data.Spikes.Waveform(:,hh(2)),'.r')
xlabel('sample #')
ylabel('Voltage (uV)')
legend({'NS6' 'NEV'})
title('Spike #2')

figure;hold on
plot(1:48,datwin(8,:),'b')
plot(1:48,NEV.Data.Spikes.Waveform(:,hh(8)),'r')
scatter(1:48,datwin(8,:),'.b')
scatter(1:48,NEV.Data.Spikes.Waveform(:,hh(8)),'.r')
xlabel('sample #')
ylabel('Voltage (uV)')
legend({'NS6' 'NEV'})
title('Spike #8')