BRnam = 'JN140825011';
BRDir = 'C:\Data\VR\Blackrock';
% BRDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\Data';

NS6 = openNSx(fullfile(BRDir,[BRnam '.ns6']),'read','c:1');
NEV = openNEV(fullfile(BRDir,[BRnam '.nev']),'read');

NS6fs = 1/NS6.MetaTags.SamplingFreq;

NS6ts = NS6fs:NS6fs:size(NS6.Data,2)*NS6fs;

% get raw data aligned to spike time, channel 1, unit "1"
chn = 1;
    
spkind = NEV.Data.Spikes.Electrode==chn*NEV.Data.Spikes.Unit==1;
spktim = NEV.Data.Spikes.TimeStamp(spkind);

raw = int16(nan(length(spktim),2000));
ft_progress('init', 'etf',     'Please wait...');
for k = 1:length(spktim)
    ft_progress(k/length(spktim), 'Processing event %d from %d', k, length(spktim));
    dum = NS6.Data(1,(spktim(k)-1102):(spktim(k)+897));
    raw(k,:) = dum-mean(dum);
end
ft_progress('close')

spkwav = NEV.Data.Spikes.Waveform(:,spkind);
    

figure
subplot(211)
plot((-1000:999)/30000,mean(raw,1))
axis tight
% xh=xlim;
% xlim([0 xh(2)])
subplot(212)
plot((1:48)/30000,mean(spkwav,2))
axis tight

