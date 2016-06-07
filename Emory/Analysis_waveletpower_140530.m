%% ns3 file

fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';

openNSx(fullfile(fildir,'JN140520002.ns3'),'read')

% samplerate=NS6.MetaTags.SamplingFreq;
samplerate=NS2.MetaTags.SamplingFreq;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

wavpowgrp=[];
for chnlop = 1:size(NS2.Data,1)
% for chnlop = 1:size(NS3.Data,1)
    
%     lfp = decimate(double(NS6.Data(chnlop,:)),30);
    lfp = double(NS2.Data(chnlop,400000:500000));

    B=single(multienergyvec(lfp,freqs,samplerate,width));

    wavpowgrp=[wavpowgrp; mean(B,2)'];
    
end

lfpavg = mean(NS2.Data(1:36,400000:500000),1);

wavpowgrp_subavg=[];
for chnlop = 1:size(NS2.Data,1)
% for chnlop = 1:size(NS3.Data,1)
    
%     lfp = decimate(double(NS6.Data(chnlop,:)),30);
    lfp = double(NS2.Data(chnlop,400000:500000))-lfpavg;

    B=single(multienergyvec(lfp,freqs,samplerate,width));

    wavpowgrp_subavg=[wavpowgrp_subavg; mean(B,2)'];
    
end

figure;plot(freqs,wavpowgrp(1:12,:)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('A-1','A-2','A-3','A-4','A-5','A-6','A-7','A-8','A-9','A-10','A-11','A-12')

set(gca,'Position',[1     1   457   382])

figure;plot(freqs,wavpowgrp(13:24,:)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('B-1','B-2','B-3','B-4','B-5','B-6','B-7','B-8','B-9','B-10','B-11','B-12')

figure;plot(freqs,wavpowgrp(25:36,:)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('C-1','C-2','C-3','C-4','C-5','C-6','C-7','C-8','C-9','C-10','C-11','C-12')


%% ns6 file

fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';

openNSx(fullfile(fildir,'JN140619002.ns6'),'read','precision','double','skipfactor',30)

% samplerate=NS6.MetaTags.SamplingFreq;
samplerate=NS6.MetaTags.SamplingFreq/30; % skip factor
% samplerate=NS3.MetaTags.SamplingFreq;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

wavpowgrp=[];
for chnlop = 1:size(NS6.Data,1)
% for chnlop = 1:size(NS3.Data,1)
    
%     lfp = decimate(double(NS6.Data(chnlop,:)),30);
%     lfp = double(NS3.Data(chnlop,:));
    lfp = double(NS6.Data(chnlop,:));

    B=single(multienergyvec(lfp,freqs,samplerate,width));

    wavpowgrp=[wavpowgrp; mean(B,2)'];
    
end