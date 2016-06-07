fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';

%% Array A

openNSx(fullfile(fildir,'JN140821002.ns6'),'read','c:1:12','t:10500000:10620000','s:30')
% this opens the file starting at sample # 10500000 and progressing for 120
% seconds (only when implementing the skip factor of 30)

datofs = []; % data offset
for k=1:size(NS6.Data,1)
    datofs(k,:) = NS6.Data(k,1:5000)+(k-1)*1000;
end
figure;plot(datofs')
legend('A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12')

datofs = []; % data offset
for k=1:size(NS6.Data,1)
    datofs(k,:) = NS6.Data(k,1:1000)+(k-1)*1000;
end
figure;plot(datofs')
legend('A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12')

samplerate=1000;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

wavpowgrp=[];
for chnlop = 1:size(NS6.Data,1)

    lfp = double(NS6.Data(chnlop,:));

    B=single(multienergyvec(lfp,freqs,samplerate,width));

    wavpowgrp=[wavpowgrp; mean(B,2)'];
    
end

figure;plot(freqs,wavpowgrp(1:12,:)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12')

%% Array B

openNSx(fullfile(fildir,'JN140821002.ns6'),'read','c:13:24','t:10500000:10620000','s:30')
% this opens the file starting at sample # 10500000 and progressing for 120
% seconds (only when implementing the skip factor of 30)

datofs = []; % data offset
for k=1:size(NS6.Data,1)
    datofs(k,:) = NS6.Data(k,1:5000)+(k-1)*1000;
end
figure;plot(datofs')
legend('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12')

datofs = []; % data offset
for k=1:size(NS6.Data,1)
    datofs(k,:) = NS6.Data(k,1:1000)+(k-1)*1000;
end
figure;plot(datofs')
legend('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12')

samplerate=1000;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

wavpowgrp=[];
for chnlop = 1:size(NS6.Data,1)

    lfp = double(NS6.Data(chnlop,:));

    B=single(multienergyvec(lfp,freqs,samplerate,width));

    wavpowgrp=[wavpowgrp; mean(B,2)'];
    
end

figure;plot(freqs,wavpowgrp(1:12,:)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12')

%% Array C

openNSx(fullfile(fildir,'JN140821002.ns6'),'read','c:25:36','t:10500000:10620000','s:30')
% this opens the file starting at sample # 10500000 and progressing for 120
% seconds (only when implementing the skip factor of 30)

datofs = []; % data offset
for k=1:size(NS6.Data,1)
    datofs(k,:) = NS6.Data(k,1:5000)+(k-1)*1000;
end
figure;plot(datofs')
legend('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12')

datofs = []; % data offset
for k=1:size(NS6.Data,1)
    datofs(k,:) = NS6.Data(k,1:1000)+(k-1)*1000;
end
figure;plot(datofs')
legend('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12')

samplerate=1000;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

wavpowgrp=[];
for chnlop = 1:size(NS6.Data,1)

    lfp = double(NS6.Data(chnlop,:));

    B=single(multienergyvec(lfp,freqs,samplerate,width));

    wavpowgrp=[wavpowgrp; mean(B,2)'];
    
end

figure;plot(freqs,wavpowgrp(1:12,:)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12')

