BRnam = 'JN140825010';
datdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\NS6decimated\';

load(fullfile(datdir,[BRnam '_NS6_SF30.mat']))

%% run pEpisode on one channel

chnlop = 25;

samplerate = 1000;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

% get energy vector for each list
fprintf('Calc. Power...');
t0 = clock;

B=multienergyvec(double(NS6.Data(chnlop,:)),freqs,samplerate,width);
fprintf('%g\n',etime(clock,t0));

% calc the mean fit
Blog = log10(double(B));
Pm = mean(Blog,2);
clear Blog

% get the fit
fprintf('Calc. fit...');
t0 = clock;
% IMPORTANS: chi_squarefit assumes that frequencies are
% logarithmically spaced!
[all,R2] = chi_squarefit(freqs,Pm);
all = all';
fprintf('%g\n',etime(clock,t0));

% test fit of regression line
pv=polyfit(log10(freqs),Pm',1); % linear regression
%     figure;plot(freqs,Pm)
%     hold on;plot(freqs,pv(2)+pv(1)*log10(freqs),'r')

% set the threshold
thresh = all(:,951);

% loop through each frequency and save the unionvector
unionvector = int8(zeros(size(freqs,2),size(B,2)));
for f = 1:size(freqs,2)
    % get the pepisode
    unionvector(f,:) = int8(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0)); % 3-cycle threshold
end


%% or, run for multiple cycle thresholds if you have a lot of memory

clear unionvector
UVdum1 = int8(zeros(size(freqs,2),size(B,2)));
UVdum2 = int8(zeros(size(freqs,2),size(B,2)));
UVdum3 = int8(zeros(size(freqs,2),size(B,2)));
UVdum4 = int8(zeros(size(freqs,2),size(B,2)));
UVdum5 = int8(zeros(size(freqs,2),size(B,2)));

for f = 1:size(freqs,2)
    % get the pepisode
    UVdum1(f,:) = int8(episodeid(B(f,:),thresh(f),1*samplerate/freqs(f),shoulder,0)); % 1-cycle threshold
    UVdum2(f,:) = int8(episodeid(B(f,:),thresh(f),2*samplerate/freqs(f),shoulder,0)); % 2-cycle threshold
    UVdum3(f,:) = int8(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0)); % 3-cycle threshold
    UVdum4(f,:) = int8(episodeid(B(f,:),thresh(f),4*samplerate/freqs(f),shoulder,0)); % 4-cycle threshold
    UVdum5(f,:) = int8(episodeid(B(f,:),thresh(f),5*samplerate/freqs(f),shoulder,0)); % 5-cycle threshold
end

unionvector = UVdum1 + UVdum2 + UVdum3 + UVdum4 + UVdum5;
clear UVdum*

% convert unionvector from single to int8 to conserve memory
unionvector=int8(unionvector);

