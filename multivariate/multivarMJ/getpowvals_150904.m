%% load variables

if strcmp(license,'375399')
    WSdir = 'C:\Data\VR';
elseif strcmp(license,'613743')
    WSdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\workspace';
end

load(fullfile(WSdir,'anaFruit_neuraldum_trainingses_newperffact_150706.mat'))


%% preprocess data

% exclude these channels
badchn = {'A02'; 'A05'; 'A08'; 'A10'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};

% exclude eye and position data for now
badchn = [badchn; {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'}];

cfg=[];
cfg.continuous    = 'no';
cfg.channel       = setxor(dum.label,badchn);
data0 = ft_preprocessing(cfg,dum);
clear dum


%% find trials with at least 5 seconds preceding banana
% this should be every trial that was included in nrlsel

bantim = nan(length(nrlsel)*2,1);
timsel = nan(length(nrlsel)*2,2);
nrltimsel = [];
trlgrp = [];
for invlop = 1:length(nrlsel)
    bantim(invlop*2-1:invlop*2) = bantimindall(nrlsel(invlop)-2:nrlsel(invlop)-1);
    if isempty(find((bantim(invlop*2-1:invlop*2))<5000,1))
        timsel(invlop*2-1:invlop*2,:) = [bantim(invlop*2-1:invlop*2)-4999 bantim(invlop*2-1:invlop*2)];
        nrltimsel = [nrltimsel; nrlsel(invlop)];
        trlgrp = [trlgrp; ones(2,1)+length(trlgrp)/2];
    else
        timsel(invlop*2-1:invlop*2,:) = repmat([nan nan],2,1);
    end
end


%% select last 5 seconds of approach path to banana in neural data
% also include 1 second after obtaining banana

data1 = data0;
for trllop = 1:length(data0.trial)
    if isempty(find(isnan(timsel(trllop,:)),1,'first')) % make sure trial wasn't deleted because of 5-second cutoff
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)+999];
        data1.time{trllop} = -4.999:0.001:1;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2)+1000);
    end
end

% remove marked trials
cfg=[];
% here choose only trials with at least 5 seconds pre-banana 
cfg.trials        = setxor(1:length(data1.trial),find(isnan(timsel(:,1))));
data1 = ft_preprocessing(cfg,data1);
clear data0


%% mark bad trials (CHANGE IF ADDING NEW DATA)
% identified with ft_rejectvisual, ft_databrowser

badtrldum = [16, 25, 30, 58, 92]; % if not re-referencing

badtrl = [];
invselexc = [];
for badlop = 1:length(badtrldum)
    badtrl = [badtrl; find(trlgrp==trlgrp(badtrldum(badlop)))];
    invselexc = [invselexc; find(trlgrp==trlgrp(badtrldum(badlop)),1,'last')/2];
end

nrlbadsel = nrltimsel(setxor(1:length(nrltimsel),invselexc));

% exclude bad trials, filter out line noise here
cfg=[];
cfg.trials = setxor(1:length(data1.time),badtrl);
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
data2 = ft_preprocessing(cfg,data1);
clear data1


%% create new data matrices for each 1-sec sample

cfg=[];
cfg.begsample = 1;
cfg.endsample = 1000;
data_1 = ft_redefinetrial(cfg,data2);

cfg.begsample = 1001;
cfg.endsample = 2000;
data_2 = ft_redefinetrial(cfg,data2);

cfg.begsample = 2001;
cfg.endsample = 3000;
data_3 = ft_redefinetrial(cfg,data2);

cfg.begsample = 3001;
cfg.endsample = 4000;
data_4 = ft_redefinetrial(cfg,data2);

cfg.begsample = 4001;
cfg.endsample = 5000;
data_5 = ft_redefinetrial(cfg,data2);

cfg.begsample = 5001;
cfg.endsample = 6000;
data_6 = ft_redefinetrial(cfg,data2);

cfg.begsample = 1;
cfg.endsample = 5000;
data_7 = ft_redefinetrial(cfg,data2);

clear data2


%% spectral analysis, all trials

% do the spectral analysis - time-averaged
cfg=[];
cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.foilim      = [1 30];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqL1 = ft_freqanalysis(cfg, data_1); % second 1
freqL2 = ft_freqanalysis(cfg, data_2); % second 2
freqL3 = ft_freqanalysis(cfg, data_3); % second 3
freqL4 = ft_freqanalysis(cfg, data_4); % second 4
freqL5 = ft_freqanalysis(cfg, data_5); % second 5
freqL6 = ft_freqanalysis(cfg, data_6); % second 6 (first second post-reward)
freqL7 = ft_freqanalysis(cfg, data_7); % whole 5-sec encoding period

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqH1 = ft_freqanalysis(cfg, data_1);
freqH2 = ft_freqanalysis(cfg, data_2);
freqH3 = ft_freqanalysis(cfg, data_3);
freqH4 = ft_freqanalysis(cfg, data_4);
freqH5 = ft_freqanalysis(cfg, data_5);
freqH6 = ft_freqanalysis(cfg, data_6);
freqH7 = ft_freqanalysis(cfg, data_7);

cfg=[];
cfg.keeptrials   = 'yes';
ftL1 = ft_freqdescriptives(cfg, freqL1);
ftL2 = ft_freqdescriptives(cfg, freqL2);
ftL3 = ft_freqdescriptives(cfg, freqL3);
ftL4 = ft_freqdescriptives(cfg, freqL4);
ftL5 = ft_freqdescriptives(cfg, freqL5);
ftL6 = ft_freqdescriptives(cfg, freqL6);
ftL7 = ft_freqdescriptives(cfg, freqL7);
ftH1 = ft_freqdescriptives(cfg, freqH1);
ftH2 = ft_freqdescriptives(cfg, freqH2);
ftH3 = ft_freqdescriptives(cfg, freqH3);
ftH4 = ft_freqdescriptives(cfg, freqH4);
ftH5 = ft_freqdescriptives(cfg, freqH5);
ftH6 = ft_freqdescriptives(cfg, freqH6);
ftH7 = ft_freqdescriptives(cfg, freqH7);

clear freq*


%% apply the correction factor based on the "excess path length" on visible trials

load('R:\Buffalo Lab\Mike\VirtualNav\MAT files\workspace\rat2d_visban_24ses.mat');

x = linspace(-11.2,11.2,81);
ctrs = x(1:end-1)+mean(gradient(x))/2;

memsel = nan(length(nrlbadsel),1);
for trllop = 1:length(nrlbadsel)
    banpossel = banposall(nrlbadsel(trllop),:);
    [~,pos1] = min(abs(ctrs-banpossel(1)));
    [~,pos2] = min(abs(ctrs-banpossel(2)));
    memsel(trllop) = excpthnrm(nrlbadsel(trllop))/rat2d(pos1,pos2);    
end

% normalize to maximum, and subtract from 1 to reverse the sign (high
% number = high memory performance)
memsel = 1-memsel/max(memsel);


%%

save('R:\Buffalo Lab\Mike\VirtualNav\multivariate\varsformultivarGiz150904.mat', ...
    'ftL1','ftL2','ftL3','ftL4','ftL5','ftL6','ftL7', ...
    'ftH1','ftH2','ftH3','ftH4','ftH5','ftH6','ftH7', ...
    'memsel');
