%% 8/25/14

load('R:\Buffalo Lab\Mike\VirtualNav\MAT files\trial data\JN140825011_navdat.mat')
trldat = data;
load('R:\Buffalo Lab\Mike\VirtualNav\MAT files\trial data\JN140825011_trldat.mat')

for k=1:length(data.time)
    data.time{k}=(data.time{k}-data.time{k}(1))/1000;
end

% bad channels based on ft_rejectartifact
badchn = {'A02'; 'A05'; 'A10'; 'A12'};

% exclude eye and position data
badchn = [badchn; {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'}];

% preprocess: remove bad channels
cfg=[];
cfg.continuous    = 'no';
cfg.channel       = setxor(data.label,badchn);
data0 = ft_preprocessing(cfg,data);

% limit to first 3 seconds
cfg=[];
cfg.toilim    = [0 3];
data0 = ft_redefinetrial(cfg, data0);

% remove bad trials (after visual inspection)
badtrl = [21, 22, 29];
cfg=[];
cfg.trials = setxor(1:length(data0.time),badtrl);
data0 = ft_preprocessing(cfg,data0);

% data1: subtract common average of each probe from all channels on that probe
data1 = data0;
clear prblaball
for k=1:length(data0.label)
    prblaball(k,1)=data0.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);
for trllop = 1:length(data0.trial)
    for prblop = 1:length(prblab)
        prbavg = mean(data0.trial{trllop}(locb==prblop,:),1);
        repavg = repmat(prbavg,length(find(locb==prblop)),1);
        mat1 = reshape(data0.trial{trllop}(locb==prblop,:),1,length(find(locb==prblop)),size(data0.trial{trllop},2));
        mat2 = reshape(repavg,1,size(repavg,1),size(repavg,2));
        data1.trial{trllop}(locb==prblop,:)=diff(cat(1,mat1,mat2),1,1);
    end
end

% filter out line noise here, demean & detrend
cfg=[];
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
cfg.demean        = 'yes';
cfg.detrend       = 'yes';
data0 = ft_preprocessing(cfg,data0);
data1 = ft_preprocessing(cfg,data1);

% do the spectral analysis - time-averaged
cfg=[];
cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.foilim      = [1 30];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqL0 = ft_freqanalysis(cfg, data0);
freqL1 = ft_freqanalysis(cfg, data1);

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqH0 = ft_freqanalysis(cfg, data0);
freqH1 = ft_freqanalysis(cfg, data1);

cfg=[];
cfg.jackknife   = 'yes';
fdL0 = ft_freqdescriptives(cfg, freqL0);
fdH0 = ft_freqdescriptives(cfg, freqH0);
fdL1 = ft_freqdescriptives(cfg, freqL1);
fdH1 = ft_freqdescriptives(cfg, freqH1);


