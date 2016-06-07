%% run after getSaccadeAlignedLFPs160414.m

%%

dataD = [];
for trllop = 1:size(sacdatD,2)
    dum = squeeze(sacdatD(:,trllop,:));
    dum = dum(:,1:length(find(~isnan(dum(1,:)))));
    dataD.trial{trllop} = dum;
%     dataD.time{trllop} = -0.2:0.001:(size(dum,2)-201)/1000;
    dataD.time{trllop} = -1:0.001:1;
end

dataF = [];
for trllop = 1:size(sacdatF,2)
    dum = squeeze(sacdatF(:,trllop,:));
    dum = dum(:,1:length(find(~isnan(dum(1,:)))));
    dataF.trial{trllop} = dum;
%     dataF.time{trllop} = -0.2:0.001:(size(dum,2)-201)/1000;
    dataF.time{trllop} = -1:0.001:1;
end

dataH = [];
for trllop = 1:size(sacdatH,2)
    dum = squeeze(sacdatH(:,trllop,:));
    dum = dum(:,1:length(find(~isnan(dum(1,:)))));
    dataH.trial{trllop} = dum;
%     dataH.time{trllop} = -0.2:0.001:(size(dum,2)-201)/1000;
    dataH.time{trllop} = -1:0.001:1;
end

dataI = [];
for trllop = 1:size(sacdatI,2)
    dum = squeeze(sacdatI(:,trllop,:));
    dum = dum(:,1:length(find(~isnan(dum(1,:)))));
    dataI.trial{trllop} = dum;
%     dataI.time{trllop} = -0.2:0.001:(size(dum,2)-201)/1000;
    dataI.time{trllop} = -1:0.001:1;
end

dataD.sampleinfo = sampleinfo;
dataF.sampleinfo = sampleinfo;
dataH.sampleinfo = sampleinfo;
dataI.sampleinfo = sampleinfo;

dataD.fsample = 1000;
dataF.fsample = 1000;
dataH.fsample = 1000;
dataI.fsample = 1000;

dataD.label = {'D01'; 'D02'; 'D03'; 'D04'; 'D05'; 'D06'; 'D07'; 'D08'; 'D09'; 'D10'; 'D11'; 'D12'; 'D13'; 'D14'; 'D15'; 'D16'};
dataF.label = {'F01'; 'F02'; 'F03'; 'F04'; 'F05'; 'F06'; 'F07'; 'F08'; 'F09'; 'F10'; 'F11'; 'F12'; 'F13'; 'F14'; 'F15'; 'F16'};
dataH.label = {'H01'; 'H02'; 'H03'; 'H04'; 'H05'; 'H06'; 'H07'; 'H08'; 'H09'; 'H10'; 'H11'; 'H12'; 'H13'; 'H14'; 'H15'; 'H16'};
dataI.label = {'I01'; 'I02'; 'I03'; 'I04'; 'I05'; 'I06'; 'I07'; 'I08'; 'I09'; 'I10'; 'I11'; 'I12'; 'I13'; 'I14'; 'I15'; 'I16'};

% clear sacdat* srp*

data = ft_appenddata([], dataD, dataF, dataH, dataI);



%%

% now run ft_rejectvisual and ft_databrowser to clean bad trials

% badtrl = [73, 104, 158]; % if re-referencing to common probe average
% 
% filindchr = filindchrdum(setxor(1:length(filindchrdum),badtrl));

% exclude bad trials, filter out line noise here
cfg=[];
% cfg.trials = setxor(1:length(datac2.time),badtrl);
cfg.continuous    = 'no';
% cfg.dftfilter     = 'yes';
% cfg.dftfreq       = [60 120];
cfg.bsfilter     = 'yes';
cfg.bsfreq       = [58 62];
data2 = ft_preprocessing(cfg,data);

cfg=[];
% cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
cfg.output      = 'powandcsd';
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'no';
cfg.foilim      = [1 30];
cfg.taper       = 'hanning';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqlowD = ft_freqanalysis(cfg, data2);

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqhighD = ft_freqanalysis(cfg, data2);

%% look at time-resolved power (time-frequency analysis)

% do the spectral analysis - time-frequency
clear cfg
cfg.output      = 'powandcsd';
cfg.method      = 'mtmconvol';

cfg.foi         = 0:2:30;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.5 * ones(1,numfoi);
cfg.tapsmofrq   = 2   * ones(1,numfoi);
cfg.toi         = -0.2:0.01:0.5;
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'no';
tfL = ft_freqanalysis(cfg, data2);

cfg.foi         = 30:4:100;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.25 * ones(1,numfoi);
cfg.tapsmofrq   = 8   * ones(1,numfoi);
cfg.toi         = -4.875:0.01:0;
tfH = ft_freqanalysis(cfg, data2);




