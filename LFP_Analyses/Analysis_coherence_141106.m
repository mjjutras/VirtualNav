%% specify filenames

% specify the directory containing log files
BRnam = 'JN140825010'; sessionDir = 'JN_14_08_25_13_12';
% BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';


logDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
trlDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\trial data';

load(fullfile(trlDir,[BRnam '_trldat.mat']))

%% remove common signal from each, and clear bad channels

% these channels are bad in the recordings on 14_08_25 (and maybe others
% around that time)
badchn = [2 5 7 10 12];

comsub = cell(1);
for trllop = 1:length(data.trial)
    
    trlavg = nanmean(data.trial{trllop}(setxor(1:length(data.label),badchn),:),1);
    
    for timlop = 1:length(data.time{trllop})
        comsub{trllop}(:,timlop) = data.trial{trllop}(:,timlop)-trlavg(timlop);
    end
    
    comsub{trllop}(badchn,:) = nan(length(badchn),size(comsub{trllop},2));
    
end

data.trial = comsub;
clear comsub trlavg

%% fix the time field

for trllop = 1:length(data.time)
    
    data.time{trllop} = 0:0.001:length(data.time{trllop})/1000-0.001;
    
end
    
%%


% do the spectral analysis - time-frequency
clear cfg
cfg.output      = 'powandcsd';
cfg.method      = 'mtmfft';
cfg.foilim      = [0 20];
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 4;
freq = ft_freqanalysis(cfg, data);

cfg.method  = 'coh';
stat = ft_connectivityanalysis(cfg, freq);

% for k=1:size(stat.cohspctrm,1)
%     figure(k)
%     imagesc(stat.time,stat.freq,squeeze(stat.cohspctrm(k,:,:)))
%     axis xy
%     title([stat.labelcmb{k,1} ' x ' stat.labelcmb{k,2}])
%     colorbar
%     pause
%     close
% end

save(fullfile('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Fieldtrip\',[BRnam '_freq_statTA.mat']),'freq','stat','cfg')

