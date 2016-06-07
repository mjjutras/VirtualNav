
fid='JN130531.2'

dataset=strcat('R:\Buffalo lab\Virtual Navigation\Free Foraging\NEX Files\',fid,'.nex');

ft_defaults

% get the header info
hdr = ft_read_header(dataset);

% set the trial structure to start and end each trial with the appearance
% of the bananas (event code 100)
cfg=[];
cfg.trialfun      = 'trialfun_nav';
cfg.dataset       = dataset;
cfg = ft_definetrial(cfg);


%code to find spike/LFP/eye channels
[lfpind,spkind,eyeind,namarr]=getchannels(dataset);
labarr=[namarr(lfpind,:); namarr(eyeind,:)]; % look at only LFPs and eye data for now

% preprocess the data, removing 60 Hz noise artifacts (and 120 Hz harmonic)
cfg.channel       = cellstr(labarr)';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
cfg.padding       = 10;
cfg.continuous    = 'yes';
cfg.detrend       = 'no';

data = ft_preprocessing(cfg);

% plot the first trial eye data
figure;plot(data.time{1},data.trial{1}(5:6,:)')
xlim(minmax(data.time{1}))
title('Plexon eye data')

%% find the trial in the Bananarchy log that matches the first trial in the Plexon data
[times, id, x, y] = textread('R:\Buffalo Lab\Virtual Navigation\Free Foraging\parsed log files\session_1191\bzflag2.par', '%f%s%f%f');

% read the events
event = ft_read_event(dataset);

numevt = length(event);
for evtlop = 1:numevt
  mrk.val(evtlop) = event(evtlop).value;
  mrk.tim(evtlop) = event(evtlop).sample;
end

log_trl_start = find(strcmp(id,'100'));
log_trl_length = diff(log_trl_start)*(1000/240);

trl_start_tim = mrk.tim(mrk.val==100);
trl_length = diff(trl_start_tim);

trlind = [];
for k=1:length(trl_length)
    trlind(k) = ft_nearest(log_trl_length,trl_length(k));
end

% check trlind to make sure the first trial is part of the sequence of all
% trials (some of the trials do not match the sequence)

time_eye = times(strcmp(id,'eye'));
x_real = x(strcmp(id,'eye'));
y_real = y(strcmp(id,'eye'));

first_trial_ind = ft_nearest(time_eye,times(log_trl_start(trlind(1)))):ft_nearest(time_eye,times(log_trl_start(trlind(2))));


eyedat_first_trial = [x_real(first_trial_ind-1) y_real(first_trial_ind-1)];

% plot the first trial eye data from Bananarchy log and match it to thelexon eye data
figure;plot(time_eye(first_trial_ind-1), eyedat_first_trial)
xlim(minmax(time_eye(first_trial_ind-1)'))
title('Bananarchy eye data')

%% quick spectral analysis

% do the spectral analysis - time-frequency
clear cfg
cfg.channel     = data.label(1:2);
cfg.trials      = 1:10;
cfg.output      = 'pow';
cfg.method      = 'mtmconvol';
cfg.foi         = 0:2:20;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.5 * ones(1,numfoi);
cfg.tapsmofrq   = 4   * ones(1,numfoi);
cfg.toi         = 0:0.1:10;
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
freq = ft_freqanalysis(cfg, data);


% plot channel 1 raw LFPsignal, trial 1, and power spectrum below
figure;subplot(2,1,1);plot(data.time{1},data.trial{1}(1,:));xlim([0 10]);set(gca,'Position',[0.13 0.583837209302326 0.651190476190476 0.341162790697674]);ylabel('voltage')
subplot(2,1,2);imagesc(freq.time,freq.freq,squeeze(freq.powspctrm(1,1,:,:)));axis xy;xlabel('time');ylabel('freq');colorbar

%% load the avatar data

[times, id, x, y] = textread('R:\Buffalo Lab\Virtual Navigation\Free Foraging\parsed log files\session_1191\avatar.par', '%f%s%f%f');

% unique(id)
% 
% ans = 
% 
%     'accel'
%     'dir'
%     'pos'
%     'speed'
%     't_accel'
%     't_speed'
    
posind = find(strcmp(id,'pos'));
xpos = x(posind);
ypos = y(posind);

dirind = find(strcmp(id,'dir'));
dir = x(dirind);

% calculate movement velocity
velx = diff(xpos);
vely = diff(ypos);
vel = sqrt(velx.^2+vely.^2);


%plot velocity in blue, direction in red
figure;plot(vel);
hold on;plot(dir(1:end-1)/1000,'r')
xlim([0 2000])

%% cross-frequency coupling

srate=data.fsample;

for chnlop = 1:4

    datbuf=[];
    % for k=1:length(data.trial)
    for k=1:5
        datbuf=[datbuf data.trial{k}(chnlop,:)];
    end

    Vlo=eegfiltMJ(datbuf,srate,3,12);
    Vhi=eegfiltMJ(datbuf,srate,90,140);

    figure
    [r, r_CI, nCtlPts] = GLM_CFC_for_paper(Vlo, Vhi, 10);
    title([fid '; r = ' num2str(r)])

    % Canolty method
    figure
    % [mod2d, flow, fhigh] = modulation_index_2d(datbuf, srate);
    [mod2d, flow, fhigh] = modulation_index_2d(datbuf, srate);
    title(fid)

end
