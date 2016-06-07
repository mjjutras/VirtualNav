
fid='JN130604.1'

dataset=strcat('R:\Buffalo Lab\Virtual Navigation\Recording Data\NEX Files\',fid,'.nex');

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

% preprocess the data, removing 60 Hz noise artifacts (and 120 Hz harmonic)
cfg.channel       = cellstr(namarr(lfpind,:));
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
cfg.padding       = 10;
cfg.continuous    = 'yes';
cfg.detrend       = 'no';
data_LFP = ft_preprocessing(cfg);

% cfg.channel       = cellstr(namarr(spkind,:));
% cfg.dftfilter     = 'no';
% data_spike = ft_preprocessing(cfg);

cfg.channel       = cellstr(namarr(eyeind,:));
cfg.dftfilter     = 'no';
data_eye = ft_preprocessing(cfg);

% data = ft_appenddata([], data_spike, data_LFP, data_eye);
data = ft_appenddata([], data_LFP, data_eye);
clear data_*

% % plot the first trial eye data
% figure;plot(data.time{1},data.trial{1}([strmatch('X',data.label) strmatch('Y',data.label)],:)')
% xlim(minmax(data.time{1}))
% title('Plexon eye data')

%% make a new information matrix with timestamps matching the NEX file, and containing eye position,
% spatial position, and heading direction

% panda samples position/direction at ~83 Hz; eye data @240
% interpolate these data points to match NEX sampling rate of 1000 Hz
% (position data will be ok, panda eye data will be used only to ensure how
% precisely in time the data matches between the two data files)

% timestamp values in panda log are already in ms

% import panda log file info
[times_eye, id_eye, x_eye, y_eye] = textread('R:\Buffalo Lab\Virtual Navigation\Free Foraging\Giz\session_1196\eye.par', '%f%s%f%f');
[times_ava, id_ava, x_ava, y_ava] = textread('R:\Buffalo Lab\Virtual Navigation\Free Foraging\Giz\session_1196\avatar.par', '%f%s%f%f');

time_eye = times_eye(strcmp(id_eye,'eye'));
xeye_real = x_eye(strcmp(id_eye,'eye'));
yeye_real = y_eye(strcmp(id_eye,'eye'));

time_pos = times_ava(strcmp(id_ava,'pos'));
xpos_real = x_ava(strcmp(id_ava,'pos'));
ypos_real = y_ava(strcmp(id_ava,'pos'));

time_dir = times_ava(strcmp(id_ava,'dir'));
dir_real = x_ava(strcmp(id_ava,'dir'));

time_vel = times_ava(strcmp(id_ava,'speed'));
vel_real = x_ava(strcmp(id_ava,'speed'));

% read the events from the nex file
event = ft_read_event(dataset);

numevt = length(event);
for evtlop = 1:numevt
  mrk.val(evtlop) = event(evtlop).value;
  mrk.tim(evtlop) = event(evtlop).sample;
end

log_trl_start = find(strcmp(id_eye,'100'));
log_trl_length = diff(times_eye(log_trl_start));

trl_start_tim = mrk.tim(mrk.val==100);
trl_length = diff(trl_start_tim);

trlind = [];
for k=1:length(trl_length)
    trlind(k) = ft_nearest(log_trl_length,trl_length(k));
end

% look in trlind to see where the matching trials fall in sequence

diftrlind = diff(trlind);
diftrlind_findones = diftrlind==1;
ll=zeros(size(diftrlind_findones));
for k=1:length(diftrlind_findones)-1
    if diftrlind_findones(k)==1
        for j=k:length(diftrlind_findones)-1
            if diftrlind_findones(j)==1 && diftrlind_findones(j+1)==1
                ll(k)=ll(k)+1;
            elseif diftrlind_findones(j+1)==0
                ll(k)=ll(k)+1;
                break
            end
        end
    end
end

% this tells you the first trial in the longest sequence of good matches
besttrlind_pand = trlind(ll==max(ll)); % best trial to line up, in panda log file
besttrlind_nex = find(ll==max(ll)); % best trial to line up, in nex log

% plot in blue, pulses where the timestamps in the NEX file mark trials
trlstartones_nex = logical(zeros(1,hdr.nSamples));
trlstartones_nex(trl_start_tim) = 1;
figure;hold on
plot(1:hdr.nSamples,trlstartones_nex)
% plot in red, pulses where the timestamps in the panda file mark trials
trlstartones_pand = logical(zeros(1,length(time_eye)));
trlstartones_pand(log_trl_start) = 1;
figure;hold on
plot(time_eye,trlstartones_pand,'r')


reftrial_pandtime = times_eye(log_trl_start(besttrlind_pand)); % reference trial start in "panda time"
reftrial_nextime = trl_start_tim(besttrlind_nex); % reference trial start in "NEX time"
newpandtime_dum = time_eye+(reftrial_nextime-reftrial_pandtime);
newpandtime_pos_dum = time_pos+(reftrial_nextime-reftrial_pandtime);
newpandtime_dir_dum = time_dir+(reftrial_nextime-reftrial_pandtime);
newpandtime_vel_dum = time_vel+(reftrial_nextime-reftrial_pandtime);


% figure;hold on
% plot(newpandtime_dum,trlstartones_pand,'r')
% plot(1:hdr.nSamples,trlstartones_nex)

% adjust based on the delay between the trial starts in the new panda time array and the nex file
trlstart_newpand_dum = times_eye(log_trl_start)+(reftrial_nextime-reftrial_pandtime);
delay_arr = [];
for k=1:length(trl_start_tim)
    delay_arr(k)=trlstart_newpand_dum(ft_nearest(trlstart_newpand_dum,trl_start_tim(k)))-trl_start_tim(k);
end

% % plot x eye data for whole NEX recording, both files, on same graph
% figure;hold on
% for trllop = 1:length(data.time)
%     plot(data.cfg.previous{1}.trl(trllop,1):data.cfg.previous{1}.trl(trllop,2), ...
%         data.trial{trllop}(strmatch('X',data.label),:)')
% end
% plot(newpandtime_dum(ft_nearest(newpandtime_dum,data.cfg.previous{1}.trl(1,1)):ft_nearest(newpandtime_dum,data.cfg.previous{1}.trl(end,2))), ...
%     xeye_real(ft_nearest(newpandtime_dum,data.cfg.previous{1}.trl(1,1)):ft_nearest(newpandtime_dum,data.cfg.previous{1}.trl(end,2)))*3.6,'r')
% 
% for k=1:length(trl_start_tim)
%     % black: plexon trial start time
%     line([trl_start_tim(k) trl_start_tim(k)],ylim,'Color','k');
%     % red: panda trial start time
%     line([trlstart_newpand_dum(ft_nearest(trlstart_newpand_dum,trl_start_tim(k))) trlstart_newpand_dum(ft_nearest(trlstart_newpand_dum,trl_start_tim(k)))],ylim,'Color','r');
% end

delay_var = 9; % ascertained this from visual inspection of eye data
newpandtime = newpandtime_dum + delay_var;

% plot(newpandtime(ft_nearest(newpandtime,data.cfg.previous{1}.trl(1,1)):ft_nearest(newpandtime,data.cfg.previous{1}.trl(end,2))), ...
%     xeye_real(ft_nearest(newpandtime,data.cfg.previous{1}.trl(1,1)):ft_nearest(newpandtime,data.cfg.previous{1}.trl(end,2)))*3.6,'g')

% delay_arr(2) = 0; % for JN130531.1; also ascertained from visual inspection

% need to resync every trial - base this on the plexon data trials, not the
% panda trials (don't always match)
trlinf = [];
for trllop = 2:length(data.trial) % start on trial 2 (first trial has sync errors)
    
    plex_trl = data.cfg.previous{1}.trl(trllop,1:2);
    
    panda_trial_sel = ft_nearest(trl_start_tim,plex_trl(1));
    
%     if delay_arr(panda_trial_sel)>-100 && delay_arr(panda_trial_sel)<100
        
        newpandtime = newpandtime_dum + delay_var - delay_arr(panda_trial_sel); % do this step for every trial loop
        newpandtime_pos = newpandtime_pos_dum + delay_var - delay_arr(panda_trial_sel);
        newpandtime_dir = newpandtime_dir_dum + delay_var - delay_arr(panda_trial_sel);
        newpandtime_vel = newpandtime_vel_dum + delay_var - delay_arr(panda_trial_sel);

        % buffer the start and end times for each trial 100 samples on each
        % side (50 for position and direction)
        trlstart_eye = ft_nearest(newpandtime,plex_trl(1))-100;
        trlend_eye = ft_nearest(newpandtime,plex_trl(2))+100;

        trlstart_pos = ft_nearest(newpandtime_pos,plex_trl(1))-50;
        trlend_pos = ft_nearest(newpandtime_pos,plex_trl(2))+50;

        trlstart_dir = ft_nearest(newpandtime_dir,plex_trl(1))-50;
        trlend_dir = ft_nearest(newpandtime_dir,plex_trl(2))+50;

        trlstart_vel = ft_nearest(newpandtime_vel,plex_trl(1))-50;
        trlend_vel = ft_nearest(newpandtime_vel,plex_trl(2))+50;

        trlinf{trllop}.eyetime = newpandtime(trlstart_eye:trlend_eye);
        trlinf{trllop}.postime = newpandtime_pos(trlstart_pos:trlend_pos);
        trlinf{trllop}.dirtime = newpandtime_dir(trlstart_dir:trlend_dir);
        trlinf{trllop}.veltime = newpandtime_vel(trlstart_vel:trlend_vel);
        trlinf{trllop}.eyedat = [xeye_real(trlstart_eye:trlend_eye)'; yeye_real(trlstart_eye:trlend_eye)'];
        trlinf{trllop}.posdat = [xpos_real(trlstart_pos:trlend_pos)'; ypos_real(trlstart_pos:trlend_pos)'];
        trlinf{trllop}.dirdat = dir_real(trlstart_dir:trlend_dir)';
        trlinf{trllop}.veldat = vel_real(trlstart_vel:trlend_vel)';
        
%     end
    
end
    
% for trllop=1:length(trlinf)
for trllop=[3 50 120 170]
    if ~isempty(trlinf{trllop})
        plex_trl = data.cfg.previous{1}.trl(trllop,1:2);

        figure;hold on
        plot(plex_trl(1):plex_trl(2),data.trial{trllop}(8,:))
        plot(trlinf{trllop}.eyetime,trlinf{trllop}.eyedat(1,:)*3.6,'r')
        plot(trlinf{trllop}.postime,trlinf{trllop}.posdat(1,:)*100,'g')
        plot(trlinf{trllop}.postime,trlinf{trllop}.posdat(2,:)*100,'m')
        plot(trlinf{trllop}.dirtime,trlinf{trllop}.dirdat*10-mean(trlinf{trllop}.dirdat*10),'k')
        plot(trlinf{trllop}.veltime,trlinf{trllop}.veldat*1000,'c')
    end
end


% % first_trial_ind = ft_nearest(time_eye,times(log_trl_start(trlind(1)))):ft_nearest(time_eye,times(log_trl_start(trlind(2))));
% first_trial_ind = ft_nearest(time_eye,time_eye(log_trl_start(1))):ft_nearest(time_eye,time_eye(log_trl_start(2)));
% 
% eyedat_first_trial = [xeye_real(first_trial_ind-1) yeye_real(first_trial_ind-1)];
% 
% % plot the first trial eye data from Bananarchy log and match it to the Plexon eye data
% figure;plot(newpandtime_eye(first_trial_ind-1), eyedat_first_trial)
% xlim(minmax(newpandtime_eye(first_trial_ind-1)'))
% title('Bananarchy eye data')

%%
% add the position and direction data into the fieldtrip data structure
% first make new data structure including only the trials with position
% information
cfg=[];
cfg.trialfun      = 'trialfun_nav';
cfg.dataset       = dataset;
cfg = ft_definetrial(cfg);

trldum = [];
goodtrial = [];
for trllop = 1:size(cfg.trl,1)
    if ~isempty(trlinf{trllop})
        trldum = [trldum; cfg.trl(trllop,:)];
        goodtrial=[goodtrial; trllop];
    end
end
cfg.trl = trldum;
trlinf = trlinf(goodtrial);

%code to find spike/LFP/eye channels
[lfpind,spkind,eyeind,namarr]=getchannels(dataset);

% for some reason extra LFP channels show up
lfpind = lfpind(1:4);

% preprocess the data, removing 60 Hz noise artifacts (and 120 Hz harmonic)
cfg.channel       = cellstr(namarr(lfpind,:));
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
cfg.padding       = 10;
cfg.continuous    = 'yes';
cfg.detrend       = 'no';
data_LFP = ft_preprocessing(cfg);

% cfg.channel       = cellstr(namarr(spkind,:));
% cfg.dftfilter     = 'no';
% data_spike = ft_preprocessing(cfg);

cfg.channel       = cellstr(namarr(eyeind,:));
cfg.dftfilter     = 'no';
data_eye = ft_preprocessing(cfg);

% data = ft_appenddata([], data_spike, data_LFP, data_eye);
data = ft_appenddata([], data_LFP, data_eye);
clear data_*

data.posdat=[];
data.dirdat=[];
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(data.trial)

    ft_progress(trllop/length(data.trial), 'Processing event %d from %d', trllop, length(data.trial));

    % interpolate position data
    postime_interp = trlinf{trllop}.postime(1):trlinf{trllop}.postime(end);
    pos_interp = nan(2,length(postime_interp));
    for k=1:length(postime_interp)
        if ismember(postime_interp(k),trlinf{trllop}.postime)
            pos_interp(:,k) = trlinf{trllop}.posdat(:,find(trlinf{trllop}.postime==postime_interp(k),1,'last'));
        end
    end
    pos_interp=inpaint_nans(pos_interp,2);

    % interpolate direction data
    dirtime_interp = trlinf{trllop}.dirtime(1):trlinf{trllop}.dirtime(end);
    dir_interp = nan(1,length(dirtime_interp));
    for k=1:length(dirtime_interp)
        if ismember(dirtime_interp(k),trlinf{trllop}.dirtime)
            dir_interp(1,k) = trlinf{trllop}.dirdat(1,find(trlinf{trllop}.dirtime==dirtime_interp(k),1,'last'));
        end
    end
    dir_interp=inpaint_nans(dir_interp,2);

    % interpolate velocity data
    veltime_interp = trlinf{trllop}.veltime(1):trlinf{trllop}.veltime(end);
    vel_interp = nan(1,length(veltime_interp));
    for k=1:length(veltime_interp)
        if ismember(veltime_interp(k),trlinf{trllop}.veltime)
            vel_interp(1,k) = trlinf{trllop}.veldat(1,find(trlinf{trllop}.veltime==veltime_interp(k),1,'last'));
        end
    end
    vel_interp=inpaint_nans(vel_interp,2);

    
    plex_trl = data.cfg.previous{1}.trl(trllop,1:2);
    
    data.posdat{trllop} = pos_interp(:,find(postime_interp==plex_trl(1)):find(postime_interp==plex_trl(2)));
    data.dirdat{trllop} = dir_interp(find(dirtime_interp==plex_trl(1)):find(dirtime_interp==plex_trl(2)));
    data.veldat{trllop} = vel_interp(find(veltime_interp==plex_trl(1)):find(veltime_interp==plex_trl(2)));
    
end
ft_progress('close')

save('C:\Users\michael.jutras\Documents\MATLAB\Virtual\FT data\JN1306041_navdat140425.mat','data')

%%

spkind = find(strncmp('sig',data.label,3));
spkloc = cell(1,length(spkind));
spkdir = cell(1,length(spkind));
avgspkratepertrl = nan(length(spkind),length(data.trial));
for spklop = 1:length(spkind)
    for trllop = 1:length(data.trial)
        spkloc{spklop} = [spkloc{spklop} data.posdat{trllop}(:,logical(data.trial{trllop}(spkind(spklop),:)))];
        spkdir{spklop} = [spkdir{spklop} data.dirdat{trllop}(logical(data.trial{trllop}(spkind(spklop),:)))];
        spkvel{spklop} = [spkvel{spklop} data.veldat{trllop}(logical(data.trial{trllop}(spkind(spklop),:)))];
        
        avgspkratepertrl(spklop,trllop) = (length(find(data.trial{trllop}(spkind(spklop),:)))/length(data.time{trllop}))*1000;
    end
end

% look at average spike rates
figure
for k=1:size(avgspkratepertrl,1)
    subplot(size(avgspkratepertrl,1),1,k)
    bar(1:size(avgspkratepertrl,2),avgspkratepertrl(k,:))
    title(data.label{k})
end

lastgoodtrial = [18 10 14 11 18];

R = [];
for spklop = 1:length(spkloc)
    
    figure; hold on
    subplot(1,2,1)
    
    S=[]; X=[]; Y=[];
    for trllop = 1:lastgoodtrial(spklop)
        S = [S data.trial{trllop}(spkind(spklop),:)];
        X = [X data.posdat{trllop}(1,:)];
        Y = [Y data.posdat{trllop}(2,:)];
    end
    
    scatter(spkloc{spklop}(1,:),spkloc{spklop}(2,:),'.r')
    line(X,Y,'Color',[0.5 0.5 0.5])

    o.gsize = [1 1];

    R{spklop} = rmapMJ(S,X,Y);
    
    
    subplot(1,2,2)
    imagesc(R{spklop}.a)
    set(gcf,'Position',[299   557   941   440])
    toptitle(data.label{spklop})
    
end


%% add the position information to the data file
% include eye information from the log file to make sure the information
% lines up across files


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

