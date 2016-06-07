if strcmp(license,'375399')
    NSDir = 'C:\Data\VR';
    trlDir = 'C:\Data\VR';
elseif strcmp(license,'613743')
    NSDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NSdat';
    trlDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trldat';
end

% % this will automatically pull files in a certain range; skip for now
% filelist=cell(1);
% c=1;
% d=dir(savDir);
% for k = 1:length(d)
%     if length(d(k).name)>=7 && strncmpi(d(k).name(end-10:end-4),'_trldat',7)
%         if strncmpi(d(k).name(4:5),'BR',2)
%             if str2num(d(k).name([7:8 10:11 13:14]))>=150410
%                 filelist{c} = d(k).name;
%                 c=c+1;
%             end
%         end
%     end
% end

% specify the sessions to use for now
seslst = {'JN_BR_15_04_10\JN_BR_15_04_10_13_03';
    'JN_BR_15_04_10\JN_BR_15_04_10_13_11';
    'JN_BR_15_04_13\JN_BR_15_04_13_13_59';
    'JN_BR_15_04_14\JN_BR_15_04_14_13_34';
%     'JN_BR_15_04_15\JN_BR_15_04_15_12_32';
    'JN_BR_15_04_16\JN_BR_15_04_16_14_22';
    'JN_BR_15_04_17\JN_BR_15_04_17_13_14';
    'JN_BR_15_04_20\JN_BR_15_04_20_14_43';
%     'JN_BR_15_04_21\JN_BR_15_04_21_12_34';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_09';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_54';
%     'JN_BR_15_04_23\JN_BR_15_04_23_13_54'; % A bad these 3 recordings
%     'JN_BR_15_04_23\JN_BR_15_04_23_14_22';
%     'JN_BR_15_04_24\JN_BR_15_04_24_14_04';
    };


% initialize dum structure
dum = [];
dum.time = cell(1);
dum.trial = cell(1);
dum.sampleinfo = [];
c=1;

alpall = []; % alpha
banfndall = []; % 1 if got banana, 0 if not
banposall = []; % banana position (first banana, includes invisible)
bantimindall = []; % latency to get banana (technically, time index)
avgbandispertrl = []; % average distance from banana
cumdispertrl = []; % cumulative distance from banana
avgwaldispertrl = []; % average distance from wall
trllngall = []; % total trial length
invprmindall = []; % invisible, prime (1st in a series of inv. bananas)
filind = []; % file index, per trial
numtrl = []; % number of trials per file
d=0;
for seslop = 1:length(seslst)
    
    [~,sesnam]=fileparts(seslst{seslop});
    disp(['Processing ' sesnam])

    load(fullfile(NSDir,[sesnam '_NSdat.mat']))
    load(fullfile(trlDir,[sesnam '_trldat.mat']))

%     sesttl=sesnam;sesttl(sesttl=='_')='-';
%     figure;plot(1:size(data.sampleinfo,1),data.sampleinfo);title(sesttl)
    
    numtrl(seslop,1) = length(trldat.time);
    
    % calculate three metrics for each trial:
    % latency to get banana (equal to time-out if banana not acquired)
    % average distance from banana
    % cumulative distance from banana
    
    invprmind = [];
    for trllop = 1:length(trldat.posdat)
        
        % find time when banana eaten, or end of trial if banana not eaten
        if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
            [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
            banfnd = 1;
        else
            bantimind = length(trldat.time{trllop});
            banfnd = 0;
        end
        
        banpos = trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),3:4);
        
        % calculate distance from banana
        xdif = trldat.posdat{trllop}(1,1:bantimind)-banpos(1);
        ydif = trldat.posdat{trllop}(2,1:bantimind)-banpos(2);
        dis = sqrt(xdif.^2+ydif.^2);
        
        % figure out distance from wall; use 11.2 as wall boundary
        waldis = nan(size(trldat.posdat{trllop},2),1);
        for timlop = 1:size(trldat.posdat{trllop},2)
            waldis(timlop) = min(abs([-11.2-trldat.posdat{trllop}(:,timlop); 11.2-trldat.posdat{trllop}(:,timlop)]));
        end
        
        alpall = [alpall; trldat.alpha{trllop}(1,2)];
        banfndall = [banfndall; banfnd];
        banposall = [banposall; banpos];
        bantimindall = [bantimindall; bantimind];
        avgbandispertrl = [avgbandispertrl; mean(dis)];
        cumdispertrl = [cumdispertrl; sum(dis)];
        avgwaldispertrl = [avgwaldispertrl; mean(waldis)];
        trllngall = [trllngall; length(trldat.time{trllop})];
        filind = [filind; seslop];
        
        % find trial with invisible banana
        if trldat.alpha{trllop}(1,2)==0
            % don't include if the banana was invisible on previous trial 
            if trllop>=2 && trldat.alpha{trllop-1}(1,2)~=0
                % don't include if previous invisible banana was in the same location
                if ~isempty(invprmind)
                    if ~prod(trldat.frtpos{invprmind(end)}(find(trldat.frtpos{invprmind(end)}(:,2)==0,1,'first'),3:4)==banpos,2)
                        invprmind = [invprmind; trllop];
                    end
                else
                    invprmind = [invprmind; trllop];
                end
            end
        end
        
    end
    
    dum.time(c:length(invprmind)+c-1) = data.time(invprmind-1);
    dum.trial(c:length(invprmind)+c-1) = data.trial(invprmind-1);
    dum.sampleinfo(c:length(invprmind)+c-1,:) = data.sampleinfo(invprmind-1,:);
    c = c+length(invprmind);

    invprmindall = [invprmindall; invprmind+d];
    d = d+length(trldat.time);
    
end

dum.fsample = data.fsample;
dum.label = data.label;


%% preprocess data

% exclude these channels
badchn = {'A02'; 'A05'; 'A08'; 'A10'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};
% badchn = {'A01'; 'A02'; 'A03'; 'A04'; 'A05'; 'A06'; 'A07'; 'A08'; 'A09'; 'A10'; 'A11'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};

% exclude eye and position data for now
badchn = [badchn; {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'}];

cfg=[];
cfg.continuous    = 'no';
cfg.channel       = setxor(data.label,badchn);
data0 = ft_preprocessing(cfg,dum);

% these data are still noisy; wait until 5-second clips are selected before
% filtering out noise


%% "banana eat" events for all invisible banana trials included in data0

alp = alpall(invprmindall-1);
bantim = bantimindall(invprmindall-1);


%% select last 5 seconds of approach path to banana in neural data

timsel = [bantim-4999 bantim+500]; % add extra 500 ms at the end to check for reward artifact
timsel(timsel(:,1)<0,:) = nan; % keep track of which trials are nans, to exclude later during behavioral analysis

data1 = data0;
for trllop = 1:length(data0.trial)
    if isempty(find(isnan(timsel(trllop,:)),1,'first'))
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1];
        data1.time{trllop} = -4.999:0.001:0.5;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2));
    end
end

% remove marked trials
cfg=[];
cfg.trials        = setxor(1:length(data1.trial),find(isnan(timsel(:,1)))); % here choose only trials with at least 5 seconds pre-banana 
data1 = ft_preprocessing(cfg,data1);


% subtract common average of each probe from all channels on that probe
data2 = data1;
clear prblaball
for k=1:length(data1.label)
    prblaball(k,1)=data1.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);
for trllop = 1:length(data1.trial)
    for prblop = 1:length(prblab)
        prbavg = mean(data1.trial{trllop}(locb==prblop,:),1);
        repavg = repmat(prbavg,length(find(locb==prblop)),1);
        mat1 = reshape(data1.trial{trllop}(locb==prblop,:),1,length(find(locb==prblop)),size(data1.trial{trllop},2));
        mat2 = reshape(repavg,1,size(repavg,1),size(repavg,2));
        data2.trial{trllop}(locb==prblop,:)=diff(cat(1,mat1,mat2),1,1);
    end
end


%% mark bad trials (CHANGE IF ADDING NEW DATA)
% identified with ft_rejectvisual, ft_databrowser

badtrl = [13 25 36 76 77];


%% do the time-locked analysis (check for reward artifact)

clear cfg
cfg.keeptrials    = 'yes';
cfg.vartrllength  =  2;
cfg.trials        = setxor(1:length(data2.time),badtrl);
timelock1 = ft_timelockanalysis(cfg,data1);
timelock2 = ft_timelockanalysis(cfg,data2);
timelock3 = ft_timelockanalysis(cfg,data3);
% timelock_dum = ft_timelockanalysis(cfg,dum); % if running ft_rejectvisual

figure;plot(timelock1.time,timelock1.avg)
legend(timelock.label);legend off


%% spectral analysis, all trials

% chop off extra time after reward
cfg=[];
cfg.trials = setxor(1:length(data2.time),badtrl);
data3 = ft_preprocessing(cfg,data2);
for trllop = 1:length(data3.time)
    ind = data3.time{trllop}<=0;
    data3.time{trllop} = data3.time{trllop}(ind);
    data3.trial{trllop} = data3.trial{trllop}(:,ind);
    data3.sampleinfo(trllop,2) = data3.sampleinfo(trllop,1)+find(ind==1,1,'last')-1;
end

% filter out noise here instead
cfg=[];
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
data3 = ft_preprocessing(cfg,data3);


% do the spectral analysis - time-averaged
cfg=[];
cfg.output      = 'fourier';
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.foilim      = [1 30];
% cfg.taper       = 'hanning';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqL = ft_freqanalysis(cfg, data3);
fdL = ft_freqdescriptives(cfg, freqL);

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqH = ft_freqanalysis(cfg, data3);
fdH = ft_freqdescriptives(cfg, freqH);


cfg=[];
cfg.method      = 'coh';
statL = ft_connectivityanalysis(cfg,freqL);
statH = ft_connectivityanalysis(cfg,freqH);


%% sample spectral analysis code for looking at noise

% do the spectral analysis - time-averaged
cfg=[];
cfg.output      = 'pow';
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.taper       = 'hanning';
cfg.foilim      = [40 130];
% cfg.channel     = data.label(find(strcmp(data.label,'C01')):find(strcmp(data.label,'C12')));
% cfg.trials      = 1:5;
% freqT1 = ft_freqanalysis(cfg, data);
% freqT2 = ft_freqanalysis(cfg, data0);
freqT3 = ft_freqanalysis(cfg, data1);
freqT4 = ft_freqanalysis(cfg, data2);
freqT5 = ft_freqanalysis(cfg, data3);
figure;hold on
% plot(freqT1.freq,squeeze(mean(freqT1.powspctrm,1)),'b')
% plot(freqT2.freq,squeeze(mean(freqT2.powspctrm,1)),'r')
plot(freqT3.freq,squeeze(mean(freqT3.powspctrm,1)),'g')
plot(freqT4.freq,squeeze(mean(freqT4.powspctrm,1)),'m')
plot(freqT5.freq,squeeze(mean(freqT5.powspctrm,1)),'c')


%% examine different behavioral measures, relationship with memory

% plot distribution of recordings represented in the data
figure;hist(filind(invprmindall),30)

figure;hold on
scatter(bantimindall(invprmindall(logical(banfndall(invprmindall)))), ...
    avgbandispertrl(invprmindall(logical(banfndall(invprmindall)))),'r')
scatter(bantimindall(invprmindall(logical(~banfndall(invprmindall)))), ...
    avgbandispertrl(invprmindall(logical(~banfndall(invprmindall)))),'b')
line([median(bantimindall(invprmindall)) median(bantimindall(invprmindall))],ylim,'Color','k')
line(xlim,[median(avgbandispertrl(invprmindall)) median(avgbandispertrl(invprmindall))],'Color','k')
xlabel('Time to acquire banana')
ylabel('Mean distance from banana')
legend('Gets banana','Doesn''t get banana')

figure;hold on
scatter(bantimindall(invprmindall(logical(banfndall(invprmindall)))), ...
    avgwaldispertrl(invprmindall(logical(banfndall(invprmindall)))),'r')
scatter(bantimindall(invprmindall(logical(~banfndall(invprmindall)))), ...
    avgwaldispertrl(invprmindall(logical(~banfndall(invprmindall)))),'b')
line([median(bantimindall(invprmindall)) median(bantimindall(invprmindall))],ylim,'Color','k')
line(xlim,[median(avgwaldispertrl(invprmindall)) median(avgwaldispertrl(invprmindall))],'Color','k')
xlabel('Time to acquire banana')
ylabel('Mean distance from wall')
legend('Gets banana','Doesn''t get banana')

figure;hold on
scatter(avgbandispertrl(invprmindall(logical(banfndall(invprmindall)))), ...
    avgwaldispertrl(invprmindall(logical(banfndall(invprmindall)))),'r')
scatter(avgbandispertrl(invprmindall(logical(~banfndall(invprmindall)))), ...
    avgwaldispertrl(invprmindall(logical(~banfndall(invprmindall)))),'b')
line([median(avgbandispertrl(invprmindall)) median(avgbandispertrl(invprmindall))],ylim,'Color','k')
line(xlim,[median(avgwaldispertrl(invprmindall)) median(avgwaldispertrl(invprmindall))],'Color','k')
xlabel('Mean distance from banana')
ylabel('Mean distance from wall')
legend('Gets banana','Doesn''t get banana')

figure
scatter3(bantimindall(invprmindall), avgbandispertrl(invprmindall), avgwaldispertrl(invprmindall),'.')
xlabel('Time to acquire banana')
ylabel('Mean distance from banana')
zlabel('Mean distance from wall')


%% categorize trials into good/bad memory
% use distance from corners of "Mean distance from banana/wall" plot
% normalize these variables first

bantimnrm = bantimindall(invprmindall)/max(bantimindall(invprmindall));
bandstnrm = avgbandispertrl(invprmindall)/max(avgbandispertrl(invprmindall));
waldstnrm = avgwaldispertrl(invprmindall)/max(avgwaldispertrl(invprmindall));


% time to acquire banana, mean dist. from banana, mean dist. from wall
mempnthi = [0 0 max(waldstnrm)];
mempntlo = [max(bantimnrm) max(bandstnrm) 0];
memdst = sqrt(sum((abs(([bantimnrm bandstnrm waldstnrm])-repmat(mempnthi,length(invprmindall),1))).^2,2));
pltttl = 'time to acquire banana, mean dist. from banana, mean dist. from wall';

% % time to acquire banana, mean dist. from banana
% mempnthi = [0 0];
% mempntlo = [max(bantimnrm) max(bandstnrm)];
% membandst = sqrt(sum((abs(([bantimnrm bandstnrm])-repmat(mempnthi,length(invprmindall),1))).^2,2));
% pltttl = 'time to acquire banana, mean dist. from banana';

% % time to acquire banana, mean dist. from wall
% mempnthi = [0 max(waldstnrm)];
% mempntlo = [max(bantimnrm) 0];
% membandst = sqrt(sum((abs(([bantimnrm waldstnrm])-repmat(mempnthi,length(invprmindall),1))).^2,2));
% pltttl = 'time to acquire banana, mean dist. from wall';

[~,memsrt] = sort(memdst);

memindhi = memsrt(1:round(length(memsrt)/2));
memindlo = memsrt(round(length(memsrt)/2)+1:end);

bantim_himem = bantimindall(invprmindall(memindhi));
bantim_lomem = bantimindall(invprmindall(memindlo));

% figure;hold on
% scatter(avgbandispertrl(invprmindall(memindhi)), avgwaldispertrl(invprmindall(memindhi)),'r')
% scatter(avgbandispertrl(invprmindall(memindlo)), avgwaldispertrl(invprmindall(memindlo)),'b')
% xlim([0 max(avgbandispertrl(invprmindall))])
% ylim([0 max(avgwaldispertrl(invprmindall))])
% line([median(avgbandispertrl(invprmindall)) median(avgbandispertrl(invprmindall))],ylim,'Color','k')
% line(xlim,[median(avgwaldispertrl(invprmindall)) median(avgwaldispertrl(invprmindall))],'Color','k')
% xlabel('Mean distance from banana')
% ylabel('Mean distance from wall')
% title(pltttl)

figure;hold on
scatter3([bantimindall(invprmindall(memindhi)); bantimindall(invprmindall(memindlo))], ...
    [avgbandispertrl(invprmindall(memindhi)); avgbandispertrl(invprmindall(memindlo))], ...
    [avgwaldispertrl(invprmindall(memindhi)); avgwaldispertrl(invprmindall(memindlo))], ...
    ones(length(invprmindall),1)*10,[repmat([1 0 0],length(memindhi),1); repmat([0 0 1],length(memindlo),1)])
xlabel('Time to acquire banana')
ylabel('Mean distance from banana')
zlabel('Mean distance from wall')
title(pltttl)


%% using memindhi and memindlo, categorize neural data trials

trlind = (1:length(timsel))'; % same size as membandst right now
trlind = trlind(~isnan(timsel(:,1))); % take out trials <5 sec pre-banana
trlind = trlind(setxor(1:length(trlind),badtrl)); % take out bad trials

memdstsel = memdst(trlind);

[~,memsrt] = sort(memdstsel);
memindhi = memsrt(1:30);
memindlo = memsrt(end-29:end);

cfg=[];
cfg.method      = 'coh';

cfg.trials = memindhi;
statL_himem = ft_connectivityanalysis(cfg,freqL);
statH_himem = ft_connectivityanalysis(cfg,freqH);
fdL_himem = ft_freqdescriptives(cfg, freqL);
fdH_himem = ft_freqdescriptives(cfg, freqH);

cfg.trials = memindlo;
statL_lomem = ft_connectivityanalysis(cfg,freqL);
statH_lomem = ft_connectivityanalysis(cfg,freqH);
fdL_lomem = ft_freqdescriptives(cfg, freqL);
fdH_lomem = ft_freqdescriptives(cfg, freqH);


%% plot

ab1=4;ab2=18;
figure;hold on;
plot(statL_himem.freq,squeeze(statL_himem.cohspctrm(ab1,ab2,:)),'g');
plot(statL_lomem.freq,squeeze(statL_lomem.cohspctrm(ab1,ab2,:)),'r');
title([statL_himem.label{ab1} ' x' statL_himem.label{ab2}]);
set(gca,'TickDir','out')
legend({'High Memory';'Low Memory'})
xlabel('Frequency (Hz)');ylabel('Coherence')
