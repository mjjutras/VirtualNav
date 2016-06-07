%% QC: CHECK LINEUP OF NEURAL DATA WITH PYTHON DATA!!!
% CHECK USING NEW BLACKROCK MATLAB FILES!!!

%%
% choose only trials where he gets the banana within 60 seconds
% check to see that there is good representation across the field
% normalize excess path length and latency by the "excess path length"
% during visible trials to account for the difference in the bubble across
% field locations

%%

if strcmp(getenv('computername'),'MIKE-PC')
    NSDir = 'C:\Data\VR\NSdat';
    trlDir = 'C:\Data\VR\trldat';
elseif strcmp(getenv('computername'),'RBU-MIKEJ')
    NSDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NSdat';
    trlDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trldat';
end

%%

% if strcmp(getenv('computername'),'MIKE-PC')
%     trlDir = 'C:\Data\VR\trldat';
% elseif strcmp(getenv('computername'),'RBU-MIKEJ')
%     trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';
% end

% specify the sessions to use for now
seslst = {'JN_BR_15_04_10\JN_BR_15_04_10_13_03';
    'JN_BR_15_04_10\JN_BR_15_04_10_13_11';
    'JN_BR_15_04_13\JN_BR_15_04_13_13_59';
    'JN_BR_15_04_14\JN_BR_15_04_14_13_34';
    'JN_BR_15_04_16\JN_BR_15_04_16_14_22';
    'JN_BR_15_04_17\JN_BR_15_04_17_13_14';
    'JN_BR_15_04_20\JN_BR_15_04_20_14_43';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_09';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_54';
    'JN_BR_15_04_28\JN_BR_15_04_28_13_07';
    'JN_BR_15_04_29\JN_BR_15_04_29_13_22';
    'JN_BR_15_04_30\JN_BR_15_04_30_14_42';
%     'JN_BR_15_05_01\JN_BR_15_05_01_13_12';
%     'JN_BR_15_05_04\JN_BR_15_05_04_13_57';
%     'JN_BR_15_05_04\JN_BR_15_05_04_14_51';
%     'JN_BR_15_05_05\JN_BR_15_05_05_13_02';
%     'JN_BR_15_05_06\JN_BR_15_05_06_13_30';
%     'JN_BR_15_05_07\JN_BR_15_05_07_12_59';
%     'JN_BR_15_05_08\JN_BR_15_05_08_14_17';
%     'JN_BR_15_05_12\JN_BR_15_05_12_12_49';
%     'JN_BR_15_05_13\JN_BR_15_05_13_13_14';
%     'JN_BR_15_05_14\JN_BR_15_05_14_13_34';
%     'JN_BR_15_05_18\JN_BR_15_05_18_13_03';
%     'JN_BR_15_05_19\JN_BR_15_05_19_13_00'; % 150624: added 15_05_01 through 15_05_19
%     'JN_BR_15_05_20\JN_BR_15_05_20_13_04';
%     'JN_BR_15_05_22\JN_BR_15_05_22_13_20';
%     'JN_BR_15_05_26\JN_BR_15_05_26_13_42';
%     'JN_BR_15_05_27\JN_BR_15_05_27_12_53';
%     'JN_BR_15_05_28\JN_BR_15_05_28_13_59';
%     'JN_BR_15_05_28\JN_BR_15_05_28_14_16';
%     'JN_BR_15_05_29\JN_BR_15_05_29_14_19';
%     'JN_BR_15_06_01\JN_BR_15_06_01_14_35';
%     'JN_BR_15_06_02\JN_BR_15_06_02_13_15';
%     'JN_BR_15_06_03\JN_BR_15_06_03_14_32';
%     'JN_BR_15_06_04\JN_BR_15_06_04_13_04';
%     'JN_BR_15_06_05\JN_BR_15_06_05_13_35';
    };

% initialize dum structure
dum = [];
dum.time = cell(1);
dum.trial = cell(1);
dum.sampleinfo = [];
b=1;

begposall = []; % begin position (trial start)
alpall = []; % alpha
banfndall = []; % 1 if got banana, 0 if not
banposall = []; % banana position (first banana, includes invisible)
bantimindall = []; % latency to get banana (technically, time index)
chrtimindall = []; % latency (time index) to get cherry/cherries
avgbandispertrl = []; % average distance from banana
cumdispertrl = []; % cumulative distance from banana
avgwaldispertrl = []; % average distance from wall
trllngall = []; % total trial length
invprmindall = []; % invisible, prime (1st in a series of inv. bananas)
nrlsel = []; % trials used to select neural activity (dum)
filind = []; % file index, per trial
numtrl = []; % number of trials per file
excpthdif = []; % excess path length, using difference b/t path length and shortest path
excpthnrm = []; % excess path length, using path length divided by shortest path
excpthind = []; % excess path length, index
d=0;
for seslop = 1:length(seslst)
    
    [~,sesnam]=fileparts(seslst{seslop});
    disp(['Processing ' sesnam])

    load(fullfile(NSDir,[sesnam '_NSdat.mat']))
    load(fullfile(trlDir,[sesnam '_trldat.mat']))

    numtrl(seslop,1) = length(trldat.time);
    
    % calculate three metrics for each trial:
    % latency to get banana (equal to time-out if banana not acquired)
    % average distance from banana
    % cumulative distance from banana
    
    invprmind = []; % invisible "prime" (first invisible banana in a sequence) index
    prvbanalp = cell(1); % previous banana alpha
    excpthdifdum = zeros(length(trldat.posdat),1);
    excpthnrmdum = zeros(length(trldat.posdat),1);
    excpthinddum = zeros(length(trldat.posdat),1);
    bantimindses = zeros(length(trldat.posdat),1);
    for trllop = 1:length(trldat.posdat)
        
        % find time when banana eaten, or end of trial if banana not eaten
        if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
            [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
            banfnd = 1;
        else
            % NOTE: if he doesn't find the banana, the actual time here
            % should be the "TIMEOUT", which is either 90 seconds or some
            % other longer amount (before 5/8/15)
            bantimind = length(trldat.time{trllop});
            banfnd = 0;
        end
        
        % get time when cherry eaten (if cherry appears)
        if ~isempty(find(trldat.frttim{trllop}(:,2)==1,1,'first'))
            chrtimind = [];
            chrind = find(ismember(trldat.frttim{trllop}(:,2),[1 2]));
            for frtlop = 1:length(chrind)
                [~,chrtimind(frtlop)] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(chrind(frtlop),1)));
            end
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
        
        begposall = [begposall; trldat.posdat{trllop}(:,1)'];
        alpall = [alpall; trldat.alpha{trllop}(1,2)];
        banfndall = [banfndall; banfnd];
        banposall = [banposall; banpos];
        bantimindall = [bantimindall; bantimind];
        chrtimindall = [chrtimindall; {chrtimind}];
        bantimindses(trllop) = bantimind;
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
                else % if this if the first invisible banana of the session
                    invprmind = [invprmind; trllop];
                end
                
                % determine alpha of each visible bananas presented in the
                % same location in the preceding trials
                if invprmind(end)==trllop
                    prvbanalp{length(invprmind),1} = [];
                    for revlop=trllop-1:-1:1
                        revpos = trldat.frtpos{revlop}(find(trldat.frtpos{revlop}(:,2)==0,1,'first'),3:4);
                        if revpos(1)~=banpos(1) || revpos(2)~=banpos(2)
                            break
                        else
                            prvbanalp{length(invprmind),1} = [prvbanalp{length(invprmind),1}; revlop+d trldat.alpha{revlop}(1,2)];
                        end
                    end
                    prvbanalp{length(invprmind),1} = flipud(prvbanalp{length(invprmind),1});                        
                end

%                 % Yoni's method: need to calculate velocity directly from
%                 % position (doesn't work for veldat from Panda log)
%                 excpthdifdum(trllop,1) = sum(conv(trldat.veldat{trllop},normpdf(-1000:1000,0,200),'same'))-dis(1);
%                 excpthnrmdum(trllop,1) = sum(conv(trldat.veldat{trllop},normpdf(-1000:1000,0,200),'same'))/dis(1);
                % this way works
                excpthdifdum(trllop,1) = sum(abs(diff(dis)))-dis(1);
                excpthnrmdum(trllop,1) = sum(abs(diff(dis)))/dis(1);
                excpthinddum(trllop,1) = (sum(abs(diff(dis)))-dis(1))/(sum(abs(diff(dis)))+dis(1));
                
            end
        end
        
    end

    % select trials for neural analysis (nrlsel)
    for invlop = 1:length(invprmind)
        % choose only "invisible banana trials" with at least 2 preceding visible bananas 
        if length(prvbanalp{invlop})>=2
            % select only trials with banana search time of 60 seconds or less
            if bantimindses(invprmind(invlop))<=60000
                % select only trials where there were at least 5 seconds of
                % search time for both of the preceding visible bananas
                if bantimindses(invprmind(invlop)-2)>=5000 && bantimindses(invprmind(invlop)-1)>=5000
                    nrlsel = [nrlsel; invprmind(invlop)+d];
                    dum.time(b:b+1) = data.time(invprmind(invlop)-2:invprmind(invlop)-1);
                    dum.trial(b:b+1) = data.trial(invprmind(invlop)-2:invprmind(invlop)-1);
                    dum.sampleinfo(b:b+1,:) = data.sampleinfo(invprmind(invlop)-2:invprmind(invlop)-1,:);
                    b = b+2;
                end
            end
        end
    end
    
    % add invprmind trials to invprmindall
    invprmindall = [invprmindall; invprmind+d];
    d = d+length(trldat.time);
    
    excpthdif = [excpthdif; excpthdifdum];
    excpthnrm = [excpthnrm; excpthnrmdum];
    excpthind = [excpthind; excpthinddum];
    
end

dum.fsample = data.fsample;
dum.label = data.label;

clear data trldat


%% save variables to network 
if strcmp(getenv('computername'),'MIKE-PC')
    WSdir = 'C:\Data\VR';
elseif strcmp(getenv('computername'),'RBU-MIKEJ')
    WSdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\workspace';
end
% save(fullfile(WSdir,'anaFruit_neuraldum_trainingses_cherry_151014.mat'),'-v7.3')

% load variables
load(fullfile(WSdir,'anaFruit_neuraldum_trainingses_cherry_151014.mat'))


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

% these data are still noisy; wait until 5-second clips are selected before
% filtering out noise

% noisy trials:
% 11, 12, 71, 72, 197, 338, 361, 413, 418, 420, 431, 432, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 549, 553, 600


%% find trials with at least 5 seconds preceding banana
% this should be every trial that was included in nrlsel

bantim = nan(length(nrlsel)*2,1);
chrtim = [];
timsel = nan(length(nrlsel)*2,2);
nrltimsel = [];
trlgrp = [];
timselchr = [];
filindnrlsel = [];
for invlop = 1:length(nrlsel)
    bantim(invlop*2-1:invlop*2) = bantimindall(nrlsel(invlop)-2:nrlsel(invlop)-1);
    chrtim = [chrtim; chrtimindall(nrlsel(invlop)-2); chrtimindall(nrlsel(invlop)-1)];
    filindnrlsel = [filindnrlsel; filind(nrlsel(invlop)-2); filind(nrlsel(invlop)-1)];
    if isempty(find((bantim(invlop*2-1:invlop*2))<5000,1))
        timsel(invlop*2-1:invlop*2,:) = [bantim(invlop*2-1:invlop*2)-4999 bantim(invlop*2-1:invlop*2)];
        nrltimsel = [nrltimsel; nrlsel(invlop)];
        trlgrp = [trlgrp; ones(2,1)+length(trlgrp)/2];
    else
        timsel(invlop*2-1:invlop*2,:) = repmat([nan nan],2,1);
    end
    for frtlop = 1:length(chrtimindall{nrlsel(invlop)-2})
        timselchr{invlop*2-1,1}(frtlop,:) = [chrtimindall{nrlsel(invlop)-2}(frtlop)-4999 chrtimindall{nrlsel(invlop)-2}(frtlop)];
    end
    for frtlop = 1:length(chrtimindall{nrlsel(invlop)-1})
        timselchr{invlop*2,1}(frtlop,:) = [chrtimindall{nrlsel(invlop)-1}(frtlop)-4999 chrtimindall{nrlsel(invlop)-1}(frtlop)];
    end
end


%% select last 5 seconds of approach path to banana in neural data

data1 = data0;
for trllop = 1:length(data0.trial)
    if isempty(find(isnan(timsel(trllop,:)),1,'first'))
%         data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1];
%         data1.time{trllop} = -4.999:0.001:0;
%         data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2));

        % use this code when doing time-frequency analysis (padding at the end)
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1+500];
        data1.time{trllop} = -4.999:0.001:0.5;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2)+500);
    end
end

% remove marked trials
cfg=[];
% here choose only trials with at least 5 seconds pre-banana 
cfg.trials        = setxor(1:length(data1.trial),find(isnan(timsel(:,1))));
data1 = ft_preprocessing(cfg,data1);
% clear data0

% subtract common average of each probe from all channels on that probe
data2 = data1;
% comment-out the following if choosing not to re-reference to common probe average
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
clear data1

% now run ft_rejectvisual and ft_databrowser to clean bad trials


%% mark bad trials (CHANGE IF ADDING NEW DATA)
% identified with ft_rejectvisual, ft_databrowser

badtrldum = [16, 25, 28, 30, 37, 94, 106, 113, 148, 215]; % if re-referencing to common probe average
% badtrldum = [16, 25, 30, 58, 92]; % if not re-referencing

badtrl = [];
invselexc = [];
for badlop = 1:length(badtrldum)
    badtrl = [badtrl; find(trlgrp==trlgrp(badtrldum(badlop)))];
    invselexc = [invselexc; find(trlgrp==trlgrp(badtrldum(badlop)),1,'last')/2];
end

nrlbadsel = nrltimsel(setxor(1:length(nrltimsel),invselexc));
filindbadsel = filindnrlsel(setxor(1:length(filindnrlsel),badtrl));

% exclude bad trials, filter out line noise here
cfg=[];
cfg.trials = setxor(1:length(data2.time),badtrl);
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
data3 = ft_preprocessing(cfg,data2);
clear data2


%% spectral analysis, all trials

% % do the spectral analysis - time-averaged
% cfg=[];
% % cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
% cfg.output      = 'powandcsd';
% cfg.method      = 'mtmfft';
% cfg.pad         = 'maxperlen';
% cfg.keeptrials  = 'yes';
% cfg.foilim      = [1 30];
% % cfg.taper       = 'hanning';
% cfg.taper       = 'dpss';
% cfg.tapsmofrq   = 2;
% freqL = ft_freqanalysis(cfg, data3);
% 
% cfg.foilim      = [30 100];
% cfg.taper       = 'dpss';
% cfg.tapsmofrq   = 8;
% freqH = ft_freqanalysis(cfg, data3);
% 
% 
% cfg=[];
% cfg.method      = 'coh';
% % cfg.jackknife   = 'yes';
% statL = ft_connectivityanalysis(cfg,freqL);
% statH = ft_connectivityanalysis(cfg,freqH);


%% get pre-cherry data

datac1 = data0;
datac1 = rmfield(datac1,{'sampleinfo' 'trial' 'time'});

filindchrdum = [];

c=1;
for trllop = 1:length(timselchr)
    
    for frtlop = 1:size(timselchr{trllop},1)
        
        datac1.sampleinfo(c,:) = [data0.sampleinfo(trllop,1)+timselchr{trllop}(frtlop,1)-1 data0.sampleinfo(trllop,1)+timselchr{trllop}(frtlop,2)-1];
        datac1.time{c} = -4.999:0.001:0;
        datac1.trial{c} = data0.trial{trllop}(:,timselchr{trllop}(frtlop,1):timselchr{trllop}(frtlop,2));

%         % use this code when doing time-frequency analysis (padding at the end)
%         datac1.sampleinfo(c,:) = [data0.sampleinfo(trllop,1)+timselchr{trllop}(frtlop,1)-1 data0.sampleinfo(trllop,1)+timselchr{trllop}(frtlop,2)-1+500];
%         datac1.time{c} = -4.999:0.001:0.5;
%         datac1.trial{c} = data0.trial{trllop}(:,timselchr{trllop}(frtlop,1):timselchr{trllop}(frtlop,2)+500);

        c=c+1;
        
        filindchrdum = [filindchrdum; filindnrlsel(trllop)];
        
    end
    
end

clear data0

% subtract common average of each probe from all channels on that probe
datac2 = datac1;
% comment-out the following if choosing not to re-reference to common probe average
clear prblaball
for k=1:length(datac1.label)
    prblaball(k,1)=datac1.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);
for trllop = 1:length(datac1.trial)
    for prblop = 1:length(prblab)
        prbavg = mean(datac1.trial{trllop}(locb==prblop,:),1);
        repavg = repmat(prbavg,length(find(locb==prblop)),1);
        mat1 = reshape(datac1.trial{trllop}(locb==prblop,:),1,length(find(locb==prblop)),size(datac1.trial{trllop},2));
        mat2 = reshape(repavg,1,size(repavg,1),size(repavg,2));
        datac2.trial{trllop}(locb==prblop,:)=diff(cat(1,mat1,mat2),1,1);
    end
end
clear datac1

% now run ft_rejectvisual and ft_databrowser to clean bad trials

badtrl = [73, 104, 158]; % if re-referencing to common probe average

filindchr = filindchrdum(setxor(1:length(filindchrdum),badtrl));

% exclude bad trials, filter out line noise here
cfg=[];
cfg.trials = setxor(1:length(datac2.time),badtrl);
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
datac3 = ft_preprocessing(cfg,datac2);
clear datac2

cfg=[];
% cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
cfg.output      = 'powandcsd';
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'no';
cfg.foilim      = [1 30];
% cfg.taper       = 'hanning';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqcL = ft_freqanalysis(cfg, datac3);

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqcH = ft_freqanalysis(cfg, datac3);


%% look at time-resolved power (time-frequency analysis)

% do the spectral analysis - time-frequency
clear cfg
cfg.output      = 'powandcsd';
cfg.method      = 'mtmconvol';
cfg.foi         = 0:2:30;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.5 * ones(1,numfoi);
cfg.tapsmofrq   = 2   * ones(1,numfoi);
cfg.toi         = -4.75:0.01:0;
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'no';
freqTFL = ft_freqanalysis(cfg, data3);
cfg.foi         = 30:4:100;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.25 * ones(1,numfoi);
cfg.tapsmofrq   = 8   * ones(1,numfoi);
cfg.toi         = -4.875:0.01:0;
freqTFH = ft_freqanalysis(cfg, data3);


% normalize by cherry power
TFLnorm = nan(size(freqTFL.powspctrm));
for timlop = 1:length(freqTFL.time)
    for frqlop = 1:length(freqTFL.freq)
        [~,frqind] = min(abs(freqcL.freq-freqTFL.freq(frqlop)));
        TFLnorm(:,frqlop,timlop) = squeeze(freqTFL.powspctrm(:,frqlop,timlop))./freqcL.powspctrm(:,frqind);
    end
end
TFHnorm = nan(size(freqTFH.powspctrm));
for timlop = 1:length(freqTFH.time)
    for frqlop = 1:length(freqTFH.freq)
        [~,frqind] = min(abs(freqcH.freq-freqTFH.freq(frqlop)));
        TFHnorm(:,frqlop,timlop) = squeeze(freqTFH.powspctrm(:,frqlop,timlop))./freqcH.powspctrm(:,frqind);
    end
end

for chnlop = 1:size(TFLnorm,1)
    
    figure

    subplot(2,1,1)
    imagesc(freqTFH.time,freqTFH.freq,squeeze(TFHnorm(chnlop,:,:)));axis xy;colormap jet
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title([freqTFH.label{chnlop} ' (30-100 Hz)'])
    set(gca,'TickDir','out')
    colorbar
    
    subplot(2,1,2)
    imagesc(freqTFL.time,freqTFL.freq(2:end),squeeze(TFLnorm(chnlop,2:end,:)));axis xy;colormap jet
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title([freqTFL.label{chnlop} ' (4-30 Hz)'])
    set(gca,'TickDir','out')
    colorbar
    
%     figdir = 'R:\Buffalo Lab\Virtual Navigation\Figures\virtual water maze, neural analysis\power';
%     saveas(gcf,fullfile(figdir,['TFpow_banana-cherry_' freqTFL.label{chnlop} '.png']),'png')
    
end


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


%% using memindhi and memindlo, categorize neural data trials

% categorize trials into good/bad memory
[~,memsrt] = sort(memsel);
% when low memsel = high memory performance
memindhi = memsrt(1:round(length(memsel)/3));
memindlo = memsrt(end-(round(length(memsel)/3)-1):end);

nrlmemindhi = sort([memindhi*2-1; memindhi*2]);
nrlmemindlo = sort([memindlo*2-1; memindlo*2]);


% do the spectral analysis - time-frequency
clear cfg
cfg.output      = 'powandcsd';
cfg.method      = 'mtmconvol';
cfg.foi         = 0:2:30;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.5 * ones(1,numfoi);
cfg.tapsmofrq   = 2   * ones(1,numfoi);
cfg.toi         = -4.75:0.01:0;
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'no';
cfg.trials = nrlmemindhi;
freqTFLh = ft_freqanalysis(cfg, data3);
cfg.trials = nrlmemindlo;
freqTFLl = ft_freqanalysis(cfg, data3);
cfg.foi         = 30:4:100;
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = 0.25 * ones(1,numfoi);
cfg.tapsmofrq   = 8   * ones(1,numfoi);
cfg.toi         = -4.875:0.01:0;
cfg.trials = nrlmemindhi;
freqTFHh = ft_freqanalysis(cfg, data3);
cfg.trials = nrlmemindlo;
freqTFHl = ft_freqanalysis(cfg, data3);

% normalize by cherry power
TFLhnorm = nan(size(freqTFLh.powspctrm));
TFLlnorm = nan(size(freqTFLl.powspctrm));
for timlop = 1:length(freqTFLh.time)
    for frqlop = 1:length(freqTFLh.freq)
        [~,frqind] = min(abs(freqcL.freq-freqTFLh.freq(frqlop)));
        TFLhnorm(:,frqlop,timlop) = squeeze(freqTFLh.powspctrm(:,frqlop,timlop))./freqcL.powspctrm(:,frqind);
        TFLlnorm(:,frqlop,timlop) = squeeze(freqTFLl.powspctrm(:,frqlop,timlop))./freqcL.powspctrm(:,frqind);
    end
end
TFHhnorm = nan(size(freqTFHh.powspctrm));
TFHlnorm = nan(size(freqTFHl.powspctrm));
for timlop = 1:length(freqTFHh.time)
    for frqlop = 1:length(freqTFHh.freq)
        [~,frqind] = min(abs(freqcH.freq-freqTFHh.freq(frqlop)));
        TFHhnorm(:,frqlop,timlop) = squeeze(freqTFHh.powspctrm(:,frqlop,timlop))./freqcH.powspctrm(:,frqind);
        TFHlnorm(:,frqlop,timlop) = squeeze(freqTFHl.powspctrm(:,frqlop,timlop))./freqcH.powspctrm(:,frqind);
    end
end

% plot high and low memory spectrograms separately
for chnlop = 1:size(TFLnorm,1)
    
    figure

    subplot(2,2,1)
    imagesc(freqTFHh.time,freqTFHh.freq,squeeze(TFHhnorm(chnlop,:,:)));axis xy;colormap jet
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title([freqTFH.label{chnlop} ' (30-100 Hz), High Memory'])
    set(gca,'TickDir','out')
    colorbar
    clim1 = get(gca,'clim');

    subplot(2,2,2)
    imagesc(freqTFHl.time,freqTFHl.freq,squeeze(TFHlnorm(chnlop,:,:)));axis xy;colormap jet
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title([freqTFH.label{chnlop} ' (30-100 Hz), Low Memory'])
    set(gca,'TickDir','out')
    colorbar
    set(gca,'clim',clim1)
    
    subplot(2,2,3)
    imagesc(freqTFLh.time,freqTFLh.freq(2:end),squeeze(TFLhnorm(chnlop,2:end,:)));axis xy;colormap jet
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title([freqTFLh.label{chnlop} ' (4-30 Hz), High Memory'])
    set(gca,'TickDir','out')
    colorbar
    clim2 = get(gca,'clim');

    subplot(2,2,4)
    imagesc(freqTFLl.time,freqTFLl.freq(2:end),squeeze(TFLlnorm(chnlop,2:end,:)));axis xy;colormap jet
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title([freqTFLl.label{chnlop} ' (4-30 Hz), Low Memory'])
    set(gca,'TickDir','out')
    colorbar
    set(gca,'clim',clim2)
    
    set(gcf,'Position',[361 184 1199 794])
 
%     figdir = 'R:\Buffalo Lab\Virtual Navigation\Figures\virtual water maze, neural analysis\power';
%     saveas(gcf,fullfile(figdir,['TFpow_banana-cherry_' freqTFL.label{chnlop} 'hilomem.png']),'png')
    
end

% plot high and low memory difference spectrogram
for chnlop = 1:size(TFLnorm,1)
    
    figure

    subplot(2,1,1)
    imagesc(freqTFHh.time,freqTFHh.freq,squeeze(TFHhnorm(chnlop,:,:))-squeeze(TFHlnorm(chnlop,:,:)));axis xy;colormap jet
    clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['High-low memory difference, ' freqTFH.label{chnlop} ' (30-100 Hz)'])
    set(gca,'TickDir','out')
    colorbar
    clim1 = get(gca,'clim');

    subplot(2,1,2)
    imagesc(freqTFLh.time,freqTFLh.freq(2:end),squeeze(TFLhnorm(chnlop,2:end,:))-squeeze(TFLlnorm(chnlop,2:end,:)));axis xy;colormap jet
    clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
    xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['High-low memory difference, ' freqTFH.label{chnlop} ' (4-30 Hz)'])
    set(gca,'TickDir','out')
    colorbar
    set(gca,'clim',clim1)
    
%     set(gcf,'Position',[361 184 1199 794])
 
%     figdir = 'R:\Buffalo Lab\Virtual Navigation\Figures\virtual water maze, neural analysis\power';
%     saveas(gcf,fullfile(figdir,['TFpow_banana-cherry_' freqTFL.label{chnlop} 'hilomem.png']),'png')
    
end


% cfg=[];
% cfg.method      = 'coh';
% cfg.jackknife   = 'yes';
% 
% cfg.trials = nrlmemindhi;
% statLh = ft_connectivityanalysis(cfg,freqL);
% statHh = ft_connectivityanalysis(cfg,freqH);
% fdLh = ft_freqdescriptives(cfg, freqL);
% fdHh = ft_freqdescriptives(cfg, freqH);
% 
% cfg.trials = nrlmemindlo;
% statLl = ft_connectivityanalysis(cfg,freqL);
% statHl = ft_connectivityanalysis(cfg,freqH);
% fdLl = ft_freqdescriptives(cfg, freqL);
% fdHl = ft_freqdescriptives(cfg, freqH);


%% time-resolved power averaged within bands, across contacts on each array
% theta 4-8
% beta 12-22
% high gamma 60-100
[~,fiT1] = min(abs(freqTFL.freq-4));
[~,fiT2] = min(abs(freqTFL.freq-8));
[~,fiB1] = min(abs(freqTFL.freq-12));
[~,fiB2] = min(abs(freqTFL.freq-22));
[~,fiG1] = min(abs(freqTFH.freq-60));
[~,fiG2] = min(abs(freqTFH.freq-100));

clear prblaball
for k=1:length(freqTFL.label)
    prblaball(k,1)=freqTFL.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);

powByProbeT = [];
powByProbeB = [];
powByProbeG = [];
for prblop=1:3
    powByProbeT{prblop}(1,:,:) = squeeze(mean(TFLhnorm(locb==prblop,fiT1:fiT2,:),2));
    powByProbeB{prblop}(1,:,:) = squeeze(mean(TFLhnorm(locb==prblop,fiB1:fiB2,:),2));
    powByProbeG{prblop}(1,:,:) = squeeze(mean(TFHhnorm(locb==prblop,fiG1:fiG2,:),2));
    powByProbeT{prblop}(2,:,:) = squeeze(mean(TFLlnorm(locb==prblop,fiT1:fiT2,:),2));
    powByProbeB{prblop}(2,:,:) = squeeze(mean(TFLlnorm(locb==prblop,fiB1:fiB2,:),2));
    powByProbeG{prblop}(2,:,:) = squeeze(mean(TFHlnorm(locb==prblop,fiG1:fiG2,:),2));
end

figure
subplot(3,1,1)
plot(freqTFL.time,squeeze(mean(powByProbeT{1}(1,:,:),2)),'r')
hold on;plot(freqTFL.time,squeeze(mean(powByProbeT{1}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeT{1}(1,:,:)),1,freqTFL.time,'r')
errorshade(squeeze(powByProbeT{1}(2,:,:)),1,freqTFL.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('Theta power (4-8 Hz), Array A')
grid on; set(gca,'GridLineStyle','--')
subplot(3,1,2)
plot(freqTFL.time,squeeze(mean(powByProbeB{1}(1,:,:),2)),'r')
hold on;plot(freqTFL.time,squeeze(mean(powByProbeB{1}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeB{1}(1,:,:)),1,freqTFL.time,'r')
errorshade(squeeze(powByProbeB{1}(2,:,:)),1,freqTFL.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('Beta power (12-22 Hz), Array A')
grid on; set(gca,'GridLineStyle','--')
subplot(3,1,3)
plot(freqTFH.time,squeeze(mean(powByProbeG{1}(1,:,:),2)),'r')
hold on;plot(freqTFH.time,squeeze(mean(powByProbeG{1}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeG{1}(1,:,:)),1,freqTFH.time,'r')
errorshade(squeeze(powByProbeG{1}(2,:,:)),1,freqTFH.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('High gamma power (60-100 Hz), Array A')
grid on; set(gca,'GridLineStyle','--')

figure
subplot(3,1,1)
plot(freqTFL.time,squeeze(mean(powByProbeT{2}(1,:,:),2)),'r')
hold on;plot(freqTFL.time,squeeze(mean(powByProbeT{2}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeT{2}(1,:,:)),1,freqTFL.time,'r')
errorshade(squeeze(powByProbeT{2}(2,:,:)),1,freqTFL.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('Theta power (4-8 Hz), Array B')
grid on; set(gca,'GridLineStyle','--')
subplot(3,1,2)
plot(freqTFL.time,squeeze(mean(powByProbeB{2}(1,:,:),2)),'r')
hold on;plot(freqTFL.time,squeeze(mean(powByProbeB{2}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeB{2}(1,:,:)),1,freqTFL.time,'r')
errorshade(squeeze(powByProbeB{2}(2,:,:)),1,freqTFL.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('Beta power (12-22 Hz), Array B')
grid on; set(gca,'GridLineStyle','--')
subplot(3,1,3)
plot(freqTFH.time,squeeze(mean(powByProbeG{2}(1,:,:),2)),'r')
hold on;plot(freqTFH.time,squeeze(mean(powByProbeG{2}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeG{2}(1,:,:)),1,freqTFH.time,'r')
errorshade(squeeze(powByProbeG{2}(2,:,:)),1,freqTFH.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('High gamma power (60-100 Hz), Array B')
grid on; set(gca,'GridLineStyle','--')

figure
subplot(3,1,1)
plot(freqTFL.time,squeeze(mean(powByProbeT{3}(1,:,:),2)),'r')
hold on;plot(freqTFL.time,squeeze(mean(powByProbeT{3}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeT{3}(1,:,:)),1,freqTFL.time,'r')
errorshade(squeeze(powByProbeT{3}(2,:,:)),1,freqTFL.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('Theta power (4-8 Hz), Array C')
grid on; set(gca,'GridLineStyle','--')
subplot(3,1,2)
plot(freqTFL.time,squeeze(mean(powByProbeB{3}(1,:,:),2)),'r')
hold on;plot(freqTFL.time,squeeze(mean(powByProbeB{3}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeB{3}(1,:,:)),1,freqTFL.time,'r')
errorshade(squeeze(powByProbeB{3}(2,:,:)),1,freqTFL.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('Beta power (12-22 Hz), Array C')
grid on; set(gca,'GridLineStyle','--')
subplot(3,1,3)
plot(freqTFH.time,squeeze(mean(powByProbeG{3}(1,:,:),2)),'r')
hold on;plot(freqTFH.time,squeeze(mean(powByProbeG{3}(2,:,:),2)),'b')
errorshade(squeeze(powByProbeG{3}(1,:,:)),1,freqTFH.time,'r')
errorshade(squeeze(powByProbeG{3}(2,:,:)),1,freqTFH.time,'b')
xlabel('Time (sec)'); ylabel('Power'); title('High gamma power (60-100 Hz), Array C')
grid on; set(gca,'GridLineStyle','--')


%% plot average spectrogram and difference across contacts on each array

% get average spectrogram for each array
clear prblaball
for k=1:size(freqTFH.label,1)
    prblaball(k,:)=freqTFH.label{k}(1);
end
prblab = unique(prblaball,'rows');
[~,locb] = ismember(prblaball,prblab,'rows');

powL_A = TFLnorm(locb==1,:,:); powH_A = TFHnorm(locb==1,:,:);
powL_B = TFLnorm(locb==2,:,:); powH_B = TFHnorm(locb==2,:,:);
% powL_C = TFLnorm(locb==3,:,:); powH_C = TFHnorm(locb==3,:,:);
powL_C = TFLnorm([18:21 23:25 27 29],:,:); powH_C = TFHnorm([18:21 23:25 27 29],:,:); % take out C05, C09 and C11, weird 60 Hz thing

% plot average coherogram for each probe combo
%A    
figure
subplot(2,1,1)
imagesc(freqTFH.time,freqTFH.freq,squeeze(mean(powH_A,1)));axis xy;colormap jet
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['A (30-100 Hz)'])
set(gca,'TickDir','out')
colorbar
subplot(2,1,2)
imagesc(freqTFL.time,freqTFL.freq(2:end),squeeze(mean(powL_A(:,2:end,:),1)));axis xy;colormap jet
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['A (4-30 Hz)'])
set(gca,'TickDir','out')
colorbar
%B    
figure
subplot(2,1,1)
imagesc(freqTFH.time,freqTFH.freq,squeeze(mean(powH_B,1)));axis xy;colormap jet
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['B (30-100 Hz)'])
set(gca,'TickDir','out')
colorbar
subplot(2,1,2)
imagesc(freqTFL.time,freqTFL.freq(2:end),squeeze(mean(powL_B(:,2:end,:),1)));axis xy;colormap jet
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['B (4-30 Hz)'])
set(gca,'TickDir','out')
colorbar
%C    
figure
subplot(2,1,1)
imagesc(freqTFH.time,freqTFH.freq,squeeze(mean(powH_C,1)));axis xy;colormap jet
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['C (30-100 Hz)'])
set(gca,'TickDir','out')
colorbar
subplot(2,1,2)
imagesc(freqTFL.time,freqTFL.freq(2:end),squeeze(mean(powL_C(:,2:end,:),1)));axis xy;colormap jet
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['C (4-30 Hz)'])
set(gca,'TickDir','out')
colorbar



% plot high and low memory difference spectrogram
powLhl_A = squeeze(mean(TFLhnorm(locb==1,:,:)-TFLlnorm(locb==1,:,:),1)); powHhl_A = squeeze(mean(TFHhnorm(locb==1,:,:)-TFHlnorm(locb==1,:,:),1));
powLhl_B = squeeze(mean(TFLhnorm(locb==2,:,:)-TFLlnorm(locb==2,:,:),1)); powHhl_B = squeeze(mean(TFHhnorm(locb==2,:,:)-TFHlnorm(locb==2,:,:),1));
% powLhl_C = squeeze(mean(TFLhnorm(locb==3,:,:)-TFLlnorm(locb==3,:,:),1)); powHhl_C = squeeze(mean(TFHhnorm(locb==3,:,:)-TFHlnorm(locb==3,:,:),1));
powLhl_C = squeeze(mean(TFLhnorm([18:21 23:25 27 29],:,:)-TFLlnorm([18:21 23:25 27 29],:,:),1)); powHhl_C = squeeze(mean(TFHhnorm([18:21 23:25 27 29],:,:)-TFHlnorm([18:21 23:25 27 29],:,:),1)); % take out C05, C09 and C11, weird 60 Hz thing

%A    
figure
subplot(2,1,1)
imagesc(freqTFH.time,freqTFH.freq,powHhl_A);axis xy;colormap jet
clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['A (30-100 Hz)'])
set(gca,'TickDir','out')
colorbar
subplot(2,1,2)
imagesc(freqTFL.time,freqTFL.freq(2:end),powLhl_A(2:end,:));axis xy;colormap jet
clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['A (4-30 Hz)'])
set(gca,'TickDir','out')
colorbar
%B    
figure
subplot(2,1,1)
imagesc(freqTFH.time,freqTFH.freq,powHhl_B);axis xy;colormap jet
clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['B (30-100 Hz)'])
set(gca,'TickDir','out')
colorbar
subplot(2,1,2)
imagesc(freqTFL.time,freqTFL.freq(2:end),powLhl_B(2:end,:));axis xy;colormap jet
clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['B (4-30 Hz)'])
set(gca,'TickDir','out')
colorbar
%C    
figure
subplot(2,1,1)
imagesc(freqTFH.time,freqTFH.freq,powHhl_C);axis xy;colormap jet
clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['C (30-100 Hz)'])
set(gca,'TickDir','out')
colorbar
subplot(2,1,2)
imagesc(freqTFL.time,freqTFL.freq(2:end),powLhl_C(2:end,:));axis xy;colormap jet
clim1=get(gca,'clim');set(gca,'clim',[-abs(max(clim1)) abs(max(clim1))]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title(['C (4-30 Hz)'])
set(gca,'TickDir','out')
colorbar



%% get trial-by-trial coherence ("leave one out")
% 
% % create cohData for multivariate
% cohData = nan(length(data3.trial),3,size(statL.labelcmb,1));
% for trllop = 1:length(data3.trial)
%     
%     fprintf('Analyzing trial %d.\n',trllop)
%     
%     % remove the marked trial
%     cfg=[];
%     cfg.trials = setxor(1:length(data3.trial),trllop);
%     dum = ft_redefinetrial(cfg,data3);
%     
%     % do the spectral analysis - time-averaged
%     cfg=[];
%     cfg.output      = 'powandcsd';
%     cfg.method      = 'mtmfft';
%     cfg.pad         = 'maxperlen';
%     cfg.keeptrials  = 'yes';
%     cfg.foilim      = [1 30];
%     cfg.taper       = 'dpss';
%     cfg.tapsmofrq   = 2;
%     freqLdum = ft_freqanalysis(cfg, dum);
% 
%     cfg.foilim      = [30 100];
%     cfg.tapsmofrq   = 8;
%     freqHdum = ft_freqanalysis(cfg, dum);
% 
%     cfg=[];
%     cfg.method      = 'coh';
%     statLdum = ft_connectivityanalysis(cfg,freqLdum);
%     statHdum = ft_connectivityanalysis(cfg,freqHdum);
%     
%     difmatL = statL.cohspctrm-statLdum.cohspctrm;
%     difmatH = statH.cohspctrm-statHdum.cohspctrm;
% 
%     [~,f1_t] = min(abs(statLdum.freq-3));
%     [~,f2_t] = min(abs(statLdum.freq-12));
%     [~,f1_lg] = min(abs(statHdum.freq-30));
%     [~,f2_lg] = min(abs(statHdum.freq-60));
%     [~,f1_hg] = min(abs(statHdum.freq-60));
%     [~,f2_hg] = min(abs(statHdum.freq-100));
% 
%     cohData(trllop,1,:) = mean(difmatL(:,f1_t:f2_t),2);
%     cohData(trllop,2,:) = mean(difmatH(:,f1_lg:f2_lg),2);
%     cohData(trllop,3,:) = mean(difmatH(:,f1_hg:f2_hg),2);
% 
% end
% 
% 
% %% plot histograms of coherence values for each probe pairing, for all trials
% 
% clear prblaball
% for k=1:length(statL.label)
%     prblaball(k,1)=statL.label{k}(1);
% end
% prblab = unique(prblaball);
% [~,locb] = ismember(prblaball,prblab);
% 
% cmblst = unique(sort([repmat(unique(locb),length(unique(locb)),1) sort(repmat(unique(locb),length(unique(locb)),1))],2),'rows');
% 
% [~,f1_t] = min(abs(statL.freq-3));
% [~,f2_t] = min(abs(statL.freq-12));
% [~,f1_lg] = min(abs(statH.freq-30));
% [~,f2_lg] = min(abs(statH.freq-60));
% [~,f1_hg] = min(abs(statH.freq-60));
% [~,f2_hg] = min(abs(statH.freq-100));
% 
% coh_theta = cell(size(cmblst,1),1); % 3-12 Hz
% coh_logamma = cell(size(cmblst,1),1); % 30-60 Hz
% coh_higamma = cell(size(cmblst,1),1); % 60-100 Hz
% for cmblop = 1:size(cmblst,1)
%     
%     cmbind1 = find(locb==cmblst(cmblop,1));
%     cmbind2 = find(locb==cmblst(cmblop,2));
%     cmblst2 = unique(sort([repmat(unique(cmbind1),length(unique(cmbind2)),1) sort(repmat(unique(cmbind2),length(unique(cmbind1)),1))],2),'rows');
%     cmblst2 = cmblst2(cmblst2(:,1)~=cmblst2(:,2),:);
%     
%     cohmat1 = []; cohmat2 = [];
%     for indlop1 = 1:size(cmblst2,1)
%         cohmat1 = cat(1,cohmat1,(squeeze(statL.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat2 = cat(1,cohmat2,(squeeze(statH.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%     end
%         
%     coh_theta{cmblop} = mean(cohmat1(:,f1_t:f2_t),2);
%     coh_logamma{cmblop} = mean(cohmat2(:,f1_lg:f2_lg),2);
%     coh_higamma{cmblop} = mean(cohmat2(:,f1_hg:f2_hg),2);
%     
% end
% 
% % theta
% edg = 0:0.01:0.5;
% edx = edg + 0.005;
% [n1,~] = histc(coh_theta{1},edg);
% [n2,~] = histc(coh_theta{2},edg);
% [n3,~] = histc(coh_theta{3},edg);
% [n4,~] = histc(coh_theta{4},edg);
% [n5,~] = histc(coh_theta{5},edg);
% [n6,~] = histc(coh_theta{6},edg);
% figure;hold on
% % bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
% bar(edx,([n2 n3 n5]),'stacked')
% xlim([0 0.4])
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
% legend({'A x B','A x C','B x C'},'Location','northeast')
% set(gca,'TickDir','out')
% box off
% title('Theta (3-12 Hz)')
% xlabel('Coherence')
% ylabel('# of pairs')
% 
% % low gamma
% edg = 0:0.01:0.5;
% edx = edg + 0.005;
% [n1,~] = histc(coh_logamma{1},edg);
% [n2,~] = histc(coh_logamma{2},edg);
% [n3,~] = histc(coh_logamma{3},edg);
% [n4,~] = histc(coh_logamma{4},edg);
% [n5,~] = histc(coh_logamma{5},edg);
% [n6,~] = histc(coh_logamma{6},edg);
% figure;hold on
% % bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
% bar(edx,([n2 n3 n5]),'stacked')
% xlim([0 0.4])
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
% legend({'A x B','A x C','B x C'},'Location','northeast')
% set(gca,'TickDir','out')
% box off
% title('Low Gamma (30-60 Hz)')
% xlabel('Coherence')
% ylabel('# of pairs')
% 
% % high gamma
% edg = 0:0.01:0.5;
% edx = edg + 0.005;
% [n1,~] = histc(coh_higamma{1},edg);
% [n2,~] = histc(coh_higamma{2},edg);
% [n3,~] = histc(coh_higamma{3},edg);
% [n4,~] = histc(coh_higamma{4},edg);
% [n5,~] = histc(coh_higamma{5},edg);
% [n6,~] = histc(coh_higamma{6},edg);
% figure;hold on
% % bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
% bar(edx,([n2 n3 n5]),'stacked')
% xlim([0 0.4])
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
% legend({'A x B','A x C','B x C'},'Location','northeast')
% set(gca,'TickDir','out')
% box off
% title('High Gamma (60-100 Hz)')
% xlabel('Coherence')
% ylabel('# of pairs')
% 
% 
% %% get invisible banana positions for selected trials
% 
% x = linspace(-10,10,17);
% % figure
% % hist3(banposall(nrlsel,:),{x,x});
% [n,c] = hist3(banposall(nrlbadsel,:),'Edges',{x,x});
% figure;imagesc(c{1},c{2},n);xlim([-10 10]);ylim([-10 10]);colormap hot;colorbar
% hold on;scatter(banposall(nrlbadsel,2),banposall(nrlbadsel,1),'bo')
% 
% % calculate distance from banana at start position
% xdif = begposall(nrlbadsel,1)-banposall(nrlbadsel,1);
% ydif = begposall(nrlbadsel,2)-banposall(nrlbadsel,2);
% begdis = sqrt(xdif.^2+ydif.^2);
% 
% 
% %% heat map, memory performance, with disk "smoothing"
% 
% eplrat3d = []; % excess path length, ratio
% mem3d = [];
% for k=1:length(nrlbadsel)
%     
%     posn = imfilter(hist3(banposall(nrlbadsel(k),:),'Edges',{x,x}),fspecial('disk',10));
%     posn(posn~=max(max(posn))) = nan; posn(~isnan(posn))=1;
%     
%     % multiply by scaling factor to correct for smoothing (to plot correct amplitudes)
%     eplrat3d(k,:,:) = posn(1:end-1,1:end-1)*excpthnrm(nrlbadsel(k));
%     mem3d(k,:,:) = posn(1:end-1,1:end-1)*memsel(k);
%     
% %     figure
% %     subplot(2,1,1)
% %     imagesc(ctrs,ctrs,n(1:end-1,1:end-1));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet
% %     subplot(2,1,2)
% %     imagesc(ctrs,ctrs,squeeze(nanmean(mem3d,1)));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet
% % %     set(gca,'clim',[0 0.002])
% %     set(gcf,'Position',[80 67 496 911])
% %     pause
% %     close
%     
% end
% 
% eplrat2d = squeeze(nanmedian(eplrat3d,1));
% mem2d = squeeze(nanmedian(mem3d,1));
% 
% figure;imagesc(ctrs,ctrs,eplrat2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (ratio) by banana location');colorbar
% figure;imagesc(ctrs,ctrs,mem2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Memory (performance factor)');colorbar
% 
% 
% %% High vs. Low memory, scatter, each pair (color coded)
% % create a different plot for each frequency/range
% % OLD THETA AND LOW/HIGH GAMMA DESIGNATIONS; USE NEXT CELL INSTEAD
% 
% clear prblaball
% for k=1:length(statLh.label)
%     prblaball(k,1)=statLh.label{k}(1);
% end
% prblab = unique(prblaball);
% [~,locb] = ismember(prblaball,prblab);
% 
% cmblst = unique(sort([repmat(unique(locb),length(unique(locb)),1) sort(repmat(unique(locb),length(unique(locb)),1))],2),'rows');
% 
% [~,f1_t] = min(abs(statLh.freq-3));
% [~,f2_t] = min(abs(statLh.freq-12));
% [~,f1_lg] = min(abs(statHh.freq-30));
% [~,f2_lg] = min(abs(statHh.freq-60));
% [~,f1_hg] = min(abs(statHh.freq-60));
% [~,f2_hg] = min(abs(statHh.freq-100));
% 
% cohmemdif_theta = cell(size(cmblst,1),1); % 3-12 Hz
% cohmemdif_logamma = cell(size(cmblst,1),1); % 30-60 Hz
% cohmemdif_higamma = cell(size(cmblst,1),1); % 60-100 Hz
% % numprs = [];
% for cmblop = 1:size(cmblst,1)
%     
%     cmbind1 = find(locb==cmblst(cmblop,1));
%     cmbind2 = find(locb==cmblst(cmblop,2));
%     cmblst2 = unique(sort([repmat(unique(cmbind1),length(unique(cmbind2)),1) sort(repmat(unique(cmbind2),length(unique(cmbind1)),1))],2),'rows');
%     cmblst2 = cmblst2(cmblst2(:,1)~=cmblst2(:,2),:);
%     
% %     numprs(cmblop) = size(cmblst2,1);
%     
%     cohmat1 = []; cohmat2 = []; cohmat3 = []; cohmat4 = [];
%     for indlop1 = 1:size(cmblst2,1)
%         
%         cohmat1 = cat(1,cohmat1,(squeeze(statLh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat2 = cat(1,cohmat2,(squeeze(statLl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat3 = cat(1,cohmat3,(squeeze(statHh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat4 = cat(1,cohmat4,(squeeze(statHl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         
%     end
%         
%     cohmemdif_theta{cmblop} = [mean(cohmat1(:,f1_t:f2_t),2) mean(cohmat2(:,f1_t:f2_t),2)];
%     cohmemdif_logamma{cmblop} = [mean(cohmat3(:,f1_lg:f2_lg),2) mean(cohmat4(:,f1_lg:f2_lg),2)];
%     cohmemdif_higamma{cmblop} = [mean(cohmat3(:,f1_hg:f2_hg),2) mean(cohmat4(:,f1_hg:f2_hg),2)];
%     
% end
% 
% % theta
% figure; hold on
% % h1=scatter(cohmemdif_theta{1}(:,1),cohmemdif_theta{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_theta{2}(:,1),cohmemdif_theta{2}(:,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(cohmemdif_theta{3}(:,1),cohmemdif_theta{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% % h4=scatter(cohmemdif_theta{4}(:,1),cohmemdif_theta{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_theta{5}(:,1),cohmemdif_theta{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% % h6=scatter(cohmemdif_theta{6}(:,1),cohmemdif_theta{6}(:,2)); set(h6,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
% % legend({'A x C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('Theta (3-12 Hz)')
% xlabel('Coherence: High Memory')
% ylabel('Coherence: Low Memory')
% 
% % low gamma
% figure; hold on
% % h1=scatter(cohmemdif_logamma{1}(:,1),cohmemdif_logamma{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_logamma{2}(:,1),cohmemdif_logamma{2}(:,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(cohmemdif_logamma{3}(:,1),cohmemdif_logamma{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% % h4=scatter(cohmemdif_logamma{4}(:,1),cohmemdif_logamma{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_logamma{5}(:,1),cohmemdif_logamma{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% % h6=scatter(cohmemdif_logamma{6}(:,1),cohmemdif_logamma{6}(:,2)); set(h6,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
% % legend({'A x C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('Low Gamma (30-60 Hz)')
% xlabel('Coherence: High Memory')
% ylabel('Coherence: Low Memory')
% 
% % high gamma
% figure; hold on
% % h1=scatter(cohmemdif_higamma{1}(:,1),cohmemdif_higamma{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_higamma{2}(:,1),cohmemdif_higamma{2}(:,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(cohmemdif_higamma{3}(:,1),cohmemdif_higamma{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% % h4=scatter(cohmemdif_higamma{4}(:,1),cohmemdif_higamma{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_higamma{5}(:,1),cohmemdif_higamma{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% % h6=scatter(cohmemdif_higamma{6}(:,1),cohmemdif_higamma{6}(:,2)); set(h6,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
% % legend({'A x C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('High Gamma (60-100 Hz)')
% xlabel('Coherence: High Memory')
% ylabel('Coherence: Low Memory')
% 
% % % do histograms instead of scatter plots
% % edg = -0.05:0.005:0.05;
% % edx = edg + 0.0025;
% % [n1,b1] = histc(-diff(cohmemdif_theta{1},1,2),edg);
% % [n2,b2] = histc(-diff(cohmemdif_theta{2},1,2),edg);
% % [n3,b3] = histc(-diff(cohmemdif_theta{3},1,2),edg);
% % [n4,b4] = histc(-diff(cohmemdif_theta{4},1,2),edg);
% % [n5,b5] = histc(-diff(cohmemdif_theta{5},1,2),edg);
% % [n6,b6] = histc(-diff(cohmemdif_theta{6},1,2),edg);
% % figure;hold on
% % bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
% % line([0 0],ylim,'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
% % set(gca,'TickDir','out')
% % box off
% % title('Theta (3-12 Hz)')
% % xlabel('Coherence difference (High - Low memory)')
% % ylabel('# of pairs')
% % 
% % edg = -0.05:0.005:0.05;
% % edx = edg + 0.0025;
% % [n1,b1] = histc(-diff(cohmemdif_logamma{1},1,2),edg);
% % [n2,b2] = histc(-diff(cohmemdif_logamma{2},1,2),edg);
% % [n3,b3] = histc(-diff(cohmemdif_logamma{3},1,2),edg);
% % [n4,b4] = histc(-diff(cohmemdif_logamma{4},1,2),edg);
% % [n5,b5] = histc(-diff(cohmemdif_logamma{5},1,2),edg);
% % [n6,b6] = histc(-diff(cohmemdif_logamma{6},1,2),edg);
% % figure;hold on
% % bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
% % line([0 0],ylim,'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
% % set(gca,'TickDir','out')
% % box off
% % title('Low Gamma (30-60 Hz)')
% % xlabel('Coherence difference (High - Low memory)')
% % ylabel('# of pairs')
% % 
% % edg = -0.05:0.005:0.05;
% % edx = edg + 0.0025;
% % [n1,b1] = histc(-diff(cohmemdif_higamma{1},1,2),edg);
% % [n2,b2] = histc(-diff(cohmemdif_higamma{2},1,2),edg);
% % [n3,b3] = histc(-diff(cohmemdif_higamma{3},1,2),edg);
% % [n4,b4] = histc(-diff(cohmemdif_higamma{4},1,2),edg);
% % [n5,b5] = histc(-diff(cohmemdif_higamma{5},1,2),edg);
% % [n6,b6] = histc(-diff(cohmemdif_higamma{6},1,2),edg);
% % figure;hold on
% % bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
% % line([0 0],ylim,'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
% % set(gca,'TickDir','out')
% % box off
% % title('High Gamma (60-100 Hz)')
% % xlabel('Coherence difference (High - Low memory)')
% % ylabel('# of pairs')
% 
% 
% %% High vs. Low memory, scatter, each pair (color coded)
% % create a different plot for each frequency/range
% % NEW THETA AND GAMMA DESIGNATIONS
% 
% clear prblaball
% for k=1:length(statLh.label)
%     prblaball(k,1)=statLh.label{k}(1);
% end
% prblab = unique(prblaball);
% [~,locb] = ismember(prblaball,prblab);
% 
% cmblst = unique(sort([repmat(unique(locb),length(unique(locb)),1) sort(repmat(unique(locb),length(unique(locb)),1))],2),'rows');
% 
% [~,f1_t] = min(abs(statLh.freq-4));
% [~,f2_t] = min(abs(statLh.freq-9));
% [~,f1_g] = min(abs(statHh.freq-40));
% [~,f2_g] = min(abs(statHh.freq-90));
% 
% cohmemdif_theta = cell(size(cmblst,1),1); % 3-9 Hz
% cohmemdif_gamma = cell(size(cmblst,1),1); % 40-90 Hz
% 
% % numprs = [];
% for cmblop = 1:size(cmblst,1)
%     
%     cmbind1 = find(locb==cmblst(cmblop,1));
%     cmbind2 = find(locb==cmblst(cmblop,2));
%     cmblst2 = unique(sort([repmat(unique(cmbind1),length(unique(cmbind2)),1) sort(repmat(unique(cmbind2),length(unique(cmbind1)),1))],2),'rows');
%     cmblst2 = cmblst2(cmblst2(:,1)~=cmblst2(:,2),:);
%     
% %     numprs(cmblop) = size(cmblst2,1);
%     
%     cohmat1 = []; cohmat2 = []; cohmat3 = []; cohmat4 = [];
%     for indlop1 = 1:size(cmblst2,1)
%         
%         cohmat1 = cat(1,cohmat1,(squeeze(statLh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat2 = cat(1,cohmat2,(squeeze(statLl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat3 = cat(1,cohmat3,(squeeze(statHh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         cohmat4 = cat(1,cohmat4,(squeeze(statHl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
%         
%     end
%         
%     cohmemdif_theta{cmblop} = [mean(cohmat1(:,f1_t:f2_t),2) mean(cohmat2(:,f1_t:f2_t),2)];
%     cohmemdif_gamma{cmblop} = [mean(cohmat3(:,f1_g:f2_g),2) mean(cohmat4(:,f1_g:f2_g),2)];
%     
% end
% 
% % theta
% figure; hold on
% % h1=scatter(cohmemdif_theta{1}(:,1),cohmemdif_theta{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_theta{2}(:,1),cohmemdif_theta{2}(:,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(cohmemdif_theta{3}(:,1),cohmemdif_theta{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% % h4=scatter(cohmemdif_theta{4}(:,1),cohmemdif_theta{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_theta{5}(:,1),cohmemdif_theta{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% % h6=scatter(cohmemdif_theta{6}(:,1),cohmemdif_theta{6}(:,2)); set(h6,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
% % legend({'A x C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('Theta (4-9 Hz)')
% xlabel('Coherence: High Memory')
% ylabel('Coherence: Low Memory')
% 
% % gamma
% figure; hold on
% % h1=scatter(cohmemdif_gamma{1}(:,1),cohmemdif_gamma{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_gamma{2}(:,1),cohmemdif_gamma{2}(:,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(cohmemdif_gamma{3}(:,1),cohmemdif_gamma{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% % h4=scatter(cohmemdif_gamma{4}(:,1),cohmemdif_gamma{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_gamma{5}(:,1),cohmemdif_gamma{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% % h6=scatter(cohmemdif_gamma{6}(:,1),cohmemdif_gamma{6}(:,2)); set(h6,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% % legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
% % legend({'A x C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('Gamma (40-90 Hz)')
% xlabel('Coherence: High Memory')
% ylabel('Coherence: Low Memory')
% 
% 
% %% power, scatter plot
% 
% % create a different plot for each frequency/range
% clear prblaball
% for k=1:length(fdLh.label)
%     prblaball(k,1)=fdLh.label{k}(1);
% end
% prblab = unique(prblaball);
% [~,locb] = ismember(prblaball,prblab);
% 
% [~,f1_t] = min(abs(fdLh.freq-4));
% [~,f2_t] = min(abs(fdLh.freq-9));
% [~,f1_g] = min(abs(fdHh.freq-40));
% [~,f2_g] = min(abs(fdHh.freq-90));
% 
% powmemdif_theta = nan(size(fdLh.powspctrm,1),2); % 4-9 Hz
% powmemdif_gamma = nan(size(fdLh.powspctrm,1),2); % 40-90 Hz
% powmat1 = []; powmat2 = []; powmat3 = []; powmat4 = [];
% % numprs = [];
% for indlop = 1:size(fdLh.powspctrm,1)
%     
%     
%     powmat1 = cat(1,powmat1,squeeze(fdLh.powspctrm(indlop,:)));
%     powmat2 = cat(1,powmat2,squeeze(fdLl.powspctrm(indlop,:)));
%     powmat3 = cat(1,powmat3,squeeze(fdHh.powspctrm(indlop,:)));
%     powmat4 = cat(1,powmat4,squeeze(fdHl.powspctrm(indlop,:)));
%         
%     powmemdif_theta(indlop,:) = [mean(powmat1(indlop,f1_t:f2_t),2) mean(powmat2(indlop,f1_t:f2_t),2)];
%     powmemdif_gamma(indlop,:) = [mean(powmat3(indlop,f1_g:f2_g),2) mean(powmat4(indlop,f1_g:f2_g),2)];
%     
% end
% 
% % theta
% figure; hold on
% h1=scatter(powmemdif_theta(locb==1,1),powmemdif_theta(locb==1,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(powmemdif_theta(locb==2,1),powmemdif_theta(locb==2,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(powmemdif_theta(locb==3,1),powmemdif_theta(locb==3,2)); set(h3,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% legend({'A' 'B' 'C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('Theta (3-12 Hz)')
% xlabel('Power: High Memory')
% ylabel('Power: Low Memory')
% 
% % gamma
% figure; hold on
% h1=scatter(powmemdif_gamma(locb==1,1),powmemdif_gamma(locb==1,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(powmemdif_gamma(locb==2,1),powmemdif_gamma(locb==2,2)); set(h2,'MarkerFaceColor','flat')
% h3=scatter(powmemdif_gamma(locb==3,1),powmemdif_gamma(locb==3,2)); set(h3,'MarkerFaceColor','flat')
% xh=xlim; yh=ylim;
% line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% legend({'A' 'B' 'C'},'Location','southeast')
% set(gca,'TickDir','out')
% title('Gamma (40-90 Hz)')
% xlabel('Power: High Memory')
% ylabel('Power: Low Memory')
% 
% % plot histograms instead of scatter plots
% % theta
% edg = -3:0.5:3;
% edx = edg + 0.25;
% [n1,~] = histc(-diff(powmemdif_theta(locb==1,:),1,2),edg);
% [n2,~] = histc(-diff(powmemdif_theta(locb==2,:),1,2),edg);
% [n3,~] = histc(-diff(powmemdif_theta(locb==3,:),1,2),edg);
% figure;hold on
% bar(edx,([n1 n2 n3]),'stacked')
% xlim([-3 3])
% line([0 0],ylim,'Color','k')
% legend({'A','B','C'},'Location','northeast')
% set(gca,'TickDir','out')
% box off
% title('Theta (3-12 Hz)')
% xlabel('Power')
% ylabel('# of contacts')
% 
% % gamma
% edg = -0.2:0.05:0.2;
% edx = edg + 0.025;
% [n1,~] = histc(-diff(powmemdif_gamma(locb==1,:),1,2),edg);
% [n2,~] = histc(-diff(powmemdif_gamma(locb==2,:),1,2),edg);
% [n3,~] = histc(-diff(powmemdif_gamma(locb==3,:),1,2),edg);
% figure;hold on
% bar(edx,([n1 n2 n3]),'stacked')
% xlim([-0.2 0.2])
% line([0 0],ylim,'Color','k')
% legend({'A','B','C'},'Location','northeast')
% set(gca,'TickDir','out')
% box off
% title('Gamma (40-90 Hz)')
% xlabel('Power')
% ylabel('# of contacts')
% 
% 
% %% correlation
% 
% cfg=[];
% cfg.keeptrials   = 'yes';
% ftL = ft_freqdescriptives(cfg, freqL);
% ftH = ft_freqdescriptives(cfg, freqH);
% 
% save('R:\Buffalo Lab\Mike\VirtualNav\multivariate\varsformultivarGiz150904.mat', ...
%     'ftL','ftH','memsel');
% 
% [~,f1_t] = min(abs(ftL.freq-4));
% [~,f2_t] = min(abs(ftL.freq-9));
% [~,f1_g] = min(abs(ftH.freq-40));
% [~,f2_g] = min(abs(ftH.freq-90));
% 
% powforcorTht = squeeze(mean(ftL.powspctrm(:,:,f1_t:f2_t),3));
% powforcorGam = squeeze(mean(ftH.powspctrm(:,:,f1_g:f2_g),3));
% memforcor = memsel(round((1:length(memsel)*2)/2));
% 
% cormatTht = corr(powforcorTht,memforcor);
% cormatGam = corr(powforcorGam,memforcor);
% 
% figure;bar([1 2],[mean(cormatTht) mean(cormatGam)])
% hold on;errorbar([1 2],[mean(cormatTht) mean(cormatGam)], ...
%     [std(cormatTht)/sqrt(length(cormatTht)) std(cormatTht)/sqrt(length(cormatTht))], ...
%     [std(cormatTht)/sqrt(length(cormatTht)) std(cormatTht)/sqrt(length(cormatTht))], ...
%     'Marker','none','LineStyle','none','Color','k')
% xlim([0.5 2.5])
% box off
% set(gca,'XTickLabel',{'Theta (4-9 Hz)' 'Gamma (40-90 Hz)'})
% ylabel('Mean Pearson Coeff.')
% % title('Re-referenced to probe average')
% title('No re-referencing')
% 
% figure
% hist(memsel)
% hold on
% line([median(memsel) median(memsel)],ylim,'Color','r','LineWidth',2)
% xlabel('Performance factor')
% ylabel('Trial count')
% box off
% set(gca,'TickDir','out')
% 
