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
%     'JN_BR_15_04_15\JN_BR_15_04_15_12_32'; % 2 bad channels on A
    'JN_BR_15_04_16\JN_BR_15_04_16_14_22';
    'JN_BR_15_04_17\JN_BR_15_04_17_13_14';
    'JN_BR_15_04_20\JN_BR_15_04_20_14_43';
%     'JN_BR_15_04_21\JN_BR_15_04_21_12_34'; % 2 bad channels on A
    'JN_BR_15_04_22\JN_BR_15_04_22_14_09';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_54';
%     'JN_BR_15_04_23\JN_BR_15_04_23_13_54'; % A bad these 3 recordings
%     'JN_BR_15_04_23\JN_BR_15_04_23_14_22';
%     'JN_BR_15_04_24\JN_BR_15_04_24_14_04';
%     'JN_BR_15_04_27\JN_BR_15_04_27_14_35';
    'JN_BR_15_04_28\JN_BR_15_04_28_13_07';
    'JN_BR_15_04_29\JN_BR_15_04_29_13_22';
    };


% initialize dum structure
dum = [];
dum.time = cell(1);
dum.trial = cell(1);
dum.sampleinfo = [];
c=1;

begposall = []; % begin position (trial start)
alpall = []; % alpha
banfndall = []; % 1 if got banana, 0 if not
banposall = []; % banana position (first banana, includes invisible)
bantimindall = []; % latency to get banana (technically, time index)
avgbandispertrl = []; % average distance from banana
cumdispertrl = []; % cumulative distance from banana
avgwaldispertrl = []; % average distance from wall
trllngall = []; % total trial length
invprmindall = []; % invisible, prime (1st in a series of inv. bananas)
prvbanalpall = cell(1); % alpha of bananas preceding each prime inv. banana, only for inv. bananas with at least 3 preceding visible
invprmnrlsel = []; % trials used to select neural activity (dum)
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

%     sesttl=sesnam;sesttl(sesttl=='_')='-';
%     figure;plot(1:size(data.sampleinfo,1),data.sampleinfo);title(sesttl)
    
    numtrl(seslop,1) = length(trldat.time);
    
    % calculate three metrics for each trial:
    % latency to get banana (equal to time-out if banana not acquired)
    % average distance from banana
    % cumulative distance from banana
    
    invprmind = [];
    prvbanalp = cell(1);
    excpthdifdum = zeros(length(trldat.posdat),1);
    excpthnrmdum = zeros(length(trldat.posdat),1);
    excpthinddum = zeros(length(trldat.posdat),1);
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
        
        begposall = [begposall; trldat.posdat{trllop}(:,1)'];
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
                
                excpthdifdum(trllop,1) = sum(abs(diff(dis)))-dis(1);
                excpthnrmdum(trllop,1) = sum(abs(diff(dis)))/dis(1);
                excpthinddum(trllop,1) = (sum(abs(diff(dis)))-dis(1))/(sum(abs(diff(dis)))+dis(1));

            end
        end
        
    end

    for invlop = 1:length(prvbanalp)
        if length(prvbanalp{invlop})>=3
            prvbanalpall(c:c+2,1) = prvbanalp(invlop);
            invprmnrlsel = [invprmnrlsel; invprmind(invlop)+d];
            dum.time(c:c+2) = data.time(invprmind(invlop)-3:invprmind(invlop)-1);
            dum.trial(c:c+2) = data.trial(invprmind(invlop)-3:invprmind(invlop)-1);
            dum.sampleinfo(c:c+2,:) = data.sampleinfo(invprmind(invlop)-3:invprmind(invlop)-1,:);
            c = c+3;
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


%% save/load variables

% run this after running previous cell
savdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\misc';
save(fullfile(savdir,'anaFruit_fillop_dum_150505.mat'), ...
    'alpall', 'banfndall', 'banposall', 'bantimindall', 'avgbandispertrl', ...
    'cumdispertrl', 'avgwaldispertrl', 'trllngall', 'invprmindall', ...
    'prvbanalpall', 'invprmnrlsel', 'filind', 'numtrl', 'dum', '-v7.3')

% use to load pre-saved variables
savdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\misc';
load(fullfile(savdir,'anaFruit_fillop_dum_150505.mat'))


%%

% all neural data is contained in dum
%
% number of trials contained in dum is a multiple of 3: each group of 3
% trials corresponds to the 3 visible banana trials preceding an
% invisible banana
%
% only invisible bananas with at least 3 visible trials preceding the
% invisible one are included
%
% invprmindall: contains an index of all the trials with an invisible
% banana that was in a new location (i.e. not presented in the same
% location as the previous invisible banana) or that was not preceded
% immediately by an invisible banana (which happened occasionaly by
% accident)
%
% invprmnrlsel: contains only the invisible bananas that were preceded by
% at least 3 visible bananas in the same location
% these represent the trials that are selected to include in dum; the
% trials indexed in 'invprmnrlsel' are not the actual trials from which
% neural data was taken (neural data was taken from the 3 preceding trials)
%
% prvbanalpall: for each invisible banana included in the neural data, this
% contains the trial number and alpha level for all bananas preceding the
% invisible banana
% trial numbers refer to the number within each recording, not the overall 
% trial number, which can be calculated using invprmnrlsel
% alpha levels are usually either 0.15 or 1; he seems to be able to see
% 0.15 without trouble, so I'm lumping them together with the completely
% visible bananas for now)
% from this variable, you can also figure out the total number of visible
% bananas presented before each invisible banana selected for dum
%
% all other "group" variables include information pertaining to every trial



%% preprocess data

% exclude these channels
badchn = {'A02'; 'A05'; 'A08'; 'A10'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};
% badchn = {'A01'; 'A02'; 'A03'; 'A04'; 'A05'; 'A06'; 'A07'; 'A08'; 'A09'; 'A10'; 'A11'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};

% exclude eye and position data for now
badchn = [badchn; {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'}];

cfg=[];
cfg.continuous    = 'no';
cfg.channel       = setxor(dum.label,badchn);
data0 = ft_preprocessing(cfg,dum);
clear dum

% these data are still noisy; wait until 5-second clips are selected before
% filtering out noise


%% find trials with at least 5 seconds preceding banana

bantim = nan(length(invprmnrlsel)*3,1);
timsel = nan(length(invprmnrlsel)*3,2);
invprmnrltimsel = [];
trlgrp = [];
for invlop = 1:length(invprmnrlsel)
    bantim(invlop*3-2:invlop*3) = bantimindall(invprmnrlsel(invlop)-3:invprmnrlsel(invlop)-1);
    if isempty(find((bantim(invlop*3-2:invlop*3))<5000,1))
        timsel(invlop*3-2:invlop*3,:) = [bantim(invlop*3-2:invlop*3)-4999 bantim(invlop*3-2:invlop*3)];
        invprmnrltimsel = [invprmnrltimsel; invprmnrlsel(invlop)];
        trlgrp = [trlgrp; ones(3,1)+length(trlgrp)/3];
    else
        timsel(invlop*3-2:invlop*3,:) = repmat([nan nan],3,1);
    end
end


%% select last 5 seconds of approach path to banana in neural data

data1 = data0;
for trllop = 1:length(data0.trial)
    if isempty(find(isnan(timsel(trllop,:)),1,'first'))
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1];
        data1.time{trllop} = -4.999:0.001:0;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2));

%         % use this code when doing time-frequency analysis (padding at the end)
%         data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1+500];
%         data1.time{trllop} = -4.999:0.001:0.5;
%         data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2)+500);
    end
end

% remove marked trials
cfg=[];
% here choose only trials with at least 5 seconds pre-banana 
cfg.trials        = setxor(1:length(data1.trial),find(isnan(timsel(:,1))));
data1 = ft_preprocessing(cfg,data1);
clear data0

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
clear data1

% now run ft_rejectvisual and ft_databrowser to clean bad trials


%% mark bad trials (CHANGE IF ADDING NEW DATA)
% identified with ft_rejectvisual, ft_databrowser

badtrldum = [29, 35, 39, 41, 45, 66, 72, 107, 143, 151, 180, 183, 241, 283];

badtrl = [];
invselexc = [];
for badlop = 1:length(badtrldum)
    badtrl = [badtrl; find(trlgrp==trlgrp(badtrldum(badlop)))];
    invselexc = [invselexc; find(trlgrp==trlgrp(badtrldum(badlop)),1,'last')/3];
end

invprmnrlbadsel = invprmnrltimsel(setxor(1:length(invprmnrltimsel),invselexc));

% exclude bad trials, filter out line noise here
cfg=[];
cfg.trials = setxor(1:length(data2.time),badtrl);
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
data3 = ft_preprocessing(cfg,data2);
clear data2

%% spectral analysis, all trials

% do the spectral analysis - time-averaged
cfg=[];
cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
% cfg.output      = 'powandcsd';
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.foilim      = [1 30];
% cfg.taper       = 'hanning';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqL = ft_freqanalysis(cfg, data3);
% fdL = ft_freqdescriptives(cfg, freqL);

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqH = ft_freqanalysis(cfg, data3);
% fdH = ft_freqdescriptives(cfg, freqH);


cfg=[];
cfg.method      = 'coh';
cfg.jackknife   = 'yes';
statL = ft_connectivityanalysis(cfg,freqL);
statH = ft_connectivityanalysis(cfg,freqH);


%% categorize trials into good/bad memory

bantimnrm = bantimindall(invprmnrlbadsel)/max(bantimindall(invprmnrlbadsel));
bandstnrm = avgbandispertrl(invprmnrlbadsel)/max(avgbandispertrl(invprmnrlbadsel));
waldstnrm = avgwaldispertrl(invprmnrlbadsel)/max(avgwaldispertrl(invprmnrlbadsel));

% time to acquire banana, mean dist. from banana, mean dist. from wall
mempnthi = [0 0 max(waldstnrm)];
mempntlo = [max(bantimnrm) max(bandstnrm) 0];
memdst = sqrt(sum((abs(([bantimnrm bandstnrm waldstnrm])-repmat(mempnthi,length(invprmnrlbadsel),1))).^2,2));

% % use memory "distance"
% [~,memsrt] = sort(memdst);
% use excess path length
[~,memsrt] = sort(excpthind(invprmnrlbadsel));

% memindhi = memsrt(1:round(length(memsrt)/2));
% memindlo = memsrt(round(length(memsrt)/2)+1:end);
memindhi = memsrt(1:30);
memindlo = memsrt(end-29:end);


%% plot distribution of banana locations

x = linspace(-10,10,17);
[n,c] = hist3(banposall(invprmnrlbadsel,:),'Edges',{x,x});
figure;imagesc(c{1},c{2},n);xlim([-10 10]);ylim([-10 10]);colormap hot;colorbar
hold on;scatter(banposall(invprmnrlbadsel,2),banposall(invprmnrlbadsel,1),'bo')


%% heat-map plot memory performance measures

% calculate distance from banana at start position
xdif = begposall(invprmnrlbadsel,1)-banposall(invprmnrlbadsel,1);
ydif = begposall(invprmnrlbadsel,2)-banposall(invprmnrlbadsel,2);
begdis = sqrt(xdif.^2+ydif.^2);

x = linspace(-11.2,11.2,81);
ctrs = x(1:end-1)+mean(gradient(x))/2;

lat3d = []; % latency to get banana
dis3d = []; % average distance from banana
wal3d = []; % average distance from wall
mem3d = []; % memory performance
beg3d = []; % starting distance
epldis3d = []; % excess path length, difference
eplrat3d = []; % excess path length, ratio
eplind3d = []; % excess path length, index
memnrm3d = []; % memory normalized by starting distance
for k=1:length(memdst)
    
    latn = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*bantimindall(invprmnrlbadsel(k)),fspecial('gaussian',30,3));
    disn = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*avgbandispertrl(invprmnrlbadsel(k)),fspecial('gaussian',30,3));
    waln = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*avgwaldispertrl(invprmnrlbadsel(k)),fspecial('gaussian',30,3));
    memn = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*memdst(k),fspecial('gaussian',30,3));
    begn = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*begdis(k),fspecial('gaussian',30,3));
    epldis = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*excpthdif(invprmnrlbadsel(k)),fspecial('gaussian',30,3));
    eplrat = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*excpthnrm(invprmnrlbadsel(k)),fspecial('gaussian',30,3));
    eplind = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*excpthind(invprmnrlbadsel(k)),fspecial('gaussian',30,3));
    memnrmn = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x})*(memdst(k)/begdis(k)),fspecial('gaussian',30,3));
    
    nd = imfilter(hist3(banposall(invprmnrlbadsel(k),:),'Edges',{x,x}),fspecial('disk',15));
    latn(nd==0) = nan;
    disn(nd==0) = nan;
    waln(nd==0) = nan;
    memn(nd==0) = nan;
    begn(nd==0) = nan;
    epldis(nd==0) = nan;
    eplrat(nd==0) = nan;
    eplind(nd==0) = nan;
    memnrmn(nd==0) = nan;
    
    lat3d(k,:,:) = latn(1:end-1,1:end-1);
    dis3d(k,:,:) = disn(1:end-1,1:end-1);
    wal3d(k,:,:) = waln(1:end-1,1:end-1);
    mem3d(k,:,:) = memn(1:end-1,1:end-1);
    beg3d(k,:,:) = begn(1:end-1,1:end-1);
    epldis3d(k,:,:) = epldis(1:end-1,1:end-1);
    eplrat3d(k,:,:) = eplrat(1:end-1,1:end-1);
    eplind3d(k,:,:) = eplind(1:end-1,1:end-1);
    memnrm3d(k,:,:) = memnrmn(1:end-1,1:end-1);
    
%     figure
%     subplot(2,1,1)
%     imagesc(ctrs,ctrs,n(1:end-1,1:end-1));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet
%     subplot(2,1,2)
%     imagesc(ctrs,ctrs,squeeze(nanmean(mem3d,1)));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet
% %     set(gca,'clim',[0 0.002])
%     set(gcf,'Position',[80 67 496 911])
%     pause
%     close
    
end

lat2d = squeeze(nanmean(lat3d,1));
dis2d = squeeze(nanmean(dis3d,1));
wal2d = squeeze(nanmean(wal3d,1));
mem2d = squeeze(nanmean(mem3d,1));
beg2d = squeeze(nanmean(beg3d,1));
epldis2d = squeeze(nanmean(epldis3d,1));
eplrat2d = squeeze(nanmean(eplrat3d,1));
eplind2d = squeeze(nanmean(eplind3d,1));
memnrm2d = squeeze(nanmean(memnrm3d,1));

figure;imagesc(ctrs,ctrs,lat2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Latency to get banana by banana location');colorbar
figure;imagesc(ctrs,ctrs,dis2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Average distance from banana by banana location');colorbar
figure;imagesc(ctrs,ctrs,wal2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Average distance from wall by banana location');colorbar
figure;imagesc(ctrs,ctrs,mem2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Memory performance by banana location');colorbar
figure;imagesc(ctrs,ctrs,beg2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Starting distance by banana location');colorbar
figure;imagesc(ctrs,ctrs,epldis2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (difference) by banana location');colorbar
figure;imagesc(ctrs,ctrs,eplrat2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (ratio) by banana location');colorbar
figure;imagesc(ctrs,ctrs,eplind2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (index) by banana location');colorbar
figure;imagesc(ctrs,ctrs,memnrm2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Normalized memory by banana location');colorbar


%% using memindhi and memindlo, categorize neural data trials

nrlmemindhi = sort([memindhi*3-2; memindhi*3-1; memindhi*3]);
nrlmemindlo = sort([memindlo*3-2; memindlo*3-1; memindlo*3]);

cfg=[];
cfg.method      = 'coh';
cfg.jackknife   = 'yes';

cfg.trials = nrlmemindhi;
statLh = ft_connectivityanalysis(cfg,freqL);
statHh = ft_connectivityanalysis(cfg,freqH);
fdLh = ft_freqdescriptives(cfg, freqL);
fdHh = ft_freqdescriptives(cfg, freqH);

cfg.trials = nrlmemindlo;
statLl = ft_connectivityanalysis(cfg,freqL);
statHl = ft_connectivityanalysis(cfg,freqH);
fdLl = ft_freqdescriptives(cfg, freqL);
fdHl = ft_freqdescriptives(cfg, freqH);


%% plot single pair

ab1=find(strcmp(statLh.label,'A06')); ab2=find(strcmp(statLh.label,'C01'));
figure;hold on
plot(statLh.freq,squeeze(statLh.cohspctrm(ab1,ab2,:)),'r');
plot(statLl.freq,squeeze(statLl.cohspctrm(ab1,ab2,:)),'b');
plot(statHh.freq,squeeze(statHh.cohspctrm(ab1,ab2,:)),'r');
plot(statHl.freq,squeeze(statHl.cohspctrm(ab1,ab2,:)),'b');
% plot(statLh.freq,squeeze(statLh.cohspctrm(ab1,ab2,:))+squeeze(statLh.cohspctrmsem(ab1,ab2,:)),'r','LineStyle','--');
% plot(statLh.freq,squeeze(statLh.cohspctrm(ab1,ab2,:))-squeeze(statLh.cohspctrmsem(ab1,ab2,:)),'r','LineStyle','--');
% plot(statLh.freq,squeeze(statLl.cohspctrm(ab1,ab2,:))+squeeze(statLl.cohspctrmsem(ab1,ab2,:)),'b','LineStyle','--');
% plot(statLh.freq,squeeze(statLl.cohspctrm(ab1,ab2,:))-squeeze(statLl.cohspctrmsem(ab1,ab2,:)),'b','LineStyle','--');
axis tight
yh=ylim;

% shade standard error instead of dashed lines
% low frequency, high memory
sempos=squeeze(statLh.cohspctrm(ab1,ab2,:))+squeeze(statLh.cohspctrmsem(ab1,ab2,:));
semneg=squeeze(statLh.cohspctrm(ab1,ab2,:))-squeeze(statLh.cohspctrmsem(ab1,ab2,:));
dum=statLh.freq;dum1=fliplr(-statLh.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'r', 'FaceAlpha',.2, 'EdgeColor','r','EdgeAlpha',0);
% low frequency, low memory
sempos=squeeze(statLl.cohspctrm(ab1,ab2,:))+squeeze(statLl.cohspctrmsem(ab1,ab2,:));
semneg=squeeze(statLl.cohspctrm(ab1,ab2,:))-squeeze(statLl.cohspctrmsem(ab1,ab2,:));
dum=statLl.freq;dum1=fliplr(-statLl.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'b', 'FaceAlpha',.2, 'EdgeColor','b','EdgeAlpha',0);
% high frequency, high memory
sempos=squeeze(statHh.cohspctrm(ab1,ab2,:))+squeeze(statHh.cohspctrmsem(ab1,ab2,:));
semneg=squeeze(statHh.cohspctrm(ab1,ab2,:))-squeeze(statHh.cohspctrmsem(ab1,ab2,:));
dum=statHh.freq;dum1=fliplr(-statHh.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'r', 'FaceAlpha',.2, 'EdgeColor','r','EdgeAlpha',0);
% low frequency, low memory
sempos=squeeze(statHl.cohspctrm(ab1,ab2,:))+squeeze(statHl.cohspctrmsem(ab1,ab2,:));
semneg=squeeze(statHl.cohspctrm(ab1,ab2,:))-squeeze(statHl.cohspctrmsem(ab1,ab2,:));
dum=statHl.freq;dum1=fliplr(-statHl.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'b', 'FaceAlpha',.2, 'EdgeColor','b','EdgeAlpha',0);

% for frqlop=1:length(statLh.freq)
%     if (squeeze(statLh.cohspctrm(ab1,ab2,frqlop))-squeeze(statLh.cohspctrmsem(ab1,ab2,frqlop)))-(squeeze(statLl.cohspctrm(ab1,ab2,frqlop))+squeeze(statLl.cohspctrmsem(ab1,ab2,frqlop)))>0
%         fill([statLh.freq(frqlop)-0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)-0.1],[yh(1) yh(1) yh(2) yh(2)],'r','FaceAlpha',0.2,'EdgeColor','none');
%     elseif (squeeze(statLl.cohspctrm(ab1,ab2,frqlop))-squeeze(statLl.cohspctrmsem(ab1,ab2,frqlop)))-(squeeze(statLh.cohspctrm(ab1,ab2,frqlop))+squeeze(statLh.cohspctrmsem(ab1,ab2,frqlop)))>0
%         fill([statLh.freq(frqlop)-0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)-0.1],[yh(1) yh(1) yh(2) yh(2)],'b','FaceAlpha',0.2,'EdgeColor','none');
%     end
% end

title([statLh.label{ab1} ' x' statLh.label{ab2}]);
set(gca,'TickDir','out')
legend({'High Memory';'Low Memory'})
xlabel('Frequency (Hz)');ylabel('Coherence')


%% plot average intra- and inter-array coherence

clear prblaball
for k=1:length(statLh.label)
    prblaball(k,1)=statLh.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);

cmblst = unique(sort([repmat(unique(locb),length(unique(locb)),1) sort(repmat(unique(locb),length(unique(locb)),1))],2),'rows');

cohmatLh = [];
cohmatLl = [];
cohmatHh = [];
cohmatHl = [];
cohdifLh = zeros(size(cmblst,1),length(statLh.freq));
cohdifLl = zeros(size(cmblst,1),length(statLl.freq));
cohdifHh = zeros(size(cmblst,1),length(statHh.freq));
cohdifHl = zeros(size(cmblst,1),length(statHl.freq));
numprs = [];
for cmblop = 1:size(cmblst,1)
    
    cmbind1 = find(locb==cmblst(cmblop,1));
    cmbind2 = find(locb==cmblst(cmblop,2));
    cmblst2 = unique(sort([repmat(unique(cmbind1),length(unique(cmbind2)),1) sort(repmat(unique(cmbind2),length(unique(cmbind1)),1))],2),'rows');
    cmblst2 = cmblst2(cmblst2(:,1)~=cmblst2(:,2),:);
    
    numprs(cmblop) = size(cmblst2,1);
    
    cohmat1 = []; cohmat2 = []; cohmat3 = []; cohmat4 = [];
    for indlop1 = 1:size(cmblst2,1)
        
        cohmat1 = cat(1,cohmat1,(squeeze(statLh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        cohmat2 = cat(1,cohmat2,(squeeze(statLl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        cohmat3 = cat(1,cohmat3,(squeeze(statHh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        cohmat4 = cat(1,cohmat4,(squeeze(statHl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        
        for frqlop=1:length(statLh.freq)
            if (squeeze(statLh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))-squeeze(statLh.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))-(squeeze(statLl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))+squeeze(statLl.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))>0
                cohdifLh(cmblop,frqlop) = cohdifLh(cmblop,frqlop)+1;
            elseif (squeeze(statLl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))-squeeze(statLl.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))-(squeeze(statLh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))+squeeze(statLh.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))>0
                cohdifLl(cmblop,frqlop) = cohdifLl(cmblop,frqlop)+1;
            end
        end

        for frqlop=1:length(statHh.freq)
            if (squeeze(statHh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))-squeeze(statHh.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))-(squeeze(statHl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))+squeeze(statHl.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))>0
                cohdifHh(cmblop,frqlop) = cohdifHh(cmblop,frqlop)+1;
            elseif (squeeze(statHl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))-squeeze(statHl.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))-(squeeze(statHh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop))+squeeze(statHh.cohspctrmsem(cmblst2(indlop1,1),cmblst2(indlop1,2),frqlop)))>0
                cohdifHl(cmblop,frqlop) = cohdifHl(cmblop,frqlop)+1;
            end
        end
        
    end
    
    cohmatLh = cat(1,cohmatLh,mean(cohmat1,1));
    cohmatLl = cat(1,cohmatLl,mean(cohmat2,1));
    cohmatHh = cat(1,cohmatHh,mean(cohmat3,1));
    cohmatHl = cat(1,cohmatHl,mean(cohmat4,1));
    
end

% for cmblop = 1:size(cmblst,1)
%     figure; hold on
%     plot(statLh.freq,cohmatLh(cmblop,:),'r')
%     plot(statLh.freq,cohmatLl(cmblop,:),'b')
%     plot(statHh.freq,cohmatHh(cmblop,:),'r')
%     plot(statHh.freq,cohmatHl(cmblop,:),'b')
%     xlabel('Frequency (Hz)')
%     ylabel('Coherence')
%     title([prblab(cmblst(cmblop,1)) ' x ' prblab(cmblst(cmblop,2))])
% end

% figure;plot(statLh.freq,sum(cohdifLh,1)/sum(numprs),'r')
% hold on;plot(statLl.freq,sum(cohdifLl,1)/sum(numprs),'b')
% hold on;plot(statHh.freq,sum(cohdifHh,1)/sum(numprs),'r')
% hold on;plot(statHl.freq,sum(cohdifHl,1)/sum(numprs),'b')
% xlabel('Frequency (Hz)')
% ylabel('P of pairs showing sig. diff.')
% title('All array combinations')

% for cmblop = 1:size(cohdifLh,1)
%     figure; hold on
%     plot(statLh.freq,cohdifLh(cmblop,:)/numprs(cmblop),'r')
%     plot(statLh.freq,cohdifLl(cmblop,:)/numprs(cmblop),'b')
%     plot(statHh.freq,cohdifHh(cmblop,:)/numprs(cmblop),'r')
%     plot(statHh.freq,cohdifHl(cmblop,:)/numprs(cmblop),'b')
%     xlabel('Frequency (Hz)')
%     ylabel('P of pairs showing sig. diff.')
%     title([prblab(cmblst(cmblop,1)) ' x ' prblab(cmblst(cmblop,2))])
% end

figure
bar([statLh.freq statHh.freq(2:end)],([cohdifLh cohdifHh(:,2:end)]/sum(numprs))','stacked')
xlim([0 100])
legend({'A x A','A x B','A x C','B x B','B x C','C x C'})
set(gca,'TickDir','out')
title('High Memory')
xlabel('Frequency (Hz)')
ylabel('P pairs showing sig. diff.')

figure
bar([statLl.freq statHl.freq(2:end)],([cohdifLl cohdifHl(:,2:end)]/sum(numprs))','stacked')
xlim([0 100])
legend({'A x A','A x B','A x C','B x B','B x C','C x C'})
set(gca,'TickDir','out')
title('Low Memory')
xlabel('Frequency (Hz)')
ylabel('P pairs showing sig. diff.')


%% High vs. Low memory, scatter, each pair (color coded)
% create a different plot for each frequency/range

clear prblaball
for k=1:length(statLh.label)
    prblaball(k,1)=statLh.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);

cmblst = unique(sort([repmat(unique(locb),length(unique(locb)),1) sort(repmat(unique(locb),length(unique(locb)),1))],2),'rows');

[~,f1_t] = min(abs(statLh.freq-3));
[~,f2_t] = min(abs(statLh.freq-12));
[~,f1_lg] = min(abs(statHh.freq-30));
[~,f2_lg] = min(abs(statHh.freq-60));
[~,f1_hg] = min(abs(statHh.freq-60));
[~,f2_hg] = min(abs(statHh.freq-100));

cohmemdif_theta = cell(size(cmblst,1),1); % 3-12 Hz
cohmemdif_logamma = cell(size(cmblst,1),1); % 30-60 Hz
cohmemdif_higamma = cell(size(cmblst,1),1); % 60-100 Hz
% numprs = [];
for cmblop = 1:size(cmblst,1)
    
    cmbind1 = find(locb==cmblst(cmblop,1));
    cmbind2 = find(locb==cmblst(cmblop,2));
    cmblst2 = unique(sort([repmat(unique(cmbind1),length(unique(cmbind2)),1) sort(repmat(unique(cmbind2),length(unique(cmbind1)),1))],2),'rows');
    cmblst2 = cmblst2(cmblst2(:,1)~=cmblst2(:,2),:);
    
%     numprs(cmblop) = size(cmblst2,1);
    
    cohmat1 = []; cohmat2 = []; cohmat3 = []; cohmat4 = [];
    for indlop1 = 1:size(cmblst2,1)
        
        cohmat1 = cat(1,cohmat1,(squeeze(statLh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        cohmat2 = cat(1,cohmat2,(squeeze(statLl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        cohmat3 = cat(1,cohmat3,(squeeze(statHh.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        cohmat4 = cat(1,cohmat4,(squeeze(statHl.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        
    end
        
    cohmemdif_theta{cmblop} = [mean(cohmat1(:,f1_t:f2_t),2) mean(cohmat2(:,f1_t:f2_t),2)];
    cohmemdif_logamma{cmblop} = [mean(cohmat3(:,f1_lg:f2_lg),2) mean(cohmat4(:,f1_lg:f2_lg),2)];
    cohmemdif_higamma{cmblop} = [mean(cohmat3(:,f1_hg:f2_hg),2) mean(cohmat4(:,f1_hg:f2_hg),2)];
    
end

% theta
figure; hold on
% h1=scatter(cohmemdif_theta{1}(:,1),cohmemdif_theta{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_theta{2}(:,1),cohmemdif_theta{2}(:,2)); set(h2,'MarkerFaceColor','flat')
h3=scatter(cohmemdif_theta{3}(:,1),cohmemdif_theta{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% h4=scatter(cohmemdif_theta{4}(:,1),cohmemdif_theta{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_theta{5}(:,1),cohmemdif_theta{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% h6=scatter(cohmemdif_theta{6}(:,1),cohmemdif_theta{6}(:,2)); set(h6,'MarkerFaceColor','flat')
xh=xlim; yh=ylim;
line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
legend({'A x C'},'Location','southeast')
set(gca,'TickDir','out')
title('Theta (3-12 Hz)')
xlabel('Coherence: High Memory')
ylabel('Coherence: Low Memory')

% low gamma
figure; hold on
% h1=scatter(cohmemdif_logamma{1}(:,1),cohmemdif_logamma{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_logamma{2}(:,1),cohmemdif_logamma{2}(:,2)); set(h2,'MarkerFaceColor','flat')
h3=scatter(cohmemdif_logamma{3}(:,1),cohmemdif_logamma{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% h4=scatter(cohmemdif_logamma{4}(:,1),cohmemdif_logamma{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_logamma{5}(:,1),cohmemdif_logamma{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% h6=scatter(cohmemdif_logamma{6}(:,1),cohmemdif_logamma{6}(:,2)); set(h6,'MarkerFaceColor','flat')
xh=xlim; yh=ylim;
line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
legend({'A x C'},'Location','southeast')
set(gca,'TickDir','out')
title('Low Gamma (30-60 Hz)')
xlabel('Coherence: High Memory')
ylabel('Coherence: Low Memory')

% high gamma
figure; hold on
% h1=scatter(cohmemdif_higamma{1}(:,1),cohmemdif_higamma{1}(:,2)); set(h1,'MarkerFaceColor','flat')
% h2=scatter(cohmemdif_higamma{2}(:,1),cohmemdif_higamma{2}(:,2)); set(h2,'MarkerFaceColor','flat')
h3=scatter(cohmemdif_higamma{3}(:,1),cohmemdif_higamma{3}(:,2)); set(h3,'MarkerFaceColor','flat')
% h4=scatter(cohmemdif_higamma{4}(:,1),cohmemdif_higamma{4}(:,2)); set(h4,'MarkerFaceColor','flat')
% h5=scatter(cohmemdif_higamma{5}(:,1),cohmemdif_higamma{5}(:,2)); set(h5,'MarkerFaceColor','flat')
% h6=scatter(cohmemdif_higamma{6}(:,1),cohmemdif_higamma{6}(:,2)); set(h6,'MarkerFaceColor','flat')
xh=xlim; yh=ylim;
line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
% legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','southeast')
% legend({'A x B','A x C','B x C'},'Location','southeast')
legend({'A x C'},'Location','southeast')
set(gca,'TickDir','out')
title('High Gamma (60-100 Hz)')
xlabel('Coherence: High Memory')
ylabel('Coherence: Low Memory')

% do histograms instead of scatter plots
edg = -0.05:0.005:0.05;
edx = edg + 0.0025;
[n1,b1] = histc(-diff(cohmemdif_theta{1},1,2),edg);
[n2,b2] = histc(-diff(cohmemdif_theta{2},1,2),edg);
[n3,b3] = histc(-diff(cohmemdif_theta{3},1,2),edg);
[n4,b4] = histc(-diff(cohmemdif_theta{4},1,2),edg);
[n5,b5] = histc(-diff(cohmemdif_theta{5},1,2),edg);
[n6,b6] = histc(-diff(cohmemdif_theta{6},1,2),edg);
figure;hold on
bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
line([0 0],ylim,'Color','k','LineStyle','--')
legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
set(gca,'TickDir','out')
box off
title('Theta (3-12 Hz)')
xlabel('Coherence difference (High - Low memory)')
ylabel('# of pairs')

edg = -0.05:0.005:0.05;
edx = edg + 0.0025;
[n1,b1] = histc(-diff(cohmemdif_logamma{1},1,2),edg);
[n2,b2] = histc(-diff(cohmemdif_logamma{2},1,2),edg);
[n3,b3] = histc(-diff(cohmemdif_logamma{3},1,2),edg);
[n4,b4] = histc(-diff(cohmemdif_logamma{4},1,2),edg);
[n5,b5] = histc(-diff(cohmemdif_logamma{5},1,2),edg);
[n6,b6] = histc(-diff(cohmemdif_logamma{6},1,2),edg);
figure;hold on
bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
line([0 0],ylim,'Color','k','LineStyle','--')
legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
set(gca,'TickDir','out')
box off
title('Low Gamma (30-60 Hz)')
xlabel('Coherence difference (High - Low memory)')
ylabel('# of pairs')

edg = -0.05:0.005:0.05;
edx = edg + 0.0025;
[n1,b1] = histc(-diff(cohmemdif_higamma{1},1,2),edg);
[n2,b2] = histc(-diff(cohmemdif_higamma{2},1,2),edg);
[n3,b3] = histc(-diff(cohmemdif_higamma{3},1,2),edg);
[n4,b4] = histc(-diff(cohmemdif_higamma{4},1,2),edg);
[n5,b5] = histc(-diff(cohmemdif_higamma{5},1,2),edg);
[n6,b6] = histc(-diff(cohmemdif_higamma{6},1,2),edg);
figure;hold on
bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
line([0 0],ylim,'Color','k','LineStyle','--')
legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
set(gca,'TickDir','out')
box off
title('High Gamma (60-100 Hz)')
xlabel('Coherence difference (High - Low memory)')
ylabel('# of pairs')


%% power, scatter plot

% create a different plot for each frequency/range
clear prblaball
for k=1:length(fdLh.label)
    prblaball(k,1)=fdLh.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);

[~,f1_t] = min(abs(fdLh.freq-3));
[~,f2_t] = min(abs(fdLh.freq-12));
[~,f1_lg] = min(abs(fdHh.freq-30));
[~,f2_lg] = min(abs(fdHh.freq-60));
[~,f1_hg] = min(abs(fdHh.freq-60));
[~,f2_hg] = min(abs(fdHh.freq-100));

powmemdif_theta = nan(size(fdLh.powspctrm,1),2); % 3-12 Hz
powmemdif_logamma = nan(size(fdLh.powspctrm,1),2); % 30-60 Hz
powmemdif_higamma = nan(size(fdLh.powspctrm,1),2); % 60-100 Hz
powmat1 = []; powmat2 = []; powmat3 = []; powmat4 = [];
% numprs = [];
for indlop = 1:size(fdLh.powspctrm,1)
    
    
    powmat1 = cat(1,powmat1,squeeze(fdLh.powspctrm(indlop,:)));
    powmat2 = cat(1,powmat2,squeeze(fdLl.powspctrm(indlop,:)));
    powmat3 = cat(1,powmat3,squeeze(fdHh.powspctrm(indlop,:)));
    powmat4 = cat(1,powmat4,squeeze(fdHl.powspctrm(indlop,:)));
        
    powmemdif_theta(indlop,:) = [mean(powmat1(indlop,f1_t:f2_t),2) mean(powmat2(indlop,f1_t:f2_t),2)];
    powmemdif_logamma(indlop,:) = [mean(powmat3(indlop,f1_lg:f2_lg),2) mean(powmat4(indlop,f1_lg:f2_lg),2)];
    powmemdif_higamma(indlop,:) = [mean(powmat3(indlop,f1_hg:f2_hg),2) mean(powmat4(indlop,f1_hg:f2_hg),2)];
    
end

% theta
figure; hold on
h1=scatter(powmemdif_theta(locb==1,1),powmemdif_theta(locb==1,2)); set(h1,'MarkerFaceColor','flat')
h2=scatter(powmemdif_theta(locb==2,1),powmemdif_theta(locb==2,2)); set(h2,'MarkerFaceColor','flat')
h3=scatter(powmemdif_theta(locb==3,1),powmemdif_theta(locb==3,2)); set(h3,'MarkerFaceColor','flat')
xh=xlim; yh=ylim;
line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
legend({'A' 'B' 'C'},'Location','southeast')
set(gca,'TickDir','out')
title('Theta (3-12 Hz)')
xlabel('Power: High Memory')
ylabel('Power: Low Memory')

% low gamma
figure; hold on
h1=scatter(powmemdif_logamma(locb==1,1),powmemdif_logamma(locb==1,2)); set(h1,'MarkerFaceColor','flat')
h2=scatter(powmemdif_logamma(locb==2,1),powmemdif_logamma(locb==2,2)); set(h2,'MarkerFaceColor','flat')
h3=scatter(powmemdif_logamma(locb==3,1),powmemdif_logamma(locb==3,2)); set(h3,'MarkerFaceColor','flat')
xh=xlim; yh=ylim;
line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
legend({'A' 'B' 'C'},'Location','southeast')
set(gca,'TickDir','out')
title('Low Gamma (30-60 Hz)')
xlabel('Power: High Memory')
ylabel('Power: Low Memory')

% high gamma
figure; hold on
h1=scatter(powmemdif_higamma(locb==1,1),powmemdif_higamma(locb==1,2)); set(h1,'MarkerFaceColor','flat')
h2=scatter(powmemdif_higamma(locb==2,1),powmemdif_higamma(locb==2,2)); set(h2,'MarkerFaceColor','flat')
h3=scatter(powmemdif_higamma(locb==3,1),powmemdif_higamma(locb==3,2)); set(h3,'MarkerFaceColor','flat')
xh=xlim; yh=ylim;
line([0 max([xh(2) yh(2)])],[0 max([xh(2) yh(2)])],'Color','k','LineStyle','--')
legend({'A' 'B' 'C'},'Location','southeast')
set(gca,'TickDir','out')
title('High Gamma (60-100 Hz)')
xlabel('Power: High Memory')
ylabel('Power: Low Memory')


%% power single pair, power

ab1=find(strcmp(fdLh.label,'A06'));
figure;hold on
plot(fdLh.freq,fdLh.powspctrm(ab1,:),'r');
plot(fdLl.freq,fdLl.powspctrm(ab1,:),'b');
plot(fdHh.freq,fdHh.powspctrm(ab1,:),'r');
plot(fdHl.freq,fdHl.powspctrm(ab1,:),'b');
axis tight
yh=ylim;

% shade standard error instead of dashed lines
% low frequency, high memory
sempos=squeeze(fdLh.powspctrm(ab1,:))+squeeze(fdLh.powspctrmsem(ab1,:));
semneg=squeeze(fdLh.powspctrm(ab1,:))-squeeze(fdLh.powspctrmsem(ab1,:));
dum=fdLh.freq;dum1=fliplr(-fdLh.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'r', 'FaceAlpha',.2, 'EdgeColor','r','EdgeAlpha',0);
% low frequency, low memory
sempos=squeeze(fdLl.powspctrm(ab1,:))+squeeze(fdLl.powspctrmsem(ab1,:));
semneg=squeeze(fdLl.powspctrm(ab1,:))-squeeze(fdLl.powspctrmsem(ab1,:));
dum=fdLl.freq;dum1=fliplr(-fdLl.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'b', 'FaceAlpha',.2, 'EdgeColor','b','EdgeAlpha',0);
% high frequency, high memory
sempos=squeeze(fdHh.powspctrm(ab1,:))+squeeze(fdHh.powspctrmsem(ab1,:));
semneg=squeeze(fdHh.powspctrm(ab1,:))-squeeze(fdHh.powspctrmsem(ab1,:));
dum=fdHh.freq;dum1=fliplr(-fdHh.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'r', 'FaceAlpha',.2, 'EdgeColor','r','EdgeAlpha',0);
% low frequency, low memory
sempos=squeeze(fdHl.powspctrm(ab1,:))+squeeze(fdHl.powspctrmsem(ab1,:));
semneg=squeeze(fdHl.powspctrm(ab1,:))-squeeze(fdHl.powspctrmsem(ab1,:));
dum=fdHl.freq;dum1=fliplr(-fdHl.freq);dum1=dum1*-1;x=[dum dum1];
dum=flipud(semneg(:));
y=[sempos(:);dum]';
hold on;
fill(x,y,'b', 'FaceAlpha',.2, 'EdgeColor','b','EdgeAlpha',0);
set(gca,'yScale','log','xScale','log')

% for frqlop=1:length(statLh.freq)
%     if (squeeze(statLh.cohspctrm(ab1,ab2,frqlop))-squeeze(statLh.cohspctrmsem(ab1,ab2,frqlop)))-(squeeze(statLl.cohspctrm(ab1,ab2,frqlop))+squeeze(statLl.cohspctrmsem(ab1,ab2,frqlop)))>0
%         fill([statLh.freq(frqlop)-0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)-0.1],[yh(1) yh(1) yh(2) yh(2)],'r','FaceAlpha',0.2,'EdgeColor','none');
%     elseif (squeeze(statLl.cohspctrm(ab1,ab2,frqlop))-squeeze(statLl.cohspctrmsem(ab1,ab2,frqlop)))-(squeeze(statLh.cohspctrm(ab1,ab2,frqlop))+squeeze(statLh.cohspctrmsem(ab1,ab2,frqlop)))>0
%         fill([statLh.freq(frqlop)-0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)+0.1 statLh.freq(frqlop)-0.1],[yh(1) yh(1) yh(2) yh(2)],'b','FaceAlpha',0.2,'EdgeColor','none');
%     end
% end

title([fdLh.label{ab1}]);
set(gca,'TickDir','out')
set(gca,'yScale','log')
legend({'High Memory';'Low Memory'})
xlabel('Frequency (Hz)');ylabel('Power')


%% stats (doesn't work yet)

cfg = [];
% dum = ones(1, size(statLh));
% cfg.design       =  [dum 2.*dum; 1:length(cnd1) 1:length(cnd1)];
dum = ones(1, 406);
cfg.design       =  [dum 2.*dum; 1:406 1:406];
cfg.method       =  'montecarlo';
cfg.statistic    =  'indepsamplesT';
cfg.parameter    =  'cohspctrm';
cfg.numrandomization  = 1000; 
% cfg.uvar         =  2;

cfg.channel      = 'sig001a';
cfg.avgoverchan  =  'no';
cfg.avgoverfreq  =  'no';

% cfg.latency      =  [0.05 0.2];
cfg.latency      =  'all';
cfg.frequency    =  'all';
% cfg.frequency    =  [phasebins(16) phasebins(40)];
% cfg.parameter    =  'powspctrm';

cfg.correctm     =  'cluster';
cfg.alpha        =  0.05;
cfg.tail         =  1;
cfg.ivar         =  1;
cfg.feedback     =  'text';
cfg.clusterstatistic  = 'maxsum';
% cfg.clusterthreshold  = 'parametric';
cfg.clusterthreshold  = 'nonparametric_common';
cfg.clusteralpha =  0.1;
% clustercritval original value: 1.96
% cfg.clustercritval    = 1.96;
cfg.clustertail  =  1;
% cfg.neighbours{1}.label = 'sig001a';
% cfg.neighbours{1}.neighblabel = {};

stat = ft_freqstatistics(cfg, statLh, statLl);


%% try time-resolved analysis

% do the spectral analysis - time-averaged
cfg=[];
cfg.output      = 'fourier';
cfg.method      = 'mtmconvol';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'yes';
cfg.foi         = 1:30;
cfg.t_ftimwin   = 0.5 * ones(size(cfg.foi));
cfg.toi         = -4.75:0.01:-0.25;
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 2;
freqTFL = ft_freqanalysis(cfg, data3);
fdTFL = ft_freqdescriptives(cfg, freqTFL);

% cfg.foi         = 30:2:100;
% cfg.t_ftimwin   = 0.25 * ones(size(cfg.foi));
% cfg.taper       = 'dpss';
% cfg.tapsmofrq   = 8;
% freqTFH = ft_freqanalysis(cfg, data3);
% fdTFH = ft_freqdescriptives(cfg, freqTFH);

cfg=[];
cfg.method      = 'coh';

cfg.trials = nrlmemindhi;
cfg.channelcmb = [freqTFL.label(4) freqTFL.label(19)];
statTFLh = ft_connectivityanalysis(cfg,freqTFL);
% statH_himem = ft_connectivityanalysis(cfg,freqH);

cfg.trials = nrlmemindlo;
statTFLl = ft_connectivityanalysis(cfg,freqTFL);
% statH_lomem = ft_connectivityanalysis(cfg,freqH);



