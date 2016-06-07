% 'JN_BR_15_04_10\JN_BR_15_04_10_13_03' % 4/10/15, part 1
% 'JN_BR_15_04_10\JN_BR_15_04_10_13_11' % 4/10/15, part 2
% 'JN_BR_15_04_13\JN_BR_15_04_13_13_59' % 4/13/15, no log entry

sessionDir = 'JN_BR_15_04_29\JN_BR_15_04_29_13_22';
trlDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\trial data';
NSDir = 'R:\Buffalo Lab\Virtual Navigation\MATLAB\MAT files\NSdat';
% savDir = 'C:\Data\VR';

[~,sesnam]=fileparts(sessionDir);
load(fullfile(trlDir,[sesnam '_trldat.mat']))
load(fullfile(NSDir,[sesnam '_NSdat.mat']))


%% pre-process the neural data in Fieldtrip

% exclude these channels
badchn = {'A02'; 'A05'; 'A08'; 'A10'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};
% exclude eye and position data for now
badchn = [badchn; {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'}];

cfg=[];
cfg.continuous    = 'no';
cfg.channel       = setxor(data.label,badchn);
data0 = ft_preprocessing(cfg,data);

% these data are still noisy; wait until 5-second clips are selected before
% filtering out noise


%% get "banana eat" events for all trials

alp = [];
bantim = [];
for k=1:length(trldat.alpha)
    alp(k,1) = trldat.alpha{k}(1,2);
    if trldat.frttim{k}(1,2)==0
        bantim(k,1) = trldat.frttim{k}(1,1)-trldat.time{k}(1);
    end
end


%% select last 5 seconds of approach path to banana in neural data

timsel = [bantim-4999 bantim];
timsel(timsel(:,1)<0,:) = nan;

data1 = data0;
for trllop = 1:length(data0.trial)
    if isempty(find(isnan(timsel(trllop,:)),1,'first'))
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1];
        data1.time{trllop} = -4.999:0.001:0;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2));
    end
end

% now filter out noise
cfg=[];
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
cfg.trials        = setxor(1:length(data1.trial),find(isnan(timsel(:,1))));
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
freqL = ft_freqanalysis(cfg, data2);
fdL = ft_freqdescriptives(cfg, freqL);

cfg.foilim      = [30 100];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 8;
freqH = ft_freqanalysis(cfg, data2);
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
cfg.trials      = setxor(1:length(data2.trial),find(isnan(timsel(:,1))));
cfg.taper       = 'hanning';
cfg.foilim      = [40 130];
cfg.channel     = data.label(find(strcmp(data.label,'C01')):find(strcmp(data.label,'C12')));
cfg.trials      = 1:5;
freqT1 = ft_freqanalysis(cfg, data);
freqT2 = ft_freqanalysis(cfg, data0);
freqT3 = ft_freqanalysis(cfg, data1);
freqT4 = ft_freqanalysis(cfg, data2);
figure;hold on
% plot(freqT1.freq,squeeze(mean(freqT1.powspctrm,1)),'b')
% plot(freqT2.freq,squeeze(mean(freqT2.powspctrm,1)),'r')
plot(freqT3.freq,squeeze(mean(freqT3.powspctrm,1)),'g')
plot(freqT4.freq,squeeze(mean(freqT4.powspctrm,1)),'m')


%% remove excluded trials from trldat

% rather than do this now, do the behavioral analysis as normal, then
% account for the removed trials when transferring results to neural data

% trldatfld = fieldnames(trldat);
% 
% for fldlop = 1:length(trldatfld)
%     trldat = setfield(trldat,trldatfld{fldlop},eval(['trldat.' trldatfld{fldlop} '(setxor(1:length(data0.trial),find(isnan(timsel(:,1)))))']));
% end


%% find banana positions on invisible trials

invind = find(alp==0);

banposinv = [];
repind = [];
for k=1:length(invind)
    banposinv = [banposinv; trldat.frtpos{invind(k)}(find(trldat.frtpos{invind(k)}(:,2)==0,1,'first'),3:4)];
end

figure;hold on
scatter(banposinv(:,1),banposinv(:,2))


%% make a 2-d histogram of location over time

grdedg = -11.2:2.8:11.2;

grdall = [];
for trllop = 1:length(trldat.posdat)
    [nx,binx]=histc(trldat.posdat{trllop}(1,:),grdedg);
    [ny,biny]=histc(trldat.posdat{trllop}(2,:),grdedg);
    grd = zeros(length(grdedg)-1,length(grdedg)-1);
    for l=1:length(binx)
        grd(binx(l),biny(l))=grd(binx(l),biny(l))+1;
    end
    grdall(trllop,:,:) = grd;
end

    
%% get average distance from banana for each trial

disall = [];
avgdispertrl = [];
banfnd = [];
bantimindall = [];
for trllop = 1:length(trldat.posdat)
    
    if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
        [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
        banfnd(trllop,1) = 1;
    else
        bantimind = length(trldat.time{trllop});
        banfnd(trllop,1) = 0;
    end
    xdif = trldat.posdat{trllop}(1,1:bantimind)-trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),3);
    ydif = trldat.posdat{trllop}(2,1:bantimind)-trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),4);
    dis = sqrt(xdif.^2+ydif.^2);

    bantimindall = [bantimindall; bantimind];
    disall = [disall; dis'];
    
    avgdis = mean(dis);
    avgdispertrl = [avgdispertrl; avgdis];
    
    figure
    subplot(3,1,1)
    plot(trldat.time{trllop}(1:bantimind),dis)
    ylabel('Distance from banana')
    title(['Trial ' num2str(trllop) '; Alpha = ' num2str(alp(trllop))])
    subplot(3,1,2)
    hold on
    plot(trldat.posdat{trllop}(1,:),trldat.posdat{trllop}(2,:))
    scatter(trldat.frtpos{trllop}(trldat.frtpos{trllop}(:,2)==0,3),trldat.frtpos{trllop}(trldat.frtpos{trllop}(:,2)==0,4),'g') % banana(s) in green
    scatter(trldat.frtpos{trllop}(trldat.frtpos{trllop}(:,2)==1,3),trldat.frtpos{trllop}(trldat.frtpos{trllop}(:,2)==1,4),'r') % cherry in red
    scatter(trldat.posdat{trllop}(1,1),trldat.posdat{trllop}(2,1),'b') % start pos in blue
    xlim([-11.2 11.2]);ylim([-11.2 11.2])
    subplot(3,1,3)
    imagesc(grdedg(1:end-1)+mean(diff(grdedg))/2,grdedg(1:end-1)+mean(diff(grdedg))/2,squeeze(grdall(trllop,:,:))');
    axis xy; colormap hot
    title(['total travel time = ' num2str(size(trldat.posdat{trllop},2)/1000) ' sec'])
    
    set(gcf,'Position',[1000 51 343 937])
    
    figdir = 'R:\Buffalo Lab\Virtual Navigation\Figures\virtual water maze, behavioral analysis';
    export_fig(fullfile(figdir,['pos_dis_JN_BR_15_04_10_13_11_trl' num2str(trllop) '.png']))
%     pause
    close
    
end

% bantimindall is basically the same as bantim, except for trials when the
% monkey doesn't get the banana (trial times out)


%% plot eat times across trials

% trial start - monkey must get banana, then cherry
% on trials where monkey never gets banana, no frttims

trlind_frttim = [];
frttimarr = [];
for k=1:length(trldat.time)
    trlind_frttim = [trlind_frttim; ones(size(trldat.frttim{k},1),1)*k];
    frttimarr = [frttimarr; trldat.frttim{k}(:,1)];
end

figure;scatter(trlind_frttim,frttimarr)

%% scatter all positions, all trials

% get dimensions of the arena

posdatall = [];
for k=1:length(trldat.posdat)
    posdatall = [posdatall; trldat.posdat{k}'];
end

figure;plot(posdatall(:,1),posdatall(:,2))


%% scatter all fruit positions, all trials

% first fruit position is the banana

frtposall = [];
banposall = [];
for k=1:length(trldat.frtpos)
    frtposall = [frtposall; trldat.frtpos{k}(:,3:4)];
    banposall = [banposall; trldat.frtpos{k}(trldat.frtpos{k}(:,2)==0,3:4)];
end

hold on;scatter(frtposall(:,1),frtposall(:,2))

% banposall can be used to see how many times a banana position is
% repeated; this may be useful to know at some point, but be careful to
% account for the fact that when the monkey gets close to a transparent or
% invisible banana, the banana appears in a new location right in front of
% him


%%
    
    
% hold on;scatter(1.5,1.5e6,'MarkerEdgeColor',[1 1 1])