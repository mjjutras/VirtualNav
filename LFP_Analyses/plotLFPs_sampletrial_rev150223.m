BRnam = 'JN140825011';

% datdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trial data\';
% datdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\trldat\';
datdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\trldat\olderTRLDATfiles\';
load(fullfile(datdir, [BRnam '_trldat.mat']))

%% plot all the trials; pause between each figure

for trllop = 1:length(data.trial)
    
    close all
    
    for chnlop = 1:36

        if chnlop>=1 && chnlop<=12
            figure(1); set(gcf,'Position',[92 510 560 420])
            yaxplc = chnlop;
        elseif chnlop>=13 && chnlop<=24
            figure(2); set(gcf,'Position',[690 510 560 420])
            yaxplc = chnlop-12;
        elseif chnlop>=25 && chnlop<=36
            figure(3); set(gcf,'Position',[1209 510 560 420])
            yaxplc = chnlop-24;
        end

        hold on
    %     if kurmat(chnlop,trltouse)<3
            plot(data.time{trllop},data.trial{trllop}(chnlop,:)+2000*yaxplc,'Color','b')
    %     else
    %         plot(data.time{trltouse},data.trial{trltouse}(chnlop,:)+2000*yaxplc,'Color',[0.8 0.8 0.8],'LineStyle','--')
    %     end

    end
    
    pause
    
end

%% plot only trial with lowest average kurtosis
% use kurtosis to select "bad" channels, and plot these in gray

kurmat = nan(size(data.trial{1},1),length(data.trial));
for trllop = 1:length(data.trial)
    for chnlop = 1:36
        kurmat(chnlop,trllop) = kurtosis(data.trial{trllop}(chnlop,:));
    end
end

if strmatch(BRnam, 'JN140825011');
    kurmat(:,4) = nan(size(kurmat,1),1);
end

figure;hist(reshape(kurmat,1,numel(kurmat)),1000)

kurpertrl = mean(kurmat,1);
trltouse = find(kurpertrl==min(kurpertrl));

figure;bar(1:size(kurmat,1),kurmat(:,trltouse));

trlsel_time = data.time{trltouse};
trlsel_trial = data.trial{trltouse};

close all
for chnlop = 1:36

    if chnlop>=1 && chnlop<=12
        figure(1)
        yaxplc = chnlop;
    elseif chnlop>=13 && chnlop<=24
        figure(2)
        yaxplc = chnlop-12;
    elseif chnlop>=25 && chnlop<=36
        figure(3)
        yaxplc = chnlop-24;
    end
    
    hold on
    if kurmat(chnlop,trltouse)<3
        plot(trlsel_time,trlsel_trial(chnlop,:)+2000*yaxplc,'Color','b')
    else
        plot(trlsel_time,trlsel_trial(chnlop,:)+2000*yaxplc,'Color',[0.8 0.8 0.8],'LineStyle','--')
    end

end


% re-reference to the common average of each array
trlsel_trial1 = trlsel_trial(1:12,:);
trlavg1 = mean(trlsel_trial1(kurmat(1:12,trltouse)<3,:),1);
for k=1:size(trlsel_trial1,1)
    trlsel_trial1(k,:) = trlsel_trial1(k,:) - trlavg1;
end
trlsel_trial2 = trlsel_trial(13:24,:);
trlavg2 = mean(trlsel_trial2(kurmat(13:24,trltouse)<3,:),1);
for k=1:size(trlsel_trial2,1)
    trlsel_trial2(k,:) = trlsel_trial2(k,:) - trlavg2;
end
trlsel_trial3 = trlsel_trial(25:36,:);
trlavg3 = mean(trlsel_trial3(kurmat(25:36,trltouse)<3,:),1);
for k=1:size(trlsel_trial3,1)
    trlsel_trial3(k,:) = trlsel_trial3(k,:) - trlavg3;
end

trlsel_trialR = cat(1,trlsel_trial1,trlsel_trial2,trlsel_trial3);

close all
for chnlop = 1:36

    if chnlop>=1 && chnlop<=12
        figure(1)
        yaxplc = chnlop;
    elseif chnlop>=13 && chnlop<=24
        figure(2)
        yaxplc = chnlop-12;
    elseif chnlop>=25 && chnlop<=36
        figure(3)
        yaxplc = chnlop-24;
    end
    
    hold on
    
%     if chnlop>=1 && chnlop<=12
%         plot(trlsel_time,trlavg1+2000*yaxplc,'Color','g')
%     elseif chnlop>=13 && chnlop<=24
%         plot(trlsel_time,trlavg2+2000*yaxplc,'Color','g')
%     elseif chnlop>=25 && chnlop<=36
%         plot(trlsel_time,trlavg3+2000*yaxplc,'Color','g')
%     end

    if kurmat(chnlop,trltouse)<3
%         plot(trlsel_time,trlsel_trial(chnlop,:)+2000*yaxplc,'Color','r')
        plot(trlsel_time,trlsel_trialR(chnlop,:)+2000*yaxplc,'Color','b')
    else
        plot(trlsel_time,trlsel_trialR(chnlop,:)+2000*yaxplc,'Color',[0.8 0.8 0.8],'LineStyle','--')
    end
        
end


%%

dataP = data;
dataP.fsample = 1000;
dataP.sampleinfo = dataP.trl;
for trllop = 1:length(data.time)
    dataP.time{trllop}=(dataP.time{trllop}-dataP.time{trllop}(1))/1000;
end
dataP.label = {'AD01'; 'AD02'; 'AD03'; 'AD04'; 'AD05'; 'AD06'; 'AD07'; ...
    'AD08'; 'AD09'; 'AD10'; 'AD11'; 'AD12'; 'AD13'; 'AD14'; 'AD15'; ...
    'AD16'; 'AD17'; 'AD18'; 'AD19'; 'AD20'; 'AD21'; 'AD22'; 'AD23'; ...
    'AD24'; 'AD25'; 'AD26'; 'AD27'; 'AD28'; 'AD29'; 'AD30'; 'AD31'; ...
    'AD32'; 'AD33'; 'AD34'; 'AD35'; 'AD36'};

% preprocess for Fieldtrip
cfg=[];
cfg.continuous = 'yes';
dataP = ft_preprocessing(cfg,dataP);

% do the spectral analysis - time-averaged
clear cfg
cfg.output      = 'pow';
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.keeptrials  = 'no';
cfg.foi         = 0:0.5:100;
cfg.tapsmofrq   = 2;
cfg.trials      = [1:3 5:11]; % first 10 "clean" trials
freq = ft_freqanalysis(cfg, dataP);

avgpow_lowfrq = mean(10*(log10(freq.powspctrm(:,1:ft_nearest(freq.freq,20)))),2);

% get out


cfg=[];
cfg.method='summary';
dum = ft_rejectvisual(cfg, dataP);

% % only run this line after inspecting data
% data=dum;

