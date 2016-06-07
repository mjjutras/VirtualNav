load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_sacdat140503.mat')
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_navdat140425.mat')
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN130604.1_pEpisode_Ch4_140503.mat')
% load('C:\Users\michael.jutras\Documents\MATLAB\Virtual\Data (MAT files)\JN130604.1_wavpow_Ch4_140505.mat')

freqs = (2^(1/8)).^(8:0.2:30);

chnsel = 4; % JN130604.1, ch 4 has the best theta

% % use 3 Hz and 6 Hz as boundaries for theta
% frqlo = ft_nearest(freqs,3);
% frqhi = ft_nearest(freqs,6);

% use 3 Hz and 5 Hz as boundaries for theta
frqlo = ft_nearest(freqs,3.5);
frqhi = ft_nearest(freqs,4.5);

% determine when theta bouts occur in trials
boutind = max(unionvector(frqlo:frqhi,:));

% determine bout duration and interbout interval
difbout = diff(boutind);
difboutpos = difbout==1;
difboutneg = difbout==-1;
bout_start = trlind(difboutpos);
bout_end = trlind(difboutneg);

%% interbout interval & bout duration

% plot to test
figure;imagesc(trlind(1:10000),freqs(frqlo:frqhi),unionvector(frqlo:frqhi,1:10000))
hold on;for k=1:3;line([bout_start(k) bout_start(k)],ylim,'Color','w');line([bout_end(k) bout_end(k)],ylim,'Color','w');end

%%% interbout interval
interbout = bout_start(2:end)-bout_end(1:end-1);

% plot interbout interval histogram
figure;hist(interbout,500)
xlim([0 5])
hold on;line([median(interbout) median(interbout)],ylim,'Color','r','LineWidth',2)
title(['median interbout interval: ' num2str(median(interbout)) ' sec'])
xlabel('Time (sec)')
ylabel('Count')
box off
set(gca,'TickDir','out')

%%% bout duration
boutdur = bout_end-bout_start;

% plot bout duration histogram
figure;hist(boutdur,500)
xlim([0 6])
hold on;line([median(boutdur) median(boutdur)],ylim,'Color','r','LineWidth',2)
title(['median bout duration: ' num2str(median(boutdur)) ' sec'])
xlabel('Time (sec)')
ylabel('Count')
box off
set(gca,'TickDir','out')

%% autocorrelation

% first concatenate data arrays from specified channel
datcat = nan(1,data.sampleinfo(end,2)-data.sampleinfo(1,1)+1);
c=1;
for trllop = 1:length(data.trial)
    datcat(c:c+size(data.trial{trllop},2)-1)=data.trial{trllop}(chnsel,:);
    c=c+size(data.trial{trllop},2);
end
timcat = data.sampleinfo(1,1):data.sampleinfo(end,2);

acf_all=nan(size(bout_start,2),1001);
for boutlop = 1:length(bout_start)
    if bout_start(boutlop)*1000>timcat(1) && bout_end(boutlop)*1000<timcat(end)
        if boutdur(boutlop)>0.5
            timind = ft_nearest(timcat,bout_start(boutlop)*1000):ft_nearest(timcat,bout_end(boutlop)*1000);
            [acf,lags,bounds] = autocorr(datcat(timind),length(timind)-1);
            if length(acf)<1001
                acf_all = [acf_all; [acf nan(1,1001-length(acf))]];
            else
                acf_all = [acf_all; acf(1:1001)];
            end
        end
    end
end

figure;plot(-1000:1000,cat(2,fliplr(nanmean(acf_all(:,1:1000),1)),nanmean(acf_all,1)))
            
%% saccade/fixation aligned LFP

sacali_notrunc = nan(size(data_eye.artifact{1},1),601);
sacali_trunc = nan(size(data_eye.artifact{1},1),2001);
c=1;
d=1;
ft_progress('init', 'etf',     'Please wait...');
for saclop = 1:size(data_eye.artifact{1},1)
    ft_progress(saclop/size(data_eye.artifact{1},1), 'Processing event %d from %d', saclop, size(data_eye.artifact{1},1));
    if (data_eye.artifact{1}(saclop,1)*1000-200)>timcat(1) && (data_eye.artifact{1}(saclop,1)*1000+400)<timcat(end)
        timsel1 = ft_nearest(timcat,data_eye.artifact{1}(saclop,1)*1000-200);
        timsel2 = ft_nearest(timcat,data_eye.artifact{1}(saclop,1)*1000+400);
        sacali_notrunc(c,:) = datcat(timsel1:timsel2);
        c=c+1;
%     else
%         sacali_notrunc=sacali_notrunc(1:end-1,:);
    end
    if saclop~=1 && saclop~=size(data_eye.artifact{1},1) && (data_eye.artifact{1}(saclop-1,1)*1000)>timcat(1) && (data_eye.artifact{1}(saclop+1,1)*1000)<timcat(end)
        dum1 = datcat(ft_nearest(timcat,data_eye.artifact{1}(saclop-1)*1000):ft_nearest(timcat,data_eye.artifact{1}(saclop)*1000-1));
        dum2 = datcat(ft_nearest(timcat,data_eye.artifact{1}(saclop)*1000):ft_nearest(timcat,data_eye.artifact{1}(saclop+1)*1000));
        if length(dum1)>1000
            dum1=dum1(end-999:end);
        end
        if length(dum2)>1001
            dum2 = dum2(1:1001);
        end
        sacali_trunc(d,:) = [nan(1,1000-length(dum1)) dum1 dum2 nan(1,1001-length(dum2))];
        d=d+1;
%     else
%         sacali_trunc=sacali_trunc(1:end-1,:);
    end
end
ft_progress('close')

save('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_sacali_acf_141105.mat','sacali_trunc','sacali_notrunc','acf_all')

%%

load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_sacali_acf_141105.mat')

figure; plot(-1000:1000,[fliplr(nanmean(acf_all,1)) nanmean(acf_all(:,2:end),1)])
line([0 0],ylim,'Color','k','LineStyle','--')
box off
set(gca,'TickDir','out')
xlabel('Time (ms)');ylabel('Auto-correlation')
