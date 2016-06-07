% load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1305311_navdat140425.mat')
% % load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN130531.1-sorted_pEpisode140428.mat')
% load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN130531.1_wavpow_Ch1_140512.mat')
% load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1305311_Pepisode_Ch1_140428.mat')
% load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1305311_sacdat140428.mat')

load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_navdat140425.mat')
% load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN130531.1-sorted_pEpisode140428.mat')
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN130604.1_wavpow_Ch4_140505.mat')
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN130604.1_pEpisode_Ch4_140503.mat')
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_sacdat140503.mat')

data.velsmooth=cell(size(data.trial));
for trllop=1:length(data.trial)
    
    e_x = data.posdat{trllop}(1,:);
    e_y = data.posdat{trllop}(2,:);
    
    % % differentiate and multiply with sampling rate without filtering
    % % this gives the eye velocity in the horizontal and vertical domains
    x_v = diff(e_x) .* data.hdr.Fs;
    y_v = diff(e_y) .* data.hdr.Fs;
    
    % combine x- and y-velocity to get eye velocity in degrees/second
    vel = abs(complex(x_v,y_v));
    
    data.velsmooth{trllop} = density(vel,100,'gauss');
    data.velsmooth{trllop}(end+1) = data.velsmooth{trllop}(end);
    
end

trltimall=[];
velsmoothall=[];
dirdatall=[];
for k=1:length(data.trial)
    
    trltimall = [trltimall data.cfg.previous{2}.trl(k,1):data.cfg.previous{2}.trl(k,2)];
    velsmoothall = [velsmoothall data.velsmooth{k}];
    dirdatall = [dirdatall density(data.dirdat{k},100,'gauss')];
    
end

dirdatall_diff = abs([diff(dirdatall) 0]);

trlindms=round(trlind*1000);

timealign=find(ismember(trlindms,trltimall));
% trlindms(timealign(1)) the same as trltimall(1)
% trlindms(timealign(end)) the same as trltimall(end)

edges = 0:0.002:3;
[n, bin] = histc(velsmoothall,edges);
figure;bar(edges,n,'histc')

% freqs = (2^(1/8)).^(8:42);
freqs = (2^(1/8)).^(8:0.2:30);

% centers = 0.01:0.02:2.99;
% pepbin=nan(length(centers),length(freqs));
% for binlop = 1:length(centers)
%     findbinind = bin==binlop;
%     pepbin(binlop,:)=mean(B(:,timealign(findbinind)),2);
% end


%%
% find the timepoints where turning speed was above threshold and forward
% velocity was below threshold

turnthresh = 5e-4;
velthresh_low = 0.1;
velthresh_high = 2.5;

turnind = dirdatall_diff>turnthresh & velsmoothall<velthresh_low;
moveind = velsmoothall>velthresh_high;

powturn = mean(B(:,timealign(turnind)),2);
powmove = mean(B(:,timealign(moveind)),2);

figure;hold on
plot(freqs,mean(B,2),'b')
plot(freqs,powturn,'r')
plot(freqs,powmove,'g')
legend('All timepoints','Turning/Exploring','Moving/Navigating')


%% look at movement speed and saccade rate histogram

sacpersec_smooth = density(data_eye.sacpersec{1},1000,'gauss');

sacpersec_turn = sacpersec_smooth(timealign(turnind));
sacpersec_move = sacpersec_smooth(timealign(moveind));
% sacpersec_turn = data_eye.sacpersec{1}(timealign(turnind));
% sacpersec_move = data_eye.sacpersec{1}(timealign(moveind));

edges = 0:0.2:10;
centers = 0.05:0.2:9.95;
n1 = histc(sacpersec_turn,edges);
n2 = histc(sacpersec_move,edges);
n_all = histc(sacpersec_smooth(timealign),edges);

n1 = n1/length(sacpersec_turn);
n2 = n2/length(sacpersec_move);

figure;hold on
bar(centers,n1(1:end-1),0.5);
yh=ylim;
set(gca,'TickDir','out')
box off
xlabel('Saccades/sec')
ylabel('Proportion')
line([median(sacpersec_turn) median(sacpersec_turn)],ylim,'Color','r')
title(['Exploring; median = ' num2str(median(sacpersec_turn))])

figure;hold on
bar(centers,n2(1:end-1),0.5);
ylim(yh)
set(gca,'TickDir','out')
box off
xlabel('Saccades/sec')
ylabel('Proportion')
line([median(sacpersec_move) median(sacpersec_move)],ylim,'Color','r')
title(['Acquiring; median = ' num2str(median(sacpersec_move))])

%%

edges = 0:0.1:8;
centers = 0.05:0.1:7.95;


[n, bin] = histc(sacpersec_smooth,edges);
[n, bin_turn] = histc(sacpersec_turn,edges);
[n, bin_move] = histc(sacpersec_move,edges);

% unionvector_turn = unionvector(:,timealign(turnind));
% unionvector_move = unionvector(:,timealign(moveind));
unionvector_turn = B(:,timealign(turnind));
unionvector_move = B(:,timealign(moveind));

pepbin=nan(length(edges)-1,length(freqs));
pepbin_turn=nan(length(edges)-1,length(freqs));
pepbin_move=nan(length(edges)-1,length(freqs));
for binlop = 1:length(edges)-1
%     pepbin(binlop,:)=nanmean(unionvector(:,bin==binlop),2);
    pepbin(binlop,:)=nanmean(B(:,bin==binlop),2);
    pepbin_turn(binlop,:)=nanmean(unionvector_turn(:,bin_turn==binlop),2);
    pepbin_move(binlop,:)=nanmean(unionvector_move(:,bin_move==binlop),2);
end

pepbin_turn(isnan(pepbin_turn))=0;
pepbin_move(isnan(pepbin_move))=0;
  
% interpolate to plot with correct axis
frq_interp = freqs(1):0.01:freqs(end);
pepbin_interp = nan(size(pepbin,1),length(frq_interp));
pepbin_turn_interp = nan(size(pepbin_turn,1),length(frq_interp));
pepbin_move_interp = nan(size(pepbin_move,1),length(frq_interp));
for cntlop = 1:size(pepbin,1)
    for frqlop = 1:size(pepbin,2)
        pepbin_interp(cntlop,ft_nearest(frq_interp,freqs(frqlop))) = pepbin(cntlop,frqlop);
        pepbin_turn_interp(cntlop,ft_nearest(frq_interp,freqs(frqlop))) = pepbin_turn(cntlop,frqlop);
        pepbin_move_interp(cntlop,ft_nearest(frq_interp,freqs(frqlop))) = pepbin_move(cntlop,frqlop);
    end
    pepbin_interp(cntlop,:)=inpaint_nans(pepbin_interp(cntlop,:),2);
    pepbin_turn_interp(cntlop,:)=inpaint_nans(pepbin_turn_interp(cntlop,:),2);
    pepbin_move_interp(cntlop,:)=inpaint_nans(pepbin_move_interp(cntlop,:),2);
end


peakfrq = nan(1,size(pepbin_interp,1));
peakfrq_turn = nan(1,size(pepbin_turn_interp,1));
peakfrq_move = nan(1,size(pepbin_move_interp,1));
for binlop = 1:size(pepbin_interp,1)
%     peakfrq(binlop) = mean(frq_interp(pepbin_interp(binlop,:)==max(pepbin_interp(binlop,:))));
%     peakfrq_turn(binlop) = mean(frq_interp(pepbin_turn_interp(binlop,:)==max(pepbin_turn_interp(binlop,:))));
%     peakfrq_move(binlop) = mean(frq_interp(pepbin_move_interp(binlop,:)==max(pepbin_move_interp(binlop,:))));
    % following lines for max power between 3 and 6 Hz
    peakfrq(binlop) = mean(frq_interp(pepbin_interp(binlop,:)==max(pepbin_interp(binlop,101:401))));
    peakfrq_turn(binlop) = mean(frq_interp(pepbin_turn_interp(binlop,:)==max(pepbin_turn_interp(binlop,101:401))));
    peakfrq_move(binlop) = mean(frq_interp(pepbin_move_interp(binlop,:)==max(pepbin_move_interp(binlop,101:401))));
end
    

%%

figure;hold on
imagesc(frq_interp(ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)), ...
    centers(ft_nearest(centers,2):ft_nearest(centers,8)), ...
    pepbin_interp(ft_nearest(centers,2):ft_nearest(centers,8),ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)))
axis xy
colormap hot
colorbar
set(gca,'clim',[2e8 5.5e8])
plot(peakfrq,centers,'Color','b','LineWidth',2)
% xlim([2 8])
% ylim([2 8])
xlim([3 6])
ylim([3 6])
xlabel('LFP freq (Hz)')
ylabel('saccade rate (Hz)')
box on
set(gca,'TickDir','out')


figure;hold on
imagesc(frq_interp(ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)), ...
    centers(ft_nearest(centers,2):ft_nearest(centers,8)), ...
    pepbin_turn_interp(ft_nearest(centers,2):ft_nearest(centers,8),ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)))
axis xy
colormap hot
colorbar
set(gca,'clim',[2e8 5.5e8])
plot(peakfrq_turn,centers,'Color','b','LineWidth',2)
% xlim([2 8])
% ylim([2 8])
xlim([3 6])
ylim([3 6])
xlabel('LFP freq (Hz)')
ylabel('Saccade rate (Hz)')
box on
set(gca,'TickDir','out')


figure;hold on
imagesc(frq_interp(ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)), ...
    centers(ft_nearest(centers,2):ft_nearest(centers,8)), ...
    pepbin_move_interp(ft_nearest(centers,2):ft_nearest(centers,8),ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)))
axis xy
colormap hot
colorbar
set(gca,'clim',[2e8 5.5e8])
plot(peakfrq_move,centers,'Color','b','LineWidth',2)
% xlim([2 8])
% ylim([2 8])
xlim([3 6])
ylim([3 6])
xlabel('LFP freq (Hz)')
ylabel('saccade rate (Hz)')
box on
set(gca,'TickDir','out')
