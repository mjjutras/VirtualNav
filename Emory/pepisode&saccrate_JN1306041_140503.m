load('C:\Users\michael.jutras\Documents\MATLAB\Virtual\Data (MAT files)\JN1306041_navdat140425.mat')


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
for k=1:length(data.trial)
    
    trltimall = [trltimall data.cfg.previous{2}.trl(k,1):data.cfg.previous{2}.trl(k,2)];
    velsmoothall = [velsmoothall data.velsmooth{k}];
    
end


timealign=find(ismember(trlindms,trltimall));
% trlindms(timealign(1)) the same as trltimall(1)
% trlindms(timealign(end)) the same as trltimall(end)


edges = 0:0.002:3;
[n, bin] = histc(velsmoothall,edges);
figure;bar(edges,n,'histc')


centers = 0.01:0.02:2.99;

pepbin=nan(length(centers),length(freqs));
for binlop = 1:length(centers)
    findbinind = bin==binlop;
    pepbin(binlop,:)=mean(unionvector(:,timealign(findbinind)),2);
end

%% plot pepisode for each movement category: fast, slow, still

movefast = velsmoothall>2.78 & velsmoothall<2.82;
moveslow = velsmoothall>0.05 & velsmoothall<0.07;
standstill = velsmoothall<0.002;

pepfast = mean(unionvector(:,timealign(movefast)),2);
pepslow = mean(unionvector(:,timealign(moveslow)),2);
pepstill = mean(unionvector(:,timealign(standstill)),2);

figure;hold on
plot(freqs,mean(unionvector,2),'b')
plot(freqs,pepfast,'r')
plot(freqs,pepslow,'m')
plot(freqs,pepstill,'g')
legend('All timepoints','Moving fast','Moving slow','Still')

%% look at movement speed and saccade rate histogram

load('C:\Users\michael.jutras\Documents\MATLAB\Virtual\FT data\JN1305311_sacdat140428.mat')

sacpersec_smooth = density(data_eye.sacpersec{1},1000,'gauss');

sacpersec_fast = sacpersec_smooth(timealign(movefast));
sacpersec_slow = sacpersec_smooth(timealign(moveslow));
sacpersec_still = sacpersec_smooth(timealign(standstill));

edges = 0:0.1:10;
centers = 0.05:0.1:9.95;
n1 = histc(sacpersec_fast,edges);
n2 = histc(sacpersec_slow,edges);
n3 = histc(sacpersec_still,edges);
n_all = histc(sacpersec_smooth(timealign),edges);


figure;h1=axes;
bar(centers,n1(1:end-1));
h2=axes;
set(h2,'Color','none')
hold on;plot(freqs,pepfast)
set(h1,'xlim',[2 8])
set(h2,'xlim',[2 8])

%% look at only fast movement speed, separate by saccade rate
% cut off at 4 Hz

% sacpersec_smooth, unionvector and trlind are same length

unionvector_fast = unionvector(:,timealign(movefast));
% already have sacpersec_fast

sacpersec_fast_pos4 = sacpersec_fast(sacpersec_fast>=4);
sacpersec_fast_neg4 = sacpersec_fast(sacpersec_fast<4);
unionvector_fast_pos4 = unionvector_fast(:,sacpersec_fast>=4);
unionvector_fast_neg4 = unionvector_fast(:,sacpersec_fast<4);

edges = 0:0.1:10;
centers = 0.05:0.1:9.95;
n1 = histc(sacpersec_fast_pos4,edges);
n2 = histc(sacpersec_fast_neg4,edges);

figure;
h1=axes;
b1=bar(centers,n1(1:end-1));
set(h1,'YTickLabel',[])
set(h1,'YTick',[])
h2=axes;
b2=bar(centers,n2(1:end-1));
set(b2,'FaceColor','r')
set(h2,'Color','none')
set(h2,'YTickLabel',[])
set(h2,'YTick',[])
ylim2=get(h2,'ylim');
set(h1,'ylim',ylim2)
h3=axes;
b3=plot(freqs,mean(unionvector_fast_pos4,2),'b');
set(b3,'LineWidth',2)
set(h3,'xlim',[0 10])
set(h3,'Color','none')
ylabel('Pepisode')
h4=axes;
b4=plot(freqs,mean(unionvector_fast_neg4,2),'r');
set(b4,'LineWidth',2)
set(h4,'xlim',[0 10])
set(h4,'Color','none')
ylim3=get(h3,'ylim');
set(h4,'ylim',ylim3)
set(h4,'YTickLabel',[])
set(h4,'YTick',[])
xlabel('Frequency (Hz)')



