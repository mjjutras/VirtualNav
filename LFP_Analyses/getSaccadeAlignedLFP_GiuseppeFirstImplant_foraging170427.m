% fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';
fildir = 'C:\Data\Blackrock_VR';

% BRnam = 'JN140813002';
BRnam = 'JN140815002';

%% load the data and create double-precision versions of eye data

% NS6 = openNSx(fullfile(fildir,'JN140821002.ns6'),'read','c:1','s:30');
% starts with sample 1 and then skips every 30, so first sample should
% align with the timestamp displayed in NS6.MetaTags.Timestamp

NS2 = openNSx(fullfile(fildir,[BRnam '.ns2']),'read');

xchan = [];
ychan = [];
for k = 1:size(NS2.ElectrodesInfo,2)
    if strncmp(NS2.ElectrodesInfo(k).Label,'ainp1',5)
        xchan = k;
    elseif strncmp(NS2.ElectrodesInfo(k).Label,'ainp2',5)
        ychan = k;
    end
end

e_x = double(NS2.Data(xchan,:));
e_y = double(NS2.Data(ychan,:));


%% investigate different filtering parameters (can skip)

% compare two sets of filtering parameters
fltord1 = 40;
lowpasfrq1 = 40;
fltord2 = 60;
lowpasfrq2 = 18;

Fs = NS2.MetaTags.SamplingFreq; % sampling rate in Hz
nyqfrq = Fs ./ 2;

flt1 = fir2(fltord1,[0,lowpasfrq1./nyqfrq,lowpasfrq1./nyqfrq,1],[1,1,0,0]);
flt2 = fir2(fltord2,[0,lowpasfrq2./nyqfrq,lowpasfrq2./nyqfrq,1],[1,1,0,0]);

% low pass filter the eye position data_eye
e_x_lowpassfilter1 = filtfilt(flt1,1, e_x); e_y_lowpassfilter1 = filtfilt(flt1,1, e_y);
e_x_lowpassfilter2 = filtfilt(flt2,1, e_x); e_y_lowpassfilter2 = filtfilt(flt2,1, e_y);

% differentiate and multiply with sampling rate to get velocity as deg/sec
% this gives the eye velocity in the horizontal and vertical domains
x_vF1 = diff(e_x_lowpassfilter1) .* Fs; y_vF1 = diff(e_y_lowpassfilter1) .* Fs;
x_vF2 = diff(e_x_lowpassfilter2) .* Fs; y_vF2 = diff(e_y_lowpassfilter2) .* Fs;

% differentiate and multiply with sampling rate without filtering
x_v = diff(e_x) .* Fs; y_v = diff(e_y) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
vel = abs(complex(x_v,y_v));
velF1 = abs(complex(x_vF1,y_vF1));
velF2 = abs(complex(x_vF2,y_vF2));

lim1 = median(velF1/0.6745)*3;
lim2 = median(velF2/0.6745)*3;

figure
ax(1)=subplot(2,1,1);
plot(1:100000,[e_x(1:100000); e_y(1:100000)])
ax(2)=subplot(2,1,2);
plot(1:100000,[vel(1:100000); velF1(1:100000); velF2(1:100000)])
hold on
line(xlim,[lim1 lim1],'Color','r')
line(xlim,[lim2 lim2],'Color','g')
linkaxes(ax,'x')
legend({'raw' 'filt1' 'filt2'})


%% filter the eye signal and calculate velocity for the original data

fltord = 40;
lowpasfrq = 40;

artpadding = 0.01;
Fs = NS2.MetaTags.SamplingFreq; % sampling rate in Hz
nyqfrq = Fs ./ 2;

flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);

% low pass filter the eye position data_eye
e_xbuf = [e_x(100:-1:1) e_x e_x(end:-1:end-99)]; % add buffer for filtering
e_ybuf = [e_y(100:-1:1) e_y e_y(end:-1:end-99)]; % add buffer for filtering
xss = filtfilt(flt,1,e_xbuf);
yss = filtfilt(flt,1,e_ybuf);
xss = xss(101:end-100); % remove buffer for filtering
yss = yss(101:end-100); % remove buffer for filtering

% differentiate and multiply with sampling rate to get velocity as deg/sec
% this gives the eye velocity in the horizontal and vertical domains
x_vF = diff(xss) .* Fs;
y_vF = diff(yss) .* Fs;

% differentiate and multiply with sampling rate without filtering
x_v = diff(e_x) .* Fs;
y_v = diff(e_y) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
vel = abs(complex(x_v,y_v));
velF = abs(complex(x_vF,y_vF));

lim = median(velF/0.6745)*5; % threshold taken from Quiroga (http://www.scholarpedia.org/article/Spike_sorting)

%detect saccade starts and saccade ends
sacbeg = find(diff(velF > lim) > 0);
sacend = find(diff(velF > lim) < 0);
if velF(end)>lim
    sacbeg = sacbeg(1:end-1);
end
if velF(1)>lim
    sacend = sacend(2:end);
end

if size(sacbeg,1)
    sacbeg = sacbeg';
    sacend = sacend';
end
    
% artifact = round([sacbeg(:) - artpadding*Fs sacend(:) + artpadding*Fs]);
artifact = round([sacbeg(:) sacend(:)]);

% figure
% ax(1) = subplot(2,1,1);
% hold on
% scatter(1:100000,e_x(1:100000),'.b')
% scatter(1:100000,e_y(1:100000),'.r')
% for k=1:find(artifact(:,2)>100000,1,'first')
%     line([artifact(k,1) artifact(k,1)],ylim,'Color','g')
%     line([artifact(k,2) artifact(k,2)],ylim,'Color','r')
% end
% ax(2) = subplot(2,1,2);
% hold on
% scatter(1:100000,velF(1:100000),'.')
% for k=1:find(artifact(:,2)>100000,1,'first')
%     line([artifact(k,1) artifact(k,1)],ylim,'Color','g')
%     line([artifact(k,2) artifact(k,2)],ylim,'Color','r')
% end
% linkaxes(ax,'x')

% h=[];
% for k=1:size(artifact,1)
%     h=[h abs(diff(e_x(artifact(k,1):artifact(k,2))))];
% end
    

% % find saccade start/end in task time
% sacdum = [NS2_timestamp(artifact(:,1))' NS2_timestamp(artifact(:,2))'];

% merge artifacts when inter-artifact interval is less than 20 ms
sacarr = [];
iai=(artifact(2:end,1)-artifact(1:end-1,2))>20;
sacdumcol1 = artifact(2:end,1);
sacdumcol2 = artifact(1:end-1,2);
sacarr(:,1) = [artifact(1,1); sacdumcol1(iai,1)];
sacarr(:,2) = [sacdumcol2(iai,1); artifact(end,2)];

% figure;hold on
% scatter(1:100000,e_x(1:100000),'.b')
% scatter(1:100000,e_y(1:100000),'.r')
% for k=1:find(sacarr(:,2)>100000,1,'first')
%     line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
%     line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
% end

% clear vel* x_* y_* e_*


%% plot x and y eye data

figure
ax(1) = subplot(3,1,1);
plot(e_x,'b')
ax(2) = subplot(3,1,2);
plot(e_y,'r')
ax(3) = subplot(3,1,3);
% plot(vel,'g')
plot(velF,'g')
linkaxes(ax,'x')

%% mark epochs where gaze hits the edges of the eyetracker range

% define "borders"
xthresh = [diff(minmax(e_x))*0.001+min(e_x) diff(minmax(e_x))*0.999+min(e_x)];
ythresh = [diff(minmax(e_y))*0.001+min(e_y) diff(minmax(e_y))*0.999+min(e_y)];

% mark 
x_bound = e_x<xthresh(1) | e_x>xthresh(2);
y_bound = e_y<ythresh(1) | e_y>ythresh(2);
eyebound = x_bound | y_bound;

% figure
% ax(1) = subplot(3,1,1);
% plot(e_x,'b')
% ax(2) = subplot(3,1,2);
% plot(e_y,'r')
% ax(3) = subplot(3,1,3);
% plot(eyebound,'g')
% ylim([-0.25 1.25])
% linkaxes(ax,'x')

% figure
% ax(1) = subplot(3,1,1);hold on
% plot(1:100000,e_x(1:100000),'.b')
% for k=1:find(sacarr(:,2)>100000,1,'first')
%     line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
%     line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
% end
% ax(2) = subplot(3,1,2);hold on
% plot(1:100000,e_y(1:100000),'.r')
% for k=1:find(sacarr(:,2)>100000,1,'first')
%     line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
%     line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
% end
% ax(3) = subplot(3,1,3);hold on
% plot(1:100000,eyebound(1:100000),'k')
% ylim([-0.25 1.25])
% for k=1:find(sacarr(:,2)>100000,1,'first')
%     line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
%     line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
% end
% linkaxes(ax,'x')

%% choose saccades that don't precede fixations at edges(borders)

sacsel = zeros(size(sacarr,1),1);
for saclop = 1:size(sacarr,1)
    if saclop~=size(sacarr,1)
        if isempty(find(eyebound(sacarr(saclop,2)+1:sacarr(saclop+1,1)-1),1))
            sacsel(saclop) = 1;
        end
    else
        if isempty(find(eyebound(sacarr(saclop,2)+1:end),1))
            sacsel(saclop) = 1;
        end
    end
end

sacsel = find(sacsel);

%% choose saccades that don't precede fixations OR FOLLOW at edges(borders)
% commented this out in preference for the previous cell, since the ERP
% obtained when running this cell is greatly dampened compared to that
% obtained when only excluding saccade that follow gaze at borders

% sacsel = zeros(size(sacarr,1),1);
% for saclop = 1:size(sacarr,1)
%     if saclop~=size(sacarr,1)
%         if saclop==1
%             eyeboundPre = eyebound(1:sacarr(saclop,1)-1);
%         else
%             eyeboundPre = eyebound(sacarr(saclop-1,2)+1:sacarr(saclop,1)-1);
%         end
%         eyeboundPost = eyebound(sacarr(saclop,2)+1:sacarr(saclop+1,1)-1);
%     else
%         eyeboundPre = eyebound(sacarr(saclop-1,2)+1:sacarr(saclop,1)-1);
%         eyeboundPost = eyebound(sacarr(saclop,2)+1:end);
%     end
%     if isempty(find(eyeboundPre,1)) && isempty(find(eyeboundPost,1))
%         sacsel(saclop) = 1;
%     end
% end
% 
% sacsel = find(sacsel);

%%

decdir = 'C:\Data\MAT\NS6 - decimated\';

% load the decimated ("skipfactor") NS6 file
load(fullfile(decdir,[BRnam '_NS6_SF30.mat']))

% NS6 file has 102 leading zeros, remove these and the signal matches the
% NS2 data without the phase delay caused by the filtering in NS2

%% get NS2 and NS6 saccade-aligned data
% note that this is hard-coded to take 36 channels of neural data
% apparently after JN140813002, stopped saving neural data to NS2 file

% sacdatNS2 = nan(length(sacsel),36,501);
sacdatNS6 = nan(length(sacsel),36,501);

ft_defaults
ft_progress('init', 'etf', 'Please wait...');
for saclop = 1:length(sacsel)
    
    ft_progress(saclop/length(sacsel), 'Processing event %d from %d', saclop, length(sacsel));

    if sacsel(saclop)~=size(sacarr,1) % don't use the very last saccade of the recording
        
        ind1 = sacarr(sacsel(saclop),1)-200; % saccade onset minus 200 ms
        ind2 = sacarr(sacsel(saclop)+1,1)-1; % onset of next saccade - 1 (end of fixation period)
        
        if ind1 > 0 && ind1+602 <= size(NS6.Data,2)
            
            if length(ind1:ind2) <= 501
                % take whole period if start of pre-saccade interval to end of
                % fixation period is less than 501 ms
%                 datbufNS2 = [double(NS2.Data(1:36,ind1:ind2)) nan(36,501-length(ind1:ind2))];
                datbufNS6 = [double(NS6.Data(1:36,(ind1:ind2)+102)) nan(36,501-length(ind1:ind2))];
            else
                % otherwise, truncate period at 300 ms after saccade onset
%                 datbufNS2 = double(NS2.Data(1:36,ind1:ind1+500));
                datbufNS6 = double(NS6.Data(1:36,(ind1:ind1+500)+102));
            end
        end
    
%         sacdatNS2(saclop,:,:) = datbufNS2;
        sacdatNS6(saclop,:,:) = datbufNS6;
        
    end
    
end
ft_progress('close')

% remove unused trial holders
% sacdatNS2 = sacdatNS2(~isnan(squeeze(sacdatNS2(:,1,1))),:,:);
sacdatNS6 = sacdatNS6(~isnan(squeeze(sacdatNS6(:,1,1))),:,:);

% trlcnt = nan(1,size(sacdatNS2,3));
% for timlop=1:size(squeeze(sacdatNS2(1,:,:)),2)
%     trlcnt(timlop) = length(find(~isnan(squeeze(sacdatNS2(:,1,timlop)))));
% end

% save(fullfile('C:\Data\MAT\SaccadeAligned\',[BRnam '_sacdat.mat']),'sacdat','NS2_timestamp','NS6_timestamp','artifact','NS2sacarr','trlcnt')


%% de-mean each LFP segment

for k=1:size(sacdatNS6,1)
    for l=1:size(sacdatNS6,2)
%         sacdatNS2(k,l,:) = sacdatNS2(k,l,:)-nanmean(sacdatNS2(k,l,:));
        sacdatNS6(k,l,:) = sacdatNS6(k,l,:)-nanmean(sacdatNS6(k,l,:));
    end
end

%% plot average/SEM saccade-aligned LFPs, one example channel

chnsel = 13;

% choose saccades where the post-saccade fixation fills the whole window
% (i.e. isn't cut off by the next saccade)
sacdatNS2_fullwin = squeeze(sacdatNS2(~isnan(sacdatNS2(:,chnsel,end)),chnsel,:));
sacdatNS6_fullwin = squeeze(sacdatNS6(~isnan(sacdatNS6(:,chnsel,end)),chnsel,:));

figure; hold on
plot(-200:300,[nanmean(sacdatNS2_fullwin,1); nanmean(sacdatNS6_fullwin,1)])
errorshade(sacdatNS2_fullwin,1,-200:300,'b')
errorshade(sacdatNS6_fullwin,1,-200:300,'r')
line([0 0],ylim,'Color','k','LineStyle','--')
legend({'NS2','NS6'})

%% save sacdat

savdir = 'C:\Data\MAT\SaccadeAligned';
% save(fullfile(savdir,[BRnam '_sacdatNS2_NS6.mat']), 'sacdatNS2', 'sacdatNS6')
save(fullfile(savdir,[BRnam '_sacdatNS6.mat']), 'sacdatNS6')

%% plot average/SEM saccade-aligned LFPs, all channels

% add offset for plotting
ofs = 30;
ofv = (ofs:ofs:12*ofs)';

figure

subplot(3,1,1)
hold on

trlsel = ~isnan(squeeze(sacdatNS2(:,1,end))); % choose trials containing data to the end of the selection window
numtrlsel = length(find(trlsel));

xaxis = -0.2:0.001:0.3;
matrix = squeeze(sacdatNS2(trlsel,1:12,:));
matavg = squeeze(mean(matrix,1));
matste = squeeze(std(matrix,0,1))/sqrt(size(matrix,1));

matavgofs = bsxfun(@plus,matavg,ofv);
plot(xaxis,matavgofs)
ColorOrder = get(gca,'ColorOrder');
ColorOrder = repmat(ColorOrder,2,1);
for k=1:12
    sempos=matavgofs(k,:)+matste(k,:);
    semneg=matavgofs(k,:)-matste(k,:);
    dum=xaxis;dum1=fliplr(-xaxis);dum1=dum1*-1;x=[dum dum1];
    dum=flipud(semneg(:));
    y=[sempos(:);dum]';

    face = ColorOrder(k,:);
    fill(x,y,face, 'FaceAlpha',.2, 'EdgeColor',face,'EdgeAlpha',0);
end

title([BRnam '; Array A; ' num2str(numtrlsel) ' fixpers'])

subplot(3,1,2)
hold on

trlsel = ~isnan(squeeze(sacdatNS2(:,1,end))); % choose trials containing data to the end of the selection window

xaxis = -0.2:0.001:0.3;
matrix = squeeze(sacdatNS2(trlsel,13:24,:));
matavg = squeeze(mean(matrix,1));
matste = squeeze(std(matrix,0,1))/sqrt(size(matrix,1));

matavgofs = bsxfun(@plus,matavg,ofv);
plot(xaxis,matavgofs)
ColorOrder = get(gca,'ColorOrder');
ColorOrder = repmat(ColorOrder,2,1);
for k=1:12
    sempos=matavgofs(k,:)+matste(k,:);
    semneg=matavgofs(k,:)-matste(k,:);
    dum=xaxis;dum1=fliplr(-xaxis);dum1=dum1*-1;x=[dum dum1];
    dum=flipud(semneg(:));
    y=[sempos(:);dum]';

    face = ColorOrder(k,:);
    fill(x,y,face, 'FaceAlpha',.2, 'EdgeColor',face,'EdgeAlpha',0);
end

title([BRnam '; Array B'])

subplot(3,1,3)
hold on

trlsel = ~isnan(squeeze(sacdatNS2(:,1,end))); % choose trials containing data to the end of the selection window

xaxis = -0.2:0.001:0.3;
matrix = squeeze(sacdatNS2(trlsel,25:36,:));
matavg = squeeze(mean(matrix,1));
matste = squeeze(std(matrix,0,1))/sqrt(size(matrix,1));

matavgofs = bsxfun(@plus,matavg,ofv);
plot(xaxis,matavgofs)
ColorOrder = get(gca,'ColorOrder');
ColorOrder = repmat(ColorOrder,2,1);
for k=1:12
    sempos=matavgofs(k,:)+matste(k,:);
    semneg=matavgofs(k,:)-matste(k,:);
    dum=xaxis;dum1=fliplr(-xaxis);dum1=dum1*-1;x=[dum dum1];
    dum=flipud(semneg(:));
    y=[sempos(:);dum]';

    face = ColorOrder(k,:);
    fill(x,y,face, 'FaceAlpha',.2, 'EdgeColor',face,'EdgeAlpha',0);
end

title([BRnam '; Array C'])
set(gcf,'Position',[777 39 560 937])

figdir = 'R:\Mike\VirtualNavigationProject\Figures\Flexshaft_firstGizImplant2014\SaccadeAligned_SigToNoise';

%% calculate signal to noise (for Chris Lewis paper)

datdir = 'C:\Data\MAT\SaccadeAligned';
chnsel = 13; % B1

f1 = load(fullfile(datdir,'JN140813002_sacdatNS2_NS6.mat'));
f2 = load(fullfile(datdir,'JN140815002_sacdatNS6.mat'));

% choose only data segments where the post-saccade fixation fills the whole window
% (i.e. isn't cut off by the next saccade)
fullwinF1 = squeeze(f1.sacdatNS6(~isnan(f1.sacdatNS6(:,chnsel,end)),chnsel,:));
fullwinF2 = squeeze(f2.sacdatNS6(~isnan(f2.sacdatNS6(:,chnsel,end)),chnsel,:));

% plot average of sacnum data segments, selected randomly; repeat
sacnum = 200;
figure; hold on
for k=1:20
    plot(-200:300,nanmean(fullwinF1(randi(size(fullwinF1,1),1,sacnum),:),1));
%     plot(-200:300,nanmean(fullwinF2(randi(size(fullwinF2,1),1,sacnum),:),1));
end

% signal window: 56-155
% noise window: -155:-56
sigwin = 56:155;
noisewin = -155:-56;
sacnum = 200;
SNR = nan(100,1);
for k=1:100 % 100 iterations
    
    trlsel = randi(size(fullwinF1,1),1,sacnum); % randomly select sacnum trials
    
    % average ABS within sigwin, for each trial
    sigavg = mean(abs(fullwinF1(trlsel,ismember(-200:300,sigwin))),2);
    % root meat square of ABS across trials
    Ps = sqrt(mean(sigavg.^2));
    
    % average ABS within noisewin, for each trial
    noiseavg = mean(abs(fullwinF1(trlsel,ismember(-200:300,noisewin))),2);
    % root meat square of ABS across trials
    Pn = sqrt(mean(noiseavg.^2));
    
    SNR(k) = 10*log10(Ps/Pn);

end


