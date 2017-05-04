% fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';
fildir = 'C:\Data\Blackrock_VR';
% BRnam = 'JN140815002';
BRnam = 'JN140813002';

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


%% investigate different filtering parameters

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

lim = median(velF/0.6745)*5;

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

figure
ax(1) = subplot(2,1,1);
hold on
scatter(1:100000,e_x(1:100000),'.b')
scatter(1:100000,e_y(1:100000),'.r')
for k=1:find(artifact(:,2)>100000,1,'first')
    line([artifact(k,1) artifact(k,1)],ylim,'Color','g')
    line([artifact(k,2) artifact(k,2)],ylim,'Color','r')
end
ax(2) = subplot(2,1,2);
hold on
scatter(1:100000,velF(1:100000),'.')
for k=1:find(artifact(:,2)>100000,1,'first')
    line([artifact(k,1) artifact(k,1)],ylim,'Color','g')
    line([artifact(k,2) artifact(k,2)],ylim,'Color','r')
end
linkaxes(ax,'x')

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

figure;hold on
scatter(1:100000,e_x(1:100000),'.b')
scatter(1:100000,e_y(1:100000),'.r')
for k=1:find(sacarr(:,2)>100000,1,'first')
    line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
    line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
end

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
x_bound = xor(e_x<xthresh(1), e_x>xthresh(2));
y_bound = xor(e_y<ythresh(1), e_y>ythresh(2));
eyebound = xor(x_bound, y_bound);

figure
ax(1) = subplot(3,1,1);
plot(e_x,'b')
ax(2) = subplot(3,1,2);
plot(e_y,'r')
ax(3) = subplot(3,1,3);
plot(eyebound,'g')
ylim([-0.25 1.25])
linkaxes(ax,'x')

%%

decdir = 'C:\Data\MAT\NS6 - decimated\';

% load the decimated ("skipfactor") NS6 file
load(fullfile(decdir,'JN140815002_NS6_SF30.mat'))

NEV = openNEV('noread',fullfile(fildir,[BRnam '.nev']));
% for some reason, NS6h.MetaTags.DataPoints not an integer when you do it
% this way:
% NS6r = openNSx('noread',fullfile(fildir,[BRnam '.ns6'])); % r is for "raw"
NS6r = openNSx(fullfile(fildir,[BRnam '.ns6']),'c:1'); % r is for "raw"

ns2DTR = NS2.MetaTags.DateTimeRaw;
ns6DTR = NS6r.MetaTags.DateTimeRaw; % (NS6 same as NS2, but Timestamp indicates slightly delayed start time)
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6r.MetaTags.Timestamp/NS6r.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

fs2 = 1/NS2.MetaTags.SamplingFreq;
fs6 = 1/NS6r.MetaTags.SamplingFreq;

% get the timestamps to match for NS2 and NS6 with 30-sample skipfactor applied:
% Blackrock code (openNSx) behaves oddly here.
% size(NS6h.Data,2) == 30259938, number of samples in original 30 kS/s file
% first nonzero value occurs at 103 in both files:
% find(NS6r.Data(1,:)~=0,1,'first') == 103
% find(NS6.Data(1,:)~=0,1,'first') == 103
% number of nonzero samples in NS6r: 30259938-102 == 30259836
% size(NS6.Data,2) == 1008764
% number of nonzero samples in NS6: 1008764-102 == 1008662
% size(1:30:30259836,2) == 1008662, meaning that the "skipfactor" file
% starts at the first nonzero value in the original NS6 file, then takes
% every SFth sample (SF = 30 in our case) until the end
% to deal with this, we need to take the timestamps for the original NS6
% file, then to get the timestamps for the "skipfactor" file, start at the
% first nonzero sample (here, sample #103)

NS2_timestamp = datenum(ns2datevec)*86400:fs2:datenum(ns2datevec)*86400+NS2.MetaTags.DataDurationSec-fs2;

% I'm choosing to strip all the zero values off of the NS6.Data matrix
% another approach would be to include the zero values in the estimation of
% the NS6 timestamps, however this means that the first n timestamps (n =
% number of zero values) would be at a different temporal resolution than
% the rest of the timestamps, due to the way Blackrock does the skipfactor
% transformation
NS6_timestamp_dum = datenum(ns6datevec)*86400:fs6:datenum(ns6datevec)*86400+NS6r.MetaTags.DataDurationSec-fs6;
% length(NS6_timestamp_dum) == size(NS6r.Data,2)
NS6_timestamp = NS6_timestamp_dum(find(NS6r.Data(1,:)~=0,1,'first'):30:end);
NS6.Data = NS6.Data(:,find(NS6r.Data(1,:)~=0,1,'first'):end);
% length(NS6_timestamp) == size(NS6.Data,2)

%% get NS6 saccade-aligned data

% find saccade start/end in task time
NS2sacarr = [NS2_timestamp(sacarr(:,1))' NS2_timestamp(sacarr(:,2))'];

sacdat = nan(size(sacarr,1),36,501);

c=1;

t0 = clock;
ft_progress('init', 'etf',     'Please wait...');
for saclop=1:size(NS2sacarr,1)-1
    
    ft_progress(saclop/(size(NS2sacarr,1)), 'Processing event %d from %d', saclop, size(NS2sacarr,1));

    if NS2sacarr(saclop,1)>NS6_timestamp(1)+0.2 && NS2sacarr(saclop,1)<NS6_timestamp(end)
        ts1 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)-0.2); % saccade onset minus 200 ms
        ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop+1,1)); % onset of next saccade (end of fixation period)
        if length(ts1:ts2)<=501
            % take whole period if start of pre-saccade interval to end of
            % fixation period is less than 501 ms
            sacdat(c,:,1:length(ts1:ts2)) = NS6.Data(:,ts1:ts2);
        else
            % otherwise, truncate period at 300 ms after saccade onset
            ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)+0.3);
            sacdat(c,:,:) = NS6.Data(:,ts1:ts2);
        end
        c=c+1;
    end
    
end
ft_progress('close')
fprintf('total analysis time: %g\n',etime(clock,t0));

% remove unused trial holders
sacdat = sacdat(~isnan(squeeze(sacdat(:,1,1))),:,:);

trlcnt = nan(1,size(sacdat,3));
for timlop=1:size(squeeze(sacdat(1,:,:)),2)
    trlcnt(timlop) = length(find(~isnan(squeeze(sacdat(:,1,timlop)))));
end

% save(fullfile('C:\Data\MAT\sacdat\',[BRnam '_sacdat.mat']),'sacdat','NS2_timestamp','NS6_timestamp','artifact','NS2sacarr','trlcnt')


%% de-mean each LFP segment
for k=1:size(sacdat,1)
    for l=1:size(sacdat,2)
        sacdat(k,l,:) = sacdat(k,l,:)-nanmean(sacdat(k,l,:));
    end
end

%% plot some saccade-aligned LFPs

% add offset for plotting
ofs = 30;
ofv = (ofs:ofs:12*ofs)';

figure

subplot(3,1,1)
hold on

trlsel = ~isnan(squeeze(sacdat(:,1,end))); % choose trials containing data to the end of the selection window
numtrlsel = length(find(trlsel));

xaxis = -0.2:0.001:0.3;
matrix = squeeze(sacdat(trlsel,1:12,:));
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

trlsel = ~isnan(squeeze(sacdat(:,1,end))); % choose trials containing data to the end of the selection window

xaxis = -0.2:0.001:0.3;
matrix = squeeze(sacdat(trlsel,13:24,:));
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

trlsel = ~isnan(squeeze(sacdat(:,1,end))); % choose trials containing data to the end of the selection window

xaxis = -0.2:0.001:0.3;
matrix = squeeze(sacdat(trlsel,25:36,:));
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
