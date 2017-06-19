% fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';
fildir = 'C:\Data\Blackrock_VR';
BRnam = 'JN140815002';

%% Array A

% NS6 = openNSx(fullfile(fildir,'JN140821002.ns6'),'read','c:1','s:30');
% starts with sample 1 and then skips every 30, so first sample should
% align with the timestamp displayed in NS6.MetaTags.Timestamp

NS2 = openNSx(fullfile(fildir,[BRnam '.ns2']),'read');

%%

e_x = double(NS2.Data(1,:));
e_y = double(NS2.Data(2,:));

% stdmat=[];
% for k=1:length(e_x)-3
%     stdmat=[stdmat std(e_x(k:k+3))];
% end

% redraw the eyetracker data, to replace data points that repeat (due to
% higher sampling rate in Blackrock compared to eyetracker) with
% interpolated data points
% *this may not be strictly necessary since the signal gets filtered for
% saccade detection
e_x_new = nan(size(e_x));
e_y_new = nan(size(e_y));
e_x_new(1) = e_x(1);
e_y_new(1) = e_y(1);
for k=2:length(e_x)
    
    % repeated data points seem to stay within 150 of the initial data point
    % however, actual deviations in measurements can stay within 150, so
    % this won't perform perfectly when the eye is relatively still
    % but this shouldn't matter since what we care about is saccade
    % detection
    if abs(e_x(k)-e_x(k-1))>150
        e_x_new(k)=e_x(k);
    end
    if k>=4 && isempty(find(~isnan(e_x_new(k-3:k)),1))
        e_x_new(k)=e_x(k);
    end
    
    if abs(e_y(k)-e_y(k-1))>150
        e_y_new(k)=e_y(k);
    end
    if k>=4 && isempty(find(~isnan(e_y_new(k-3:k)),1))
        e_y_new(k)=e_y(k);
    end

end

e_x_new(end) = e_x(end);
e_y_new(end) = e_y(end);

% fill in nan values with interpolation
e_x_new = inpaint_nans(e_x_new,2);
e_y_new = inpaint_nans(e_y_new,2);


figure;hold on
scatter(1:10000,e_x(1:10000),'.b')
scatter(1:10000,e_x_new(1:10000),'ob')
% scatter(1:10000,e_y(1:10000),'.r')
% scatter(1:10000,e_y_new(1:10000),'or')

%%
fltord = 40;
lowpasfrq = 80;
lim = 120000;
artpadding = 0.01;
Fs = NS2.MetaTags.SamplingFreq; % sampling rate in Hz


%low pass filter the eye position data_eye
nyqfrq = Fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
e_x_lowpassfilter=filtfilt(flt,1, e_x);
e_y_lowpassfilter=filtfilt(flt,1, e_y);
e_x_lowpassfilter_new=filtfilt(flt,1, e_x_new);
e_y_lowpassfilter_new=filtfilt(flt,1, e_y_new);

%differentiate and multiply with sampling rate to get velocity as deg/sec
x_vF = diff(e_x_lowpassfilter) .* Fs;
y_vF = diff(e_y_lowpassfilter) .* Fs;
x_vF_new = diff(e_x_lowpassfilter_new) .* Fs;
y_vF_new = diff(e_y_lowpassfilter_new) .* Fs;

% % differentiate and multiply with sampling rate without filtering
% % this gives the eye velocity in the horizontal and vertical domains
x_v = diff(e_x) .* Fs;
y_v = diff(e_y) .* Fs;
x_v_new = diff(e_x_new) .* Fs;
y_v_new = diff(e_y_new) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
vel = abs(complex(x_v,y_v));
velF = abs(complex(x_vF,y_vF));
vel_new = abs(complex(x_v_new,y_v_new));
velF_new = abs(complex(x_vF_new,y_vF_new));

% figure;hold on
% plot(vel(1:100000))
% plot(velF(1:100000),'r')
% plot(vel_new(1:100000),'m')
% plot(velF_new(1:100000),'g')
% line(xlim,[lim lim],'Color','r')

%detect saccade begins and saccade ends
sacbeg = find(diff(velF_new > lim) > 0);
sacend = find(diff(velF_new > lim) < 0);
if velF_new(end)>lim
    sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
end
if velF_new(1)>lim
    sacend = sacend(2:end);
end

if size(sacbeg,1)
    sacbeg=sacbeg';
    sacend=sacend';
end
    
% artifact = round([sacbeg(:) - artpadding*Fs sacend(:) + artpadding*Fs]);
artifact = round([sacbeg(:) sacend(:)]);

% figure;hold on
% scatter(1:10000,e_x(1:10000),'.b')
% scatter(1:10000,e_x_new(1:10000),'ob')
% scatter(1:10000,e_y(1:10000),'.r')
% scatter(1:10000,e_y_new(1:10000),'or')
% for k=1:find(artifact(:,2)>10000,1,'first')
%     line([artifact(k,1) artifact(k,1)],ylim,'Color','g')
%     line([artifact(k,2) artifact(k,2)],ylim,'Color','r')
% end

% h=[];
% for k=1:size(artifact,1)
%     h=[h abs(diff(e_x(artifact(k,1):artifact(k,2))))];
% end
    

% % find saccade start/end in task time
% sacdum = [NS2_timestamp(artifact(:,1))' NS2_timestamp(artifact(:,2))'];

% merge artifacts when inter-artifact interval is less than 10 ms
sacarr = [];
iai=(artifact(2:end,1)-artifact(1:end-1,2))>10;
sacdumcol1 = artifact(2:end,1);
sacdumcol2 = artifact(1:end-1,2);
sacarr(:,1) = [artifact(1,1); sacdumcol1(iai,1)];
sacarr(:,2) = [sacdumcol2(iai,1); artifact(end,2)];

figure;hold on
scatter(1:10000,e_x(1:10000),'.b')
scatter(1:10000,e_x_new(1:10000),'ob')
scatter(1:10000,e_y(1:10000),'.r')
scatter(1:10000,e_y_new(1:10000),'or')
for k=1:find(sacarr(:,2)>10000,1,'first')
    line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
    line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
end

clear vel* x_* y_* e_*

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

%% load pre-saved variables

BRnam = 'JN140815002';

sacdat = load(fullfile('C:\Data\MAT\SaccadeAligned\', [BRnam '_sacdatNS6.mat']));
sacdat = sacdat.sacdatNS6;

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

% figdir = 'R:\Mike\VirtualNavigationProject\Figures\Flexshaft_firstGizImplant2014\SaccadeAligned_SigToNoise';
