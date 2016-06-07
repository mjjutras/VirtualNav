fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';
BRnam = 'JN140821002';

%% Array A

% NS6 = openNSx(fullfile(fildir,'JN140821002.ns6'),'read','c:1','s:30');
% starts with sample 1 and then skips every 30, so first sample should
% align with the timestamp displayed in NS6.MetaTags.Timestamp

NS2 = openNSx(fullfile(datdir,[BRnam '.ns2']),'read');

%%
% last sample: 5790147 (after this signal changes for unknown reason)

e_x = double(NS2.Data(37,1:5790147));
e_y = double(NS2.Data(38,1:5790147));

% stdmat=[];
% for k=1:length(e_x)-3
%     stdmat=[stdmat std(e_x(k:k+3))];
% end

e_x_new = nan(size(e_x));
e_y_new = nan(size(e_y));
e_x_new(1) = e_x(1);
e_y_new(1) = e_y(1);
for k=2:length(e_x)
    
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

e_x_new = inpaint_nans(e_x_new,2);
e_y_new = inpaint_nans(e_y_new,2);


% figure;hold on
% scatter(1:10000,e_x(1:10000),'.b')
% scatter(1:10000,e_x_new(1:10000),'ob')
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

datdir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';

load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\JN140821\ChA1-12_skipfactor30.mat');
NS6a = NS6;
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\JN140821\ChB1-12_skipfactor30.mat');
NS6b = NS6;
load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\JN140821\ChC1-12_skipfactor30.mat');
NS6c = NS6;

NEV = openNEV(fullfile(datdir,[BRnam '.nev']));
NS6 = openNSx(fullfile(datdir,[BRnam '.ns6']));

ns2DTR = NS2.MetaTags.DateTimeRaw;
ns6DTR = NS6.MetaTags.DateTimeRaw; % (NS6 same as NS2, but Timestamp indicates slightly delayed start time)
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

fs = 1/NS2.MetaTags.SamplingFreq;

NS2_timestamp = datenum(ns2datevec)*86400:fs:datenum(ns2datevec)*86400+NS2.MetaTags.DataDurationSec-fs;
NS6_timestamp = datenum(ns6datevec)*86400:fs:datenum(ns6datevec)*86400+NS6.MetaTags.DataDurationSec-fs;

%% get NS6 saccade-aligned data

% find saccade start/end in task time
NS2sacarr = [NS2_timestamp(sacarr(:,1))' NS2_timestamp(sacarr(:,2))'];

sacdatA = nan(12,size(sacarr,1),501);
sacdatB = nan(12,size(sacarr,1),501);
sacdatC = nan(12,size(sacarr,1),501);

c=1;

t0 = clock;
ft_progress('init', 'etf',     'Please wait...');
for saclop=1:size(NS2sacarr,1)-1
    
    ft_progress(saclop/(size(NS2sacarr,1)), 'Processing event %d from %d', saclop, size(NS2sacarr,1));

    if NS2sacarr(saclop,1)>NS6_timestamp(1)+0.2
        ts1 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)-0.2);
        ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop+1,1));
        if length(ts1:ts2)<=501
            sacdatA(:,c,1:length(ts1:ts2)) = NS6a.Data(:,ts1:ts2);
            sacdatB(:,c,1:length(ts1:ts2)) = NS6b.Data(:,ts1:ts2);
            sacdatC(:,c,1:length(ts1:ts2)) = NS6c.Data(:,ts1:ts2);
        else
            ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)+0.3);
            sacdatA(:,c,:) = NS6a.Data(:,ts1:ts2);
            sacdatB(:,c,:) = NS6b.Data(:,ts1:ts2);
            sacdatC(:,c,:) = NS6c.Data(:,ts1:ts2);
        end
        c=c+1;
    end
    
end
ft_progress('close')
fprintf('total analysis time: %g\n',etime(clock,t0));

save('C:\Users\michael.jutras\Documents\Virtual Navigation Study\JN140821\JN140821002_sacdat141009.mat','sacdatA','sacdatB','sacdatC','NS2_timestamp','NS6_timestamp','artifact','NS2sacarr')

%% get NS2 saccade-aligned data

% find saccade start/end in task time
NS2sacarr = [NS2_timestamp(sacarr(:,1))' NS2_timestamp(sacarr(:,2))'];

sacdatA_NS2 = nan(12,size(sacarr,1),501);
sacdatB_NS2 = nan(12,size(sacarr,1),501);
sacdatC_NS2 = nan(12,size(sacarr,1),501);

c=1;

t0 = clock;
ft_progress('init', 'etf',     'Please wait...');
for saclop=1:size(NS2sacarr,1)-1
    
    ft_progress(saclop/(size(NS2sacarr,1)), 'Processing event %d from %d', saclop, size(NS2sacarr,1));

    if NS2sacarr(saclop,1)>NS2_timestamp(1)+0.2
        ts1 = ft_nearest(NS2_timestamp,NS2sacarr(saclop,1)-0.2);
        ts2 = ft_nearest(NS2_timestamp,NS2sacarr(saclop+1,1));
        if length(ts1:ts2)<=501
            sacdatA_NS2(:,c,1:length(ts1:ts2)) = NS2.Data(1:12,ts1:ts2);
            sacdatB_NS2(:,c,1:length(ts1:ts2)) = NS2.Data(13:24,ts1:ts2);
            sacdatC_NS2(:,c,1:length(ts1:ts2)) = NS2.Data(25:36,ts1:ts2);
        else
            ts2 = ft_nearest(NS2_timestamp,NS2sacarr(saclop,1)+0.3);
            sacdatA_NS2(:,c,:) = NS2.Data(1:12,ts1:ts2);
            sacdatB_NS2(:,c,:) = NS2.Data(13:24,ts1:ts2);
            sacdatC_NS2(:,c,:) = NS2.Data(25:36,ts1:ts2);
        end
        c=c+1;
    end
    
end
ft_progress('close')
fprintf('total analysis time: %g\n',etime(clock,t0));

save('C:\Users\michael.jutras\Documents\Virtual Navigation Study\JN140821\JN140821002_sacdatNS2_141009.mat','sacdatA_NS2','sacdatB_NS2','sacdatC_NS2','NS2_timestamp','artifact','NS2sacarr')



