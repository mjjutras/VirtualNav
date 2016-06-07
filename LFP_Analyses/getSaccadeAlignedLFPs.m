BRnam = 'JN140620002';
pandata = '14_06_20_13_58';

%%
datdir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';

NS2 = openNSx(fullfile(datdir,[BRnam '.ns2']),'read');
% NS6 = openNSx(fullfile(datdir,[BRnam '.ns6']),'read','skipfactor',30);
NEV = openNEV(fullfile(datdir,[BRnam '.nev']),'read');

ns2DTR = NS2.MetaTags.DateTimeRaw;
% ns6DTR = NS6.MetaTags.DateTimeRaw;
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
% ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

fs = 1/NS2.MetaTags.SamplingFreq;

NS2_timestamp = datenum(ns2datevec)*86400:fs:datenum(ns2datevec)*86400+NS2.MetaTags.DataDurationSec-fs;

%%
fltord = 40;
lowpasfrq = 40;
lim = 70000;
artpadding = 0.01;
Fs = 1/fs; % sampling rate in Hz

% data_eye.artifact = cell(size(data_eye.time));
e_x = double(NS2.Data(1,:));
e_y = double(NS2.Data(2,:));

%low pass filter the eye position data_eye
nyqfrq = Fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
e_x_lowpassfilter=filtfilt(flt,1, e_x);
e_y_lowpassfilter=filtfilt(flt,1, e_y);

%differentiate and multiply with sampling rate to get velocity as deg/sec
x_vF = diff(e_x_lowpassfilter) .* Fs;
y_vF = diff(e_y_lowpassfilter) .* Fs;

% % differentiate and multiply with sampling rate without filtering
% % this gives the eye velocity in the horizontal and vertical domains
x_v = diff(e_x) .* Fs;
y_v = diff(e_y) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
vel = abs(complex(x_v,y_v));
velF = abs(complex(x_vF,y_vF));

%detect saccade begins and saccade ends
sacbeg = find(diff(velF > lim) > 0);
sacend = find(diff(velF > lim) < 0);
if velF(end)>lim
    sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
end
if velF(1)>lim
    sacend = sacend(2:end);
end

if size(sacbeg,1)
    sacbeg=sacbeg';
    sacend=sacend';
end
    
artifact = round([sacbeg(:) - artpadding*Fs sacend(:) + artpadding*Fs]);

% find saccade start/end in task time
sacdum = [NS2_timestamp(artifact(:,1))' NS2_timestamp(artifact(:,2))'];

% merge artifacts when inter-artifact interval is less than 10 ms
sacarr = [];
iai=(sacdum(2:end,1)-sacdum(1:end-1,2))>0.01; sacdumcol1=sacdum(2:end,1); sacdumcol2=sacdum(1:end-1,2);
sacarr(:,1) = [sacdum(1,1); sacdumcol1(iai,1)];
sacarr(:,2) = [sacdumcol2(iai,1); sacdum(end,2)];

% figure;plot(NS2_timestamp(1:10000)-min(NS2_timestamp),velF(1:10000));hold on;for k=1:40;line([sacarr(k,1)-min(NS2_timestamp) sacarr(k,1)-min(NS2_timestamp)],ylim,'Color','g');line([sacarr(k,2)-min(NS2_timestamp) sacarr(k,2)-min(NS2_timestamp)],ylim,'Color','r');end

%%
clear vel* x_* y_* e_*

parsedir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\parsed continuous data';

clear f1 f2 f3 f4

load(fullfile(parsedir, BRnam,'chan1', 'data_chan1.mat'));
load(fullfile(parsedir, BRnam,'chan1', 'timestamp_chan1.mat'));
f1.data = data; f1.timestamp = timestamp;
load(fullfile(parsedir, BRnam,'chan2', 'data_chan2.mat'));
load(fullfile(parsedir, BRnam,'chan2', 'timestamp_chan2.mat'));
f2.data = data; f2.timestamp = timestamp;
load(fullfile(parsedir, BRnam,'chan3', 'data_chan3.mat'));
load(fullfile(parsedir, BRnam,'chan3', 'timestamp_chan3.mat'));
f3.data = data; f3.timestamp = timestamp;
load(fullfile(parsedir, BRnam,'chan4', 'data_chan4.mat'));
load(fullfile(parsedir, BRnam,'chan4', 'timestamp_chan4.mat'));
f4.data = data; f4.timestamp = timestamp;

clear data timestamp

sacdat_f1 = nan(size(sacarr,1),501);
sacdat_f2 = nan(size(sacarr,1),501);
sacdat_f3 = nan(size(sacarr,1),501);
sacdat_f4 = nan(size(sacarr,1),501);

c1=1; c2=1; c3=1; c4=1;

t0 = clock;
ft_progress('init', 'etf',     'Please wait...');
for saclop=1:size(sacarr,1)-1
    
    ft_progress(saclop/(size(sacarr,1)), 'Processing event %d from %d', saclop, size(sacarr,1));

    if sacarr(saclop,1)>f1.timestamp(1)+0.2
        ts1 = ft_nearest(f1.timestamp,sacarr(saclop,1)-0.2);
        ts2 = ft_nearest(f1.timestamp,sacarr(saclop+1,1));
        if length(ts1:30:ts2)<=501
            sacdat_f1(c1,1:length(ts1:30:ts2)) = f1.data(ts1:30:ts2);
        else
            ts2 = ft_nearest(f1.timestamp,sacarr(saclop,1)+0.3);
            sacdat_f1(c1,:) = f1.data(ts1:30:ts2);
        end
        c1=c1+1;
    end
%     if sacarr(saclop,1)>f2.timestamp(1)+0.2
%         ts1 = ft_nearest(f2.timestamp,sacarr(saclop,1)-0.25);
%         ts2 = ft_nearest(f2.timestamp,sacarr(saclop,1)+0.25);
%         sacdat_f2(c2,:) = f2.data(ts1:30:ts2);
%         c2=c2+1;
%     end
%     if sacarr(saclop,1)>f3.timestamp(1)+0.2
%         ts1 = ft_nearest(f3.timestamp,sacarr(saclop,1)-0.25);
%         ts2 = ft_nearest(f3.timestamp,sacarr(saclop,1)+0.25);
%         sacdat_f3(c3,:) = f3.data(ts1:30:ts2);
%         c3=c3+1;
%     end
%     if sacarr(saclop,1)>f4.timestamp(1)+0.2
%         ts1 = ft_nearest(f4.timestamp,sacarr(saclop,1)-0.25);
%         ts2 = ft_nearest(f4.timestamp,sacarr(saclop,1)+0.25);
%         sacdat_f4(c4,:) = f4.data(ts1:30:ts2);
%         c4=c4+1;
%     end
    
end
ft_progress('close')
fprintf('total analysis time: %g\n',etime(clock,t0));

save('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\saccade-aligned\JN140620002_c1_sacdat.mat','sacdat_f1')



