%%
% This code calculates the saccade-aligned LFP out to the end of each fixation period,
% and uses this to derive a signal-to-noise ratio in decibels.
%
% Modifed from getSaccadeAlignedLFP_JN140815002.m
% Mike Jutras 6/19/2017

%%

fildir = 'C:\Data\Blackrock_VR';
% BRnam = 'JN140813002'; % didn't run; saccades are weird
% BRnam = 'JN140815002';
BRnam = 'JN140911002';
% BRnam = 'JN141013002';
% BRnam = 'JN141212003'; % don't use this one, from a calibration session?
% BRnam = 'JN141212004';

%% get eye data from NS2 and treat to address sampling rate issues
% (Blackrock samples at 1000 Hz vs. eyetracker 240 Hz sampling rate,
% and interpolates with repeated values between eyetracker samples, so we
% need to fix this)

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


% figure;hold on
% scatter(1:10000,e_x(1:10000),'.b')
% scatter(1:10000,e_x_new(1:10000),'ob')
% % scatter(1:10000,e_y(1:10000),'.r')
% % scatter(1:10000,e_y_new(1:10000),'or')

%% detect saccades and return sacarr marking start/end times of saccades 

fltord = 40;
lowpasfrq = 40;
% lim = 120000; % replaced with line below
artpadding = 0.01;
Fs = NS2.MetaTags.SamplingFreq; % sampling rate in Hz

%low pass filter the eye position data_eye
nyqfrq = Fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
e_x_lowpassfilter=filtfilt(flt,1, e_x);
e_y_lowpassfilter=filtfilt(flt,1, e_y);
e_x_lowpassfilter_new=filtfilt(flt,1, e_x_new);
e_y_lowpassfilter_new=filtfilt(flt,1, e_y_new);

% differentiate and multiply with sampling rate to get velocity as deg/sec
x_vF = diff(e_x_lowpassfilter) .* Fs;
y_vF = diff(e_y_lowpassfilter) .* Fs;
x_vF_new = diff(e_x_lowpassfilter_new) .* Fs;
y_vF_new = diff(e_y_lowpassfilter_new) .* Fs;

% differentiate and multiply with sampling rate without filtering
% this gives the eye velocity in the horizontal and vertical domains
x_v = diff(e_x) .* Fs;
y_v = diff(e_y) .* Fs;
x_v_new = diff(e_x_new) .* Fs;
y_v_new = diff(e_y_new) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
vel = abs(complex(x_v,y_v));
velF = abs(complex(x_vF,y_vF));
vel_new = abs(complex(x_v_new,y_v_new));
velF_new = abs(complex(x_vF_new,y_vF_new));

lim = median(velF/0.6745)*5; % threshold taken from Quiroga (http://www.scholarpedia.org/article/Spike_sorting)

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

% merge artifacts when inter-artifact interval is less than 20 ms
% DON'T DO THIS for this analysis, since it messes up the evoked response
% to saccade onset
% instead, later on select the saccades that are not followed or preceded
% by another saccade ("artifact") within some window
sacarr = [];
% iai=(artifact(2:end,1)-artifact(1:end-1,2))>20;
iai=(artifact(2:end,1)-artifact(1:end-1,2))>0;
sacdumcol1 = artifact(2:end,1);
sacdumcol2 = artifact(1:end-1,2);
sacarr(:,1) = [artifact(1,1); sacdumcol1(iai,1)];
sacarr(:,2) = [sacdumcol2(iai,1); artifact(end,2)];

% inspect to make sure saccades make sense
figure;hold on
scatter(1:100000,e_x(1:100000),'.b')
scatter(1:100000,e_x_new(1:100000),'ob')
scatter(1:100000,e_y(1:100000),'.r')
scatter(1:100000,e_y_new(1:100000),'or')
for k=1:find(sacarr(:,2)>100000,1,'first')
    line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g')
    line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r')
end

% clear vel* x_* y_* e_*

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
% 
% sacsel1 = zeros(size(sacarr,1),1);
% for saclop = 1:size(sacarr,1)
%     if saclop~=size(sacarr,1)
%         if isempty(find(eyebound(sacarr(saclop,2)+1:sacarr(saclop+1,1)-1),1))
%             sacsel1(saclop) = 1;
%         end
%     else
%         if isempty(find(eyebound(sacarr(saclop,2)+1:end),1))
%             sacsel1(saclop) = 1;
%         end
%     end
% end
% 
% sacsel1 = find(sacsel1);

%% choose saccades that don't precede fixations OR FOLLOW at edges(borders)
% usually comment this out in preference for the previous cell, since the ERP
% obtained when running this cell can be dampened compared to that
% obtained when only excluding saccade that follow gaze at borders
% however in the presence of lots of saccade/fixations this may be
% appropriate

sacsel1 = zeros(size(sacarr,1),1);
for saclop = 1:size(sacarr,1)
    if saclop~=size(sacarr,1)
        if saclop==1
            eyeboundPre = eyebound(1:sacarr(saclop,1)-1);
        else
            eyeboundPre = eyebound(sacarr(saclop-1,2)+1:sacarr(saclop,1)-1);
        end
        eyeboundPost = eyebound(sacarr(saclop,2)+1:sacarr(saclop+1,1)-1);
    else
        eyeboundPre = eyebound(sacarr(saclop-1,2)+1:sacarr(saclop,1)-1);
        eyeboundPost = eyebound(sacarr(saclop,2)+1:end);
    end
    if isempty(find(eyeboundPre,1)) && isempty(find(eyeboundPost,1))
        sacsel1(saclop) = 1;
    end
end

sacsel1 = find(sacsel1);

%% further select saccades that are not followed or preceded by another saccade
% within a 150(/lagwin) ms window

lagwin = 20;

sacsel2 = zeros(size(sacarr,1),1);
for saclop = 1:size(sacarr,1)
    if saclop==1
        if sacarr(saclop,1)>=lagwin
            sacsel2(saclop) = 1;
        end
    elseif saclop==size(sacarr,1)
        if size(NS2.Data,2)-sacarr(saclop,2)>=lagwin
            sacsel2(saclop) = 1;
        end
    elseif diff(sacarr(saclop-1:saclop+1,1))>=lagwin
        sacsel2(saclop) = 1;
    end
end

sacsel2 = find(sacsel2);

% combine sacsel1 & sacsel2
sacsel = intersect(sacsel1,sacsel2);

%% load decimated NS6 neural data and sync with NS2
% (grrr Blackrock)

decdir = 'C:\Data\MAT\NS6 - decimated\';

% load the decimated ("skipfactor") NS6 file
load(fullfile(decdir,[BRnam '_NS6_SF30.mat']))

NEV = openNEV('noread',fullfile(fildir,[BRnam '.nev']));
NS6r = openNSx(fullfile(fildir,[BRnam '.ns6']),'c:1'); % r is for "raw"

ns2DTR = NS2.MetaTags.DateTimeRaw;
ns6DTR = NS6r.MetaTags.DateTimeRaw;
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6r.MetaTags.Timestamp/NS6r.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

fs2 = 1/NS2.MetaTags.SamplingFreq;
fs6 = 1/NS6r.MetaTags.SamplingFreq;

% get the timestamps to match for NS2 and NS6 with 30-sample skipfactor applied
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

% this next step revised from getSaccadeAlignedLFP_JN140815002
% since updating Blackrock NPMK some time between 11/8/2016 and 2/6/2017,
% openNSx with a skipfactor performs differently... now there are no longer
% the same number of zeros at the beginning of NS6 (decimated) and NS6r
% (raw), but the first non-zero value is the same for each file.
NS6_timestamp_dum = datenum(ns6datevec)*86400:fs6:datenum(ns6datevec)*86400+NS6r.MetaTags.DataDurationSec-fs6;
NS6_timestamp = NS6_timestamp_dum(find(NS6r.Data(1,:)~=0,1,'first'):30:end);
% simple enough fix to the above: change 'NS6r' in the following line to 'NS6'
% NS6.Data = NS6.Data(:,find(NS6r.Data(1,:)~=0,1,'first'):end);
NS6.Data = NS6.Data(:,find(NS6.Data(1,:)~=0,1,'first'):end);

% fix JN141212004, recording is a mess at the very end
if strmatch('JN141212004',BRnam)
    NS6.Data = NS6.Data(:,NS6_timestamp<=63585703122);
    NS6_timestamp = NS6_timestamp(NS6_timestamp<=63585703122);
end

%% get NS6 saccade-aligned data
% revised this cell: changed sacdat to a cell array and changed loop to a
% parfor loop, decreasing analysis time significantly and including all
% data (not cutoff at 300 ms after each saccade)

% 0 to include all saccades; 1 to exclude saccades to outside of eyetracker range (border)
noBorder = 1;

% find saccade start/end in task time
NS2sacarr = [NS2_timestamp(sacarr(:,1))' NS2_timestamp(sacarr(:,2))'];

t0 = clock;
if noBorder==0

    % sacdat as a cell array:
    % don't include the last fixation period of the session; only include
    % fixations that terminate in an observed saccade
    sacdat = cell(size(NS2sacarr,1)-1,1);

    parfor saclop=1:size(NS2sacarr,1)-1

        % due to this if statement, it is possible there will be empty cells at
        % the beginning and the end of sacdat
        if NS2sacarr(saclop,1)>NS6_timestamp(1)+0.2 && NS2sacarr(saclop,1)<NS6_timestamp(end)

            % start segment at saccade onset minus 200 ms
            ts1 = ft_nearest(NS6_timestamp,NS2sacarr(saclop,1)-0.2);

            % end segment at onset of next saccade (end of fixation period)
            ts2 = ft_nearest(NS6_timestamp,NS2sacarr(saclop+1,1));

            sacdat{saclop} = NS6.Data(:,ts1:ts2);

            % demean the data in this step
            sacdat{saclop} = sacdat{saclop}-int16(repmat(mean(sacdat{saclop},2),1,size(sacdat{saclop},2)));

        end

    end
    
elseif noBorder==1
    
    sacdat = cell(length(sacsel),1);
    
    parfor saclop=1:length(sacsel)
        
        if sacsel(saclop)~=size(sacarr,1) % don't use the very last saccade of the recording
            if NS2sacarr(sacsel(saclop),1)>NS6_timestamp(1)+0.2 && NS2sacarr(sacsel(saclop),1)<NS6_timestamp(end)
                
                % start segment at saccade onset minus 200 ms
                ts1 = ft_nearest(NS6_timestamp,NS2sacarr(sacsel(saclop),1)-0.2);
                
                % end segment at onset of next saccade (end of fixation period)
                ts2 = ft_nearest(NS6_timestamp,NS2sacarr(sacsel(saclop)+1,1));

                
                sacdat{saclop} = NS6.Data(:,ts1:ts2);
                
                % demean the data in this step
                sacdat{saclop} = sacdat{saclop}-int16(repmat(mean(sacdat{saclop},2),1,size(sacdat{saclop},2)));
                
            end
        end
        
    end
    
end
fprintf('total analysis time: %g\n',etime(clock,t0));


%% look at data

chnsel = 13;

trlsiz = nan(size(sacdat));
for k=1:length(sacdat)
    trlsiz(k) = size(sacdat{k},2);
end

datchnsel = nan(length(sacdat),max(trlsiz));
for k=1:length(sacdat)
    if ~isempty(sacdat{k})
        datchnsel(k,1:size(sacdat{k},2)) = sacdat{k}(chnsel,:);
    end
end

figure;plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel,1)))
hold on;errorshade(datchnsel,1,(1:size(datchnsel,2))-200,'b')

figure; hold on
plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(isnan(datchnsel(:,350)),:),1)))
plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(~isnan(datchnsel(:,350))&isnan(datchnsel(:,400)),:),1)),'r')
plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(~isnan(datchnsel(:,400))&isnan(datchnsel(:,450)),:),1)),'g')
plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(~isnan(datchnsel(:,450))&isnan(datchnsel(:,500)),:),1)),'m')
plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(~isnan(datchnsel(:,500))&isnan(datchnsel(:,600)),:),1)),'c')
plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(~isnan(datchnsel(:,600))&isnan(datchnsel(:,700)),:),1)),'k')

% figure;hold on
% plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(1:10,:),1)))
% plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(11:20,:),1)),'r')
% plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(21:30,:),1)),'g')
% plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(31:40,:),1)),'m')
% plot((1:size(datchnsel,2))-200,squeeze(nanmean(datchnsel(41:50,:),1)),'k')


%% calculate signal to noise

% designate "signal" and "noise" windows in units of sample after saccade
% start (samples correspond to ms with 1 kHz sampling rate)
% following are "signal" and "noise" windows with calculated SNR (JN140815002) in comment
% the following all taken from all saccades (not excluding saccades to/from
% outside eyetracker window)
% winS = [56 155]; winN = [-155 -56]; % -0.501
% winS = [56 155]; winN = [400 499]; % -0.300
% winS = [70 139]; winN = [-89 -20]; % 0.267
% winS = [80 139]; winN = [-79 -20]; % 0.400
% winS = [85 134]; winN = [-80 -31]; % 0.358
% winS = [85 134]; winN = [-64 -15]; % 0.609

% after excluding saccades to/from outside eyetracker window:
% JN140815002 (n = 1027)
% winS = [85 134]; winN = [-64 -15]; % 0.234
% winS = [85 134]; winN = [-94 -45]; % -0.319

% JN140911002
% winS = [60 109]; winN = [-64 -15]; % 0.327
% winS = [60 109]; winN = [-94 -45]; % 0.035
winS = [70 119]; winN = [-84 -35]; % 0.035

% JN141013002
% winS = [75 124]; winN = [-94 -45]; % 0.536

% JN141212004
% winS = [65 114]; winN = [-94 -45]; % 0.246
% winS = [65 114]; winN = [-64 -15]; % 0.602

trltim = (1:max(trlsiz))-200;
selS = find(trltim==winS(1)):find(trltim==winS(2));
selN = find(trltim==winN(1)):find(trltim==winN(2));

dattimsel = datchnsel(~isnan(datchnsel(:,max([selS selN]))),:);

close all
figure;hold on
plot((1:max([selS selN]))-200,squeeze(nanmean(dattimsel(:,1:max([selS selN])),1)))
errorshade(dattimsel(:,1:max([selS selN])),1,(1:max([selS selN]))-200,'b')
% shade "signal" and "noise" windows
yh=ylim;
fill([winS(1):winS(2) fliplr(winS(1):winS(2))],[ones(size(winS(1):winS(2)))*yh(1) ones(size(winS(1):winS(2)))*yh(2)],'g', 'FaceAlpha',.2, 'EdgeColor','g','EdgeAlpha',0);
fill([winN(1):winN(2) fliplr(winN(1):winN(2))],[ones(size(winN(1):winN(2)))*yh(1) ones(size(winN(1):winN(2)))*yh(2)],'r', 'FaceAlpha',.2, 'EdgeColor','r','EdgeAlpha',0);

% store these figures in R:\Mike\VirtualNavigationProject\Figures\Flexshaft_firstGizImplant2014\SaccadeAligned_SigToNoise

% make sure length(selS) and length(selN) are equal
if length(selS) == length(selN)
    
    valS = max(abs(dattimsel(:,selS)),[],2);
    valN = max(abs(dattimsel(:,selN)),[],2);
    
    Ps = mean(valS.^2);
    Pn = mean(valN.^2);
    
    10*log10(Ps/Pn)
    
else
    
    fprintf('selS and selN are not the same length\n')
    
end

   