% BRnam = 'JN160405002';
% BRnam = '160406001';
% BRnam = '160406002';
% BRnam = 'jn160407002';
% BRnam = 'JN160412001';
% BRnam = 'JN160412002';
% BRnam = 'JN160414001';
% BRnam = 'JN160414002';
% BRnam = 'JN160414003';
% BRnam = 'JN160418003';
% BRnam = 'JN160421003';


%%
if str2num(BRnam(3:4))==14
    datdir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
elseif str2num(BRnam(3:4))==16
    datdir = 'R:\Buffalo Lab\Virtual Navigation\UnityVR\NeuralData\';
end

NS2 = openNSx(fullfile(datdir,[BRnam '.ns2']),'read');

fs = 1/NS2.MetaTags.SamplingFreq;

% JN160405002: include one of the following lines
% NS2.Data = NS2.Data{1}(:,1:1813522); % task 1
% NS2.Data = NS2.Data{1}(:,1999352:3126110); % task 2
% NS2.Data = NS2.Data{1}(:,3349671:end); % task 3

% jn160407002: include one of the following lines
% NS2.Data = NS2.Data(:,1:2016756); % task 1
% NS2.Data = NS2.Data(:,2431085:4564847); % task 2
% NS2.Data = NS2.Data(:,4625860:end); % task 3


%%
fltord = 40;
lowpasfrq = 40;
lim = 70000;
artpadding = 0.01;
Fs = 1/fs; % sampling rate in Hz

clear xind yind
for k=1:length(NS2.ElectrodesInfo)
    if strmatch('Eye1',NS2.ElectrodesInfo(k).Label)
        xind = k;
    elseif strmatch('Eye2',NS2.ElectrodesInfo(k).Label)
        yind = k;
    end
end

% data_eye.artifact = cell(size(data_eye.time));
e_x = double(NS2.Data(xind,:));
e_y = double(NS2.Data(yind,:));
% e_x = double(NS2.Data{1}(xind,:));
% e_y = double(NS2.Data{1}(yind,:));

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

% % find saccade start/end in task time
% sacdum = [NS2_timestamp(artifact(:,1))' NS2_timestamp(artifact(:,2))'];

% merge artifacts when inter-artifact interval is less than 10 ms
sacarr = [];
iai=(artifact(2:end,1)-artifact(1:end-1,2))>0.01; sacdumcol1=artifact(2:end,1); sacdumcol2=artifact(1:end-1,2);
sacarr(:,1) = [artifact(1,1); sacdumcol1(iai,1)];
sacarr(:,2) = [sacdumcol2(iai,1); artifact(end,2)];

% figure;plot(NS2_timestamp(1:10000)-min(NS2_timestamp),velF(1:10000));hold on;for k=1:40;line([sacarr(k,1)-min(NS2_timestamp) sacarr(k,1)-min(NS2_timestamp)],ylim,'Color','g');line([sacarr(k,2)-min(NS2_timestamp) sacarr(k,2)-min(NS2_timestamp)],ylim,'Color','r');end


%%

indD = [];
indE = [];
indF = [];
indG = [];
indH = [];
indI = [];
for k=1:length(NS2.ElectrodesInfo)
    
    if ~isempty(strmatch('D',NS2.ElectrodesInfo(k).Label))
        indD = [indD; k str2num(NS2.ElectrodesInfo(k).Label(2:end))];
    elseif ~isempty(strmatch('E',NS2.ElectrodesInfo(k).Label)) && isempty(strmatch('Eye',NS2.ElectrodesInfo(k).Label))
        indE = [indE; k str2num(NS2.ElectrodesInfo(k).Label(2:end))];
    elseif ~isempty(strmatch('F',NS2.ElectrodesInfo(k).Label))
        indF = [indF; k str2num(NS2.ElectrodesInfo(k).Label(2:end))];
    elseif ~isempty(strmatch('G',NS2.ElectrodesInfo(k).Label))
        indG = [indG; k str2num(NS2.ElectrodesInfo(k).Label(2:end))];
    elseif ~isempty(strmatch('H',NS2.ElectrodesInfo(k).Label))
        indH = [indH; k str2num(NS2.ElectrodesInfo(k).Label(2:end))];
    elseif ~isempty(strmatch('I',NS2.ElectrodesInfo(k).Label))
        indI = [indI; k str2num(NS2.ElectrodesInfo(k).Label(2:end))];
    end
    
end

[~,sD] = sort(indD(:,2));
[~,sE] = sort(indE(:,2));
[~,sF] = sort(indF(:,2));
[~,sG] = sort(indG(:,2));
[~,sH] = sort(indH(:,2));
[~,sI] = sort(indI(:,2));

indD = indD(sD,:);
indE = indE(sE,:);
indF = indF(sF,:);
indG = indG(sG,:);
indH = indH(sH,:);
indI = indI(sI,:);

%%
clear e_* vel* x_* y_*
pack

fixint = diff(sacarr(:,1));

fixThresh = 300;

% sacdat = nan(96,size(sacarr,1),501); % all fixation periods
sacdat = nan(96,length(find(fixint>fixThresh)),701); % OPTION 1: only >fixThresh (ms) fixation periods, -200:500 (nans when fixation period ends before 500ms)
% sacdat = nan(96,length(find(fixint>fixThresh)),2001); % OPTION 2: only >fixThresh (ms) fixation periods, -1000:1000 (keep saccade info in trialinfo, otherwise data are padded with continuous signal on either side for TF analysis)

c=1;

sampleinfo = [];
trialinfo = [];

for saclop=1:size(sacarr,1)-1
    
    if sacarr(saclop+1,1)-sacarr(saclop,1)>fixThresh && sacarr(saclop,1)>fixThresh % OPTION 1
%     if sacarr(saclop+1,1)-sacarr(saclop,1)>fixThresh && sacarr(saclop,1)>1000 && sacarr(saclop)+1000<size(NS2.Data,2) % OPTION 2

        trialinfo(c,:) = [sacarr(saclop,1) sacarr(saclop+1,1)];
        
        timsel = sacarr(saclop,1)-200:min([sacarr(saclop,1)+500 sacarr(saclop+1,1)-1]); % OPTION 1
%         timsel = sacarr(saclop,1)-1000:sacarr(saclop,1)+1000; % OPTION 2

        sacdat(:,c,1:length(timsel)) = double(NS2.Data(1:96,timsel))-repmat(nanmean(NS2.Data(1:96,timsel),2),1,length(timsel));
        sampleinfo(c,:) = [timsel(1) timsel(end)];
        c = c+1;
        
    end
    
end

sacdat = sacdat(:,~isnan(squeeze(sacdat(1,:,1))),:);

% save('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\saccade-aligned\JN140620002_c1_sacdat.mat','sacdat_f1')

% calculate average saccade-related potential
% remove outlying events (based on variance in 60 Hz bandpass filtered signal)
Fbp=[58 62];
N=4;
Fn=1000/2;
[B,A]=butter(N, [min(Fbp)/Fn max(Fbp)/Fn]);

varall = nan(size(sacdat,1),size(sacdat,2));
for trllop = 1:size(sacdat,2)
    sacdatflt = (filtfilt(B,A,squeeze(sacdat(:,trllop,1:find(~isnan(sacdat(1,trllop,:)),1,'last')))'))';
    varall(:,trllop) = var(sacdatflt,0,2,'omitnan');
end
w_t = 4; % weight for trials
% find outlying trials using max variance for each trial (across channels)
trlthresh1 = quantile(max(varall),0.75) + w_t*(quantile(max(varall),0.75) - quantile(max(varall),0.25));
trlselvardum = find(max(varall)>trlthresh1);
trlsel = setxor(1:size(varall,2),trlselvardum);

sampleinfo = sampleinfo(trlsel,:);
trialinfo = trialinfo(trlsel,:);


% srpD = squeeze(nanmean(sacdat(indD(:,1),:,:),2));
% srpE = squeeze(nanmean(sacdat(indE(:,1),:,:),2));
% srpF = squeeze(nanmean(sacdat(indF(:,1),:,:),2));
% srpG = squeeze(nanmean(sacdat(indG(:,1),:,:),2));
% srpH = squeeze(nanmean(sacdat(indH(:,1),:,:),2));
% srpI = squeeze(nanmean(sacdat(indI(:,1),:,:),2));

sacdatD = sacdat(indD(:,1),trlsel,:);
srpD = squeeze(nanmean(sacdatD,2));
sacdatE = sacdat(indE(:,1),trlsel,:);
srpE = squeeze(nanmean(sacdatE,2));
sacdatF = sacdat(indF(:,1),trlsel,:);
srpF = squeeze(nanmean(sacdatF,2));
sacdatG = sacdat(indG(:,1),trlsel,:);
srpG = squeeze(nanmean(sacdatG,2));
sacdatH = sacdat(indH(:,1),trlsel,:);
srpH = squeeze(nanmean(sacdatH,2));
sacdatI = sacdat(indI(:,1),trlsel,:);
srpI = squeeze(nanmean(sacdatI,2));


%%

figure
subplot(6,1,1); plot(-200:500,srpD); xlim([-200 250]); hold on; line([0 0],ylim); title('D')
subplot(6,1,2); plot(-200:500,srpE); xlim([-200 250]); hold on; line([0 0],ylim); title('E')
subplot(6,1,3); plot(-200:500,srpF); xlim([-200 250]); hold on; line([0 0],ylim); title('F')
subplot(6,1,4); plot(-200:500,srpG); xlim([-200 250]); hold on; line([0 0],ylim); title('G')
subplot(6,1,5); plot(-200:500,srpH); xlim([-200 250]); hold on; line([0 0],ylim); title('H')
subplot(6,1,6); plot(-200:500,srpI); xlim([-200 250]); hold on; line([0 0],ylim); title('I')

toptitle(BRnam)

toptitle([BRnam ' W maze'])
toptitle([BRnam ' scene manip'])
toptitle([BRnam ' foraging'])


%% CSDs

doTheta = 0;
doGamma = 0;

keepsplines = 1;

bad = [];

CSDmethod = 'iCSD_splines';
% CSDmethod = 'standard';

avgForBad = 1; % should definitely do this instead of leaving them out!

dttrl = 'yes'; % detrend the trials, should probably do this too
% dttrl = 'no';


goodset = 1:16;
good = setdiff(goodset,bad);
labels = {'AD01','AD02','AD03','AD04','AD05','AD06','AD07','AD08','AD09','AD10','AD11','AD12','AD13','AD14','AD15','AD16'};

data= [];
c=1;
for trllop = 1:size(sacdatG,2)
    if isempty(find(isnan(squeeze(sacdatG(1,trllop,1:451))), 1))
        data.trial{c} = squeeze(sacdatG(good,trllop,1:451));
        data.time{c} = -0.200:0.001:0.25;
        c=c+1;
    end
end

lfpind = 1:16;

data.label = labels(good);

if doTheta
    %     fltLFP = 1;fltfreqs = [1 20];fs = 1e3;
    %     fltLFP = 1;fltfreqs = [6 12];fs = 1e3;
    fltLFP = 1;fltfreqs = [3 12];fs = 1e3;
    
elseif doGamma
    fltLFP = 1;fltfreqs = [50 150];fs = 1e3;
else
    fltLFP = 0;
end

if ~exist('keepsplines','var'), keepsplines = 1;end

dttrl = 0;


LFPdatatemp = [];

shift = 100;

locstmp        = ((0:200:3000)+shift)/1000;%mm

% for missing channels, a dummy variable of all zeros needs to be created
% so that a replacement average can be taken
% this also assumes that these 'bad' channnels are in the spreadsheet
possibleLFPs = {'AD01','AD02','AD03','AD04','AD05','AD06','AD07','AD08','AD09','AD10','AD11','AD12','AD13','AD14','AD15','AD16'};
for k = 1:length(possibleLFPs)
    if ~any(strcmp(possibleLFPs{k},data.label))
        data.label = [data.label possibleLFPs{k}];
        [data.label,sortind] = sort(data.label);
        for kk = 1:length(data.trial)
            data.trial{kk} = [data.trial{kk};  zeros(1,size(data.trial{kk},2))];
            data.trial{kk} = data.trial{kk}(sortind,:);
        end
    end
end

if strmatch(dttrl,'yes') %remove temporal average from each channel
    for k = 1:length(data.trial)
        data.trial{k} = data.trial{k}-repmat(nanmean(data.trial{k},2),1,size(data.trial{k},2));
    end
end
if fltLFP
    for k = 1:length(data.trial)
        for kk = 1:size(data.trial{1},1)
            data.trial{k}(kk,:) = bpf( data.trial{k}(kk,:),fltfreqs,4,fs);
        end
    end
end


% Assumes that voltage is in Volts to start with, around 0.4 uA/mm3 is the max CSD you'd expect to see after hundreds of
% trials on the average, plexon is in microvolts, need to do *1e-6
%                 cond = 0.3;%S/m
%                 cond_top = 0.3;
%                 diam = 0.5*1e-3;%mm->m
%                 gauss_sigma = 0.15*1e-3;%mm->m
%                 filter_range = 5*gauss_sigma;
% the CSD can be thought of as:
% CSD (of positive charge) for an external/arbitrary location   =   V''/R
% CSD (of negative charge wrt the outside of the membrane) or of positive charge
%     wrt the internal wrt membrane         = - V''/R
%       this is the 'CSD', so the sign of the LFP is the same as the sign of the CSD
%       meaning you should negate the CSD to put the sinks as being more positive
% basically, a locally negative LFP should mean local excitation because the internal
% part of the cell is becoming positive
%  V'' has -2Vo, so a neg. LFP leads to a pos. CSD (ext. source/internal sink)
% -V'' means neg. LFP -> neg. CSD
% standard method as coded here is defined correctly, so that it needs to be negated for sinks (excitation) to be positive
%
% so here more positive values mean excitation
% SET PARAMETERS FOR THE CONDUCTANCE AND THE SMOOTHING:
%---------------------------------------------------------
cond = 0.3;%S/m S = 1 A/V, this is 0.3 mhos/meter = 3.33... ohms/meter
cond_top = 0.3;
diam = 0.50*1e-3;%mm->m
% IMPORTANT NOTE:
% standard CSD assumes infinite disc sources so that you can reduce to 1-dimension
% as disc diameter of the presumed sources approaches infinity, this becomes the standard method
% as the diameter becomes very small (tiny sources), the CSD is then proportional to the LFP!

gauss_sigma = 0.1*1e-3;%mm->m
filter_range = 5*gauss_sigma;
%                 num_zs = 200;
num_zs = 50;

% there's a little bit of sub-optimal stuff here leftover from a previous version
labels = data.label;
clear CL
for kk = 1:length(labels)
    fc1 = labels{kk}(1);
    if fc1 == 'A';CL(kk,1) = str2num(labels{kk}(3:4))+100;%100 for AD
    elseif fc1 == 's'
        CL(kk,1) = str2num(labels{kk}(4:6))+200;         %200 for spks
    else
        CL(kk,1) = nan;         % other type
    end
end
rows = [];chi = 1;
if avgForBad,goodkeep = good;good = goodset;else goodkeep = good;end
for Cvs = good
    ACv = Cvs+100;
    row = find(CL(:,1) == ACv);
    if ~isempty(row)
        rows(chi) = row;
    else
        rows(chi) = nan;
    end
    chi = chi+1;
end

el_pos = locstmp(good)*1e-3;%mm->m
%         if ~keeptip, el_pos = el_pos-0.0009;end
%         if ~keeptip, el_pos = el_pos-0.0009;end

[el_pos si] = sort(el_pos);
chselind = rows(si);
oldchannelorder = good;
newchannelorder = good(si);
newgoodorder = newchannelorder;newgoodorder(ismember(newgoodorder,bad)) = nan;
Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);
ft_progress('init', 'none', 'computing CSD')
Nrpt = length(data.trial);
for trli = 1:Nrpt %should do for each trial if want to use the trials later, e.g. for SFC
    %             progress(trli/length(data.trial));
    dtmp = data.trial{trli}(rows(si),:);
    
    if avgForBad,
        % average potentials for bad channels to create surrogate data there
        % to improve the estimate
        for bi = 1:length(bad)
            bind = find(newchannelorder==bad(bi));
            avgind1 = nearest(newgoodorder,bad(bi)-1);
            avgind2 = nearest(newgoodorder,bad(bi)+1);
            if avgind1~=avgind2
                dtmp(bind,:) = nanmean(dtmp([avgind1 avgind2],:),1);
            else
                dtmp(bind,:) = dtmp(avgind1,:);
            end
            data.trial{trli}(rows( oldchannelorder == bad(bi)), : ) = dtmp(bind,:);
        end
        
    end
    
    % do the CSD calculation
    [zs,CSD_cs] = make_cubic_splines(el_pos,dtmp,Fcs,num_zs);
    if gauss_sigma~=0 %filter iCSD
        [zs,CSD_cs] = gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
    end;
    unit_scale  = 1e-3 *1e-6; % 1e-6 is to correct for the data being in microvolts, A/m^3 -> muA/mm^3, schroeder uses mV/mm^2, CSD = mV/mm^2 * g, g = 0.3 1/ohm-m
    CSD_cs      = CSD_cs*unit_scale;
    
    if keepsplines
        timetmp = data.time{trli};
        LFPdatatemp.trial{trli} = [];
        LFPdatatemp.trial{trli} = CSD_cs;
        LFPdatatemp = var2field(LFPdatatemp,el_pos,zs,gauss_sigma,filter_range,cond,cond_top,diam,0);
        if trli == 1
            labels2 = cell(length(zs),1);
            for k = 1:length(zs)
                labels2{k,1} = num2str(rsig(zs(k)*1e6,0));%microns
            end
            LFPdatatemp.label = labels2;
        end
    else
        % choose nearest points
        %                     if trli == 1
        %                         newrows = zeros(size(el_pos))';
        %                         for pi = 1:length(el_pos)
        %                             newrows(pi) = nearest(zs,el_pos(pi));
        %                         end
        %                     end
        %do interpolation -better than nearest
        timetmp =   data.time{trli};% this allows for variable trial lengths!!
        [X Y]   =   meshgrid(timetmp,zs);%current interpolated data points
        [XI YI] =   meshgrid(timetmp,el_pos);%points to extract
        CSD_new =   interp2(X,Y,CSD_cs,XI,YI);
        LFPdatatemp.trial{trli} = [];
        %     LFPdatatemp.trial{trli} = CSD_cs(newrows,:);%for nearest
        LFPdatatemp.trial{trli} = CSD_new;
        if trli == 1,LFPdatatemp.label = labels(rows(si));end
    end
    
    
    
end
ft_progress('close');

h=[];
for k=1:length(LFPdatatemp.trial)
    h(k,:,:)=LFPdatatemp.trial{k};
end
figure;imagesc(-0.2:0.001:0.25,LFPdatatemp.zs,squeeze(nanmean(h,1))); axis xy
set(gca,'YTick',el_pos)
hold on;line([0 0],ylim,'Color','k')
