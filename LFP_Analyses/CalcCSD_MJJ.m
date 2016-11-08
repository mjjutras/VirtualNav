load('R:\Buffalo Lab\Virtual Navigation\JN140825\JN140825011_sacdat141013.mat')

%%

doTheta = 0;
doGamma = 0;

keepsplines = 1;

bad = [2 5 10 12];
% bad = [];

CSDmethod = 'iCSD_splines';
% CSDmethod = 'standard';

avgForBad = 1; % should definitely do this instead of leaving them out!

dttrl = 'yes'; % detrend the trials, should probably do this too
% dttrl = 'no';


goodset = 1:12;
good = setdiff(goodset,bad);
labels = {'AD01','AD02','AD03','AD04','AD05','AD06','AD07','AD08','AD09','AD10','AD11','AD12'};

data= [];
c=1;
for trllop = 1:size(sacdatB,1)
    if isempty(find(isnan(squeeze(sacdatA(trllop,1,:))), 1))
        data.trial{c} = squeeze(sacdatA(trllop,good,:));
        data.time{c} = -0.2:0.001:0.3;
        c=c+1;
    end
end

lfpind = 1:12;

data.label = labels(good);


%%

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

locstmp        = ((0:200:2200)+shift)/1000;%mm

% for missing channels, a dummy variable of all zeros needs to be created
% so that a replacement average can be taken
% this also assumes that these 'bad' channnels are in the spreadsheet
possibleLFPs = {'AD01','AD02','AD03','AD04','AD05','AD06','AD07','AD08','AD09','AD10','AD11','AD12'};
for k =1:length(possibleLFPs)
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
figure;imagesc(-0.2:0.001:0.3,LFPdatatemp.zs,squeeze(nanmean(h,1))); axis xy
set(gca,'YTick',el_pos)
