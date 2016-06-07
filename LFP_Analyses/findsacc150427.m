BRnam = 'JN150422001';
sesDir = 'JN_BR_15_04_22\JN_BR_15_04_22_14_54';

NSDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NSdat';
BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';



[~,sesnam]=fileparts(sesDir);

disp(['Processing ' sesnam])
load(fullfile(NSDir,[sesnam '_NSdat.mat']))


NS2 = openNSx(fullfile(BRDir,[BRnam '.NS2']),'read','precision','double');

% create label arrays for NS2 and NS6 channels
NS2lab = cell(1,length(NS2.ElectrodesInfo));
for lablop = 1:length(NS2.ElectrodesInfo)
    NS2lab{lablop} = NS2.ElectrodesInfo(lablop).Label;
end

% line up behavioral data with raw LFPs (adapted from quickLineup)

% find eye & position signals in NS2
eyeindx = find(strncmpi(NS2lab,'ainp1',5));
eyeindy = find(strncmpi(NS2lab,'ainp2',5));

%%

fltord = 40;
lowpasfrq = 80;
lim = 100000;
artpadding = 0.01;
Fs = NS2.MetaTags.SamplingFreq; % sampling rate in Hz

e_x = double(NS2.Data(eyeindx,:));
e_y = double(NS2.Data(eyeindy,:));

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

figure;hold on
plot(vel(1:100000))
% plot(velF(1:100000),'r')
line(xlim,[lim lim],'Color','k')
for k=1:20;line([artifact(k,1) artifact(k,1)],ylim,'Color','g');line([artifact(k,2) artifact(k,2)],ylim,'Color','r');end

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
% artifact = round([sacbeg(:) sacend(:)]);

%% get NS6 saccade-aligned data

clear sacdat
c=1;
ft_progress('init', 'etf',     'Please wait...');
for trllop = 1:length(data.time)
    ft_progress(trllop/length(data.time), 'Processing event %d from %d', trllop, length(data.time));
   
    artind = find(artifact(:,1)>=data.sampleinfo(trllop,1) & artifact(:,2)<=data.sampleinfo(trllop,2));
    
%     figure;hold on
%     plot(data.time{trllop},data.trial{trllop}(ismember(data.label,{'eyeX_B' 'eyeY_B'}),:))
%     for saclop=1:length(artind)
%         line([data.time{trllop}(artifact(artind(saclop),1)-data.sampleinfo(trllop,1)) ...
%             data.time{trllop}(artifact(artind(saclop),1)-data.sampleinfo(trllop,1))],ylim,'Color','g')
%         line([data.time{trllop}(artifact(artind(saclop),2)-data.sampleinfo(trllop,1))...
%             data.time{trllop}(artifact(artind(saclop),2)-data.sampleinfo(trllop,1))],ylim,'Color','r')
%     end

    for saclop = 1:length(artind)
        trlsacind = artifact(artind(saclop),1)-data.sampleinfo(trllop,1);
        if trlsacind>=200 && trlsacind+200<=length(data.time{trllop})
            sacdat(c,:,:) = data.trial{trllop}(1:36,trlsacind-200:trlsacind+200);
            c=c+1;
        end
    end
    
end
ft_progress('close')



