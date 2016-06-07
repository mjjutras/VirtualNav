%% save variables to network 
if strcmp(getenv('computername'),'MIKE-PC')
    WSdir = 'C:\Data\VR';
elseif strcmp(getenv('computername'),'RBU-MIKEJ')
    WSdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\workspace';
end
% save(fullfile(WSdir,'anaFruit_neuraldum_trainingses_newperffact_150706.mat'),'-v7.3')

% load variables
load(fullfile(WSdir,'anaFruit_neuraldum_trainingses_newperffact_150706.mat'))


%% preprocess data

% exclude these channels
badchn = {'A02'; 'A05'; 'A08'; 'A10'; 'A12'; 'B02'; 'B12'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};

% exclude eye and position data for now
badchn = [badchn; {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'}];

cfg=[];
cfg.continuous    = 'no';
cfg.channel       = setxor(dum.label,badchn);
data0 = ft_preprocessing(cfg,dum);
clear dum

% these data are still noisy; wait until 5-second clips are selected before
% filtering out noise

% noisy trials:
% 11, 12, 71, 72, 197, 338, 361, 413, 418, 420, 431, 432, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 549, 553, 600


%% find trials with at least 5 seconds preceding banana
% this should be every trial that was included in nrlsel

bantim = nan(length(nrlsel)*2,1);
timsel = nan(length(nrlsel)*2,2);
nrltimsel = [];
trlgrp = [];
for invlop = 1:length(nrlsel)
    bantim(invlop*2-1:invlop*2) = bantimindall(nrlsel(invlop)-2:nrlsel(invlop)-1);
    if isempty(find((bantim(invlop*2-1:invlop*2))<5000,1))
        timsel(invlop*2-1:invlop*2,:) = [bantim(invlop*2-1:invlop*2)-4999 bantim(invlop*2-1:invlop*2)];
        nrltimsel = [nrltimsel; nrlsel(invlop)];
        trlgrp = [trlgrp; ones(2,1)+length(trlgrp)/2];
    else
        timsel(invlop*2-1:invlop*2,:) = repmat([nan nan],2,1);
    end
end


%% select last 5 seconds of approach path to banana in neural data

data1 = data0;
for trllop = 1:length(data0.trial)
    if isempty(find(isnan(timsel(trllop,:)),1,'first'))

        % use this code for modulation index (1-sec padding at the end)
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1+1000];
        data1.time{trllop} = -4.999:0.001:1.0;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2)+1000);
        
    end
end

% remove marked trials
cfg=[];
% here choose only trials with at least 5 seconds pre-banana 
cfg.trials        = setxor(1:length(data1.trial),find(isnan(timsel(:,1))));
data1 = ft_preprocessing(cfg,data1);
clear data0

% subtract common average of each probe from all channels on that probe
data2 = data1;
% comment-out the following if choosing not to re-reference to common probe average
clear prblaball
for k=1:length(data1.label)
    prblaball(k,1)=data1.label{k}(1);
end
prblab = unique(prblaball);
[~,locb] = ismember(prblaball,prblab);
for trllop = 1:length(data1.trial)
    for prblop = 1:length(prblab)
        prbavg = mean(data1.trial{trllop}(locb==prblop,:),1);
        repavg = repmat(prbavg,length(find(locb==prblop)),1);
        mat1 = reshape(data1.trial{trllop}(locb==prblop,:),1,length(find(locb==prblop)),size(data1.trial{trllop},2));
        mat2 = reshape(repavg,1,size(repavg,1),size(repavg,2));
        data2.trial{trllop}(locb==prblop,:)=diff(cat(1,mat1,mat2),1,1);
    end
end
clear data1

% now run ft_rejectvisual and ft_databrowser to clean bad trials


%% mark bad trials (CHANGE IF ADDING NEW DATA)
% identified with ft_rejectvisual, ft_databrowser

badtrldum = [16, 25, 28, 30, 37, 106, 113, 148]; % if re-referencing to common probe average

badtrl = [];
invselexc = [];
for badlop = 1:length(badtrldum)
    badtrl = [badtrl; find(trlgrp==trlgrp(badtrldum(badlop)))];
    invselexc = [invselexc; find(trlgrp==trlgrp(badtrldum(badlop)),1,'last')/2];
end

nrlbadsel = nrltimsel(setxor(1:length(nrltimsel),invselexc));

% exclude bad trials, filter out line noise here
cfg=[];
cfg.trials = setxor(1:length(data2.time),badtrl);
cfg.continuous    = 'no';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [60 120];
data3 = ft_preprocessing(cfg,data2);
clear data2


%%
% construct the time series by selecting 3-second epochs from pre-stimulus
% fixation period, which is the period required for looking down to 3 Hz
% (we're missing the first second due to the need to buffer, probably ok
% since theta doesn't seem to show up until a little later)

datarr = [];
c=1;
for trllop = 1:length(data3.trial)
    datarr(c,:,:) = data3.trial{trllop}(:,1:3000); % 1001:2000
    datarr(c+1,:,:) = data3.trial{trllop}(:,1001:4000); % 2001:3000
    datarr(c+2,:,:) = data3.trial{trllop}(:,2001:5000); % 3001:4000
    datarr(c+3,:,:) = data3.trial{trllop}(:,3001:6000); % 4001:5000
    c=c+4;
end
    

%% get the modulation index

flow1 = (3:1:20);       %The phase frequencies.
flow2 = flow1+2.0;      %3-22 Hz with 1 Hz steps, 2 Hz bandwidth

fhigh1 = (4:2:100);     %The amp frequencies.
fhigh2 = fhigh1+4.0;    %4-104 Hz with 2 Hz steps, 4 Hz bandwidth


srate=1000;

mod4d = nan(length(data3.label),size(datarr,1),length(flow1), length(fhigh1));
for chnlop = 1:length(data3.label)
    mod3d=[];
    for trllop=1:size(datarr,1)
        fprintf(['Trial # ' num2str(trllop) '\n'])
        mod2d = zeros(length(flow1), length(fhigh1));
        for i=1:length(flow1)
            theta=eegfiltMJ(squeeze(datarr(trllop,chnlop,:))',srate,flow1(i),flow2(i));  %Compute the low freq signal.
            theta=theta(1001:2000);           %Drop the first and last second.
            phase = angle(hilbert(theta));                  %Compute the low freq phase.
            fprintf(['Loops remaining = ', num2str(length(flow1)-i) '\n'])
            for j=1:length(fhigh1)
                gamma=eegfiltMJ(squeeze(datarr(trllop,chnlop,:))',srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
                gamma=gamma(1001:2000);         %Drop the first and last second.
                amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.

                %%%%%Compute the modulation index%%%%%
                newphase=zeros(1,length(phase))+nan;
                for k=1:length(phase)
                    if phase(k)>=0
                        newphase(k)=phase(k);
                    else
                        newphase(k)=2*pi+phase(k);
                    end
                end
                newphase=circ_rad2ang(newphase);

                binedges=0:20:360;
                [nphs,binphs] = histc(newphase,binedges);
                nphs=nphs(1:end-1);
                A=zeros(1,length(nphs))+nan; %mean_amp_at_each_phase
                for k=1:length(nphs)
                    A(k)=mean(amp(binphs==k));
                end

                %entropy
                p=zeros(1,length(nphs))+nan;
                for k=1:length(nphs)
                    p(k)=A(k)/nansum(A);
                end

                H=-sum(p.*log10(p));
                Hmax=log10(length(p));
                MI=(Hmax-H)/Hmax;
                %%%%%END%%%%%

                mod2d(i,j) = MI;
            end
        end

        mod3d=cat(1,mod3d,reshape(mod2d,1,size(mod2d,1),size(mod2d,2)));
    end
    mod4d(chnlop,:,:,:) = mod3d;
    clear mod3d
end

%Plot the two-dimensional modulation index.

modavg=squeeze(mean(mod3d,1))';

% mod2d=mod2d';

figure
flow = (flow1 + flow2) / 2;
fhigh = (fhigh1 + fhigh2) / 2;
imagesc(flow, fhigh, modavg);  colorbar;
axis xy
set(gca, 'FontSize', 12);
xlabel('Phase Frequency [Hz]');  ylabel('Envelope Frequency [Hz]');

newmod=[];
for k=1:size(mod3d,1)
    newmod(k,:,:)=squeeze(mod3d(k,:,:))-surthresh;
end

newmodavg=squeeze(mean(newmod,1))';

figure
flow = (flow1 + flow2) / 2;
fhigh = (fhigh1 + fhigh2) / 2;
imagesc(flow, fhigh, newmodavg);  colorbar;
axis xy
set(gca, 'FontSize', 12);
xlabel('Phase Frequency [Hz]');  ylabel('Envelope Frequency [Hz]');



%% construct trial shuffled composite time series then get surrogate MI values

% rndind1=[];
% rndind2=[];
% for k=1:200
%     rand1=randperm(size(datbuf,1));
%     rand2=randperm(size(datbuf,1));
%     rndind1(k)=rand1(1);
%     rndind2(k)=rand2(1);
%     i=2;
%     while rndind1(k)==rndind2(k)
%         rndind2(k)=rand2(i);
%         i=i+1;
%     end
% end

rndind1=randperm(size(datbuf,1));
rndind2=randperm(size(datbuf,1));

flow1 = (1:1:20);       %The phase frequencies.
flow2 = flow1+2.0;      %3-22 Hz with 1 Hz steps, 2 Hz bandwidth

fhigh1 = (4:2:100);     %The amp frequencies.
fhigh2 = fhigh1+4.0;    %4-104 Hz with 2 Hz steps, 4 Hz bandwidth

srate=1000;

sur3d=[];
for trllop=1:length(rndind1)
    fprintf(['Trial # ' num2str(trllop) '\n'])
    sur2d = zeros(length(flow1), length(fhigh1));
    for i=1:length(flow1)
        theta=eegfiltMJ(datbuf(rndind1(trllop),:),srate,flow1(i),flow2(i));  %Compute the low freq signal.
        theta=theta(5000:5999);           %Drop the first and last 5000 ms.
        phase = angle(hilbert(theta));                  %Compute the low freq phase.
        fprintf(['Loops remaining = ', num2str(length(flow1)-i) '\n'])
        for j=1:length(fhigh1)
            gamma=eegfiltMJ(datbuf(rndind2(trllop),:),srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
            %         gamma=gamma(srate:length(s)-srate-1);         %Drop the first and last second.
            gamma=gamma(5000:5999);         %Drop the first and last 5000 ms.
            amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.

            %%%%%Compute the modulation index%%%%%
            newphase=zeros(1,length(phase))+nan;
            for k=1:length(phase)
                if phase(k)>=0
                    newphase(k)=phase(k);
                else
                    newphase(k)=2*pi+phase(k);
                end
            end
            newphase=circ_rad2ang(newphase);

            binedges=0:20:360;
            [nphs,binphs] = histc(newphase,binedges);
            nphs=nphs(1:end-1);
            A=zeros(1,length(nphs))+nan; %mean_amp_at_each_phase
            for k=1:length(nphs)
                A(k)=mean(amp(binphs==k));
            end

            %entropy
            p=zeros(1,length(nphs))+nan;
            for k=1:length(nphs)
                p(k)=A(k)/nansum(A);
            end

            H=-sum(p.*log10(p));
            Hmax=log10(length(p));
            MI=(Hmax-H)/Hmax;
            %%%%%END%%%%%

            sur2d(i,j) = MI;
        end
    end
    
    sur3d=cat(1,sur3d,reshape(sur2d,1,size(sur2d,1),size(sur2d,2)));
end

%Plot the two-dimensional surrogate modulation index.

suravg=squeeze(mean(sur3d,1))';

% mod2d=mod2d';

figure
flow = (flow1 + flow2) / 2;
fhigh = (fhigh1 + fhigh2) / 2;
imagesc(flow, fhigh, suravg);  colorbar;
axis xy
set(gca, 'FontSize', 12);
xlabel('Phase Frequency [Hz]');  ylabel('Envelope Frequency [Hz]');

% load('SUR_MP0702015c1.mat')

surthresh=zeros(size(sur3d,2),size(sur3d,3))+nan;
for k=1:size(sur3d,2)
    for l=1:size(sur3d,3)
        srtval=sort(sur3d(:,k,l));
        surthresh(k,l)=srtval(196);
    end
end

figure
flow = (flow1 + flow2) / 2;
fhigh = (fhigh1 + fhigh2) / 2;
imagesc(flow, fhigh, surthresh');  colorbar;
axis xy
set(gca, 'FontSize', 12);
xlabel('Phase Frequency [Hz]');  ylabel('Envelope Frequency [Hz]');


%% interpolate
interpnum=50;
clear interp1 interp2
for loop=1:size(mod2d,1)
    interp1(loop,:)=interp(mod2d(loop,:),interpnum);
end
for loop=1:size(interp1,2)
    interp2(:,loop)=interp(interp1(:,loop),interpnum);
end
xinterp=interp(flow,interpnum);
yinterp=interp(fhigh,interpnum);
figure;imagesc(xinterp,yinterp,interp2);colorbar;axis xy
% xlim([-0.2 0.4])
% ylim([0 20])
% title('both')
axis xy
set(gca, 'FontSize', 12);
xlabel('Phase Frequency [Hz]');  ylabel('Envelope Frequency [Hz]');
