
fid='JN130531.2'

dataset=strcat('R:\Buffalo Lab\Virtual Navigation\Recording Data\NEX Files\',fid,'.nex');

ft_defaults

% get the header info
hdr = ft_read_header(dataset);

% set the trial structure to start and end each trial with the appearance
% of the bananas (event code 100)
cfg=[];
cfg.trialfun      = 'trialfun_nav';
cfg.dataset       = dataset;
cfg = ft_definetrial(cfg);


%code to find spike/LFP/eye channels
[lfpind,spkind,eyeind,namarr]=getchannels(dataset);
labarr=[namarr(lfpind,:); namarr(eyeind,:)]; % look at only LFPs and eye data for now

% preprocess the data, removing 60 Hz noise artifacts (and 120 Hz harmonic)
cfg.channel       = cellstr(labarr)';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
cfg.padding       = 10;
cfg.continuous    = 'yes';
cfg.detrend       = 'no';

data = ft_preprocessing(cfg);

srate=data.fsample;

datbuf=[];
for k=1:length(data.trial)
    datbuf=[datbuf data.trial{k}(1,:)];
end

%%

flow1 = (1:1:20);       %The phase frequencies.
flow2 = flow1+2.0;      %3-22 Hz with 1 Hz steps, 2 Hz bandwidth

fhigh1 = (4:2:100);     %The amp frequencies.
fhigh2 = fhigh1+4.0;    %4-104 Hz with 2 Hz steps, 4 Hz bandwidth


srate=1000;

mod3d=[];
for trllop=1:size(timsel,1)
    fprintf(['Trial # ' num2str(trllop) '\n'])
    mod2d = zeros(length(flow1), length(fhigh1));
    for i=1:length(flow1)
        theta=eegfiltMJ(datbuf(timsel(trllop,1):timsel(trllop,2)),srate,flow1(i),flow2(i));  %Compute the low freq signal.
        %     theta=theta(srate:length(s)-srate-1);           %Drop the first and last second.
        theta=theta(5001:6000);           %Drop the first and last 5000 ms.
        phase = angle(hilbert(theta));                  %Compute the low freq phase.
        fprintf(['Loops remaining = ', num2str(length(flow1)-i) '\n'])
        for j=1:length(fhigh1)
            gamma=eegfiltMJ(datbuf(timsel(trllop,1):timsel(trllop,2)),srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
            %         gamma=gamma(srate:length(s)-srate-1);         %Drop the first and last second.
            gamma=gamma(5001:6000);         %Drop the first and last 5000 ms.
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