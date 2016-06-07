fid='JN130604.1'

data_eyeset=strcat('R:\Buffalo Lab\Virtual Navigation\Recording Data\NEX Files\',fid,'.nex');

% get the header info
hdr = ft_read_header(data_eyeset);

% set the trial structure to contain one trial, with all samples contained
% within that trial
cfg=[];
cfg.datafile       = data_eyeset;
cfg.trl           = [10 hdr.nSamples 10]; % start at 10th sample because the first few samples are NaNs

%code to find spike/LFP/eye channels
[lfpind,spkind,eyeind,namarr]=getchannels(data_eyeset);

cfg.channel       = cellstr(namarr(eyeind,:));
cfg.dftfilter     = 'no';
data_eye = ft_preprocessing(cfg);

clear cfg
cfg.artfctdef                     = [];
cfg.artfctdef.xysaccade           = [];
cfg.artfctdef.xysaccade.sgn       = {'X' 'Y'};
cfg.artfctdef.xysaccade.fltord    = 40;
cfg.artfctdef.xysaccade.lpfreq    = 40;
cfg.artfctdef.xysaccade.threshold = 10000;
cfg.artfctdef.xysaccade.artpadding= 0.01;

data_eye.artifact = cell(size(data_eye.time));
for trllop = 1:length(data_eye.trial)
    
    e_x = data_eye.trial{trllop}(strmatch('X',data_eye.label),:);
    e_y = data_eye.trial{trllop}(strmatch('Y',data_eye.label),:);

    %low pass filter the eye position data_eye 
    fltord = cfg.artfctdef.xysaccade.fltord;
    lowpasfrq = cfg.artfctdef.xysaccade.lpfreq;
    nyqfrq = data_eye.hdr.Fs ./ 2;
    flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
    e_x_lowpassfilter=filtfilt(flt,1, e_x);
    e_y_lowpassfilter=filtfilt(flt,1, e_y);

    %differentiate and multiply with sampling rate to get velocity as deg/sec
    x_vF = diff(e_x_lowpassfilter) .* data_eye.hdr.Fs; 
    y_vF = diff(e_y_lowpassfilter) .* data_eye.hdr.Fs;

    % % differentiate and multiply with sampling rate without filtering
    % % this gives the eye velocity in the horizontal and vertical domains
    x_v = diff(e_x) .* data_eye.hdr.Fs; 
    y_v = diff(e_y) .* data_eye.hdr.Fs;

    % combine x- and y-velocity to get eye velocity in degrees/second
    vel = abs(complex(x_v,y_v));
    velF = abs(complex(x_vF,y_vF));

    % set an eye velocity threshold: saccades are eye velocities greater than
    % this threshold
    lim = cfg.artfctdef.xysaccade.threshold; 

    %detect saccade begins and saccade ends
    sacbeg = find(diff(velF > lim) > 0);
    sacend = find(diff(velF > lim) < 0);
    % next three lines unique for JN130604.1
    if velF(1)>lim
        sacbeg = [11 sacbeg];
    end
    if velF(end)>lim
        sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
    end
    
    artifact = round([sacbeg(:) - cfg.artfctdef.xysaccade.artpadding .* ...
        data_eye.hdr.Fs sacend(:) + cfg.artfctdef.xysaccade.artpadding .* data_eye.hdr.Fs]);
    
    data_eye.artifact{trllop} = data_eye.time{trllop}(artifact);
  
end

% calculate time-windowed saccade rate
data_eye.sacpersec = cell(size(data_eye.time));
for trllop = 1:length(data_eye.time)
    
    % average saccade start & end time to get saccade time
    % (timestamp @ center of saccade)
    arttim = mean(data_eye.artifact{trllop},2);
    
    data_eye.sacpersec{trllop}=nan(size(data_eye.time{trllop}));
    ft_progress('init', 'etf',     'Please wait...');
    for timlop = 501:length(data_eye.time{trllop})-500
        ft_progress(timlop/(length(data_eye.time{trllop})-500), 'Processing event %d from %d', timlop, (length(data_eye.time{trllop})-500));
        data_eye.sacpersec{trllop}(timlop) = length(find(arttim>=data_eye.time{trllop}(timlop)-0.5 & arttim<=data_eye.time{trllop}(timlop)+0.5));
    end
    ft_progress('close')
    
end

save('C:\Users\michael.jutras\Documents\MATLAB\Virtual\FT data\JN1306041_sacdat140503.mat','data_eye')

%%
load('C:\Users\michael.jutras\Documents\MATLAB\Virtual\Data (MAT files)\JN1306041_sacdat140503.mat')
pep=load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\Emory\JN130604.1_pEpisode_Ch4_140503.mat');

% freqs = (2^(1/8)).^(8:42);
freqs = (2^(1/8)).^(8:0.2:30);

chnlop = 4;

sacratperfrq = nan(size(freqs));

sacpersec_smooth = density(data_eye.sacpersec{1},1000,'gauss');
sac_hist=[];
edges = 0:0.01:10;
for frqlop = 1:length(freqs)
    
    sacratperfrq(frqlop) = nanmedian(sacpersec_smooth(logical(pep.unionvector(frqlop,:))));

    n = histc(sacpersec_smooth(logical(pep.unionvector(frqlop,:))),edges);
    sac_hist(frqlop,:) = n;

%     figure
%     hist(sacpersec_smooth(logical(pep.unionvector_all{chnlop}(frqlop,:))),500)
%     xlim([0 15]);
%     ylim([0 40000]);
%     title(num2str(freqs(frqlop)))
    
end

% interpolate to plot with correct axis
frq_interp = freqs(1):0.01:freqs(end);
sac_hist_interp = nan(length(frq_interp),size(sac_hist,2));
for edglop = 1:size(sac_hist,2)
    for frqlop = 1:size(sac_hist,1)
        sac_hist_interp(ft_nearest(frq_interp,freqs(frqlop)),edglop) = sac_hist(frqlop,edglop);
    end
    sac_hist_interp(:,edglop)=inpaint_nans(sac_hist_interp(:,edglop),2);
end


%%

r = length(sacpersec_smooth);
edges = 0:0.1:8;
[n, bin] = histc(sacpersec_smooth,edges);
figure;bar(edges,n/r,'histc')
xlabel('saccade rate (Hz)')
ylabel('proportion')
title('proportion of time @ each saccade rate')

pepbin=nan(length(edges)-1,length(freqs));
for binlop = 1:length(edges)-1
%     pepbin(binlop,:)=mean(B(:,bin==binlop),2);
    pepbin(binlop,:)=mean(pep.unionvector(:,bin==binlop),2);
end

centers = 0.05:0.1:7.95;
  
% interpolate to plot with correct axis
frq_interp = freqs(1):0.01:freqs(end);
pepbin_interp = nan(size(pepbin,1),length(frq_interp));
for cntlop = 1:size(pepbin,1)
    for frqlop = 1:size(pepbin,2)
        pepbin_interp(cntlop,ft_nearest(frq_interp,freqs(frqlop))) = pepbin(cntlop,frqlop);
    end
    pepbin_interp(cntlop,:)=inpaint_nans(pepbin_interp(cntlop,:),2);
end

peakfrq = nan(1,size(pepbin_interp,1));
for binlop = 1:size(pepbin_interp,1)
    peakfrq(binlop) = mean(frq_interp(pepbin_interp(binlop,:)==max(pepbin_interp(binlop,:))));
end

figure;hold on
imagesc(frq_interp(ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)),centers(ft_nearest(centers,2):ft_nearest(centers,8)),pepbin_interp(ft_nearest(centers,2):ft_nearest(centers,8),ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)))
axis xy
colormap hot
plot(peakfrq,centers,'Color','b','LineWidth',2)
xlim([2 8])
ylim([2 8])
xlabel('LFP (Pepisode) freq (Hz)')
ylabel('saccade rate (Hz)')
title('JN1306041 ch4 Pepisode peak frequency for each saccade rate bin')

%%

% figure;imagesc(edges,freqs,sac_hist);axis xy
% clim = get(gca,'clim');
% close
% 
% %[1 0 0] red
% %[0 1 0] green
% %[0 0 1] blue
% c1 = [zeros(size(0:2/clim(2):1)) 0:2/clim(2):1];
% c2 = [0:2/clim(2):1 1:-2/clim(2):0];
% c3 = [1:-1/clim(2):0 zeros(size(1:-1/clim(2):0))];
% 
% figure;hold on
% % for frqlop=1:size(sac_hist,1)
% for frqlop=25:75
%     for edglop = 1:size(sac_hist,2)
%         sacval = round(sac_hist(frqlop,edglop));
%         if sacval~=0
%             h=scatter(edges(edglop),freqs(frqlop),'.');
%             set(h,'MarkerEdgeColor',[c1(sacval) c2(sacval) c3(sacval)])
%         end
%     end
% end
    

