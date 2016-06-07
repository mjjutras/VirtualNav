
fid='JN130531.1-sorted'

dataset=strcat('R:\Buffalo lab\Virtual Navigation\Free Foraging\NEX Files\',fid,'.nex');

tag='_pEpisode140428';
resfil = strcat('C:\Users\michael.jutras\Documents\MATLAB\Virtual\Pepisode & CFC\',fid,tag,'.mat');

ft_defaults % 

% get the header info
hdr = ft_read_header(dataset);

% set the trial structure to contain one trial, with all samples contained
% within that trial
cfg=[];
cfg.dataset       = dataset;
cfg.trl           = [10 hdr.nSamples 10]; % start at 10th sample because the first few samples are NaNs

%code to find spike/LFP channels
header=getnexheader(dataset);
numvar = size(header.varheader,2);
for varlop = 1:numvar
    typarr(varlop) = header.varheader(varlop).typ;
end
for varlop = 1:numvar
    namarr(varlop,1:64) = char(header.varheader(varlop).nam');
end
spkind = find(typarr == 0);
analogind = find(typarr == 5);
lfpindbzganalog = find(double(namarr(analogind,1)) == 65);
lfpind = analogind(lfpindbzganalog);
sortindx=find(double(namarr(spkind,7)) ~= 105);
spkind=spkind(sortindx);
% labarr=[namarr(spkind,:);namarr(lfpind,:)];
labarr=[namarr(lfpind,:)];

cfg.channel       = cellstr(labarr(1:4,:))';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
cfg.padding       = 10;
cfg.continuous    = 'yes';
cfg.detrend       = 'no';

data = ft_preprocessing(cfg); 

clear Pm_all pv_all eeg_all unionvector_all

% convert each data channel into one continuous signal
for chnlop=1:length(data.label)
    
    eeg=nan(1,sum(cfg.trl(:,2)-cfg.trl(:,1))+size(cfg.trl,1));
    trlind=nan(1,sum(cfg.trl(:,2)-cfg.trl(:,1))+size(cfg.trl,1));
    for l=1:length(data.trial)
        ind=find(isnan(eeg),1,'first'):(find(isnan(eeg),1,'first')+length(data.trial{l})-1);
        eeg(ind)=data.trial{l}(chnlop,:);
        trlind(ind)=data.time{l};
    end
    
    samplerate=data.fsample;
    
    freqs = (2^(1/8)).^(8:42);
    
    width=7;

    shoulderMS = 500;
    shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

    t_total = clock;

    % get energy vector for each list
    fprintf('Calc. Power...');
    t0 = clock;

    B=single(multienergyvec(eeg,freqs,samplerate,width));
    fprintf('%g\n',etime(clock,t0));

    % calc the mean fit
    Blog = log10(double(B));
%     Pm = log10(mean(B,2));
    Pm = mean(Blog,2);

    % get the fit
    fprintf('Calc. fit...');
    t0 = clock;
    % IMPORTANS: chi_squarefit assumes that frequencies are
    % logarithmically spaced!
    [all,R2] = chi_squarefit(freqs,Pm);
    all = all';
    fprintf('%g\n',etime(clock,t0));

    % test fit of regression line
    pv=polyfit(log10(freqs),Pm',1); % linear regression
%     figure;plot(freqs,Pm)
%     hold on;plot(freqs,pv(2)+pv(1)*log10(freqs),'r')
    
    % set the threshold
    thresh = all(:,951);

    % loop through each frequency and save the unionvector
    fprintf('Calc. union...');
    t0 = clock;
    unionvector = single(zeros(size(freqs,2),size(B,2)));

    for f = 1:size(freqs,2)
        % get the pepisode
        unionvector(f,:) = single(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0));
    end

    fprintf('%g\n',etime(clock,t0));

    % convert unionvector from single to int8 to conserve memory
    unionvector=int8(unionvector);

    Pm_all{chnlop}=Pm;
    pv_all{chnlop}=pv;
    unionvector_all{chnlop}=unionvector;

end

% save the results
save(resfil, 'Pm_all', 'pv_all', 'unionvector_all', 'trlind');
disp(strcat('Generated:',resfil))

%%

fid='JN130531.2'
tag='_pEpisode';
resfil = strcat('C:\Users\michael.jutras\Documents\MATLAB\Virtual\',fid,tag,'.mat');

load(resfil)

figure;plot(freqs,mean(unionvector_all{1},2))
hold on;plot(freqs,mean(unionvector_all{2},2),'r')
plot(freqs,mean(unionvector_all{3},2),'g')
plot(freqs,mean(unionvector_all{4},2),'m')
xlabel('freq'); ylabel('Pepisode')
legend('Ch. 1','Ch. 2','Ch. 3','Ch. 4')
title(fid)