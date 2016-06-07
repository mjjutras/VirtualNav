dataset=strcat('R:\Buffalo Lab\Virtual Navigation\Recording Data\NEX Files\JN130531.2.nex');

chnsel = 1;

lofrq1 = 3;
lofrq2 = 12;
hifrq1 = 20;
hifrq2 = 40;

%%
% get the data
cfg=[];
cfg.trialfun      = 'trialfun_nav';
cfg.dataset       = dataset;
cfg = ft_definetrial(cfg);

[lfpind,spkind,eyeind,namarr]=getchannels(dataset);

labarr=[namarr(lfpind,:)];

cfg.channel       = cellstr(labarr)';
cfg.dftfilter     = 'yes';
cfg.dftfreq       = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2];
cfg.padding       = 10;
cfg.continuous    = 'yes';
cfg.detrend       = 'no';

data = ft_preprocessing(cfg);

srate=data.fsample;

%%

datbuf=[];
for k=1:length(data.trial)
    datbuf=[datbuf data.trial{k}(chnsel,:)];
end

%%
Vlo=eegfiltMJ(datbuf,srate,lofrq1,lofrq2);
Vhi=eegfiltMJ(datbuf,srate,hifrq1,hifrq2);

[r, r_CI, nCtlPts] = GLM_CFC_for_paper(Vlo, Vhi, 10);
title([dataset(find(dataset=='\',1,'last')+1:end) '; chn ' num2str(chnsel) '; r = ' num2str(r) '; ' num2str(lofrq1) '-' num2str(lofrq2) '; ' num2str(hifrq1) '-' num2str(hifrq2)])
box off
set(gca,'TickDir','out')


