load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\Emory\JN1306041_sacali_acf_141105.mat')

fixlng = []
for k=1:size(sacali_trunc,1)
    if ~isempty(find(isnan(sacali_trunc(k,1000:end)),1,'first'))
        fixlng(k) = find(isnan(sacali_trunc(k,1000:end)),1,'first');
    end
end

medfixlng = median(fixlng);

lfp = nan(1,size(sacali_trunc,2));
c=1;
for k=1:size(sacali_trunc,1)
    if find(isnan(sacali_trunc(k,1000:end)),1,'first') >= medfixlng
        lfp(c,:) = sacali_trunc(k,:);
        c=c+1;
    end
end

% calculate the phase using the hilbert transform
angraw=nan(size(lfp));
for l=1:size(lfp,1)
    h=hilbert(lfp(l,~isnan(lfp(l,:))));
    ang1=mod(angle(h)+2*pi,2*pi);
    angraw(l,~isnan(lfp(l,:)))=ang1;
end

% calculate rayleigh statistic
phasebins = 0:pi/12:2*pi;

% for all trials
anglp_raw = histc(angraw,phasebins);
anglp_raw = anglp_raw ./ repmat(sum(anglp_raw),size(anglp_raw,1),1);

unitv_raw = exp(1i*angraw);
sumv_raw = squeeze(nanmean(unitv_raw,1));

rayleigh_raw=abs(sumv_raw);
critval_raw=sqrt(-log(0.01)/size(angraw,1));

% group variables
anglp_raw_all(c,:,:)=anglp_raw;
rayleigh_raw_all(c,:)=rayleigh_raw;
critval_raw_all(c)=critval_raw;

%% plot with interpolation/smoothing

phasebins = 0:pi/12:4*pi;
phasebins=phasebins(1:end-1);
timex = -1000:1000;

% pad phases on upper and lower limits to eliminate edge distortion
phasebins=[phasebins phasebins+4*pi phasebins+8*pi];

phaseinterp=interp(phasebins,8);
timeinterp=interp(timex,10);

anglp_raw_cat=cat(1,anglp_raw(1:end-1,:),anglp_raw(1:end-1,:));
anglp_raw_cat=cat(1,anglp_raw_cat,anglp_raw_cat,anglp_raw_cat);

% % limit the time period
anglp_raw_cat=anglp_raw_cat(:,ft_nearest(timex,-400):ft_nearest(timex,800));
timex=timex(ft_nearest(timex,-400):ft_nearest(timex,800));
timeinterp=interp(timex,10);

% interpolate
clear anglp_raw_cat_interp1
for phslop=1:size(anglp_raw_cat,1)
    anglp_raw_cat_interp1(phslop,:)=interp(anglp_raw_cat(phslop,:),10);
end
clear anglp_raw_cat_interp2
for timlop=1:size(anglp_raw_cat_interp1,2)
    anglp_raw_cat_interp2(:,timlop)=interp(anglp_raw_cat_interp1(:,timlop),8);
end

% smooth
clear anglp_raw_cat_smooth1
for k=1:size(anglp_raw_cat_interp2,1)
    anglp_raw_cat_smooth1(k,:)=smooth(anglp_raw_cat_interp2(k,:),15);
end
phasesmooth=smooth(phaseinterp,15);
clear anglp_raw_cat_smooth2
for k=1:size(anglp_raw_cat_smooth1,2)
    anglp_raw_cat_smooth2(:,k)=smooth(anglp_raw_cat_smooth1(:,k),15);
end
timesmooth=smooth(timeinterp,15);

% remove padding
anglp_raw_cat_smooth2=anglp_raw_cat_smooth2(ft_nearest(phasesmooth,4*pi):ft_nearest(phasesmooth,8*pi),:);
phasesmooth=phasesmooth(ft_nearest(phasesmooth,0):ft_nearest(phasesmooth,4*pi));

%%

% change so 0 is the trough instead of the peak
anglp_raw_cat_smooth2_shift=[anglp_raw_cat_smooth2(ft_nearest(phasesmooth,pi):end,:); ...
    anglp_raw_cat_smooth2(1:ft_nearest(phasesmooth,pi)-1,:)];


figure;imagesc(timesmooth,phasesmooth,anglp_raw_cat_smooth2_shift);axis xy
colorbar
set(gca,'TickDir','out')
ylim([0 pi*4])
set(gca,'YTick',[0 pi pi*2 pi*3 pi*4])
line([0 0],ylim,'Color','w')
