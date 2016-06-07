datdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated\';
BRnam = 'JN140825010'

load(fullfile(datdir,[BRnam '_NS6_SF30.mat']))


%%

samplerate = 1000;

freqs = (2^(1/8)).^(8:42);

width=7;

shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

Pm_all = cell(1,size(NS6.Data,1));
pv_all = cell(1,size(NS6.Data,1));
unionvector_all = cell(1,size(NS6.Data,1));

ft_progress('init', 'etf',     'Please wait...');
for chnlop = 1:size(NS6.Data,1)
    
    ft_progress(chnlop/size(NS6.Data,1), 'Processing event %d from %d', chnlop, size(NS6.Data,1));
    
    t_total = clock;

    % get energy vector for each list
    fprintf('Calc. Power...');
    t0 = clock;

    B=multienergyvec(double(NS6.Data(chnlop,:)),freqs,samplerate,width);
    fprintf('%g\n',etime(clock,t0));

    % calc the mean fit
    Pm = log10(mean(B,2));

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
ft_progress('close')

save(['C:\Users\michael.jutras\Documents\Virtual Navigation Study\Analysis\pEpisode\' BRnam '_pEpisode_141031.mat'],'unionvector_all','Pm_all','pv_all','freqs','-v7.3')

%%

pepdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\pEpisode';

freqs = (2^(1/8)).^(8:42);

p1=load(fullfile(pepdir,'JN140825010_pEpisode_141031'));
p2=load(fullfile(pepdir,'JN140825011_pEpisode_141016'));


%%

load('C:\Users\michael.jutras\Documents\Virtual Navigation Study\JN140825\JN140825011_pEpisode_141016.mat')

trlind = (round(NS6_timestamp*1000)-round(NS6_timestamp(1)*1000))/1000;

freqs = (2^(1/8)).^(8:42);

% 2:4 Hz = freqs(1:9)
% 2.59:4 Hz = freqs(4:9)
% 2.59:5.66 Hz = freqs(4:13)

% first designate LFPs recorded from channels with SUA to avoid wasting
% time

boutdur=cell(1);
interbout=cell(1);
c=1;
for chnlop = 1:length(unionvector_all)
    
    % find theta bouts throughout VPLT blocks
    [~,j]=find(unionvector_all{chnlop}(ft_nearest(freqs,3):ft_nearest(freqs,6),:));
    
    j=unique(j)'; % j represents indices of bouts in unionvector
    
    boutind=trlind(j)';
    
    % find breaks between bouts
    boutbreak=find(diff(boutind)>0.002);
    
    % bout definitions: timestamps (according to trlind) designating the
    % occurrence of theta bouts throughout VPLT blocks (timestamps are
    % taken from entire session, as in cfg.trl)
    % (converted to ms)
    boutdef=([boutind([1; boutbreak+1]) boutind([boutbreak; length(boutind)])])*1000;
    
    boutdef=round(boutdef);
    
    boutdur{c}=diff(boutdef,1,2);
    
    interbout{c}=[];    
    for k=1:size(boutdef,1)-1
        interbout{c}=[interbout{c} boutdef(k+1,1)-boutdef(k,2)];
    end
    
    c=c+1;
    
end


%% bout duration & interbout interval

boutdurA=[];
interboutA=[];
boutdurB=[];
interboutB=[];
boutdurC=[];
interboutC=[];
for k=1:length(boutdur)
    if ismember(k,1:12)
        if ismember(k,[2 5 10 12])
            boutdurA=[boutdurA; boutdur{k}];
            interboutA=[interboutA interbout{k}];
        end
    elseif ismember(k,13:24)
        boutdurB=[boutdurB; boutdur{k}];
        interboutB=[interboutB interbout{k}];
    elseif ismember(k,25:36)
        boutdurC=[boutdurC; boutdur{k}];
        interboutC=[interboutC interbout{k}];
    end
end
    
% edges_boutdur = 0:250:10000;
% edges_interbout = 0:0.25:5;
centers_boutdur = (0:250:5000)+125;
centers_interbout = (0:1000:20000)+500;

figure
subplot(2,1,1)
h=hist(boutdurA,centers_boutdur);
h=h/length(boutdurA);
bar(centers_boutdur,h)
xlim([0 5250])
box off
set(gca,'TickDir','out')
% set(gca,'YTick',0:0.05:0.3)
set(gca,'XTick',0:1000:5000)
xlabel('Time [ms]');ylabel('Proportion of bouts')
line([median(boutdurA) median(boutdurA)],ylim,'Color','k','LineStyle','--')
title(['Array A: bout duration, median = ' num2str(median(boutdurA)) ' ms'])
subplot(2,1,2)
h=hist(interboutA,centers_interbout);
h=h/length(interboutA);
bar(centers_interbout,h)
xlim([0 21000])
box off
set(gca,'TickDir','out')
% set(gca,'YTick',0:0.05:0.35)
set(gca,'XTick',0:5000:20000)
xlabel('Time [ms]');ylabel('Proportion of bouts')
line([median(interboutA) median(interboutA)],ylim,'Color','k','LineStyle','--')
title(['Array A: interbout interval, median = ' num2str(median(interboutA)) ' ms'])

figure
subplot(2,1,1)
h=hist(boutdurB,centers_boutdur);
h=h/length(boutdurB);
bar(centers_boutdur,h)
xlim([0 5250])
box off
set(gca,'TickDir','out')
% set(gca,'YTick',0:0.05:0.3)
set(gca,'XTick',0:1000:5000)
xlabel('Time [ms]');ylabel('Proportion of bouts')
line([median(boutdurB) median(boutdurB)],ylim,'Color','k','LineStyle','--')
title(['Array B: bout duration, median = ' num2str(median(boutdurB)) ' ms'])
subplot(2,1,2)
h=hist(interboutB,centers_interbout);
h=h/length(interboutB);
bar(centers_interbout,h)
xlim([0 21000])
box off
set(gca,'TickDir','out')
set(gca,'YTick',0:0.05:0.35)
set(gca,'XTick',0:5000:20000)
xlabel('Time [ms]');ylabel('Proportion of bouts')
line([median(interboutB) median(interboutB)],ylim,'Color','k','LineStyle','--')
title(['Array B: interbout interval, median = ' num2str(median(interboutB)) ' ms'])

figure
subplot(2,1,1)
h=hist(boutdurC,centers_boutdur);
h=h/length(boutdurC);
bar(centers_boutdur,h)
xlim([0 5250])
box off
set(gca,'TickDir','out')
% set(gca,'YTick',0:0.05:0.3)
set(gca,'XTick',0:1000:5000)
xlabel('Time [ms]');ylabel('Proportion of bouts')
line([median(boutdurC) median(boutdurC)],ylim,'Color','k','LineStyle','--')
title(['Array C: bout duration, median = ' num2str(median(boutdurC)) ' ms'])
subplot(2,1,2)
h=hist(interboutC,centers_interbout);
h=h/length(interboutC);
bar(centers_interbout,h)
xlim([0 21000])
box off
set(gca,'TickDir','out')
set(gca,'YTick',0:0.05:0.35)
set(gca,'XTick',0:5000:20000)
xlabel('Time [ms]');ylabel('Proportion of bouts')
line([median(interboutC) median(interboutC)],ylim,'Color','k','LineStyle','--')
title(['Array C: interbout interval, median = ' num2str(median(interboutC)) ' ms'])

