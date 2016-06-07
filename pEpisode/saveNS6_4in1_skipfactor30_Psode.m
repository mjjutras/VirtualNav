fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
savdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\NS6 - decimated';

d=dir([fildir '\*.ns6']);
clear filelist
for k=1:25
    
    filelist(k,:)=d(k).name;
    
    BRnam = filelist(k,1:11);
    
    disp(['Loading ' BRnam])
    NS6 = openNSx(fullfile(fildir,filelist(k,:)),'read','skipfactor',30);

    save(fullfile(savdir,[BRnam '_NS6_SF30.mat']),'NS6')
    disp(['Saved ' BRnam '_NS6_SF30.mat'])

end

%% p'sode

samplerate = 1000;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

decdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\NS6 - decimated';
pepdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\pEpisode';

d=dir([decdir '\*.mat']);

ft_progress('init', 'etf',     'Please wait...');
for fillop = 1:25

    ft_progress(fillop/25, 'Processing file %d of %d', fillop, 25);
    
    BRnam = d(fillop).name(1:11);
    
    load(fullfile(decdir,d(fillop).name));

%     disp(['Analyzing ' BRnam])
    
    Pm_all = cell(1,size(NS6.Data,1));
    pv_all = cell(1,size(NS6.Data,1));
%     B_all = cell(1,size(NS6.Data,1));
    unionvector_all = cell(1,size(NS6.Data,1));
    
    for chnlop = 1:size(NS6.Data,1)
        
%         disp(['Channel ' num2str(chnlop)])

%         % filter the LFP
%         raw = double(NS6.Data(chnlop,:));
%         Fbp=[1.7 100];
%         N=4;
%         Fn=1000/2;
%         [B,A]=butter(N, [min(Fbp)/Fn max(Fbp)/Fn]);
%         flt = filtfilt(B,A,raw);
        
        % get energy vector for each list

        B = multienergyvec(double(NS6.Data(chnlop,:)),freqs,samplerate,width);
%         B = multienergyvec(flt,freqs,samplerate,width);

        % calc the mean fit
        Pm = log10(mean(B,2));

        % get the fit
        % IMPORTANS: chi_squarefit assumes that frequencies are
        % logarithmically spaced!
        [all,R2] = chi_squarefit(freqs,Pm);
        all = all';

        % test fit of regression line
        pv=polyfit(log10(freqs),Pm',1); % linear regression
    %     figure;plot(freqs,Pm)
    %     hold on;plot(freqs,pv(2)+pv(1)*log10(freqs),'r')

        % set the threshold
        thresh = all(:,951);

        % loop through each frequency and save the unionvector
        unionvector = single(zeros(size(freqs,2),size(B,2)));

        for f = 1:size(freqs,2)
            % get the pepisode
            unionvector(f,:) = single(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0));
        end

        % convert unionvector from single to int8 to conserve memory
        unionvector=int8(unionvector);

        Pm_all{chnlop}=Pm;
        pv_all{chnlop}=pv;
%         B_all{chnlop}=B;
        unionvector_all{chnlop}=unionvector;
        
    end
    
    save(fullfile(pepdir,[BRnam '_pEpisode.mat']),'unionvector_all','Pm_all','pv_all','freqs','-v7.3')
    
end
ft_progress('close')

%% all array

pepdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MAT files\pEpisode';

d=dir([pepdir '\*.mat']);

namarr = [];
Pmarr = [];
pvarr = [];
UVarr = [];

for fillop = 1:25

    load(fullfile(pepdir,d(fillop).name));

    BRnam = d(fillop).name(1:11);
    
    for chnlop = 1:length(Pm_all)
    
        namarr = [namarr; BRnam '.' num2str(chnlop)];
        Pmarr = [Pmarr; Pm_all{chnlop}'];
        pvarr = [pvarr; pv_all{chnlop}];
        UVarr = [UVarr; mean(unionvector_all{chnlop},2)'];
        
    end
    
end


    
