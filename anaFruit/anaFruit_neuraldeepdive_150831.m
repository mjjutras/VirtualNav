%% compare neural signals across recordings to see whether any changes occur that might
%% impact the analysis, such as changes in power or coherence
%% focus on the 5-second period preceding the acquisition of visible bananas (encoding
%% phase), [i.e. during the approach path] since this is the period we have been focusing
%% on so far; may revise this in the future

%%

NSDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NSdat';
trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';


%%

seslst = {'JN_BR_15_04_10\JN_BR_15_04_10_13_03';
    'JN_BR_15_04_10\JN_BR_15_04_10_13_11';
    'JN_BR_15_04_13\JN_BR_15_04_13_13_59';
    'JN_BR_15_04_14\JN_BR_15_04_14_13_34';
    'JN_BR_15_04_15\JN_BR_15_04_15_12_32';
    'JN_BR_15_04_16\JN_BR_15_04_16_14_22';
    'JN_BR_15_04_17\JN_BR_15_04_17_13_14';
    'JN_BR_15_04_20\JN_BR_15_04_20_14_43';
    'JN_BR_15_04_21\JN_BR_15_04_21_12_34';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_09';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_54';
    'JN_BR_15_04_23\JN_BR_15_04_23_13_54';
    'JN_BR_15_04_23\JN_BR_15_04_23_14_22';
    'JN_BR_15_04_24\JN_BR_15_04_24_14_04';
    'JN_BR_15_04_28\JN_BR_15_04_28_13_07';
    'JN_BR_15_04_29\JN_BR_15_04_29_13_22';
    'JN_BR_15_04_30\JN_BR_15_04_30_14_42';
    'JN_BR_15_05_01\JN_BR_15_05_01_13_12';
    'JN_BR_15_05_04\JN_BR_15_05_04_13_57';
    'JN_BR_15_05_04\JN_BR_15_05_04_14_51';
    'JN_BR_15_05_05\JN_BR_15_05_05_13_02';
    'JN_BR_15_05_06\JN_BR_15_05_06_13_30';
    'JN_BR_15_05_07\JN_BR_15_05_07_12_59';
    'JN_BR_15_05_08\JN_BR_15_05_08_14_17';
    'JN_BR_15_05_12\JN_BR_15_05_12_12_49';
    'JN_BR_15_05_13\JN_BR_15_05_13_13_14';
    'JN_BR_15_05_14\JN_BR_15_05_14_13_34';
    'JN_BR_15_05_18\JN_BR_15_05_18_13_03';
    'JN_BR_15_05_19\JN_BR_15_05_19_13_00';
    'JN_BR_15_05_20\JN_BR_15_05_20_13_04';
    'JN_BR_15_05_22\JN_BR_15_05_22_13_20';
    'JN_BR_15_05_26\JN_BR_15_05_26_13_42';
    'JN_BR_15_05_27\JN_BR_15_05_27_12_53';
    'JN_BR_15_05_28\JN_BR_15_05_28_13_59';
    'JN_BR_15_05_28\JN_BR_15_05_28_14_16';
    'JN_BR_15_05_29\JN_BR_15_05_29_14_19';
    'JN_BR_15_06_01\JN_BR_15_06_01_14_35';
    'JN_BR_15_06_02\JN_BR_15_06_02_13_15';
    'JN_BR_15_06_03\JN_BR_15_06_03_14_32';
    'JN_BR_15_06_04\JN_BR_15_06_04_13_04';
    'JN_BR_15_06_05\JN_BR_15_06_05_13_35'};


%% define variables for neural analysis

badtrlall = cell(size(seslst));
badchnall = cell(size(seslst));

numtrl = []; % total number of trials per session
numnrltrl = []; % number of trials included in neural analysis per session

powLmat = [];
powHmat = [];

cohLmat = [];
cohHmat = [];

%%

for seslop = 1:length(seslst)
    
    numtrl(seslop,1) = length(trldat.time);
    
    [~,sesnam]=fileparts(seslst{seslop});
    disp(['Processing ' sesnam])

    load(fullfile(NSDir,[sesnam '_NSdat.mat']))
    load(fullfile(trlDir,[sesnam '_trldat.mat']))
    
    bantimindses = zeros(length(trldat.posdat),1);
    banfndses =  zeros(length(trldat.posdat),1);
    trlalp = nan(length(trldat.posdat),1); % banana alpha per trial
    for trllop = 1:length(trldat.posdat)

        trlalp(trllop,1) = trldat.alpha{trllop}(1,2);

        % find time when banana eaten, or end of trial if banana not eaten
        if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
            [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
            banfndses(trllop) = 1;
        else
            % NOTE: if he doesn't find the banana, the actual time here
            % should be the "TIMEOUT", which is either 90 seconds or some
            % other longer amount (before 5/8/15)
            bantimind = length(trldat.time{trllop});
            banfndses(trllop) = 0;
        end
        
        bantimindses(trllop) = bantimind;

    end

    % select trials for neural analysis (nrlsel)
    nrlsel = find(trlalp~=0 & banfndses==1 & bantimindses>5000);
    bantimsel = bantimindses(nrlsel);
    
    % initialize dum structure
    clear dum
    dum.time = data.time(nrlsel);
    dum.trial = data.trial(nrlsel);
    dum.sampleinfo = data.sampleinfo(nrlsel,:);
    dum.fsample = data.fsample;
    dum.label = data.label;
    
    % exclude eye and position data
    excchn = {'eyeX_B'; 'eyeY_B'; 'posX_B'; 'posY_B'; 'eyeX_P'; 'eyeY_P'; 'posX_P'; 'posY_P'};
    
    clear cfg
    cfg.channel       = setxor(dum.label,excchn);
    data = ft_preprocessing(cfg,dum);
    clear dum

    % run automated data cleaning
    
    % concatenate variance and kurtosis for both data structures
    varall = nan(length(data.label),length(data.trial));
    for trllop = 1:length(data.trial)
        varall(:,trllop) = var(data.trial{trllop},0,2);
    end
    
    % set weight for determining outliers
    % default for boxplot is 1.5
    % adjust upward to be more lenient, downward to be more aggressive
    % the two variables will balance each other; in general, adjusting one will
    % cause the other to effectively compensate in the other direction
    w_t = 2.3; % weight for trials
    w_c = 8.0; % weight for channels
    
    % find outlying trials using max variance for each trial (across channels)
    trlthresh1 = quantile(max(varall),0.75) + w_t*(quantile(max(varall),0.75) - quantile(max(varall),0.25));
    trlselvardum = find(max(varall)>trlthresh1);
    
    % find outlying channels using max variance for each channel (across trials)
    chnthresh1 = quantile(max(varall,[],2),0.75) + w_c*(quantile(max(varall,[],2),0.75) - quantile(max(varall,[],2),0.25));
    chnselvardum = find(max(varall,[],2)>chnthresh1);
    
    trlselvar = [];
    chnselvar = [];
    varalldum = varall;
    while(1)
        while(~isempty(trlselvardum))
            % remove the furthest outlying trial in each iteration
            [~,maxvartrlind] = max(max(varalldum(:,trlselvardum)));
            varalldum(:,trlselvardum(maxvartrlind)) = nan(size(varalldum,1),1);
            chnthreshnew = quantile(max(varalldum,[],2),0.75) + w_c*(quantile(max(varalldum,[],2),0.75) - quantile(max(varalldum,[],2),0.25));
            [~,chnrej] = setdiff(chnselvardum,find(max(varalldum,[],2)>chnthreshnew));
            % if removing this trial changed the number of outlying channels, then
            % set this trial aside for possible permanent exclusion
            if ~isempty(chnrej)
                chnselvardum = chnselvardum(setxor(1:length(chnselvardum),chnrej));
                trlselvar = [trlselvar trlselvardum(maxvartrlind)];
                disp(['Adding trial #' num2str(trlselvardum(maxvartrlind)) ' to trlselvar'])
            end
            trlselvardum = trlselvardum(setxor(1:length(trlselvardum),maxvartrlind));
        end
        
        % after removing all outlying trials, check to see if any channels are still
        % outlying; if so, set the furthest outlying channel aside for exclusion,
        % then clear the temporary trial exclusion list
        if ~isempty(chnselvardum)
            [~,maxvarchnind] = max(max(varalldum(chnselvardum,:),[],2));
            chnselvar = [chnselvar; chnselvardum(maxvarchnind)];
            disp(['Adding channel #' num2str(chnselvardum(maxvarchnind)) ' to chnselvar'])
            disp('Clearing trlselvar')
            trlselvar = [];
        end
        
        % reset the variance measures; keep permanent excluded trials and channels out
        varalldum = varall;
        varalldum(:,trlselvar) = nan(size(varalldum,1),length(trlselvar));
        varalldum(chnselvar,:) = nan(length(chnselvar),size(varalldum,2));
        
        % check again for outlying channels; if there are none, then break
        chnthreshdum = quantile(max(varalldum,[],2),0.75) + w_c*(quantile(max(varalldum,[],2),0.75) - quantile(max(varalldum,[],2),0.25));
        chnselvardum = find(max(varalldum,[],2)>chnthreshdum);
        if isempty(chnselvardum)
            break
        end
        % with outlying channels excluded, repopulate temporary outlying trial list
        trlthreshdum = quantile(max(varalldum),0.75) + w_t*(quantile(max(varalldum),0.75) - quantile(max(varalldum),0.25));
        trlselvardum = find(max(varalldum)>trlthreshdum);
    end
    trlselvar = sort(trlselvar); % trials to exclude based on variance
    chnselvar = sort(chnselvar); % channels to exclude based on variance
    
    badtrlall{seslop} = trlselvar;
    badchnall{seslop} = chnselvar;
    
    % preprocess in Fieldtrip to remove the bad trials/channels
    cfg=[];
    cfg.channel = data.label(setxor(1:length(data.label),chnselvar));
    cfg.trials = setxor(1:length(data.trial),trlselvar);
    data0 = ft_preprocessing(cfg,data);
    clear data
    bantimsel = bantimsel(cfg.trials);
    
    timsel = [bantimsel-4999 bantimsel];

    data1 = data0;
    data1TF = data0;
    for trllop = 1:length(data0.trial)
        % data1: time-averaged
        data1.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1];
        data1.time{trllop} = -4.999:0.001:0;
        data1.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2));
        
        % data1TF: time-frequency analysis (padding at the beginning and end)
        data1TF.sampleinfo(trllop,:) = [data0.sampleinfo(trllop,1)+timsel(trllop,1)-1 data0.sampleinfo(trllop,1)+timsel(trllop,2)-1+1500];
        data1TF.time{trllop} = -4.999:0.001:1.5;
        data1TF.trial{trllop} = data0.trial{trllop}(:,timsel(trllop,1):timsel(trllop,2)+1500);
    end
    clear data0

    % filter out line noise
    cfg=[];
    cfg.dftfilter     = 'yes';
    cfg.dftfreq       = 60;
    data2 = ft_preprocessing(cfg,data1);
    data2TF = ft_preprocessing(cfg,data1TF);
    clear data1*

    % spectral analysis

    % do the spectral analysis - time-averaged
    cfg=[];
    cfg.output      = 'fourier'; % specify 'fourier' to get chan_chan_freq in coherence
    % cfg.output      = 'powandcsd';
    cfg.method      = 'mtmfft';
    cfg.pad         = 'maxperlen';
    cfg.keeptrials  = 'yes';
    cfg.foilim      = [1 30];
    % cfg.taper       = 'hanning';
    cfg.taper       = 'dpss';
    cfg.tapsmofrq   = 2;
    freqL = ft_freqanalysis(cfg, data2);
    fdL = ft_freqdescriptives(cfg, freqL);

    cfg.foilim      = [30 200];
    cfg.taper       = 'dpss';
    cfg.tapsmofrq   = 8;
    freqH = ft_freqanalysis(cfg, data2);
    fdH = ft_freqdescriptives(cfg, freqH);

    cfg=[];
    cfg.method      = 'coh';
    cfg.jackknife   = 'yes';
    statL = ft_connectivityanalysis(cfg,freqL);
    statH = ft_connectivityanalysis(cfg,freqH);

    clear prblaball
    for k=1:length(statL.label)
        prblaball(k,1)=statL.label{k}(1);
    end
    prblab = unique(prblaball);
    [~,locb] = ismember(prblaball,prblab);

    cmblst = unique(sort([repmat(unique(locb),length(unique(locb)),1) sort(repmat(unique(locb),length(unique(locb)),1))],2),'rows');

    powL = [];
    powH = [];
    for prblop = 1:length(prblab)
        ind = find(locb==prblop);
        powL(prblop,:) = squeeze(mean(mean(fdL.powspctrm(:,ind,:),1),2))';
        powH(prblop,:) = squeeze(mean(mean(fdH.powspctrm(:,ind,:),1),2))';
    end

    cohL = [];
    cohH = [];
    for cmblop = 1:size(cmblst,1)

        cmbind1 = find(locb==cmblst(cmblop,1));
        cmbind2 = find(locb==cmblst(cmblop,2));
        cmblst2 = unique(sort([repmat(unique(cmbind1),length(unique(cmbind2)),1) sort(repmat(unique(cmbind2),length(unique(cmbind1)),1))],2),'rows');
        cmblst2 = cmblst2(cmblst2(:,1)~=cmblst2(:,2),:);

        cohmat1 = []; cohmat2 = [];
        for indlop1 = 1:size(cmblst2,1)
            cohmat1 = cat(1,cohmat1,(squeeze(statL.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
            cohmat2 = cat(1,cohmat2,(squeeze(statH.cohspctrm(cmblst2(indlop1,1),cmblst2(indlop1,2),:)))');
        end
        
        cohL(cmblop,:) = mean(cohmat1,1);
        cohH(cmblop,:) = mean(cohmat2,1);

    end

    powLmat(seslop,:,:) = powL;
    powHmat(seslop,:,:) = powH;
    
    cohLmat(seslop,:,:) = cohL;
    cohHmat(seslop,:,:) = cohH;

    
    % do the spectral analysis - time-frequency
    clear cfg
    cfg.output      = 'fourier';
    cfg.method      = 'mtmconvol';
    cfg.foi         = 1:2:30;
    numfoi          = length(cfg.foi);
    cfg.t_ftimwin   = 0.5 * ones(1,numfoi);
    cfg.tapsmofrq   = 4   * ones(1,numfoi);
    cfg.toi         = -4.75:.01:1.25;
    cfg.pad         = 'maxperlen';
    cfg.keeptrials  = 'yes';
    freqTFL = ft_freqanalysis(cfg, data2TFL);
    fdTFL = ft_freqdescriptives(cfg, freqTFL);
    
    cfg=[];
    cfg.method      = 'coh';
    statTFL = ft_connectivityanalysis(cfg,freqTFL);

    numnrltrl(seslop,1) = length(data2.trial);
    
end

    
