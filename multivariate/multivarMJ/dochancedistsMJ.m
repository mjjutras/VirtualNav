% build powerData

load('varsformultivarGiz150904.mat')
freqBins = [1 3; 3 9; 40 70; 70 100];

[~,f1_1] = min(abs(ftL1.freq-freqBins(1,1)));
[~,f1_2] = min(abs(ftL1.freq-freqBins(1,2)));
[~,f2_1] = min(abs(ftL1.freq-freqBins(2,1)));
[~,f2_2] = min(abs(ftL1.freq-freqBins(2,2)));

[~,f3_1] = min(abs(ftH1.freq-freqBins(3,1)));
[~,f3_2] = min(abs(ftH1.freq-freqBins(3,2)));
[~,f4_1] = min(abs(ftH1.freq-freqBins(4,1)));
[~,f4_2] = min(abs(ftH1.freq-freqBins(4,2)));

clear powerData

powerData(:,1,1,:) = reshape(mean(ftL1.powspctrm(:,:,f1_1:f1_2),3),size(ftL1.powspctrm,1),1,1,size(ftL1.powspctrm,2));
powerData(:,2,1,:) = reshape(mean(ftL1.powspctrm(:,:,f2_1:f2_2),3),size(ftL1.powspctrm,1),1,1,size(ftL1.powspctrm,2));
powerData(:,3,1,:) = reshape(mean(ftH1.powspctrm(:,:,f3_1:f3_2),3),size(ftH1.powspctrm,1),1,1,size(ftH1.powspctrm,2));
powerData(:,4,1,:) = reshape(mean(ftH1.powspctrm(:,:,f4_1:f4_2),3),size(ftH1.powspctrm,1),1,1,size(ftH1.powspctrm,2));

powerData(:,1,2,:) = reshape(mean(ftL2.powspctrm(:,:,f1_1:f1_2),3),size(ftL2.powspctrm,1),1,1,size(ftL2.powspctrm,2));
powerData(:,2,2,:) = reshape(mean(ftL2.powspctrm(:,:,f2_1:f2_2),3),size(ftL2.powspctrm,1),1,1,size(ftL2.powspctrm,2));
powerData(:,3,2,:) = reshape(mean(ftH2.powspctrm(:,:,f3_1:f3_2),3),size(ftH2.powspctrm,1),1,1,size(ftH2.powspctrm,2));
powerData(:,4,2,:) = reshape(mean(ftH2.powspctrm(:,:,f4_1:f4_2),3),size(ftH2.powspctrm,1),1,1,size(ftH2.powspctrm,2));

powerData(:,1,3,:) = reshape(mean(ftL3.powspctrm(:,:,f1_1:f1_2),3),size(ftL3.powspctrm,1),1,1,size(ftL3.powspctrm,2));
powerData(:,2,3,:) = reshape(mean(ftL3.powspctrm(:,:,f2_1:f2_2),3),size(ftL3.powspctrm,1),1,1,size(ftL3.powspctrm,2));
powerData(:,3,3,:) = reshape(mean(ftH3.powspctrm(:,:,f3_1:f3_2),3),size(ftH3.powspctrm,1),1,1,size(ftH3.powspctrm,2));
powerData(:,4,3,:) = reshape(mean(ftH3.powspctrm(:,:,f4_1:f4_2),3),size(ftH3.powspctrm,1),1,1,size(ftH3.powspctrm,2));

powerData(:,1,4,:) = reshape(mean(ftL4.powspctrm(:,:,f1_1:f1_2),3),size(ftL4.powspctrm,1),1,1,size(ftL4.powspctrm,2));
powerData(:,2,4,:) = reshape(mean(ftL4.powspctrm(:,:,f2_1:f2_2),3),size(ftL4.powspctrm,1),1,1,size(ftL4.powspctrm,2));
powerData(:,3,4,:) = reshape(mean(ftH4.powspctrm(:,:,f3_1:f3_2),3),size(ftH4.powspctrm,1),1,1,size(ftH4.powspctrm,2));
powerData(:,4,4,:) = reshape(mean(ftH4.powspctrm(:,:,f4_1:f4_2),3),size(ftH4.powspctrm,1),1,1,size(ftH4.powspctrm,2));

powerData(:,1,5,:) = reshape(mean(ftL5.powspctrm(:,:,f1_1:f1_2),3),size(ftL5.powspctrm,1),1,1,size(ftL5.powspctrm,2));
powerData(:,2,5,:) = reshape(mean(ftL5.powspctrm(:,:,f2_1:f2_2),3),size(ftL5.powspctrm,1),1,1,size(ftL5.powspctrm,2));
powerData(:,3,5,:) = reshape(mean(ftH5.powspctrm(:,:,f3_1:f3_2),3),size(ftH5.powspctrm,1),1,1,size(ftH5.powspctrm,2));
powerData(:,4,5,:) = reshape(mean(ftH5.powspctrm(:,:,f4_1:f4_2),3),size(ftH5.powspctrm,1),1,1,size(ftH5.powspctrm,2));

powerData(:,1,6,:) = reshape(mean(ftL6.powspctrm(:,:,f1_1:f1_2),3),size(ftL6.powspctrm,1),1,1,size(ftL6.powspctrm,2));
powerData(:,2,6,:) = reshape(mean(ftL6.powspctrm(:,:,f2_1:f2_2),3),size(ftL6.powspctrm,1),1,1,size(ftL6.powspctrm,2));
powerData(:,3,6,:) = reshape(mean(ftH6.powspctrm(:,:,f3_1:f3_2),3),size(ftH6.powspctrm,1),1,1,size(ftH6.powspctrm,2));
powerData(:,4,6,:) = reshape(mean(ftH6.powspctrm(:,:,f4_1:f4_2),3),size(ftH6.powspctrm,1),1,1,size(ftH6.powspctrm,2));

[~,f1_1] = min(abs(ftL7.freq-freqBins(1,1)));
[~,f1_2] = min(abs(ftL7.freq-freqBins(1,2)));
[~,f2_1] = min(abs(ftL7.freq-freqBins(2,1)));
[~,f2_2] = min(abs(ftL7.freq-freqBins(2,2)));

[~,f3_1] = min(abs(ftH7.freq-freqBins(3,1)));
[~,f3_2] = min(abs(ftH7.freq-freqBins(3,2)));
[~,f4_1] = min(abs(ftH7.freq-freqBins(4,1)));
[~,f4_2] = min(abs(ftH7.freq-freqBins(4,2)));

powerData(:,1,7,:) = reshape(mean(ftL7.powspctrm(:,:,f1_1:f1_2),3),size(ftL7.powspctrm,1),1,1,size(ftL7.powspctrm,2));
powerData(:,2,7,:) = reshape(mean(ftL7.powspctrm(:,:,f2_1:f2_2),3),size(ftL7.powspctrm,1),1,1,size(ftL7.powspctrm,2));
powerData(:,3,7,:) = reshape(mean(ftH7.powspctrm(:,:,f3_1:f3_2),3),size(ftH7.powspctrm,1),1,1,size(ftH7.powspctrm,2));
powerData(:,4,7,:) = reshape(mean(ftH7.powspctrm(:,:,f4_1:f4_2),3),size(ftH7.powspctrm,1),1,1,size(ftH7.powspctrm,2));


% build other variables

% size of feature matrix
nElecs = size(powerData,4);
nTimes = size(powerData,3);
nFreqs = size(powerData,2);

params = multiParams();
params.freqBins = freqBins;
params.timeBins = params.timeBins(2:end,:);
params.timeBinLabels = {'Drive1' 'Drive2' 'Drive3' 'Drive4' 'Drive5' 'Post' 'Enc'};
    
trialInds = round((1:length(memsel)*2)/2);
nFolds = size(memsel,1);
folds = false(nFolds,size(powerData,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end

doBinary = params.doBinary;
nCV = params.nCV;

t = nTimes;
newpowdat = squeeze(powerData(:,:,t,:));

numIters = 100;

perf_all = nan(numIters,1);
auc_all  = nan(numIters,1);
parfor i = 1:numIters

    % response data
    trialInds = round((1:length(memsel)*2)/2);
    Y = memsel(trialInds);
    if doBinary
        Y  = Y > median(Y);
    end
    
    % permute the responses if desired
    randOrder = randperm(size(memsel,1));
    Ytmp = reshape(Y,2,[]);
    Y = reshape(Ytmp(:,randOrder),[],1);

%     perf = NaN(1,nTimes);
%     lambdas = NaN(1,nTimes);
%     AUC = NaN(1,nTimes);
    
        
    % reshape into # trials x # features
    X = reshape(newpowdat,size(powerData,1),nFreqs*nElecs);
    
    [~,lambda] = calcLambda(X,Y,doBinary,nCV);
    
    res = [];
    % will hold results from each fold
    [res.yPred,res.yTest,res.A,res.intercept,res.err] = deal(cell(nFolds,1));
    
    % run for each fold
    for iFold = 1:nFolds
%         fprintf('Time %d of %d, Fold %d of %d.\n',t,nTimes,iFold,nFolds)
        [res.yPred{iFold},...
            res.yTest{iFold},...
            res.A{iFold},...
            res.intercept{iFold},...
            res.err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda,nCV);
    end
    perf = mean(vertcat(res.err{:}));
    if doBinary
        yPred = vertcat(res.yPred{:});
        [~,~,~,AUC,~] = perfcurve(Y,yPred,true);
    end

    perf_all(i) = perf;
    auc_all(i) = AUC;
end

save('chancedists150905.mat','perf_all','auc_all')

