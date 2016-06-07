%%

% powerData: in example R1061T, 96 x 4 x 8 x 108
% 96: number of trials
% 4: number of frequency bins (nFreqs)
% 8: number of time bins (nTimes)
% 108: number of channels (nElecs)

% Y: response data (1 or 0 when binary)
% always same for same block ("learning pair")


%% build powerData

load('R:\Buffalo Lab\Mike\VirtualNav\multivariate\varsformultivarGiz150904.mat')
freqBins = [1 3; 3 9; 30 60; 60 100];

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

% % take out this group, doesn't work
% powerData(:,1,6,:) = reshape(mean(ftL6.powspctrm(:,:,f1_1:f1_2),3),size(ftL6.powspctrm,1),1,1,size(ftL6.powspctrm,2));
% powerData(:,2,6,:) = reshape(mean(ftL6.powspctrm(:,:,f2_1:f2_2),3),size(ftL6.powspctrm,1),1,1,size(ftL6.powspctrm,2));
% powerData(:,3,6,:) = reshape(mean(ftH6.powspctrm(:,:,f3_1:f3_2),3),size(ftH6.powspctrm,1),1,1,size(ftH6.powspctrm,2));
% powerData(:,4,6,:) = reshape(mean(ftH6.powspctrm(:,:,f4_1:f4_2),3),size(ftH6.powspctrm,1),1,1,size(ftH6.powspctrm,2));

[~,f1_1] = min(abs(ftL7.freq-freqBins(1,1)));
[~,f1_2] = min(abs(ftL7.freq-freqBins(1,2)));
[~,f2_1] = min(abs(ftL7.freq-freqBins(2,1)));
[~,f2_2] = min(abs(ftL7.freq-freqBins(2,2)));

[~,f3_1] = min(abs(ftH7.freq-freqBins(3,1)));
[~,f3_2] = min(abs(ftH7.freq-freqBins(3,2)));
[~,f4_1] = min(abs(ftH7.freq-freqBins(4,1)));
[~,f4_2] = min(abs(ftH7.freq-freqBins(4,2)));

powerData(:,1,6,:) = reshape(mean(ftL7.powspctrm(:,:,f1_1:f1_2),3),size(ftL7.powspctrm,1),1,1,size(ftL7.powspctrm,2));
powerData(:,2,6,:) = reshape(mean(ftL7.powspctrm(:,:,f2_1:f2_2),3),size(ftL7.powspctrm,1),1,1,size(ftL7.powspctrm,2));
powerData(:,3,6,:) = reshape(mean(ftH7.powspctrm(:,:,f3_1:f3_2),3),size(ftH7.powspctrm,1),1,1,size(ftH7.powspctrm,2));
powerData(:,4,6,:) = reshape(mean(ftH7.powspctrm(:,:,f4_1:f4_2),3),size(ftH7.powspctrm,1),1,1,size(ftH7.powspctrm,2));


%% build other variables

% size of feature matrix
nElecs = size(powerData,4);
nTimes = size(powerData,3);
nFreqs = size(powerData,2);

params = multiParams();
params.freqBins = freqBins;
params.timeBins = params.timeBins(2:end,:);
% params.timeBinLabels = {'Drive1' 'Drive2' 'Drive3' 'Drive4' 'Drive5' 'Post' 'Enc'};
params.timeBinLabels = {'Drive1' 'Drive2' 'Drive3' 'Drive4' 'Drive5' 'Enc'};
   
% response data
trialInds = round((1:length(memsel)*2)/2);
Y = memsel(trialInds);
if params.doBinary
    Y  = Y > median(Y);
end

nFolds = size(memsel,1);
folds = false(nFolds,size(powerData,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end

%%

% We can model time points seperately, so # features = # freqs x # elecs,
% or we can model it all together, so # features = # times x # freqs x #
% elecs.
res = [];
perf = NaN(1,nTimes);
lambdas = NaN(1,nTimes);
for t = 1:nTimes
    
    % reshape into # trials x # features
    X = reshape(squeeze(powerData(:,:,t,:)),size(powerData,1),nFreqs*nElecs);
    
    % see if we are precomputing lambda based on all the data
    lambda = [];
    if params.crossValStrictness == 0
        if ~isempty(params.lambda)
            lambda = params.lambda(t);
        else
            fprintf('Computing optimal lambda.\n')
            t0 = clock;
            [stats,lambda] = calcLambda(X,Y,params.doBinary,params.nCV);
            fprintf('%g seconds to run calcLambda\n',etime(clock,t0));
        end
        lambdas(t) = lambda;
    end
    
    % will hold results from each fold
    [res(t).yPred,res(t).yTest,res(t).A,res(t).intercept,res(t).err] = deal(cell(nFolds,1));
    
    % run for each fold
    for iFold = 1:nFolds
        fprintf('Time %d of %d, Fold %d of %d.\n',t,nTimes,iFold,nFolds)
        [res(t).yPred{iFold},...
            res(t).yTest{iFold},...
            res(t).A{iFold},...
            res(t).intercept{iFold},...
            res(t).err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda,params.nCV);
    end
    perf(t) = mean(vertcat(res(t).err{:}));
    res(t).perf = perf(t);
    if params.doBinary
        yPred = vertcat(res(t).yPred{:});
        [res(t).xAUC,res(t).yAUC,~,res(t).AUC,res(t).optPoint] = perfcurve(Y,yPred,true);
        AUC(t) = res(t).AUC;
    end
end
lambda = lambdas;
