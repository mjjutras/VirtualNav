% load cohData

% datdir = 'C:\Data\VR';
datdir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\workspace';

load(fullfile(datdir,'cohData150908.mat'))
load(fullfile(datdir,'memsel150908.mat'))


%%

% remove array-array pairs
% A-A
remind = intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2)));
% A-B
remind = union(intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2))), ...
    intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2))));
% A-C
remind = union(intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2))), ...
    intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2))));
% B-B
remind = intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2)));
% B-C
remind = union(intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2))), ...
    intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2))));
% C-C
remind = intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2)));

cohData = cohData(:,:,setxor(1:size(cohData,3),remind));

% select only array-array pairs with certain C channels
% spacing of 1
selind = union(find(ismember(statL.labelcmb(:,1),{'C03' 'C04' 'C05' 'C06'})), ...
    find(ismember(statL.labelcmb(:,2),{'C03' 'C04' 'C05' 'C06'})));
selind = union(find(ismember(statL.labelcmb(:,1),{'C07' 'C08' 'C09' 'C10'})), ...
    find(ismember(statL.labelcmb(:,2),{'C07' 'C08' 'C09' 'C10'})));
selind = union(find(ismember(statL.labelcmb(:,1),{'C01' 'C02' 'C03' 'C04'})), ...
    find(ismember(statL.labelcmb(:,2),{'C01' 'C02' 'C03' 'C04'})));
% spacing of 2
selind = union(find(ismember(statL.labelcmb(:,1),{'C01' 'C03' 'C05' 'C07'})), ...
    find(ismember(statL.labelcmb(:,2),{'C01' 'C03' 'C05' 'C07'})));
selind = union(find(ismember(statL.labelcmb(:,1),{'C02' 'C04' 'C06' 'C08'})), ...
    find(ismember(statL.labelcmb(:,2),{'C02' 'C04' 'C06' 'C08'})));
selind = union(find(ismember(statL.labelcmb(:,1),{'C05' 'C07' 'C09' 'C11'})), ...
    find(ismember(statL.labelcmb(:,2),{'C05' 'C07' 'C09' 'C11'})));
% spacing of 3
selind = union(find(ismember(statL.labelcmb(:,1),{'C01' 'C04' 'C07' 'C10'})), ...
    find(ismember(statL.labelcmb(:,2),{'C01' 'C04' 'C07' 'C10'})));
selind = union(find(ismember(statL.labelcmb(:,1),{'C02' 'C05' 'C08' 'C11'})), ...
    find(ismember(statL.labelcmb(:,2),{'C02' 'C05' 'C08' 'C11'})));
selind = union(find(ismember(statL.labelcmb(:,1),{'C03' 'C06' 'C09' 'C12'})), ...
    find(ismember(statL.labelcmb(:,2),{'C03' 'C06' 'C09' 'C12'})));

cohData = cohData(:,:,selind);


%% build other variables

% size of feature matrix
nElecs = size(cohData,3);
nTimes = 1;
nFreqs = size(cohData,2);

params = multiParams();
params.freqBins = [3 12; 30 60; 60 100];
params.timeBins = params.timeBins(end,:);
params.timeBinLabels = {'Enc'};
   
% response data
trialInds = round((1:length(memsel)*2)/2);
Y = memsel(trialInds);
if params.doBinary
    Y  = Y > median(Y);
end

nFolds = size(memsel,1);
folds = false(nFolds,size(cohData,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end

%% add fourth series to "freq band" dimension containing array-pair categories

load(fullfile(datdir,'statL_150916.mat'))

% A-A
indAA = intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2)));
% A-B
indAB = union(intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2))), ...
    intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2))));
% A-C
indAC = union(intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2))), ...
    intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2))));
% B-B
indBB = intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2)));
% B-C
indBC = union(intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2))), ...
    intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2))));
% C-C
indCC = intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2)));

cohData(:,4,:) = nan(size(cohData,1),1,size(cohData,3));
cohData(:,4,indAA) = 1;
cohData(:,4,indAB) = 2;
cohData(:,4,indAC) = 3;
cohData(:,4,indBB) = 4;
cohData(:,4,indBC) = 5;
cohData(:,4,indCC) = 6;

nFreqs=nFreqs+1;


%%

% We can model time points seperately, so # features = # freqs x # elecs,
% or we can model it all together, so # features = # times x # freqs x #
% elecs.
res = [];
perf = NaN(1,nTimes);
lambdas = NaN(1,nTimes);
for t = 1:nTimes
    
    % reshape into # trials x # features
    X = reshape(cohData,size(cohData,1),nFreqs*nElecs);
    
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

% save('C:\Data\VR\multivarres_coh1_150908.mat','perf','lambda','res','perf','AUC')


%% get chance distributions

numIters = 1000;

perf_all = nan(numIters,1);
auc_all  = nan(numIters,1);
for i = 1:numIters
    
    fprintf('Running iteration #%d\n',i)

    % response data
    trialInds = round((1:length(memsel)*2)/2);
    Y = memsel(trialInds);
    if params.doBinary
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
    X = reshape(cohData,size(cohData,1),nFreqs*nElecs);
    
    fprintf('Computing optimal lambda.\n')
    t0 = clock;
    [~,lambda] = calcLambda(X,Y,params.doBinary,params.nCV);
    fprintf('%g seconds to run calcLambda\n',etime(clock,t0));
    
    res = [];
    % will hold results from each fold
    [res.yPred,res.yTest,res.A,res.intercept,res.err] = deal(cell(nFolds,1));
    
    % run for each fold
    for iFold = 1:nFolds
        fprintf('Fold %d of %d.\n',iFold,nFolds)
        [res.yPred{iFold},...
            res.yTest{iFold},...
            res.A{iFold},...
            res.intercept{iFold},...
            res.err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda,params.nCV);
    end
    perf = mean(vertcat(res.err{:}));
    if params.doBinary
        yPred = vertcat(res.yPred{:});
        [~,~,~,AUC,~] = perfcurve(Y,yPred,true);
    end

    perf_all(i) = perf;
    auc_all(i) = AUC;
end


%% top/bottom third

load('C:\Data\VR\cohData150908.mat')
load('C:\Data\VR\memsel150908.mat')

[~,memsrt] = sort(memsel);

memindhi = memsrt(end-(round(length(memsel)/3)-1):end);
memindlo = memsrt(1:round(length(memsel)/3));
nrlmemindhi = sort([memindhi*2-1; memindhi*2]);
nrlmemindlo = sort([memindlo*2-1; memindlo*2]);

cohData = cohData([nrlmemindhi; nrlmemindlo],:,:);
memsel = memsel([memindhi; memindlo]);


% size of feature matrix
nElecs = size(cohData,3);
nTimes = 1;
nFreqs = size(cohData,2);

params = multiParams();
params.freqBins = [3 12; 30 60; 60 100];
params.timeBins = params.timeBins(end,:);
params.timeBinLabels = {'Enc'};
   
% response data
trialInds = round((1:length(memsel)*2)/2);
Y = memsel(trialInds);
if params.doBinary
    Y  = Y > median(Y);
end

nFolds = size(memsel,1);
folds = false(nFolds,size(cohData,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end

res = [];
perf = NaN(1,nTimes);
lambdas = NaN(1,nTimes);
for t = 1:nTimes
    
    % reshape into # trials x # features
    X = reshape(cohData,size(cohData,1),nFreqs*nElecs);
    
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

