function [yPred,yTest,A,intercept,err] = lassoReg(X,Y,trainInds,lambda,nCV)
%
% This does lasso. 
% X = # trials x # features
% Y = # trials vector of responses
% trainInds = logical vector of training/test
%

if isempty(nCV)
    nCV = round(length(Y)/2);
end

doBinary = false;
if islogical(Y)
    doBinary = true;
end

% We will do binary classification Y is logical, which is currently the
% default
if doBinary
    
    % I'm under sampling the larger class so that we have equal numbers.
    % This isn't a great idea of more skewed dataset, but since we are
    % using a median threshold, it doesn't really matter here
    yTrainBool = Y(trainInds);
    
    % figure out which observations to remove from larger class
    numToRemove = sum(yTrainBool) - sum(~yTrainBool);
    toRemove = [];
    if numToRemove > 0
        toRemove = randsample(find(yTrainBool),abs(numToRemove));
    elseif numToRemove < 0
        toRemove = randsample(find(~yTrainBool),abs(numToRemove));
    end
    
    % training set x
    xTrain = X(trainInds,:)';
    xTrain(:,toRemove) = [];
    
    % training set y
    yTrainBool(toRemove) = [];
    
    % compute model
    if isempty(lambda)
        [A_lasso, stats] = lassoglm(xTrain',yTrainBool,'binomial','CV', nCV, 'NumLambda', 50);
        
        % get the best cooefficients and intercept
        A = A_lasso(:,stats.IndexMinDeviance);
        intercept = stats.Intercept(stats.IndexMinDeviance);
    else
        [A, stats] = lassoglm(xTrain',yTrainBool,'binomial','Lambda',lambda);
        intercept = stats.Intercept;
    end
    
    % testing set
    xTest = X(~trainInds,:)';
    yTest = Y(~trainInds);
    
    % predict
    B1 = [intercept;A];
    yPred = glmval(B1,xTest','logit');
    
    % see if the predictions match the actual results
    err = (yPred > mean(yTrainBool)) == Y(~trainInds);
    
% if Y is not logical, do regression
else
    
    % training set x
    xTrain = X(trainInds,:)';    
    xHat = mean(xTrain, 2);
    xTrain = xTrain - xHat * ones(1,size(xTrain,2));
    
    % training set y
    yTrain = Y(trainInds);
    intercept = mean(yTrain);
    yTrain = round(yTrain - intercept,14);
    
    % compute model
    if isempty(lambda)
        [A_lasso, stats] = lasso(xTrain', yTrain, 'CV', nCV, 'NumLambda', 50);
        A = A_lasso(:,stats.IndexMinMSE);
    else
        A = lasso(xTrain', yTrain, 'Lambda',lambda);
    end
    
    % testing set
    xTest = X(~trainInds,:)';
    yTest = Y(~trainInds);
    
    % Double check this
    yPred = (xTest - xHat*ones(1,sum(~trainInds)))' * A + intercept;
    err = mean((yTest - yPred).^2);
           
end


