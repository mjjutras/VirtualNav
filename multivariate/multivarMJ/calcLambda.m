function [stats,lambda] = calcLambda(X,Y,doBinary,nCV)
if isempty(nCV)
    nCV = round(length(Y)/2);
end

if doBinary
  
    [~,stats] = lassoglm(X,Y,'binomial','CV', nCV, 'NumLambda', 50);
    lambda    = stats.LambdaMinDeviance;
    
else
    % center x
    xHat = mean(X, 1)';
    xCentered = X' - xHat * ones(1,size(X',2));
    
    % center y
    intercept = mean(Y);
    yCentered = round(Y - intercept,14);
    
    % compute optimal lamda
    [~, stats] = lasso(xCentered', yCentered, 'CV', nCV, 'NumLambda', 50);
    lambda     = stats.LambdaMinMSE;
end


