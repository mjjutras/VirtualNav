function [errorPercentile,chanceError] = calc_YC_error(objectLocation,responseLocation,nSamples)
% function error_percentile = calc_YC_error(objectLocation,responseLocation,nSamples)
%
% INPUTS
%
%         objectLocation: [x,y] coordinates of the object location
%       responseLocation: [x,y] coordinates of the response location
%               nSamples: number random points used to calculate the
%                         percentile. Default = 100000
%
% OUTPUTS
%
%       errorPercentile: percentile for the error, relative to all
%                        possible errors. 0 = best, 1 = worst
%   
%           chanceError: 50th percentile of all possible errors

% number of samples to use
if ~exist('nSamples','var') || isempty(nSamples)
    nSamples = 100000;
end

% set of random points to sample environment
x = 32.4 + (-32.4-32.4).*rand(nSamples,1);
y = 18 + (-18-18).*rand(nSamples,1);
randomPoints = [x y];

% difference between object location and all samples
possibleErrors = sqrt(sum((repmat(objectLocation,[size(randomPoints,1) 1]) - randomPoints).^2,2));

% calculate actual error
actErr = sqrt(sum((objectLocation - responseLocation).^2));

% percentile of error
errorPercentile = mean(possibleErrors<=actErr);

% 50th percentile (chance)
chanceError = prctile(possibleErrors,50);
