load('R:\Buffalo Lab\Mike\VirtualNav\multivariate\vars150816.mat')
load('R:\Buffalo Lab\Mike\VirtualNav\multivariate\R1008J_lasso.mat')

% config.distributedParams.freqBins:
% 50 frequencies logarithmically spaced between 1 and 200 Hz

% config.distributedParams.timeBins:
% 80 100-ms bins from -999 to 7000
          
          
% powerData:  #trials x #freqbins x #times x #elecs

% trials: all encoding (2 per test trial)

% params.freqBins
% ans =
%      1     3
%      3     9
%     40    70
%     70   200

% params.timeBins
% ans =
%         -999           0
%            1        1000
%         1001        2000
%         2001        3000
%         3001        4000
%         4001        5000
%         5001        6000