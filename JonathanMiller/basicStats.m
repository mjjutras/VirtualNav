function [trainingTrials, meanTraining, seekTrials, meanSeek] = basicStats(times,type,mode)


training = strcmp(type, 'pickup') & (strcmp(mode, 'train_i') | strcmp(mode, 'train_v'));
trainingTimes = times(training);
trainingIntervals = diff(trainingTimes);
meanTraining = mean(trainingIntervals)/1000;

trainingTrials = length(trainingTimes);


seeks = strcmp(type,'pickup') & strcmp(mode, 'seek');
seekTimes = times(seeks);
seekIntervals = diff(seekTimes);
meanSeek = mean(seekIntervals)/1000;

seekTrials = length(seekTimes);
