function timestamp = seekStart(times, type, mode, path)
  training = strcmp(type,'pickup') & strcmp(mode, 'train_i');
  trainingNumerical = find(training);
  timestamp = times(trainingNumerical(end));
