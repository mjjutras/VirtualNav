function curve = latencyLearningCurve(times, type, mode)

pickups = strcmp(type,'pickup') & strcmp(mode, 'seek');
pickup_times = times(pickups);
pickup_intervals = diff([seekStart(times,type,mode);pickup_times])/1000;

curve = pickup_intervals;
