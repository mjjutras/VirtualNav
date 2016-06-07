function curve = blockedLatencyLearningCurve(times, type, mode)

pickups = strcmp(type,'pickup') & strcmp(mode, 'seek');
pickup_times = times(pickups);
pickup_intervals = [diff(pickup_times)/1000];

fm = @(block_struct) sum(block_struct.data);
blocked_pickup_intervals = blockproc(pickup_intervals, [4,1],fm);

curve = blocked_pickup_intervals;
