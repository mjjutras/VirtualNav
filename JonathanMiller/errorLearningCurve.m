function curve = errorLearningCurve(times, type, mode, x, y)

% Distances between each movement
distances = sqrt(diff(x).^2 + diff(y).^2);

pickups = strcmp(type,'pickup') & strcmp(mode, 'seek');
pickups_numerical = find(pickups);

% Initialize vector for the total distance between each pickup
pickup_distances = zeros(length(pickups_numerical)-1,1);

% Fill in the total distances between each pickup
for i=1:length(pickup_distances)
  pickup_distances(i) = sum(distances(pickups_numerical(i):pickups_numerical(i+1)));
end

% Straight line distance between each pickup
pickup_true_distances = sqrt(diff(x(pickups)).^2 + diff(y(pickups)).^2);

% Percent error of actual distance traveled from the true distance
%distance_errors = 100*(pickup_distances - pickup_true_distances)./pickup_true_distances;

% Excess distance
distance_errors = pickup_distances - pickup_true_distances;
	
curve = distance_errors;
