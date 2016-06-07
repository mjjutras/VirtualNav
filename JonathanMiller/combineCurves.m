function [curveFile, histFile] = combineCurves(curves, name, xAxis, yAxis, pathstr)

meanCurve = zeros(size(curves, 1),1);
for i=1:length(meanCurve)
  meanCurve(i) = nanmean(curves(i,:));
end
trials = transpose(1:length(meanCurve));
p = polyfit(trials,meanCurve,1);
r = polyval(p, trials);

[rho, pval] = corr(trials, meanCurve);

scatter(trials,meanCurve);

hold on;

xlabel(xAxis);
ylabel(yAxis);
plot(trials, r);
text(3*max(trials)/4, max(meanCurve),0, sprintf('r = %0.3f  p = %0.3f',rho, pval));


curveFile = sprintf('%s.eps',name);

print(sprintf('%s/%s',pathstr,curveFile));


hold off;

hist(meanCurve,25);


histFile = sprintf('%s_hist.eps',name);

print(sprintf('%s/%s',pathstr,histFile));


print(histFile);

end
