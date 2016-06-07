function [file] = occupancy(times, type, mode, x, y, speed, pathstr)

nBins=20;

dtSec=[diff(times);0]/1000;

%exclude periods of stillness
good=abs(speed)>.1 & strcmp(mode, 'seek') & strcmp(type,'refresh');
x=x(good); y=y(good);dtSec=dtSec(good);

bullshitFactor=.01;
xEdges=linspace(min(x),max(x)+bullshitFactor,nBins);
yEdges=linspace(min(y),max(y)+bullshitFactor,nBins);

[~,xBin]=histc(x,xEdges);
[~,yBin]=histc(y,yEdges);


occ=accumarray([yBin,xBin],dtSec);
imagesc(occ);
colorbar;

print(strcat(pathstr,'/occupancy.eps'));
file = 'occupancy.eps';
