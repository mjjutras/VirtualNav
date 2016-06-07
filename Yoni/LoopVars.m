%% Add the necessary data folders to path

addpath('/Users/yonibrowning/Documents/BananaRecall/MatFiles');%Data folder;
addpath('/Users/yonibrowning/Documents/BananaRecall/FeatureAnalysis'); %Scripts
cd('/Users/yonibrowning/Documents/BananaRecall/MatFiles'); %Need to be in data folder

%% Create storage matricies for each of the variables we want to look at.
% Transparency of banana
alpha = []; %-)clear
% Time to banana
tmeToBanana = [];%-)
% Time in area. Eventual matrix must be 9xnumtrls, since there are 9 areas
tmeInEachArea = [];%-)
% Area that the banana is in
bananaArea = [];%-);
% Time in correct Area;
tmeInArea = [];%-);
% Average distance from wall
wallDist = [];
% Average distnace from corner
cornerDist = [];
% Average distance from center
centerDist = [];
% Average distance from banana
bananaDist = [];
% Time near wall. 
%%%%Time in center is rest of time. 
tmeNearWall = [];
% time near banana
tmeNearBanana = [];
% time near corner
tmeNearCorner = [];
% mean Velocity
meanVel = [];
% time going forward
tmeForward = [];
% mean angular velocity
meanAngleV = [];
% time turning 
tmeTurning = [];
% Excess Path Length
excessPath = [];
% "Meander"--> amount turned per distance traveled
meander = [];

%%%% time stopped  = totalTime-timeTurning-timeForward
%% Loop through trials, loading one file at a time
c  = 0;
list = dir;
for l = 3:length(dir);
    load(list(l).name);
    numtrls = length(trldat.time);
    for trl = 1:numtrls
        % count the current trial
        c = c+1;
        
        % Calculate alpha for this trial
        try
            thisAlpha = trldat.alpha{trl}(end-1,2);
        catch
            thisAlpha = trldat.alpha{trl}(end,2);
        end
        alpha = [alpha thisAlpha];
        
        % Calculate time to banana
        tme0 = trldat.time{trl}(1);
        try
            tme1 = trldat.frttim{trl}(1);
        catch
            tme1 = trldat.time{trl}(end);
        end
        thisTmeToBanana = tme1-tme0;
        tmeToBanana = [tmeToBanana thisTmeToBanana];
        
        % find banana area
        banx = trldat.frtpos{trl}(1,3);
        bany = trldat.frtpos{trl}(1,4);
        bananaArea = [bananaArea findarea(banx,bany)];
        
        % find time in each area
        x=trldat.posdat{trl}(1,1:thisTmeToBanana);
        y=trldat.posdat{trl}(2,1:thisTmeToBanana);
        areain = findarea(x,y);
        thisTmeInEachArea = zeros(9,1);
        for a = 1:9
            thisTmeInEachArea(a) = nnz(areain == a);
        end
        tmeInEachArea = [tmeInEachArea thisTmeInEachArea];
        
        % find time in correct area
        tmeInArea =[tmeInArea thisTmeInEachArea(bananaArea(c))];
        
        % find average distance from corner
        A = sqrt([(-10-x).^2+(-10-y).^2;...
            (-10-x).^2+(10-y).^2;...
            (10-x).^2+(-10-y).^2;...
            (10-x).^2+(-10-y).^2]);
        distFromCorner = min(A);
        cornerDist = [cornerDist mean(distFromCorner)];
        CornerThresh = 5;
        tmeNearCorner = [tmeNearCorner nnz(distFromCorner<=CornerThresh)];

        % find average distance from wall
        A = abs([(10-x);(-10-x);(10-y);(-10-y)]);
        distFromWall = min(A);
        wallDist = [wallDist mean(distFromWall)];
        
        % find average time near wall
        wallThresh = 3;
        tmeNearWall = [tmeNearWall nnz(distFromWall<=wallThresh)];
        
        % find average distance from center
        centerDist = [centerDist mean(sqrt(x.^2+y.^2))];

        % find average distance from banana, time spent near banana
        distFromBanana = sqrt((banx-x).^2+(bany-y).^2);
        BananaThresh = 6;
        tmeNearBanana = [tmeNearBanana nnz(distFromBanana<=BananaThresh)];
        bananaDist = [bananaDist mean(distFromBanana)];
        
        % Find Mean Angular Velocity, time turning
        angleV  = diff(trldat.dirdat{trl});
        angleV(angleV>=pi) = angleV(angleV>=pi)-2*pi;
        meanAngleV = [meanAngleV mean(angleV(angleV~=0))];
        tmeTurning = [tmeTurning nnz(angleV~=0)];
        
        % Find Mean Velocity, time going forward
        velocity = sqrt(diff(x).^2+diff(y).^2);
        meanVel = [meanVel mean(velocity)];
        tmeForward = [tmeForward nnz(velocity~=0)];
        
        % Excess Path length
        excessPath = [excessPath sum(velocity)-(distFromBanana(1)-distFromBanana(end))];
        
        % Meander
        meander = [meander sum(abs(angleV))/sum(velocity)];
        
        % store path so that we can look at it later
        paths.pos{c} = trldat.posdat{trl}(:,1:thisTmeToBanana);
        paths.dir{c} = trldat.dirdat{trl}(1:thisTmeToBanana);
    end
end

%%
M = [tmeInArea'./tmeToBanana',wallDist',cornerDist',...
    bananaDist',tmeNearWall'./tmeToBanana',tmeNearCorner'./tmeToBanana',meanVel',...
    tmeForward'./tmeToBanana',meanAngleV',tmeTurning'./tmeToBanana',excessPath',meander'];
M(isinf(M)) = 0;
M(isnan(M)) = 0;
M = M(alpha==0,:);
[coeff,score,L] = princomp(M)

%%
numClus = 3
idx = kmeans([excessPath',tmeNearWall'./tmeToBanana',bananaDist'],numClus)
for i = 1:numClus
    plot(excessPath(idx==i),tmeNearWall(idx==i)./tmeToBanana(idx==i),'.');
    hold on;
end

%%
H = hist3([excessPath',tmeNearWall'./tmeToBanana'],{linspace(1,max(excessPath),100),linspace(1,max(tmeNearWall'./tmeToBanana'),100)});
pcolor(log(H));hold on;shading interp
plot(excessPath,tmeNearWall/tmeToBanana,'.');