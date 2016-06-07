clear all;close all;clc
load parNavDat.mat
%% import data
trl = 3;

figure
pos  = data.posdat{trl}';
vel = data.veldat{trl};
dir = data.dirdat{trl};
accel = data.acceldat{trl};
tvel = data.tveldat{trl};
taccel = data.acceldat{trl};

dd = diff(dir);
dd(abs(dd)>1) = 0;
X= [vel(1:length(vel)-1)' abs(dd)'];
%X = [vel' tvel'];


%% Cluster data using k-means
nk = 3;
ind = kmeans(X,nk,'Replicates',10);

plot(pos(1,1),pos(1,2),'go','MarkerSize',50);hold on;



color = ['r','g','b','m','k','y','c'];
for i = 1:nk
    plot(pos((ind==i),1),pos((ind==i),2),strcat('*',color(i)));hold on;
end

axis([-10 10 -10 10])

plot(data.banpos{trl}(:,2),data.banpos{trl}(:,3),strcat('*',color(i+1)),'MarkerSize',20)

%% plot the optimal path
    points = [pos(1,:);data.banpos{trl}(:,2),data.banpos{trl}(:,3)];
    distance = zeros(10);
    for i = 1:11;
        for j = 1:11;
            distance(i,j) = sqrt((points(i,1)-points(j,1))^2+(points(i,2)-points(j,2))^2);
        end
    end

    P = perms(1:10);
    P = P+1;
    minu = inf;
    optipath = 0;
    for i = 1:length(P(:,1));
        path = distance(1,P(i,1));
        for k = 1:9
            path = path+distance(P(i,k),P(i,k+1));
        end

        minu = min(minu,path);
        if minu==path;
            optipath = i;
        end
    end

    plot(points([1 P(optipath,:)],1),points([1 P(optipath,:)],2),'c','LineWidth',3)
    drawnow

%%
figure
for i = 1:nk
    plot(vel((ind==i)),abs(dd(ind==i)),strcat('*',color(i))); hold on;
end
xlabel('linear velocity')
ylabel('angular velocity')

%% Monkey position to all bananas
figure
start  = pos(1,:);
bends = data.baneat{trl}(:,2);
d=sqrt((start(1)-data.banpos{trl}(:,2)).^2+(start(2)-data.banpos{trl}(:,3)).^2);

bar(d(bends+1));
