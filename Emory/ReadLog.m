session = 1106;
Logfile = ['C:\GizVideo500\session_' num2str(session) '\log.txt'];

fid = fopen(Logfile);

%Pull log file in as one big string structure
line_count = 1;
l = 1;
while l>0
   l = fgets(fid);
   if l>0
       logtxt{line_count,1} = l;
       line_count = line_count+1;
   end
end

% EYE DATA
eyedata.timeEarly = nan(length(logtxt), 1);
eyedata.timeLate = nan(length(logtxt), 1);
eyedata.X = nan(length(logtxt), 1);
eyedata.Y = nan(length(logtxt), 1);
for i = 1:length(logtxt)
   if ~isempty(strfind(logtxt{i}, 'EyeData'))
       textscan(logtxt{i}, '%d64 %d %*s %f64 %f64')
       eyedata.timeEarly(i) = ans{1};
       eyedata.timeLate(i) = ans{2};
       eyedata.X(i) = ans{3};
       eyedata.Y(i) = ans{4};
   end  
end

eyedata.Xreal = eyedata.X(find(~isnan(eyedata.X)));
eyedata.Yreal = eyedata.Y(find(~isnan(eyedata.Y)));
eyedata.timeProlly = eyedata.timeEarly(find(~isnan(eyedata.timeEarly))) + eyedata.timeLate(find(~isnan(eyedata.timeLate)));



% AVATAR POSITION
avatarPos.timeEarly = nan(length(logtxt), 1);
avatarPos.timeLate = nan(length(logtxt), 1);
avatarPos.X = nan(length(logtxt), 1);
avatarPos.Y = nan(length(logtxt), 1);
avatarPos.Z = nan(length(logtxt), 1);
for i = 1:length(logtxt)
    if ~isempty(strfind(logtxt{i}, 'VROBJECT_POS	PandaEPL_avatar	LPoint3f('))
        textscan(logtxt{i}, '%d64 %d %*s %*s LPoint3f(%f64 %f64 %f64', 1, 'delimiter', {',', '\t'})
        avatarPos.timeEarly(i) = ans{1};
        avatarPos.timeLate(i) = ans{2};
        avatarPos.X(i) = ans{3};
        avatarPos.Y(i) = ans{4};
        avatarPos.Z(i) = ans{5};
    end
end


avatarPos.Xreal = avatarPos.X(find(~isnan(avatarPos.X)));
avatarPos.Yreal = avatarPos.Y(find(~isnan(avatarPos.Y)));
avatarPos.timeProlly = avatarPos.timeEarly(find(~isnan(avatarPos.timeEarly))); %is this right? 

% BANANA POSITIONS

% BANANA HITS
bananaHits.time = nan(length(logtxt), 1);
bananaHits.banana = cell(length(logtxt), 1);
for i = 1:length(logtxt)
    if (~isempty(strfind(logtxt{i}, 'VROBJECT_STASHED'))) && (~isempty(strfind(logtxt{i}, 'True')))
        textscan(logtxt{i}, '%d64 %d %*s %s %*s')
        bananaHits.time(i) = ans{1};
        bananaHits.banana{i} = ans{3}{1};
    end
end

for i = 1:length(bananaHits.banana)
    if ~isempty(bananaHits.banana{i})
        if bananaHits.time(i) == 1365004146698
            disp(i)
        end
    end
end


% TRIAL INFO
trials.timestamps = nan(length(logtxt), 1);
trials.trialNum = nan(length(logtxt), 1);
for i = 1:length(logtxt)
     if ~isempty(strfind(logtxt{i}, 'NewTrial'))
         textscan(logtxt{i}, '%d64 %d %*s %d')
         trials.timestamps(i) = ans{1};
         trials.trialNum(i) = ans{3};
     end
end

trials.timeReal = trials.timestamps(find(~isnan(trials.timestamps)));
trials.trialNumReal = trials.trialNum(find(~isnan(trials.trialNum)));


%scatter(avatarPos.Xreal(1508:4020), avatarPos.Yreal(1508:4020))
%plot(avatarPos.Xreal(1508:4020), avatarPos.Yreal(1508:4020))
numFrames = 113380;
M1(1:numFrames) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);

figure;
ylim([-5.5 5.5]);
xlim([-5.5 5.5]);
start = 51353;
hold on
for i = 1:numFrames
    a = plot(avatarPos.X(start), avatarPos.Y(start));  
    set(a, 'Color', 'red', 'LineWidth', 2);
    if (~isnan(avatarPos.X(start)))
        M1(i) = getframe;
    elseif i == 1
        M1(i) = getframe;
    else
        M1(i) = M1(i-1);
    end
    start=start+1;
    hold on
    i
end

outputVidPath = 'c:\gizvideo500\pathvid4pt2.avi';
pathvid = VideoWriter(outputVidPath);
%pathvid.FrameRate = numFrames/59.85;
open(pathvid);
writeVideo(pathvid, M1(20000:40000));
close(pathvid);
    

%%
    
    clear bananaHits
clear eyedata
clear logtxt
clear trials

numFrames = 113380;
M1(1:numFrames) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);

figure;
ylim([-5.5 5.5]);
xlim([-5.5 5.5]);
start = 51353;
hold on
startInd=1:(numFrames/5):numFrames;
for k=1:length(startInd);
    M1(1:numFrames) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);

   for i = startInd(k):(startInd(k)+((numFrames/5)-1))
       a = plot(avatarPos.X(start:(start + 1)), avatarPos.Y(start:(start + 1)));
       set(a, 'Color', 'red', 'LineWidth', 2);
       if ~isnan(avatarPos.X(start:(start+1)))
           M1(i) = getframe;
       else
           M1(i) = M1(i-1);
       end
       start=start+1;
       hold on
   end
   outputVidPath = ['c:\gizvideo500\pathvid5_part',num2str(k),'.avi'];
   pathvid = VideoWriter(outputVidPath);
   %pathvid.FrameRate = numFrames/59.85;
   open(pathvid);
   writeVideo(pathvid, M1);
   close(pathvid);
end


numFrames = 113380;

figure;
ylim([-5.5 5.5]);
xlim([-5.5 5.5]);
start = 51353;
hold on
startInd=1:(numFrames/5):numFrames;
for k=1:length(startss);   
   
   M1(1:(numFrames/5)) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);
   
   for i = startInd(k):(startInd(k)+((numFrames/5)-1))
       a = plot(avatarPos.X(start), avatarPos.Y(start));
       set(a, 'Color', 'red', 'LineWidth', 2);
       if ~isnan(avatarPos.X(start)
           M1(i) = getframe;
       else
           M1(i) = M1(i-1);
       end
       start=start+1;
       hold on
   end
   
   outputVidPath = ['c:\gizvideo500\pathvid5_part',num2str(k),'.avi'];
   pathvid = VideoWriter(outputVidPath);
   %pathvid.FrameRate = numFrames/59.85;
   open(pathvid);
   writeVideo(pathvid, M1);
   close(pathvid);
   clear M1
end

%%
numFrames = 113380;
%M1(1:numFrames/10) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);
pt = 0;

figure;
ylim([-5.5 5.5]);
xlim([-5.5 5.5]);
start = 51353;
hold on
for i = 1:numFrames
    if mod(i, numFrames/10) == 1
        M1(1:numFrames/10) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);
        pt = pt + 1;
    end
    a = plot(avatarPos.X(start), avatarPos.Y(start));  
    set(a, 'Color', 'red', 'LineWidth', 2);
    if (~isnan(avatarPos.X(start)))
        M1(i) = getframe;
    elseif i == 1
        M1(i) = getframe;
    else
        M1(i) = M1(i-1);
    end
    start=start+1;
    hold on
    i
    
    if mod(i, numFrames/10) == 0
        pathvid = ['c:\gizvideo500\pathvidFinal_' num2str(pt) 'of10.avi'];
        pathvid = VideoWriter(outputVidPath);
        open(pathvid);
        writeVideo(pathvid, M1);
        close(pathvid);
    end
        
end

%% This Worked
%numChunk = 5;
%numChunk = 7;
numChunk=17;
%numFrames = 113380;
%numFrames = 59913;
numFrames = ceil(59913 / 25);
chunk=numFrames/numChunk;

% figure;
% ylim([-5.5 5.5]);
% xlim([-5.5 5.5]);
% %startOrig = 51353;
 startOrig = 1;
% hold on
startInd=1:chunk:numFrames;

%Setup Figure
hf = figure;
set(hf, 'position', [5 50 1500 1400])
set(gca, 'Units', 'pixels')
%set(gca, 'Position', [20 10 1280 1000])
set(gca, 'Position', [20 10 1200 960])

% Plot Eye Position if true, plot Bubble/Searchlight if false
eyePos=false;

% Size of Marker
markPos=500;
markBub=5000;

for k=1:length(startInd);   
   k
   %M1(1:chunk) = struct('cdata', zeros(343, 435, 3, 'uint8'), 'colormap', []);
   M1(1:chunk) = struct('cdata', zeros(1024, 1280, 3, 'uint8'), 'colormap', []);
   
   start=(startInd(k)-1)+startOrig;
   
   for i = 1:chunk
       %a = plot(avatarPos.X(start), avatarPos.Y(start));
       %a = plot(avatarPos.ourSection(start, 2), avatarPos.ourSection(start, 3));
       if eyePos
           markStr='eyePos';
           sizeStr=num2str(markPos);
           a = scatter(eyedata.ourSection(start:start+24, 2), eyedata.ourSection(start:start+24, 3), markPos, 'or', 'MarkerFaceColor', 'red');
       else
           markStr='bubble';
           sizeStr=num2str(markBub);
           a = scatter(eyedata.ourSection(start:start+24, 2), eyedata.ourSection(start:start+24, 3), markBub, 'w', 'filled', 'marker', 'o');
       end       
       set(gca, 'Color', 'k')
       %xlim([-640 640])
       %ylim([-512 512])
       xlim([-700 700]) %effectively adjusted the gain.
       ylim([-600 600]) %effectively adjusted the gain.
       
       %set(a, 'Color', 'red', 'LineWidth', 2);
       %if ~isnan(avatarPos.X(start))
       %if ~isnan(avatarPos.ourSection(start, 2))
       %if (~isnan(eyedata.ourSection(start, 2))) && (~isnan(eyedata.ourSection(start+1, 2))) && (~isnan(eyedata.ourSection(start+2, 2))) && (~isnan(eyedata.ourSection(start+3, 2))) && (~isnan(eyedata.ourSection(start+4, 2)))
           M1(i) = getframe;
       %elseif i == 1
       %     M1(i) = getframe;
       %else
       %    M1(i) = M1(i-1);
       %end
       start=start+25;
       %hold on
       i
   end
   
   outputVidPath = ['c:\gizvideo500\pathvid_',markStr,sizeStr,'_part',num2str(k),'.avi'];
   pathvid = VideoWriter(outputVidPath);
   %pathvid.FrameRate = numFrames/59.85;
   pathvid.FrameRate = 40;
   open(pathvid);
   writeVideo(pathvid, M1);
   close(pathvid);
   clear M1
end

%%
avatarPos.ourSection = nan(59913, 1);
avatarPos.ourSection(:, 1) = 1365004146698:1365004206610;
avatarPos.ourSection(:, 2) = nan(59913, 1);
avatarPos.ourSection(:, 3) = nan(59913, 1);

for i = 1:length(avatarPos.timeProlly)
    find(avatarPos.ourSection(:, 1) == avatarPos.timeProlly(i))
    if ~isempty(ans)
        avatarPos.ourSection(ans, 2) = avatarPos.Xreal(i);
        avatarPos.ourSection(ans, 3) = avatarPos.Yreal(i);
    end    
end
    
plot(avatarPos.ourSection(:, 2), avatarPos.ourSection(:, 3))


%%
eyedata.ourSection = nan(59913, 1);
eyedata.ourSection(:, 1) = 1365004146698:1365004206610;
eyedata.ourSection(:, 2) = nan(59913, 1);
eyedata.ourSection(:, 3) = nan(59913, 1);

for i = 1:length(eyedata.timeProlly)
    find(eyedata.ourSection(:, 1) == eyedata.timeProlly(i))
    if ~isempty(ans)
        eyedata.ourSection(ans, 2) = eyedata.Xreal(i);
        eyedata.ourSection(ans, 3) = eyedata.Yreal(i);
    end    
end


