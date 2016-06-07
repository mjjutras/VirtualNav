session = 1107;
dirDat= 'R:\Buffalo Lab\eblab\Virtual Navigation\Free Foraging\Giz\';
Logfile = [dirDat 'session_' num2str(session) '\log.txt'];

% Logfile = ['S:\Kiril\GizVideo500\session_' num2str(session) '\log.txt'];

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

%Clear unnecessary variables
fclose(fid); clear fid line_count Logfile session l ans;

% TRIAL INFO
%trials.timestamps = nan(length(logtxt), 1);
%trials.trialNum = nan(length(logtxt), 1);
trials.startPoints = [];

for i = 1:length(logtxt)
     if ~isempty(strfind(logtxt{i}, 'NewTrial'))
         textscan(logtxt{i}, '%d64 %d64 %*s %d')
         trials.startPoints(length(trials.startPoints) + 1, 1) = ans{1};
         %trials.timestamps(i) = ans{1};
         %trials.trialNum(i) = ans{3};
     % TRIAL ONE IS ACTUALLY DIFFERENT BECAUSE OF THE INSTRUCTION DISPLAY!!!
     elseif ~isempty(strfind(logtxt{i}, 'INSTRUCT_OFF'))
         textscan(logtxt{i}, '%d64', 1)
         trials.startPoints(1) = ans{1};
     end
end

endTime = textscan(logtxt{i}, '%d64', 1);
endTime = endTime{1};

eyedata(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
eyedata(:, 2) = nan(length(eyedata), 1);
eyedata(:, 3) = nan(length(eyedata), 1);
eyedata(:, 1) = (trials.startPoints(1):endTime);

avatar.Position(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
avatar.Position(:, 2) = nan(length(avatar.Position), 1);
avatar.Position(:, 3) = nan(length(avatar.Position), 1);
avatar.Position(:, 1) = (trials.startPoints(1):endTime);

avatar.Header(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
avatar.Header(:, 2) = nan(length(avatar.Header), 1);
avatar.Header(:, 3) = nan(length(avatar.Header), 1);
avatar.Header(:, 1) = (trials.startPoints(1):endTime);

avatar.TurningSpeed(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
avatar.TurningSpeed(:, 2) = nan(length(avatar.TurningSpeed), 1);
avatar.TurningSpeed(:, 1) = (trials.startPoints(1):endTime);

avatar.LinearSpeed(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
avatar.LinearSpeed(:, 2) = nan(length(avatar.LinearSpeed), 1);
avatar.LinearSpeed(:, 1) = (trials.startPoints(1):endTime);

avatar.TurningAccel(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
avatar.TurningAccel(:, 2) = nan(length(avatar.TurningAccel), 1);
avatar.TurningAccel(:, 1) = (trials.startPoints(1):endTime);

avatar.LinearAccel(:, 1) = nan(endTime - trials.startPoints(1) + 1, 1);
avatar.LinearAccel(:, 2) = nan(length(avatar.LinearAccel), 1);
avatar.LinearAccel(:, 1) = (trials.startPoints(1):endTime);

%bananas.StashedHits = cell(endTime - trials.startPoints(1) + 1, 2);
%for i = 1:length(bananas.StashedHits)
%    bananas.StashedHits{i, 1} = trials.startPoints(1) + i - 1;
%    disp(i)
%end
%bananas.StashedHits = cell(length(trials.startPoints)*10, 2);
bananas.StashedHits = {}
for i = 1:length(logtxt)
   if ~isempty(strfind(logtxt{i}, 'EyeData'))
       textscan(logtxt{i}, '%d64 %d64 %*s %f64 %f64');
       timeIndex = find(eyedata(:, 1) == ans{1} + ans{2});
       eyedata(timeIndex, 2) = ans{3};
       eyedata(timeIndex, 3) = ans{4};
   elseif ~isempty(strfind(logtxt{i}, 'VROBJECT_POS	PandaEPL_avatar	LPoint3f('))
       textscan(logtxt{i}, '%d64 %d64 %*s %*s LPoint3f(%f64 %f64 %f64', 1, 'delimiter', {',', '\t'});
       timeIndex = find(avatar.Position(:, 1) == ans{1}); %Is this right?
       avatar.Position(timeIndex, 2) = ans{3};
       avatar.Position(timeIndex, 3) = ans{4};
   elseif ~isempty(strfind(logtxt{i}, 'VROBJECT_HEADING	PandaEPL_avatar'))
       textscan(logtxt{i}, '%d64 %d64 %*s %*s  %f64');
       timeIndex = find(avatar.Header(:, 1) == ans{1}); %Is this right?
       avatar.Header(timeIndex, 2) = ans{3};
       avatar.Header(timeIndex, 3) = mod(avatar.Header(timeIndex, 2), 360);
   elseif ~isempty(strfind(logtxt{i}, 'MOVINGOBJECT_TURNINGSPEED	PandaEPL_avatar'))
       textscan(logtxt{i}, '%d64 %d64 %*s %*s  %f64');
       timeIndex = find(avatar.TurningSpeed(:, 1) == ans{1});
       avatar.TurningSpeed(timeIndex, 2) = ans{3};
   elseif ~isempty(strfind(logtxt{i}, 'MOVINGOBJECT_TURNINGACCEL	PandaEPL_avatar'))
       textscan(logtxt{i}, '%d64 %d64 %*s %*s %f64');
       timeIndex = find(avatar.TurningAccel(:, 1) == ans{1});
       avatar.TurningAccel(timeIndex, 2) = ans{3};
   elseif ~isempty(strfind(logtxt{i}, 'MOVINGOBJECT_LINEARACCEL	PandaEPL_avatar'))
       textscan(logtxt{i}, '%d64 %d64 %*s %*s %f64');
       timeIndex = find(avatar.LinearAccel(:, 1) == ans{1});
       avatar.LinearAccel(timeIndex, 2) = ans{3};
   elseif ~isempty(strfind(logtxt{i}, 'MOVINGOBJECT_LINEARSPEED	PandaEPL_avatar'))
       textscan(logtxt{i}, '%d64 %d64 %*s %*s %f64');
       timeIndex = find(avatar.LinearSpeed(:, 1) == ans{1});
       avatar.LinearSpeed(timeIndex, 2) = ans{3};
   elseif (~isempty(strfind(logtxt{i}, 'VROBJECT_STASHED'))) && (~isempty(strfind(logtxt{i}, 'True'))) && (~isempty(strfind(logtxt{i}, 'banana')))
       textscan(logtxt{i}, '%d64 %d64 %*s %s %*s');
       %timeIndex = find(bananas.StashedHits(:, 1) == ans{1});
       %bananas.StashedHits{timeIndex, 2} = ans{3};
       bananas.StashedHits{size(bananas.StashedHits, 1) + 1, 1} = ans{1};
       bananas.StashedHits{size(bananas.StashedHits, 1), 2} = ans{3};
   end  
   disp([num2str(i/length(logtxt)*100) '% Done'])
end

%Banana Captures
%Banana Position Info
%%
%trials.timeReal = trials.timestamps(find(~isnan(trials.timestamps)));
%trials.trialNumReal = trials.trialNum(find(~isnan(trials.trialNum)));

% EYE DATA
eyedata.timeEarly = nan(length(logtxt), 1);
eyedata.timeLate = nan(length(logtxt), 1);
%eyedata.X = nan(length(logtxt), 1);
%eyedata.Y = nan(length(logtxt), 1);

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

scatter(eyedata.ourSection(:, 2), eyedata.ourSection(:, 3), 1, 'filled');
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
%  hold on
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
   
   start=(25*startInd(k)-1)+startOrig;
   
   for i = 1:chunk
       %a = plot(avatarPos.X(start), avatarPos.Y(start));
       %a = plot(avatarPos.ourSection(start, 2), avatarPos.ourSection(start, 3));
       %a = plot(avatarPos.ourSection(start:start+24, 2), avatarPos.ourSection(start:start+24, 3));
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
   
   outputVidPath = ['c:\gizvideo500\pathvidTestEYE_',markStr,sizeStr,'_part',num2str(k),'.avi'];
   pathvid = VideoWriter(outputVidPath);
   %pathvid.FrameRate = numFrames/59.85;
   pathvid.FrameRate = 40;
   open(pathvid);
   writeVideo(pathvid, M1);
   close(pathvid);
   clear M1
end



