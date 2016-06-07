%% questions

% does performance vary based on starting position (distance) relative to the banana?
% does performance vary based on the position of the banana, relative to walls/center?


%%

trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';

% specify the sessions to use for now
seslst = {
    'JN_BR_15_05_20\JN_BR_15_05_20_13_04';
    'JN_BR_15_05_22\JN_BR_15_05_22_13_20';
    'JN_BR_15_05_26\JN_BR_15_05_26_13_42';
    'JN_BR_15_05_27\JN_BR_15_05_27_12_53';
    'JN_BR_15_05_28\JN_BR_15_05_28_13_59';
    'JN_BR_15_05_28\JN_BR_15_05_28_14_16';
    'JN_BR_15_05_29\JN_BR_15_05_29_14_19';
    'JN_BR_15_06_01\JN_BR_15_06_01_14_35';
    'JN_BR_15_06_02\JN_BR_15_06_02_13_15';
    'JN_BR_15_06_03\JN_BR_15_06_03_14_32';
    'JN_BR_15_06_04\JN_BR_15_06_04_13_04';
    'JN_BR_15_06_05\JN_BR_15_06_05_13_35';
%     'JN_BR_15_06_08\JN_BR_15_06_08_13_14'; % exclude this day because behavior bad all week
    };

begposall = []; % begin position (trial start)
alpall = []; % alpha
banfndall = []; % 1 if got banana, 0 if not
banposall = []; % banana position (first banana, includes invisible)
bantimindall = []; % latency to get banana (technically, time index)
avgbandispertrl = []; % average distance from banana
cumdispertrl = []; % cumulative distance from banana
avgwaldispertrl = []; % average distance from wall
trllngall = []; % total trial length
invprmindall = []; % invisible, prime (1st in a series of inv. bananas)
invprmnrlsel = []; % trials used to select neural activity (dum)
filind = []; % file index, per trial
numtrl = []; % number of trials per file
excpthdif = []; % excess path length, using difference b/t path length and shortest path
excpthnrm = []; % excess path length, using path length divided by shortest path
excpthind = []; % excess path length, index
c=1;
d=0;
for seslop = 1:length(seslst)
    
    [~,sesnam]=fileparts(seslst{seslop});
    disp(['Processing ' sesnam])

    load(fullfile(trlDir,[sesnam '_trldat.mat']))

    numtrl(seslop,1) = length(trldat.time);
    
    % calculate three metrics for each trial:
    % latency to get banana (equal to time-out if banana not acquired)
    % average distance from banana
    % cumulative distance from banana
    
    invprmind = []; % invisible "prime" (first invisible banana in a sequence) index
    prvbanalp = cell(1); % previous banana alpha
    excpthdifdum = zeros(length(trldat.posdat),1);
    excpthnrmdum = zeros(length(trldat.posdat),1);
    excpthinddum = zeros(length(trldat.posdat),1);
    for trllop = 1:length(trldat.posdat)
        
        % find time when banana eaten, or end of trial if banana not eaten
        if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
            [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
            banfnd = 1;
        else
            bantimind = length(trldat.time{trllop});
            banfnd = 0;
        end
        
        banpos = trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),3:4);
        
        % calculate distance from banana
        xdif = trldat.posdat{trllop}(1,1:bantimind)-banpos(1);
        ydif = trldat.posdat{trllop}(2,1:bantimind)-banpos(2);
        dis = sqrt(xdif.^2+ydif.^2);
        
        % figure out distance from wall; use 11.2 as wall boundary
        waldis = nan(size(trldat.posdat{trllop},2),1);
        for timlop = 1:size(trldat.posdat{trllop},2)
            waldis(timlop) = min(abs([-11.2-trldat.posdat{trllop}(:,timlop); 11.2-trldat.posdat{trllop}(:,timlop)]));
        end
        
        begposall = [begposall; trldat.posdat{trllop}(:,1)'];
        alpall = [alpall; trldat.alpha{trllop}(1,2)];
        banfndall = [banfndall; banfnd];
        banposall = [banposall; banpos];
        bantimindall = [bantimindall; bantimind];
        avgbandispertrl = [avgbandispertrl; mean(dis)];
        cumdispertrl = [cumdispertrl; sum(dis)];
        avgwaldispertrl = [avgwaldispertrl; mean(waldis)];
        trllngall = [trllngall; length(trldat.time{trllop})];
        filind = [filind; seslop];
        
        % find trial with invisible banana
        if trldat.alpha{trllop}(1,2)==0
            % don't include if the banana was invisible on previous trial 
            if trllop>=2 && trldat.alpha{trllop-1}(1,2)~=0
                
                % don't include if previous invisible banana was in the same location
                if ~isempty(invprmind)
                    if ~prod(trldat.frtpos{invprmind(end)}(find(trldat.frtpos{invprmind(end)}(:,2)==0,1,'first'),3:4)==banpos,2)
                        invprmind = [invprmind; trllop];
                    end
                else % if this if the first invisible banana of the session
                    invprmind = [invprmind; trllop];
                end
                
                % determine alpha of each visible bananas presented in the
                % same location in the preceding trials
                if invprmind(end)==trllop
                    prvbanalp{length(invprmind),1} = [];
                    for revlop=trllop-1:-1:1
                        revpos = trldat.frtpos{revlop}(find(trldat.frtpos{revlop}(:,2)==0,1,'first'),3:4);
                        if revpos(1)~=banpos(1) || revpos(2)~=banpos(2)
                            break
                        else
                            prvbanalp{length(invprmind),1} = [prvbanalp{length(invprmind),1}; revlop+d trldat.alpha{revlop}(1,2)];
                        end
                    end
                    prvbanalp{length(invprmind),1} = flipud(prvbanalp{length(invprmind),1});                        
                end

%                 % Yoni's method: need to calculate velocity directly from
%                 % position (doesn't work for veldat from Panda log)
%                 excpthdifdum(trllop,1) = sum(conv(trldat.veldat{trllop},normpdf(-1000:1000,0,200),'same'))-dis(1);
%                 excpthnrmdum(trllop,1) = sum(conv(trldat.veldat{trllop},normpdf(-1000:1000,0,200),'same'))/dis(1);
                % this way works
                excpthdifdum(trllop,1) = sum(abs(diff(dis)))-dis(1);
                excpthnrmdum(trllop,1) = sum(abs(diff(dis)))/dis(1);
                excpthinddum(trllop,1) = (sum(abs(diff(dis)))-dis(1))/(sum(abs(diff(dis)))+dis(1));
                
            end
        end
        
    end

    for invlop = 1:length(prvbanalp)
        if length(prvbanalp{invlop})>=2
            invprmnrlsel = [invprmnrlsel; invprmind(invlop)+d];
        end
    end
    
    % add invprmind trials to invprmindall
    invprmindall = [invprmindall; invprmind+d];
    d = d+length(trldat.time);
    
    excpthdif = [excpthdif; excpthdifdum];
    excpthnrm = [excpthnrm; excpthnrmdum];
    excpthind = [excpthind; excpthinddum];
    
end


%%

% banposall = []; % banana position (first banana, includes invisible)
% [begposall(invprmnrlsel) banfndall(invprmnrlsel) bantimindall(invprmnrlsel) avgbandispertrl(invprmnrlsel) avgwaldispertrl(invprmnrlsel)]

x = linspace(-10,10,17);
% figure
% hist3(banposall(invprmnrlsel,:),{x,x});
[n,c] = hist3(banposall(invprmnrlsel,:),'Edges',{x,x});
figure;imagesc(c{1},c{2},n);xlim([-10 10]);ylim([-10 10]);colormap hot;colorbar
hold on;scatter(banposall(invprmnrlsel,2),banposall(invprmnrlsel,1),'bo')

% calculate distance from banana at start position
xdif = begposall(invprmnrlsel,1)-banposall(invprmnrlsel,1);
ydif = begposall(invprmnrlsel,2)-banposall(invprmnrlsel,2);
begdis = sqrt(xdif.^2+ydif.^2);


%% categorize trials into good/bad memory
% use distance from corners of "Mean distance from banana/wall" plot
% normalize these variables first

bantimnrm = bantimindall(invprmnrlsel)/max(bantimindall(invprmnrlsel));
bandstnrm = avgbandispertrl(invprmnrlsel)/max(avgbandispertrl(invprmnrlsel));
waldstnrm = avgwaldispertrl(invprmnrlsel)/max(avgwaldispertrl(invprmnrlsel));

% time to acquire banana, mean dist. from banana, mean dist. from wall
mempnthi = [0 0 max(waldstnrm)];
mempntlo = [max(bantimnrm) max(bandstnrm) 0];
memdst = sqrt(sum((abs(([bantimnrm bandstnrm waldstnrm])-repmat(mempnthi,length(invprmnrlsel),1))).^2,2));
pltttl = 'time to acquire banana, mean dist. from banana, mean dist. from wall';

[~,memsrt] = sort(memdst);

% % memindhi = memsrt(1:round(length(memsrt)/2));
% % memindlo = memsrt(round(length(memsrt)/2)+1:end);
% memindhi = memsrt(1:30);
% memindlo = memsrt(end-29:end);
% 
% bantim_himem = bantimindall(invprmnrlsel(memindhi));
% bantim_lomem = bantimindall(invprmnrlsel(memindlo));
% 
% figure;hold on
% scatter3([bantimindall(invprmnrlsel(memindhi)); bantimindall(invprmnrlsel(memindlo))], ...
%     [avgbandispertrl(invprmnrlsel(memindhi)); avgbandispertrl(invprmnrlsel(memindlo))], ...
%     [avgwaldispertrl(invprmnrlsel(memindhi)); avgwaldispertrl(invprmnrlsel(memindlo))], ...
%     ones(length(memindhi)+length(memindlo),1)*10,[repmat([0 1 0],length(memindhi),1); repmat([1 0 0],length(memindlo),1)])
% xlabel('Time to acquire banana')
% ylabel('Mean distance from banana')
% zlabel('Mean distance from wall')
% title(pltttl)


%% 2d histogram (large bins), memory performance

x = linspace(-11.2,11.2,33);
ctrs = x(1:end-1)+mean(gradient(x))/2;

mem3d = [];
for k=1:length(memdst)
    
    n = hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*memdst(k);
    n(n==0) = nan;
    mem3d(k,:,:) = n(1:end-1,1:end-1);
    
end

mem2d = squeeze(nanmean(mem3d,1));

figure;imagesc(ctrs,ctrs,mem2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet



%% heat map, memory performance

x = linspace(-11.2,11.2,81);
ctrs = x(1:end-1)+mean(gradient(x))/2;

lat3d = []; % latency to get banana
dis3d = []; % average distance from banana
wal3d = []; % average distance from wall
mem3d = []; % memory performance
beg3d = []; % starting distance
epldis3d = []; % excess path length, difference
eplrat3d = []; % excess path length, ratio
memnrm3d = []; % memory normalized by starting distance
for k=1:length(memdst)
    
    latn = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*bantimindall(invprmnrlsel(k)),fspecial('gaussian',30,3));
    disn = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*avgbandispertrl(invprmnrlsel(k)),fspecial('gaussian',30,3));
    waln = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*avgwaldispertrl(invprmnrlsel(k)),fspecial('gaussian',30,3));
    memn = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*memdst(k),fspecial('gaussian',30,3));
    begn = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*begdis(k),fspecial('gaussian',30,3));
    epldis = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*excpthdif(invprmnrlsel(k)),fspecial('gaussian',30,3));
    eplrat = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*excpthnrm(invprmnrlsel(k)),fspecial('gaussian',30,3));
    eplind = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*excpthind(invprmnrlsel(k)),fspecial('gaussian',30,3));
    memnrmn = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x})*(memdst(k)/begdis(k)),fspecial('gaussian',30,3));
    
    nd = imfilter(hist3(banposall(invprmnrlsel(k),:),'Edges',{x,x}),fspecial('disk',15));
    latn(nd==0) = nan;
    disn(nd==0) = nan;
    waln(nd==0) = nan;
    memn(nd==0) = nan;
    begn(nd==0) = nan;
    epldis(nd==0) = nan;
    eplrat(nd==0) = nan;
    eplind(nd==0) = nan;
    memnrmn(nd==0) = nan;
    
    % multiply by scaling factor to correct for smoothing (to plot correct amplitudes)
    lat3d(k,:,:) = latn(1:end-1,1:end-1)*58.141424919673;
    dis3d(k,:,:) = disn(1:end-1,1:end-1)*58.141424919673;
    wal3d(k,:,:) = waln(1:end-1,1:end-1)*58.141424919673;
    mem3d(k,:,:) = memn(1:end-1,1:end-1)*58.141424919673;
    beg3d(k,:,:) = begn(1:end-1,1:end-1)*58.141424919673;
    epldis3d(k,:,:) = epldis(1:end-1,1:end-1)*58.141424919673;
    eplrat3d(k,:,:) = eplrat(1:end-1,1:end-1)*58.141424919673;
    eplind3d(k,:,:) = eplind(1:end-1,1:end-1)*58.141424919673;
    memnrm3d(k,:,:) = memnrmn(1:end-1,1:end-1)*58.141424919673;
    
%     figure
%     subplot(2,1,1)
%     imagesc(ctrs,ctrs,n(1:end-1,1:end-1));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet
%     subplot(2,1,2)
%     imagesc(ctrs,ctrs,squeeze(nanmean(mem3d,1)));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet
% %     set(gca,'clim',[0 0.002])
%     set(gcf,'Position',[80 67 496 911])
%     pause
%     close
    
end

lat2d = squeeze(nanmean(lat3d,1));
dis2d = squeeze(nanmean(dis3d,1));
wal2d = squeeze(nanmean(wal3d,1));
mem2d = squeeze(nanmean(mem3d,1));
beg2d = squeeze(nanmean(beg3d,1));
epldis2d = squeeze(nanmean(epldis3d,1));
eplrat2d = squeeze(nanmean(eplrat3d,1));
eplind2d = squeeze(nanmean(eplind3d,1));
memnrm2d = squeeze(nanmean(memnrm3d,1));

figure;imagesc(ctrs,ctrs,lat2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Latency to get banana by banana location');colorbar
figure;imagesc(ctrs,ctrs,dis2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Average distance from banana by banana location');colorbar
figure;imagesc(ctrs,ctrs,wal2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Average distance from wall by banana location');colorbar
figure;imagesc(ctrs,ctrs,mem2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Memory performance by banana location');colorbar
figure;imagesc(ctrs,ctrs,beg2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Starting distance by banana location');colorbar
figure;imagesc(ctrs,ctrs,epldis2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (difference) by banana location');colorbar
figure;imagesc(ctrs,ctrs,eplrat2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (ratio) by banana location');colorbar
figure;imagesc(ctrs,ctrs,eplind2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Excess path length (index) by banana location');colorbar
figure;imagesc(ctrs,ctrs,memnrm2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('Normalized memory by banana location');colorbar

% figure;imagesc(ctrs,ctrs,imfilter(mem2d,fspecial('gaussian',3,1)));xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet

