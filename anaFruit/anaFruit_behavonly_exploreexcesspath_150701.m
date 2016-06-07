%% questions

% look at the path distance as a function of target (banana) location when
% the bananas are visible:
% can we normalize the performance metric (excess path length) by this
% distance?
% should we still take into account the time spent in proximity to the
% invisible target somehow, since excess path length does not address this?


%% info

% trial always starts with banana, then two cherries
%
% on invisible banana trials, the two cherries appear in trldat.frtpos, but
% they are not presented in the trial; the next trial starts immediately
% when the monkey gets the invisible banana


%%

trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';

% specify the sessions to use for now
seslst = {
    'JN_BR_15_05_01\JN_BR_15_05_01_13_12';
    'JN_BR_15_05_04\JN_BR_15_05_04_13_57';
    'JN_BR_15_05_04\JN_BR_15_05_04_14_51';
    'JN_BR_15_05_05\JN_BR_15_05_05_13_02';
    'JN_BR_15_05_06\JN_BR_15_05_06_13_30';
    'JN_BR_15_05_07\JN_BR_15_05_07_12_59';
    'JN_BR_15_05_08\JN_BR_15_05_08_14_17';
    'JN_BR_15_05_12\JN_BR_15_05_12_12_49';
    'JN_BR_15_05_13\JN_BR_15_05_13_13_14';
    'JN_BR_15_05_14\JN_BR_15_05_14_13_34';
    'JN_BR_15_05_18\JN_BR_15_05_18_13_03';
    'JN_BR_15_05_19\JN_BR_15_05_19_13_00'; % 150624: added 15_05_01 through 15_05_19
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

banposall = []; % banana position
begposall = []; % begin position (trial start)
alpall = []; % alpha
pthlngall = []; % absolute path length
pthratall = []; % path length divided by distance from starting point (shortest path)
banfndall = []; % 1 if got banana, 0 if not
bantimall = []; % time to get banana (or time out)
begdisall = []; % distance from start position to banana
filind = []; % file index, per trial
for seslop = 1:length(seslst)
    
    [~,sesnam]=fileparts(seslst{seslop});
    disp(['Processing ' sesnam])

    load(fullfile(trlDir,[sesnam '_trldat.mat']))

    for trllop = 1:length(trldat.posdat)
        
        % find time when banana eaten, or end of trial if banana not eaten
        bantimind = find(trldat.frttim{trllop}(:,2)==0,1,'first');
        if ~isempty(bantimind)
            [~,bantim] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(bantimind,1)));
            banfnd = 1;
        else
            bantim = length(trldat.time{trllop});
            banfnd = 0;
        end
        
        banpos = trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),3:4);
        
        % calculate distance from banana
        xdif = trldat.posdat{trllop}(1,1:bantim)-banpos(1);
        ydif = trldat.posdat{trllop}(2,1:bantim)-banpos(2);
        dis = sqrt(xdif.^2+ydif.^2);

        % calculate path length (from start position)            
        pthlng = sum(abs(diff(dis)));
        pthrat = sum(abs(diff(dis)))/dis(1);
        
        pthlngall = [pthlngall; pthlng];
        pthratall = [pthratall; pthrat];
        banposall = [banposall; banpos];
        begposall = [begposall; trldat.posdat{trllop}(:,1)'];
        alpall = [alpall; trldat.alpha{trllop}(1,2)];
        banfndall = [banfndall; banfnd];
        bantimall = [bantimall; bantim];
        begdisall = [begdisall; dis(1)];
        filind = [filind; seslop];

    end

end
            
%% heat maps

x = linspace(-11.2,11.2,81);
ctrs = x(1:end-1)+mean(gradient(x))/2;

banposvis = banposall(logical(alpall),:);
pthlngvis = pthlngall(logical(alpall));
pthratvis = pthratall(logical(alpall));
bantimvis = bantimall(logical(alpall));
begdisvis = begdisall(logical(alpall));

lng3d = []; % path length (distance travelled)
rat3d = []; % path length divided by distance from starting point (shortest path)
tim3d = []; % time to get banana
beg3d = []; % distance from starting position to banana
for k=1:size(banposvis,1)
    
    posn = imfilter(hist3(banposvis(k,:),'Edges',{x,x}),fspecial('disk',10));
    posn(posn~=max(max(posn))) = nan; posn(~isnan(posn))=1;
    
    lng3d(k,:,:) = posn(1:end-1,1:end-1)*pthlngvis(k);
    rat3d(k,:,:) = posn(1:end-1,1:end-1)*pthratvis(k);
    tim3d(k,:,:) = posn(1:end-1,1:end-1)*bantimvis(k);
    beg3d(k,:,:) = posn(1:end-1,1:end-1)*begdisall(k);
    
end

lng2d = squeeze(nanmedian(lng3d,1));
rat2d = squeeze(nanmedian(rat3d,1));
tim2d = squeeze(nanmedian(tim3d,1));
beg2d = squeeze(nanmedian(beg3d,1));

figure;imagesc(ctrs,ctrs,lng2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('path length (distance travelled)');colorbar
figure;imagesc(ctrs,ctrs,rat2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('path length divided by distance from starting point (shortest path)');colorbar
figure;imagesc(ctrs,ctrs,tim2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('time to get banana');colorbar
figure;imagesc(ctrs,ctrs,beg2d);xlim([-11.2 11.2]);ylim([-11.2 11.2]);colormap jet;title('distance from starting position to banana');colorbar

save('R:\Buffalo Lab\Mike\VirtualNav\MAT files\workspace\rat2d_visban_24ses.mat','rat2d')
