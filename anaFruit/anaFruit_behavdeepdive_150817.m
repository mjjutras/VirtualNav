%%

NSDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NSdat';
trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';


%%

seslst = {'JN_BR_15_04_10\JN_BR_15_04_10_13_03';
    'JN_BR_15_04_10\JN_BR_15_04_10_13_11';
    'JN_BR_15_04_13\JN_BR_15_04_13_13_59';
    'JN_BR_15_04_14\JN_BR_15_04_14_13_34';
    'JN_BR_15_04_15\JN_BR_15_04_15_12_32';
    'JN_BR_15_04_16\JN_BR_15_04_16_14_22';
    'JN_BR_15_04_17\JN_BR_15_04_17_13_14';
    'JN_BR_15_04_20\JN_BR_15_04_20_14_43';
    'JN_BR_15_04_21\JN_BR_15_04_21_12_34';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_09';
    'JN_BR_15_04_22\JN_BR_15_04_22_14_54';
    'JN_BR_15_04_23\JN_BR_15_04_23_13_54';
    'JN_BR_15_04_23\JN_BR_15_04_23_14_22';
    'JN_BR_15_04_24\JN_BR_15_04_24_14_04';
    'JN_BR_15_04_28\JN_BR_15_04_28_13_07';
    'JN_BR_15_04_29\JN_BR_15_04_29_13_22';
    'JN_BR_15_04_30\JN_BR_15_04_30_14_42';
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
    'JN_BR_15_05_19\JN_BR_15_05_19_13_00';
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
    'JN_BR_15_06_05\JN_BR_15_06_05_13_35'};


%% define variables for behavioral analysis

numtrl = []; % total number of trials per session
numtrlV = []; % number of VISIBLE banana trials per session
numtrlT = []; % number of TRANSLUCENT banana trials per session
numtrlI = []; % number of INVISIBLE banana trials per session
alpTlevs = []; % translucent alpha range used per session
% two columns, first is lowest alpha range, second is highest

frtbubavg = []; % fruit "bubble" (distance between avatar and fruit to trigger acquisition), average
frtbubste = []; % fruit "bubble", standard error
invbubavg = []; % bubble that triggers invisible banana to become visible, average
invbubste = []; % invisible banana bubble, standard error

numchrall = []; % number of cherries per visible banana trial, average per session
invrepall = []; % prop. invisible bananas presented in same location as previous invisible banana

begposall = []; % begin position (trial start)
alpall = []; % alpha
banfndall = []; % 1 if got banana, 0 if not
prpbanfnd = []; % proportion of invisible bananas found
banposall = []; % banana position (first banana, includes invisible)

bantimindall = []; % latency to get banana (technically, time index)
avginvbantimall = []; % mean latency to get invisible bananas actually acquired
steinvbantimall = []; % standard error of latency to get invisible bananas actually acquired

avgbandispertrl = []; % average distance from banana
cumdispertrl = []; % cumulative distance from banana
avgwaldispertrl = []; % average distance from wall
trllngall = []; % total trial length
invprmindall = []; % invisible, prime (1st in a series of inv. bananas)
nrlsel = []; % trials used to select neural activity (dum)
filind = []; % file index, per trial

% excess path length: median value per session (invisible banana trials)
excpthdif = []; % excess path length, using difference b/t path length and shortest path
excpthnrm = []; % excess path length, using path length divided by shortest path

angpertrlall = []; % cumulative angle (direction) traversed during trial (invisible bananas)

numvisall = []; % average number of visible bananas per invisible banana, per session
numtrnall = []; % average number of translucent bananas per invisible banana, per session

x = linspace(-11.2,11.2,9); % grid squares to estimate coverage of environment
numgrdall = []; % number of grid squares covered by invisible bananas
prpgrdall = []; % number of grid squares divided by number of invisible bananas

% d=0;

% % dum structure contains neural data for analysis
% % initialize dum structure
% dum = [];
% dum.time = cell(1);
% dum.trial = cell(1);
% dum.sampleinfo = [];
% c=1;

%% loop sessions for behavioral analysis

for seslop = 1:length(seslst)
    
    [~,sesnam]=fileparts(seslst{seslop});
    disp(['Processing ' sesnam])

    load(fullfile(NSDir,[sesnam '_NSdat.mat']))
    load(fullfile(trlDir,[sesnam '_trldat.mat']))

    numtrl(seslop,1) = length(trldat.time);
    
    invrepind = zeros(length(trldat.posdat),1); % index of repeated invisible bananas:
    % 1 when banana is presented in the same location as the previous invisible banana
    prvbanalp = cell(1); % previous banana alpha
    pthlng = zeros(length(trldat.posdat),1);
    bantimindses = zeros(length(trldat.posdat),1);
    frtdisses = [];
    invbub = zeros(length(trldat.posdat),1);
    trlalp = nan(length(trldat.posdat),1); % banana alpha per trial
    numvisses = zeros(length(trldat.posdat),1);
    numtrnses = zeros(length(trldat.posdat),1);
    numchrses = zeros(length(trldat.posdat),1); % # cherries
    banposses = zeros(length(trldat.posdat),2);
    banfndses = zeros(length(trldat.posdat),1);
    begdisses = zeros(length(trldat.posdat),1); % distance from banana at trial start
    angpertrl = zeros(length(trldat.posdat),1);
    for trllop = 1:length(trldat.posdat)
        
        trlalp(trllop,1) = trldat.alpha{trllop}(1,2);
        
        angpertrl(trllop,1) = sum(abs(diff(unwrap(trldat.dirdat{trllop}))));
        
        numchrses(trllop,1) = length(find(trldat.frtpos{trllop}(:,2)>0));
        
        % find time when banana eaten, or end of trial if banana not eaten
        if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
            [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
            banfnd = 1;
        else
            % NOTE: if he doesn't find the banana, the actual time here
            % should be the "TIMEOUT", which is either 90 seconds or some
            % other longer amount (before 5/8/15)
            bantimind = length(trldat.time{trllop});
            banfnd = 0;
        end
        
        bantimindses(trllop) = bantimind;
        
        % position of banana
        banpos = trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),3:4);
        
        % calculate distance from banana across trial (until getting banana)
        xdif = trldat.posdat{trllop}(1,1:bantimind)-banpos(1);
        ydif = trldat.posdat{trllop}(2,1:bantimind)-banpos(2);
        dis = sqrt(xdif.^2+ydif.^2);

        % path length to banana (calculate excess path length later)
        pthlng(trllop,1) = sum(sqrt(diff(xdif).^2+diff(ydif).^2));
        
        % figure out distance from wall; use 11.2 as wall boundary
        waldis = nan(size(trldat.posdat{trllop},2),1);
        for timlop = 1:size(trldat.posdat{trllop},2)
            waldis(timlop) = min(abs([-11.2-trldat.posdat{trllop}(:,timlop); 11.2-trldat.posdat{trllop}(:,timlop)]));
        end
        
        % figure out bubble size: distance from fruit target when acquired
        % this will be useful when using shortest distance to target when
        % calculating excess path length
        for frtlop = 1:size(trldat.frttim{trllop},1)
            getfrtpos = trldat.posdat{trllop}(:,trldat.time{trllop}==trldat.frttim{trllop}(frtlop,1)); % position when got fruit
            frtloc = trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==trldat.frttim{trllop}(frtlop,2),1,'last'),3:4); % fruit location
            % takes into account the fact that the invisible banana may move when it becomes visible
            frtdisses = [frtdisses; sqrt((getfrtpos(1)-frtloc(1))^2+(getfrtpos(2)-frtloc(2))^2)];
        end

        % find trial with invisible banana
        if trlalp(trllop)==0
                
            % if this if the first invisible banana of the session
            if length(find(trlalp==0))==1
                invrepind(trllop,1) = 0;
            else % don't include if previous invisible banana was in the same location
                prvinv = find(trlalp(1:trllop-1)==0,1,'last'); % trial with last invisible banana
                if prod(trldat.frtpos{prvinv}(1,3:4)==banpos,2)
                    invrepind(trllop,1) = 1;
                else
                    invrepind(trllop,1) = 0;
                end
            end
            
            % determine alpha of each visible bananas presented in the
            % same location in the preceding trials
            %%% finish this; add contingency so it stops at the last
            %%% invisible banana
            numvis = 0;
            numtrn = 0;
            for revlop=trllop-1:-1:1
                revpos = trldat.frtpos{revlop}(1,3:4);
                if revpos(1)~=banpos(1) || revpos(2)~=banpos(2)
                    break
                elseif trlalp(revlop)==0
                    break
                else
                    if trldat.alpha{revlop}(1,2)==1
                        numvis = numvis+1;
                    elseif trldat.alpha{revlop}(1,2)>0 && trldat.alpha{revlop}(1,2)<1
                        numtrn = numtrn+1;
                    end
                end
            end
            numvisses(trllop) = numvis;
            numtrnses(trllop) = numtrn;
            
            if banfnd
                trnvistim = trldat.alpha{trllop}(find(trldat.alpha{trllop}(:,2)==1,1,'last'),1);
                vispos = trldat.posdat{trllop}(:,trldat.time{trllop}==trnvistim);
                invpos = trldat.frtpos{trllop}(find(trldat.frtpos{trllop}(:,2)==0,1,'first'),3:4);
                invbub(trllop) = sqrt((vispos(1)-invpos(1))^2+(vispos(2)-invpos(2))^2);
            end
            
        else
            invrepind(trllop,1) = 0;
        end
        
        begposall = [begposall; trldat.posdat{trllop}(:,1)'];
        banfndall = [banfndall; banfnd];
        banfndses(trllop,1) = banfnd;
        begdisses(trllop,1) = dis(1);
        
        banposall = [banposall; banpos];
        banposses(trllop,:) = banpos;
        
        bantimindall = [bantimindall; bantimind];
        avgbandispertrl = [avgbandispertrl; mean(dis)];
        cumdispertrl = [cumdispertrl; sum(dis)];
        avgwaldispertrl = [avgwaldispertrl; mean(waldis)];
        trllngall = [trllngall; length(trldat.time{trllop})];
        filind = [filind; seslop];
        alpall = [alpall; trlalp];
        
    end
    
    numtrlV(seslop,1) = length(find(trlalp==1));
    numtrlT(seslop,1) = length(find(trlalp>0 & trlalp<1));
    numtrlI(seslop,1) = length(find(trlalp==0));
    
    if numtrlT(seslop,1)>0
        alpTlevs(seslop,:) = [min(trlalp(trlalp>0 & trlalp<1)) max(trlalp(trlalp>0 & trlalp<1))];
    else
        alpTlevs(seslop,:) = [nan nan];
    end

    % calculate excess path length, taking into account the buffer (bubble)
    % around the banana when getting reward
    excpthdifses = pthlng-(begdisses-mean(frtdisses));
    excpthnrmses = pthlng./(begdisses-mean(frtdisses));

%     % this section used to select trials for neural analysis
%     % select trials for neural analysis (nrlsel)
%     for invlop = 1:length(invprmind)
%         % choose only "invisible banana trials" with at least 2 preceding visible bananas 
%         if length(prvbanalp{invlop})>=2
%             % select only trials with banana search time of 60 seconds or less
%             if bantimindses(invprmind(invlop))<=60000;
%                 % select only trials where there were at least 5 seconds of
%                 % search time for both of the preceding visible bananas
%                 if bantimindses(invprmind(invlop)-2)>=5000 && bantimindses(invprmind(invlop)-1)>=5000
%                     nrlsel = [nrlsel; invprmind(invlop)+d];
%                     dum.time(c:c+1) = data.time(invprmind(invlop)-2:invprmind(invlop)-1);
%                     dum.trial(c:c+1) = data.trial(invprmind(invlop)-2:invprmind(invlop)-1);
%                     dum.sampleinfo(c:c+1,:) = data.sampleinfo(invprmind(invlop)-2:invprmind(invlop)-1,:);
%                     c = c+2;
%                 end
%             end
%         end
%     end
    
%     % add invprmind trials to invprmindall
%     invprmindall = [invprmindall; invprmind+d];
%     d = d+length(trldat.time);
    
    [n,~] = hist3(banposses(trlalp==0,:),'Edges',{x,x});
    numgrdall = [numgrdall; length(find(n(1:8,1:8)))];
    prpgrdall = [prpgrdall; length(find(n(1:8,1:8)))/length(find(trlalp==0))];
    
    prpbanfnd = [prpbanfnd; sum(banfndses(trlalp==0))/length(find(trlalp==0))];

    excpthdif = [excpthdif; median(excpthdifses(trlalp==0))];
    excpthnrm = [excpthnrm; median(excpthnrmses(trlalp==0))];
    
    frtbubavg = [frtbubavg; mean(frtdisses)];
    frtbubste = [frtbubste; std(frtdisses)/sqrt(length(frtdisses))];
    invbubavg = [invbubavg; mean(invbub(invbub>0))];
    invbubste = [invbubste; std(invbub(invbub>0))/sqrt(length(invbub>0))];
    
    numvisall = [numvisall; mean(numvisses(trlalp==0))];
    numtrnall = [numtrnall; mean(numtrnses(trlalp==0))];

    numchrall = [numchrall; mean(numchrses(trlalp~=0))];
    
    invrepall = [invrepall; length(find(invrepind(trlalp==0)))/length(find(trlalp==0))];
    
    avginvbantimall = [avginvbantimall; mean(bantimindses(logical((trlalp==0).*banfndses)))];
    steinvbantimall = [steinvbantimall; std(bantimindses(logical((trlalp==0).*banfndses)))/sqrt(length(find((trlalp==0).*banfndses)))];
    
    angpertrlall = [angpertrlall; mean(angpertrl(trlalp==0))];
    
end

% dum.fsample = data.fsample;
% dum.label = data.label;
% 
% clear data trldat

