%%

if strcmp(license,'375399')
    NSDir = 'C:\Data\VR';
    trlDir = 'C:\Data\VR';
elseif strcmp(license,'613743')
    NSDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NSdat';
    trlDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trldat';
end


%%

trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';

% specify the session to use
ses = 'JN_BR_15_06_05\JN_BR_15_06_05_13_35';

[~,sesnam]=fileparts(ses);
disp(['Processing ' sesnam])

% get the date from the sesnam
c = textscan(sesnam,'%s','Delimiter','_');
yer = str2num(['20' c{1}{3}]);
mon = str2num(c{1}{4});
day = str2num(c{1}{5});

% load(fullfile(NSDir,[sesnam '_NSdat.mat']))
load(fullfile(trlDir,[sesnam '_trldat.mat']))

numtrl = length(trldat.time);


%% make a 2-d histogram of location over time, for each trial

grdedg = -11.2:2.8:11.2;

grdall = [];
for trllop = 1:length(trldat.posdat)
    [nx,binx]=histc(trldat.posdat{trllop}(1,:),grdedg);
    [ny,biny]=histc(trldat.posdat{trllop}(2,:),grdedg);
    grd = zeros(length(grdedg)-1,length(grdedg)-1);
    for l=1:length(binx)
        grd(binx(l),biny(l))=grd(binx(l),biny(l))+1;
    end
    grdall(trllop,:,:) = grd;
end


%%

% calculate three metrics for each trial:
% latency to get banana (equal to time-out if banana not acquired)
% average distance from banana
% cumulative distance from banana

begposall = []; % begin position (trial start)
alpall = []; % alpha
banfndall = []; % 1 if got banana, 0 if not
banposall = []; % banana position (first banana, includes invisible)
bantimindall = []; % latency to get banana (technically, time index)
avgbandispertrl = []; % average distance from banana
cumdispertrl = []; % cumulative distance from banana
avgwaldispertrl = []; % average distance from wall
trllngall = []; % total trial length
% invprmindall = []; % invisible, prime (1st in a series of inv. bananas)
nrlsel = []; % trials used to select neural activity (dum)
% filind = []; % file index, per trial
% numtrl = []; % number of trials per file
invprmind = []; % invisible "prime" (first invisible banana in a sequence) index
% prvbanalp = cell(1); % previous banana alpha
excpthdif = zeros(length(trldat.posdat),1); % excess path length, using difference b/t path length and shortest path
excpthnrm = zeros(length(trldat.posdat),1); % excess path length, using path length divided by shortest path
excpthind = zeros(length(trldat.posdat),1); % excess path length, index
% bantimindses = zeros(length(trldat.posdat),1);
d=0;
for trllop = 1:length(trldat.posdat)
    
    % find time when banana eaten, or end of trial if banana not eaten
    if ~isempty(find(trldat.frttim{trllop}(:,2)==0,1,'first'))
        [~,bantimind] = min(abs(trldat.time{trllop}-trldat.frttim{trllop}(find(trldat.frttim{trllop}(:,2)==0,1,'first'),1)));
        banfnd = 1;
    else
        % NOTE: if he doesn't find the banana, the actual time here
        % should be the "TIMEOUT", which is either 90 seconds or some
        % other longer amount (before 5/8/15)
%         bantimind = length(trldat.time{trllop});
        if datenum(yer,mon,day) < 736092
            % don't know the timeout before 5/8/15
        else
            bantimind = 90000;
        end
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
    bantimindses(trllop) = bantimind;
    avgbandispertrl = [avgbandispertrl; mean(dis)];
    cumdispertrl = [cumdispertrl; sum(dis)];
    avgwaldispertrl = [avgwaldispertrl; mean(waldis)];
    trllngall = [trllngall; length(trldat.time{trllop})];
    
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

            % calculate excess path length
            excpthdif(trllop,1) = sum(abs(diff(dis)))-dis(1);
            excpthnrm(trllop,1) = sum(abs(diff(dis)))/dis(1);
            excpthind(trllop,1) = (sum(abs(diff(dis)))-dis(1))/(sum(abs(diff(dis)))+dis(1));

        end
    end
    
    figure
    subplot(3,1,1)
    plot(trldat.time{trllop}(1:bantimind),dis)
    ylabel('Distance from banana')
    title(['Trial ' num2str(trllop) '; Alpha = ' num2str(trldat.alpha{trllop}(1,2))])
    subplot(3,1,2)
    hold on
    plot(trldat.posdat{trllop}(1,:),trldat.posdat{trllop}(2,:))
    scatter(trldat.frtpos{trllop}(trldat.frtpos{trllop}(:,2)==0,3),trldat.frtpos{trllop}(trldat.frtpos{trllop}(:,2)==0,4),'g') % banana(s) in green
    scatter(trldat.frtpos{trllop}(ismember(trldat.frtpos{trllop}(:,2),[1 2]),3),trldat.frtpos{trllop}(ismember(trldat.frtpos{trllop}(:,2),[1 2]),4),'r') % cherries in red
    scatter(trldat.posdat{trllop}(1,1),trldat.posdat{trllop}(2,1),'b') % start pos in blue
    xlim([-11.2 11.2]);ylim([-11.2 11.2])
    subplot(3,1,3)
    imagesc(grdedg(1:end-1)+mean(diff(grdedg))/2,grdedg(1:end-1)+mean(diff(grdedg))/2,squeeze(grdall(trllop,:,:))');
    axis xy; colormap hot
    title(['total travel time = ' num2str(size(trldat.posdat{trllop},2)/1000) ' sec'])
    set(gcf,'Position',[1000 51 343 937])
    
end

% select trials for neural analysis (nrlsel)
for invlop = 1:length(invprmind)
    % choose only "invisible banana trials" with at least 2 preceding visible bananas
    if length(prvbanalp{invlop})>=2
        % select only trials with banana search time of 60 seconds or less
        if bantimindses(invprmind(invlop))<=60000;
            % select only trials where there were at least 5 seconds of
            % search time for both of the preceding visible bananas
            if bantimindses(invprmind(invlop)-2)>=5000 && bantimindses(invprmind(invlop)-1)>=5000
                nrlsel = [nrlsel; invprmind(invlop)+d];
            end
        end
    end
end



