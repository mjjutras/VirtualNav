% This code calculates pEpisode for a select few VR sessions and uses three
% different cycle lengths, giving pEpisode measures for each cycle length
% Also contains code to parse virtual movement based on forward momentum
% and turning, marking "bouts" of each, and calculates pEpisode within
% bouts of each behavioral condition 

% datdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated\';
datdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\NS6decimated\';
pepdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\pEpisode\pEpisode_mult_cycles\';

samplerate = 1000;
freqs = (2^(1/8)).^(8:42);
width=7;
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

seslst = [...
    'JN140825010';
    'JN150209001';
    'JN150410003';
    ];

for seslop = 1:size(seslst,1)
    
    BRnam = seslst(seslop,:);

    load(fullfile(datdir,[BRnam '_NS6_SF30.mat']))

    Pm_all = cell(1,size(NS6.Data,1));
    pv_all = cell(1,size(NS6.Data,1));
    unionvector_all = cell(1,size(NS6.Data,1));

    ft_progress('init', 'etf',     'Please wait...');
    for chnlop = 1:size(NS6.Data,1)

        ft_progress(chnlop/size(NS6.Data,1), 'Processing event %d from %d', chnlop, size(NS6.Data,1));

        t_total = clock;

        % get energy vector for each list
    %     fprintf('Calc. Power...');
    %     t0 = clock;

        B=multienergyvec(double(NS6.Data(chnlop,:)),freqs,samplerate,width);
    %     fprintf('%g\n',etime(clock,t0));

        % calc the mean fit
        Blog = log10(double(B));
        Pm = mean(Blog,2);
        clear Blog

        % get the fit
    %     fprintf('Calc. fit...');
    %     t0 = clock;
        % IMPORTANS: chi_squarefit assumes that frequencies are
        % logarithmically spaced!
        [all,R2] = chi_squarefit(freqs,Pm);
        all = all';
    %     fprintf('%g\n',etime(clock,t0));

        % test fit of regression line
        pv=polyfit(log10(freqs),Pm',1); % linear regression
    %     figure;plot(freqs,Pm)
    %     hold on;plot(freqs,pv(2)+pv(1)*log10(freqs),'r')

        % set the threshold
        thresh = all(:,951);

        % loop through each frequency and save the unionvector
        UVdum1 = int8(zeros(size(freqs,2),size(B,2)));
        UVdum2 = int8(zeros(size(freqs,2),size(B,2)));
        UVdum3 = int8(zeros(size(freqs,2),size(B,2)));
        UVdum4 = int8(zeros(size(freqs,2),size(B,2)));
        UVdum5 = int8(zeros(size(freqs,2),size(B,2)));

        for f = 1:size(freqs,2)
            % get the pepisode
            UVdum1(f,:) = int8(episodeid(B(f,:),thresh(f),1*samplerate/freqs(f),shoulder,0)); % 1-cycle threshold
            UVdum2(f,:) = int8(episodeid(B(f,:),thresh(f),2*samplerate/freqs(f),shoulder,0)); % 2-cycle threshold
            UVdum3(f,:) = int8(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0)); % 3-cycle threshold
            UVdum4(f,:) = int8(episodeid(B(f,:),thresh(f),4*samplerate/freqs(f),shoulder,0)); % 4-cycle threshold
            UVdum5(f,:) = int8(episodeid(B(f,:),thresh(f),5*samplerate/freqs(f),shoulder,0)); % 5-cycle threshold
        end
        clear B

        unionvector = UVdum1 + UVdum2 + UVdum3 + UVdum4 + UVdum5;
        clear UVdum*

        % convert unionvector from single to int8 to conserve memory
        unionvector=int8(unionvector);

        Pm_all{chnlop}=Pm;
        pv_all{chnlop}=pv;
        unionvector_all{chnlop}=unionvector;

    end
    ft_progress('close')

    save(fullfile(pepdir,[BRnam '_pEpisode.mat']),'unionvector_all','Pm_all','pv_all','freqs','-v7.3')

end


%% plot average pepisode across session, for each array

sesnam = 'JN140825010';
% sesnam = 'JN150209001';
% sesnam = 'JN150410003';

pepdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\pEpisode\pEpisode_mult_cycles\';
load(fullfile(pepdir,[sesnam '_pEpisode.mat']))

% average pepisode across session
clear UV*
for chnlop = 1:length(unionvector_all)
    UV1(chnlop,:) = mean(unionvector_all{chnlop}>=1,2);
    UV2(chnlop,:) = mean(unionvector_all{chnlop}>=2,2);
    UV3(chnlop,:) = mean(unionvector_all{chnlop}>=3,2);
    UV4(chnlop,:) = mean(unionvector_all{chnlop}>=4,2);
    UV5(chnlop,:) = mean(unionvector_all{chnlop}>=5,2);
end

% add offset for plotting
ofs = 0.1;
ofv = (ofs:ofs:12*ofs)';

clear UVofs*
UVofs1{1} = bsxfun(@plus,UV1(1:12,:),ofv);
UVofs1{2} = bsxfun(@plus,UV1(13:24,:),ofv);
UVofs1{3} = bsxfun(@plus,UV1(25:36,:),ofv);

UVofs2{1} = bsxfun(@plus,UV2(1:12,:),ofv);
UVofs2{2} = bsxfun(@plus,UV2(13:24,:),ofv);
UVofs2{3} = bsxfun(@plus,UV2(25:36,:),ofv);

UVofs3{1} = bsxfun(@plus,UV3(1:12,:),ofv);
UVofs3{2} = bsxfun(@plus,UV3(13:24,:),ofv);
UVofs3{3} = bsxfun(@plus,UV3(25:36,:),ofv);

UVofs4{1} = bsxfun(@plus,UV4(1:12,:),ofv);
UVofs4{2} = bsxfun(@plus,UV4(13:24,:),ofv);
UVofs4{3} = bsxfun(@plus,UV4(25:36,:),ofv);

UVofs5{1} = bsxfun(@plus,UV5(1:12,:),ofv);
UVofs5{2} = bsxfun(@plus,UV5(13:24,:),ofv);
UVofs5{3} = bsxfun(@plus,UV5(25:36,:),ofv);

figure
subplot(1,3,1); plot(freqs,UVofs1{1})
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Pepisode')
subplot(1,3,2); plot(freqs,UVofs1{2})
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); plot(freqs,UVofs1{3})
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', 1 cycle'],'t')

figure
subplot(1,3,1); plot(freqs,UVofs3{1})
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Pepisode')
subplot(1,3,2); plot(freqs,UVofs3{2})
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); plot(freqs,UVofs3{3})
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', 3 cycles'],'t')

figure
subplot(1,3,1); plot(freqs,UVofs5{1})
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Pepisode')
subplot(1,3,2); plot(freqs,UVofs5{2})
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); plot(freqs,UVofs5{3})
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', 5 cycles'],'t')

% plot power for comparison
clear Pm
for k=1:length(Pm_all)
    Pm(k,:) = Pm_all{k};
end

clear Pmofs
Pmofs{1} = bsxfun(@plus,Pm(1:12,:),ofv);
Pmofs{2} = bsxfun(@plus,Pm(13:24,:),ofv);
Pmofs{3} = bsxfun(@plus,Pm(25:36,:),ofv);

figure
subplot(1,3,1); plot(freqs,Pmofs{1})
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Power')
subplot(1,3,2); plot(freqs,Pmofs{2})
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); plot(freqs,Pmofs{3})
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', Power'],'t')


%% examine relationship between theta and movement
% moving forward vs. turning vs. standing still

% load trial data
trlDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\trldat';
load(fullfile(trlDir,'JN_14_08_25_13_12_trldat.mat'))
NSDir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\NSdat';
load(fullfile(NSDir,'JN_14_08_25_13_12_NSdat.mat'))

% define movement patterns
velcat = [];
for trllop = 1:length(trldat.time)
    velcat = [velcat trldat.veldat{trllop}];
end
velcat = velcat(~isnan(velcat));


%% use threshold method to detect and mark turning bouts

% each direction sample is taken only when the direction changes, for a
% maximum resolution of 60 Hz; if direction doesn't change after 40
% ms, fill in the last sample taken 25 ms after the previous sample in
% order to get a smooth sampling of direction values across time for each
% trial
for trllop = 1:length(trldat.dirdat)
    if isnan(trldat.dirdat{trllop}(1))
        trldat.dirdat{trllop}(1) = trldat.dirdat{trllop}(find(~isnan(trldat.dirdat{trllop}),1,'first'));
    end
    fnddirval = [find(~isnan(trldat.dirdat{trllop})) length(trldat.dirdat{trllop})];
    while ~isempty(find(diff(fnddirval)>=40))
        for vallop = 2:length(fnddirval)
            if fnddirval(vallop)-fnddirval(vallop-1)>=40
                trldat.dirdat{trllop}(fnddirval(vallop-1)+25) = trldat.dirdat{trllop}(fnddirval(vallop-1));
            end
        end
        fnddirval = [find(~isnan(trldat.dirdat{trllop})) length(trldat.dirdat{trllop})];
    end
    if isnan(trldat.dirdat{trllop}(end))
        trldat.dirdat{trllop}(end) = trldat.dirdat{trllop}(find(~isnan(trldat.dirdat{trllop}),1,'last'));
    end
end

% calculate angular difference and interpolate across timepoints
% this includes smoothing which shouldn't greatly affect the timecourse of
% the angular difference but will make it easier to identify start/end of
% turning bouts
trldat.angdif = cell(size(trldat.dirdat));
for trllop = 1:length(trldat.dirdat)
    trldat.angdif{trllop} = nan(size(trldat.dirdat{trllop}));
    trldat.angdif{trllop}(~isnan(trldat.dirdat{trllop})) = [0 smooth(abs(angdiff(unwrap(trldat.dirdat{trllop}(~isnan(trldat.dirdat{trllop}))))))'];
    trldat.angdif{trllop} = inpaint_nans(trldat.angdif{trllop},2);
%     for ii = find(~isnan(trldat.angdif{trllop}))
%         I = trldat.angdif{trllop}(ii);
%         if ~isempty(find(~isnan(trldat.angdif{trllop}(ii+1:end)),1,'first'))
%             subind = ii+1:find(~isnan(trldat.angdif{trllop}(ii+1:end)),1,'first')+ii-1;
%         else
%             subind = ii+1:length(trldat.angdif{trllop});
%         end
%         trldat.angdif{trllop}(subind) = repmat(I,1,length(subind));
%     end
end

% threshold for detecting "turning bout"
anglim = 1e-3;

for trllop = 1:length(trldat.angdif)
    
    ang_below_thresh = trldat.angdif{trllop} <= anglim;

    angdifdif = [diff(trldat.angdif{trllop}) 0]; % create differential of filtered velocity signal

    % differential change marker: finds "reversals", when turning velocity
    % goes from negative to positive (speeds up from slowing down)
    difchgmrk = zeros(size(angdifdif));
    for k = 2:length(angdifdif)-1
        if angdifdif(k) == 0
            difchgmrk(k) = 1;
        elseif angdifdif(k-1) < 0 && angdifdif(k) > 0
            difchgmrk(k) = 1;
        end
    end

    % detect turn begins and turn ends, based on when ang. velocity crosses
    % threshold set above
    trnbegdum = find(diff(trldat.angdif{trllop} > anglim) > 0)+1;
    trnenddum = find(diff(trldat.angdif{trllop} > anglim) < 0)+1;
    if trldat.angdif{trllop}(1) > anglim
        trnenddum = trnenddum(2:end);
    end
    if trldat.angdif{trllop}(end) > anglim
        trnbegdum = trnbegdum(1:end-1);
    end

    trnbeg = nan(1,length(trnbegdum));
    trnend = nan(1,length(trnenddum));
    for k = 1:length(trnbegdum)

        % diff index: time periods 500 ms before "dummy" saccade start or after
        % "dummy" saccade end
        difind1 = trnbegdum(k)-500:trnbegdum(k);
        difind2 = trnenddum(k):trnenddum(k)+500;

        % make sure time periods don't include before start or after end of
        % trial
        difind1 = difind1(difind1>0 & difind1<=length(difchgmrk));
        difind2 = difind2(difind2>0 & difind2<=length(difchgmrk));

        % make sure only includes velocity "reversals" that are not part of
        % a turn (so only at beginning/end of turn)
        ang_below_thresh_dum1 = ang_below_thresh(difind1);
        ang_below_thresh_dum2 = ang_below_thresh(difind2);

        % find velocity "reversals" occurring within each time period
        difchgtmp1 = find(difchgmrk(difind1)) + difind1(1) - 1;
        difchgtmp2 = find(difchgmrk(difind2)) + difind2(1) - 1;

        difchgtmp1 = intersect(difchgtmp1, find(ang_below_thresh_dum1) + difind1(1) - 1);
        difchgtmp2 = intersect(difchgtmp2, find(ang_below_thresh_dum2) + difind2(1) - 1);

        if ~isempty(difchgtmp1) && ~isempty(difchgtmp2)
            trnbeg(k) = difchgtmp1(ft_nearest(difchgtmp1, trnbegdum(k)));
            trnend(k) = difchgtmp2(ft_nearest(difchgtmp2, trnenddum(k)));
        end

    end
    
    trnbeg = trnbeg(~isnan(trnbeg));
    trnend = trnend(~isnan(trnend));
    
    % create trndir, 1 if left turn, 2 if right
    trndir = nan(1,length(trnbeg));
    dirnew = unwrap(trldat.dirdat{trllop});
    for ii = find(~isnan(dirnew))
        I = dirnew(ii);
        if ~isempty(find(~isnan(dirnew(ii+1:end)),1,'first'))
            subind = ii+1:find(~isnan(dirnew(ii+1:end)),1,'first')+ii-1;
        else
            subind = ii+1:length(dirnew);
        end
        dirnew(subind) = repmat(I,1,length(subind));
    end
    for trnlop = 1:length(trnbeg)
        if dirnew(trnbeg(trnlop)) > dirnew(trnend(trnlop))
            trndir(trnlop) = 2;
        elseif dirnew(trnbeg(trnlop)) < dirnew(trnend(trnlop))
            trndir(trnlop) = 1;
        end
    end
    
    % fix weird glitch where turn bouts occur back to back (only when
    % adjacent turns are in the same direction)
    [~,ia,ib] = intersect(trnbeg,trnend);
    trnrem = trndir(ia)==trndir(ib);
    ia = ia(trnrem);
    ib = ib(trnrem);
    trnbeg = trnbeg(setxor(1:length(trnbeg),ia));
    trnend = trnend(setxor(1:length(trnend),ib));
    trndir = trndir(setxor(1:length(trndir),ia));
    
    trldat.trnmrk{trllop} = [trnbeg; trnend; trndir];
    
end

% % plot to test
% figure;hold on;
% subplot(3,1,1);scatter(trldat.time{trllop},trldat.veldat{trllop},'.')
% for k=1:size(trldat.trnmrk{trllop},2);line([trldat.time{trllop}(trldat.trnmrk{trllop}(1,k)) trldat.time{trllop}(trldat.trnmrk{trllop}(1,k))],ylim,'Color','g');line([trldat.time{trllop}(trldat.trnmrk{trllop}(2,k)) trldat.time{trllop}(trldat.trnmrk{trllop}(2,k))],ylim,'Color','r');end
% subplot(3,1,2);plot(trldat.time{trllop},trldat.angdif{trllop})
% for k=1:size(trldat.trnmrk{trllop},2);line([trldat.time{trllop}(trldat.trnmrk{trllop}(1,k)) trldat.time{trllop}(trldat.trnmrk{trllop}(1,k))],ylim,'Color','g');line([trldat.time{trllop}(trldat.trnmrk{trllop}(2,k)) trldat.time{trllop}(trldat.trnmrk{trllop}(2,k))],ylim,'Color','r');end
% subplot(3,1,3);scatter(trldat.time{trllop},unwrap(trldat.dirdat{trllop}),'.')
% for k=1:size(trldat.trnmrk{trllop},2);line([trldat.time{trllop}(trldat.trnmrk{trllop}(1,k)) trldat.time{trllop}(trldat.trnmrk{trllop}(1,k))],ylim,'Color','g');line([trldat.time{trllop}(trldat.trnmrk{trllop}(2,k)) trldat.time{trllop}(trldat.trnmrk{trllop}(2,k))],ylim,'Color','r');end


%% use threshold method to detect and mark bouts of forward momentum (trldat.veldat)

% interpolate velocity data with last non-nan value at each timepoint
for trllop = 1:length(trldat.veldat)
    if isnan(trldat.veldat{trllop}(1))
        trldat.veldat{trllop}(1) = 0;
    end
    for ii = find(~isnan(trldat.veldat{trllop}))
        I = trldat.veldat{trllop}(ii);
        if ~isempty(find(~isnan(trldat.veldat{trllop}(ii+1:end)),1,'first'))
            subind = ii+1:find(~isnan(trldat.veldat{trllop}(ii+1:end)),1,'first')+ii-1;
        else
            subind = ii+1:length(trldat.veldat{trllop});
        end
        trldat.veldat{trllop}(subind) = repmat(I,1,length(subind));
    end
end

% threshold for detecting "movement bout"
movlim = 1;

for trllop = 1:length(trldat.veldat)
    
    vel_below_thresh = trldat.veldat{trllop} <= movlim;

    veldif = [diff(trldat.veldat{trllop}) 0]; % create differential of filtered velocity signal

    % differential change marker: finds "reversals", when velocity goes
    % down to zero
    difchgmrk = trldat.veldat{trllop}==0;

    %detect movement starts and ends, based on when velocity crosses
    %threshold set above
    movbegdum = find(diff(trldat.veldat{trllop} > movlim) > 0)+1;
    movenddum = find(diff(trldat.veldat{trllop} > movlim) < 0)+1;
    if trldat.veldat{trllop}(1) > movlim
        movenddum = movenddum(2:end);
    end
    if trldat.veldat{trllop}(end) > movlim
        movbegdum = movbegdum(1:end-1);
    end

    movbeg = nan(1,length(movbegdum));
    movend = nan(1,length(movenddum));
    for k = 1:length(movbegdum)

        % diff index: time periods 800 ms before "dummy" saccade start or after
        % "dummy" saccade end
        difind1 = movbegdum(k)-800:movbegdum(k);
        difind2 = movenddum(k):movenddum(k)+800;

        % make sure time periods don't include before start or after end of
        % trial
        difind1 = difind1(difind1>0 & difind1<=length(difchgmrk));
        difind2 = difind2(difind2>0 & difind2<=length(difchgmrk));

        % make sure only includes velocity "reversals" that are not part of a
        % saccade (so only at beginning/end of saccade)
        vel_below_thresh_dum1 = vel_below_thresh(difind1);
        vel_below_thresh_dum2 = vel_below_thresh(difind2);

        % find velocity "reversals" occurring within each time period
        difchgtmp1 = find(difchgmrk(difind1)) + difind1(1) - 1;
        difchgtmp2 = find(difchgmrk(difind2)) + difind2(1) - 1;

        difchgtmp1 = intersect(difchgtmp1, find(vel_below_thresh_dum1) + difind1(1) - 1);
        difchgtmp2 = intersect(difchgtmp2, find(vel_below_thresh_dum2) + difind2(1) - 1);

        if ~isempty(difchgtmp1) && ~isempty(difchgtmp2)
            movbeg(k) = difchgtmp1(ft_nearest(difchgtmp1, movbegdum(k)));
            movend(k) = difchgtmp2(ft_nearest(difchgtmp2, movenddum(k)));
        end

    end
    
    movbeg = movbeg(~isnan(movbeg));
    movend = movend(~isnan(movend));
    
    trldat.movmrk{trllop} = [movbeg; movend];
    
end

% % plot to test
% figure;hold on;
% plot(trldat.time{trllop},trldat.veldat{trllop})
% for k=1:size(trldat.movmrk{trllop},2);line([trldat.time{trllop}(trldat.movmrk{trllop}(1,k)) trldat.time{trllop}(trldat.movmrk{trllop}(1,k))],ylim,'Color','g');line([trldat.time{trllop}(trldat.movmrk{trllop}(2,k)) trldat.time{trllop}(trldat.movmrk{trllop}(2,k))],ylim,'Color','r');end


%% characterize turn bouts: duration, velocity, angle

trndur = [];
trndir = [];
trnang = [];
trnvel = [];
trnvelmed = [];
trnsmp = [];
c = 1;
for trllop = 1:length(trldat.time)

    datsmp = data.sampleinfo(trllop,1):data.sampleinfo(trllop,2);
    
    dirnew = unwrap(trldat.dirdat{trllop});
    for ii = find(~isnan(dirnew))
        I = dirnew(ii);
        if ~isempty(find(~isnan(dirnew(ii+1:end)),1,'first'))
            subind = ii+1:find(~isnan(dirnew(ii+1:end)),1,'first')+ii-1;
        else
            subind = ii+1:length(dirnew);
        end
        dirnew(subind) = repmat(I,1,length(subind));
    end
    
    for trnlop = 1:size(trldat.trnmrk{trllop},2)
        
        trndur(c,1) = diff(trldat.trnmrk{trllop}(1:2,trnlop));
        trndir(c,1) = trldat.trnmrk{trllop}(3,trnlop);
        trnang(c,1) = diff(dirnew(trldat.trnmrk{trllop}(1:2,trnlop)));
        trnvel(c,1) = mean(trldat.veldat{trllop}(trldat.trnmrk{trllop}(1,trnlop):trldat.trnmrk{trllop}(2,trnlop)));
        trnvelmed(c,1) = median(trldat.veldat{trllop}(trldat.trnmrk{trllop}(1,trnlop):trldat.trnmrk{trllop}(2,trnlop)));
        trnsmp(c,:) = datsmp(trldat.trnmrk{trllop}(1:2,trnlop));
        
        c=c+1;
        
    end
    
end


%% characterize movement bouts

movdur = [];
movang = [];
movsmp = [];
c = 1;
for trllop = 1:length(trldat.time)

    datsmp = data.sampleinfo(trllop,1):data.sampleinfo(trllop,2);
    
    dirnew = unwrap(trldat.dirdat{trllop});
    for ii = find(~isnan(dirnew))
        I = dirnew(ii);
        if ~isempty(find(~isnan(dirnew(ii+1:end)),1,'first'))
            subind = ii+1:find(~isnan(dirnew(ii+1:end)),1,'first')+ii-1;
        else
            subind = ii+1:length(dirnew);
        end
        dirnew(subind) = repmat(I,1,length(subind));
    end
    
    for movlop = 1:size(trldat.movmrk{trllop},2)

        movdur(c,1) = diff(trldat.movmrk{trllop}(1:2,movlop));
        movang(c,1) = diff(dirnew(trldat.movmrk{trllop}(1:2,movlop)));
        movsmp(c,:) = datsmp(trldat.movmrk{trllop}(1:2,movlop));
        
        c=c+1;
        
    end
    
end

        
%% select turn bouts where median velocity is .25 or less, look at pepisode

sesnam = 'JN140825010';

pepdir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\pEpisode_mult_cycles\';
load(fullfile(pepdir,[sesnam '_pEpisode.mat']))

selsmp = trnsmp(trnvelmed<0.25,:);
UVt = cell(size(unionvector_all));
for sellop = 1:size(selsmp)
    for chnlop = 1:length(unionvector_all)
        UVt{chnlop} = [UVt{chnlop} unionvector_all{chnlop}(:,selsmp(sellop,1):selsmp(sellop,2))];
    end
end


%% select movement bouts, look at pepisode

% sesnam = 'JN140825010';
% 
% pepdir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\pEpisode_mult_cycles\';
% load(fullfile(pepdir,[sesnam '_pEpisode.mat']))

UVm = cell(size(unionvector_all));
for sellop = 1:size(movsmp)
    for chnlop = 1:length(unionvector_all)
        UVm{chnlop} = [UVm{chnlop} unionvector_all{chnlop}(:,movsmp(sellop,1):movsmp(sellop,2))];
    end
end

% average pepisode across session
for chnlop = 1:length(unionvector_all)
    UV1(chnlop,:) = mean(unionvector_all{chnlop}>=1,2); UVt1(chnlop,:) = mean(UVt{chnlop}>=1,2); UVm1(chnlop,:) = mean(UVm{chnlop}>=1,2);
    UV3(chnlop,:) = mean(unionvector_all{chnlop}>=3,2); UVt3(chnlop,:) = mean(UVt{chnlop}>=3,2); UVm3(chnlop,:) = mean(UVm{chnlop}>=3,2);
    UV5(chnlop,:) = mean(unionvector_all{chnlop}>=5,2); UVt5(chnlop,:) = mean(UVt{chnlop}>=5,2); UVm5(chnlop,:) = mean(UVm{chnlop}>=5,2);
end

% add offset for plotting
ofs = 0.1;
ofv = (ofs:ofs:12*ofs)';

clear UVofs*
UVofs1{1} = bsxfun(@plus,UV1(1:12,:),ofv); UVofst1{1} = bsxfun(@plus,UVt1(1:12,:),ofv); UVofsm1{1} = bsxfun(@plus,UVm1(1:12,:),ofv);
UVofs1{2} = bsxfun(@plus,UV1(13:24,:),ofv); UVofst1{2} = bsxfun(@plus,UVt1(13:24,:),ofv); UVofsm1{2} = bsxfun(@plus,UVm1(13:24,:),ofv);
UVofs1{3} = bsxfun(@plus,UV1(25:36,:),ofv); UVofst1{3} = bsxfun(@plus,UVt1(25:36,:),ofv); UVofsm1{3} = bsxfun(@plus,UVm1(25:36,:),ofv);

UVofs3{1} = bsxfun(@plus,UV3(1:12,:),ofv); UVofst3{1} = bsxfun(@plus,UVt3(1:12,:),ofv); UVofsm3{1} = bsxfun(@plus,UVm3(1:12,:),ofv);
UVofs3{2} = bsxfun(@plus,UV3(13:24,:),ofv); UVofst3{2} = bsxfun(@plus,UVt3(13:24,:),ofv); UVofsm3{2} = bsxfun(@plus,UVm3(13:24,:),ofv);
UVofs3{3} = bsxfun(@plus,UV3(25:36,:),ofv); UVofst3{3} = bsxfun(@plus,UVt3(25:36,:),ofv); UVofsm3{3} = bsxfun(@plus,UVm3(25:36,:),ofv);

UVofs5{1} = bsxfun(@plus,UV5(1:12,:),ofv); UVofst5{1} = bsxfun(@plus,UVt5(1:12,:),ofv); UVofsm5{1} = bsxfun(@plus,UVm5(1:12,:),ofv);
UVofs5{2} = bsxfun(@plus,UV5(13:24,:),ofv); UVofst5{2} = bsxfun(@plus,UVt5(13:24,:),ofv); UVofsm5{2} = bsxfun(@plus,UVm5(13:24,:),ofv);
UVofs5{3} = bsxfun(@plus,UV5(25:36,:),ofv); UVofst5{3} = bsxfun(@plus,UVt5(25:36,:),ofv); UVofsm5{3} = bsxfun(@plus,UVm5(25:36,:),ofv);

figure
subplot(1,3,1); hold on
plot(freqs,UVofs1{1},'b')
plot(freqs,UVofst1{1},'r')
plot(freqs,UVofsm1{1},'g')
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Pepisode')
subplot(1,3,2); hold on
plot(freqs,UVofs1{2},'b')
plot(freqs,UVofst1{2},'r')
plot(freqs,UVofsm1{2},'g')
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); hold on
plot(freqs,UVofs1{3},'b')
plot(freqs,UVofst1{3},'r')
plot(freqs,UVofsm1{3},'g')
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', 1 cycle'],'t')

figure
subplot(1,3,1); hold on
plot(freqs,UVofs3{1},'b')
plot(freqs,UVofst3{1},'r')
plot(freqs,UVofsm3{1},'g')
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Pepisode')
subplot(1,3,2); hold on
plot(freqs,UVofs3{2},'b')
plot(freqs,UVofst3{2},'r')
plot(freqs,UVofsm3{2},'g')
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); hold on
plot(freqs,UVofs3{3},'b')
plot(freqs,UVofst3{3},'r')
plot(freqs,UVofsm3{3},'g')
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', 3 cycles'],'t')

figure
subplot(1,3,1); hold on
plot(freqs,UVofs5{1},'b')
plot(freqs,UVofst5{1},'r')
plot(freqs,UVofsm5{1},'g')
title('Array A'); xlabel('Frequency (Hz)'); ylabel('Pepisode')
subplot(1,3,2); hold on
plot(freqs,UVofs5{2},'b')
plot(freqs,UVofst5{2},'r')
plot(freqs,UVofsm5{2},'g')
title('Array B'); xlabel('Frequency (Hz)')
subplot(1,3,3); hold on
plot(freqs,UVofs5{3},'b')
plot(freqs,UVofst5{3},'r')
plot(freqs,UVofsm5{3},'g')
title('Array C'); xlabel('Frequency (Hz)')
suplabel([sesnam ', 5 cycles'],'t')


%% find periods of stillness (no forward movement or turning)

for trllop = 1:length(trldat.time)
    
    movind = zeros(size(trldat.time{trllop}));
    for trnlop = 1:size(trldat.trnmrk{trllop},2)
        movind(trldat.trnmrk{trllop}(1,trnlop):trldat.trnmrk{trllop}(2,trnlop)) = 1;
    end
    for movlop = 1:size(trldat.movmrk{trllop},2)
        movind(trldat.movmrk{trllop}(1,movlop):trldat.movmrk{trllop}(2,movlop)) = 1;
    end

    stlbegdum = find(diff(movind)==-1)+1;
    stlenddum = find(diff(movind)==1);
    
    stlbeg = [];
    stlend = [];
    c = 1;
    for k=1:length(stlbegdum)
        if ~isempty(find(stlenddum-stlbegdum(k)>0,1,'first'))
            stlbeg(c) = stlbegdum(k);
            stlend(c) = stlenddum(find(stlenddum-stlbegdum(k)>0,1,'first'));
            c=c+1;
        end
    end
    
    trldat.stlmrk{trllop} = [stlbeg; stlend];
    
end

stldur = [];
stlsmp = [];
for trllop = 1:length(trldat.time)
    
    datsmp = data.sampleinfo(trllop,1):data.sampleinfo(trllop,2);

    stldur = [stldur; diff(trldat.stlmrk{trllop},1,1)'];
    stlsmp = [stlsmp; datsmp(trldat.stlmrk{trllop})'];
    
end

stlsel = find(stldur>=500);

UVs = cell(size(unionvector_all));
for sellop = 1:size(stlsel)
    for chnlop = 1:length(unionvector_all)
        UVs{chnlop} = [UVs{chnlop} unionvector_all{chnlop}(:,stlsmp(stlsel(sellop),1):stlsmp(stlsel(sellop),2))];
    end
end

figure;plot(freqs,mean(unionvector_all{16}>=3,2))
hold on;plot(freqs,mean(UVs{16}>=3,2),'c')


%% create quiver plots to visualize each trial

trlsel = 1;

dirfil = trldat.dirdat{trlsel};
for ii = find(~isnan(dirfil))
    I = dirfil(ii);
    subind = ii+1:find(~isnan(dirfil(ii+1:end)),1,'first')+ii-1;
    dirfil(subind) = repmat(I,1,length(subind));
end

velfil = trldat.veldat{trlsel};
for ii = find(~isnan(velfil))
    I = velfil(ii);
    subind = ii+1:find(~isnan(velfil(ii+1:end)),1,'first')+ii-1;
    velfil(subind) = repmat(I,1,length(subind));
end


[u,v] = pol2cart(dirfil,velfil);

figure
quiver(trldat.posdat{trlsel}(1,1:100:end),trldat.posdat{trlsel}(2,1:100:end),u(1:100:end),v(1:100:end))
xlim([-11 11])
ylim([-11 11])
title(['Trial ' num2str(trlsel)])


%% make a video

trlsel = 1;

clear M1
figure
c=1;
for timlop = 1:100:length(data.time{trlsel})

    subplot(2,1,1)
    quiver(trldat.posdat{trlsel}(1,timlop),trldat.posdat{trlsel}(2,timlop),u(timlop),v(timlop))
    xlim([-11 11])
    ylim([-11 11])
    title(num2str(trldat.time{trlsel}(timlop)))
    
    subplot(2,1,2)
    delete(h)
    delete(l)
    h = scatter(trldat.time{1},unwrap(dirfil),'.b');
    hold on
    l = line([trldat.time{trlsel}(timlop) trldat.time{trlsel}(timlop)],ylim,'Color','r');
    set(gcf,'Position',[680    68   422   910])
    
    M1(c)=getframe;
    
    c=c+1;
    
end

% movie2avi(M1,fullfile(sessionDir,['trl' num2str(trlsel) '_nav.avi']),'FPS',20)


%% correlate saccade rate with theta peak frequency

clear cfg
cfg.artfctdef                     = [];
cfg.artfctdef.xysaccade           = [];
cfg.artfctdef.xysaccade.sgn       = {'eyeX_B' 'eyeY_B'};
cfg.artfctdef.xysaccade.fltord    = 40;
cfg.artfctdef.xysaccade.lpfreq    = 40;
cfg.artfctdef.xysaccade.threshold = 6e4;
cfg.artfctdef.xysaccade.artpadding= 0.01;

data.artifact = cell(size(data.time));
for trllop = 1:length(data.trial)
    
    e_x = data.trial{trllop}(strmatch(cfg.artfctdef.xysaccade.sgn{1},data.label),:);
    e_y = data.trial{trllop}(strmatch(cfg.artfctdef.xysaccade.sgn{2},data.label),:);

    %low pass filter the eye position data 
    fltord = cfg.artfctdef.xysaccade.fltord;
    lowpasfrq = cfg.artfctdef.xysaccade.lpfreq;
    nyqfrq = data.fsample ./ 2;
    flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
    e_x_lowpassfilter=filtfilt(flt,1, e_x);
    e_y_lowpassfilter=filtfilt(flt,1, e_y);

    %differentiate and multiply with sampling rate to get velocity as deg/sec
    x_vF = diff(e_x_lowpassfilter) .* data.fsample; 
    y_vF = diff(e_y_lowpassfilter) .* data.fsample;

    % % differentiate and multiply with sampling rate without filtering
    % % this gives the eye velocity in the horizontal and vertical domains
    x_v = diff(e_x) .* data.fsample; 
    y_v = diff(e_y) .* data.fsample;

    % combine x- and y-velocity to get eye velocity in degrees/second
    vel = abs(complex(x_v,y_v));
    velF = abs(complex(x_vF,y_vF));

    % set an eye velocity threshold: saccades are eye velocities greater than
    % this threshold
    lim = cfg.artfctdef.xysaccade.threshold; 

    %detect saccade begins and saccade ends
    sacbeg = find(diff(velF > lim) > 0);
    sacend = find(diff(velF > lim) < 0);
    % next three lines unique for JN130604.1
    if velF(1)>lim
        sacbeg = [11 sacbeg];
    end
    if velF(end)>lim
        sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
    end
    
%     artifact = round([sacbeg(:) - cfg.artfctdef.xysaccade.artpadding .* ...
%         data.fsample sacend(:) + cfg.artfctdef.xysaccade.artpadding .* data.fsample]);

    artifact = round([sacbeg(:) sacend(:)]);

    % fix this later; some artifacts overlap
    
    data.artifact{trllop} = data.time{trllop}(artifact);
    

end

% calculate time-windowed saccade rate
data.sacpersec = cell(size(data.time));
for trllop = 1:length(data.time)
    
    % average saccade start & end time to get saccade time
    % (timestamp @ center of saccade)
    arttim = mean(data.artifact{trllop},2);
    
    data.sacpersec{trllop}=nan(size(data.time{trllop}));
    ft_progress('init', 'etf',     'Please wait...');
    for timlop = 501:length(data.time{trllop})-500
        ft_progress(timlop/(length(data.time{trllop})-500), 'Processing event %d from %d', timlop, (length(data.time{trllop})-500));
        data.sacpersec{trllop}(timlop) = length(find(arttim>=data.time{trllop}(timlop)-0.5 & arttim<=data.time{trllop}(timlop)+0.5));
    end
    ft_progress('close')
    
end

pepdir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\pEpisode_mult_cycles\';
pep = load(fullfile(pepdir,'JN140825010_pEpisode.mat'));

freqs = pep.freqs;

chnlop = 16;

sacratperfrq = nan(size(freqs));

sacpersec_smooth = density(data.sacpersec{1},1000,'gauss');
sac_hist=[];
edges = 0:0.01:10;
for frqlop = 1:length(freqs)
    
    sacratperfrq(frqlop) = nanmedian(sacpersec_smooth(pep.unionvector_all{chnlop}(frqlop,:)>=3));

    n = histc(sacpersec_smooth(logical(pep.unionvector(frqlop,:))),edges);
    sac_hist(frqlop,:) = n;

%     figure
%     hist(sacpersec_smooth(logical(pep.unionvector_all{chnlop}(frqlop,:))),500)
%     xlim([0 15]);
%     ylim([0 40000]);
%     title(num2str(freqs(frqlop)))
    
end

% interpolate to plot with correct axis
frq_interp = freqs(1):0.01:freqs(end);
sac_hist_interp = nan(length(frq_interp),size(sac_hist,2));
for edglop = 1:size(sac_hist,2)
    for frqlop = 1:size(sac_hist,1)
        sac_hist_interp(ft_nearest(frq_interp,freqs(frqlop)),edglop) = sac_hist(frqlop,edglop);
    end
    sac_hist_interp(:,edglop)=inpaint_nans(sac_hist_interp(:,edglop),2);
end


%%

r = length(sacpersec_smooth);
edges = 0:0.1:8;
[n, bin] = histc(sacpersec_smooth,edges);
figure;bar(edges,n/r,'histc')
xlabel('saccade rate (Hz)')
ylabel('proportion')
title('proportion of time @ each saccade rate')

pepbin=nan(length(edges)-1,length(freqs));
for binlop = 1:length(edges)-1
%     pepbin(binlop,:)=mean(B(:,bin==binlop),2);
    pepbin(binlop,:)=mean(pep.unionvector(:,bin==binlop),2);
end

centers = 0.05:0.1:7.95;
  
% interpolate to plot with correct axis
frq_interp = freqs(1):0.01:freqs(end);
pepbin_interp = nan(size(pepbin,1),length(frq_interp));
for cntlop = 1:size(pepbin,1)
    for frqlop = 1:size(pepbin,2)
        pepbin_interp(cntlop,ft_nearest(frq_interp,freqs(frqlop))) = pepbin(cntlop,frqlop);
    end
    pepbin_interp(cntlop,:)=inpaint_nans(pepbin_interp(cntlop,:),2);
end

peakfrq = nan(1,size(pepbin_interp,1));
for binlop = 1:size(pepbin_interp,1)
    peakfrq(binlop) = mean(frq_interp(pepbin_interp(binlop,:)==max(pepbin_interp(binlop,:))));
end

figure;hold on
imagesc(frq_interp(ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)),centers(ft_nearest(centers,2):ft_nearest(centers,8)),pepbin_interp(ft_nearest(centers,2):ft_nearest(centers,8),ft_nearest(frq_interp,2):ft_nearest(frq_interp,8)))
axis xy
colormap hot
plot(peakfrq,centers,'Color','b','LineWidth',2)
xlim([2 8])
ylim([2 8])
xlabel('LFP (Pepisode) freq (Hz)')
ylabel('saccade rate (Hz)')
title('JN1306041 ch4 Pepisode peak frequency for each saccade rate bin')