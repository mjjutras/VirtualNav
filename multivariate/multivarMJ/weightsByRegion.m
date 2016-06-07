% run domultivarcohMJ first

% datdir = 'C:\Data\VR';
datdir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\workspace';

load(fullfile(datdir,'statL_150916.mat'))


%%

nTimes = size(params.timeBins,1);
nFreqs = size(params.freqBins,1);
spectTimeFreq        = NaN(size(params.freqBins,1),size(params.timeBins,1));
regionDataAll        = NaN(6,size(params.timeBins,1));
regionDataNonZeroAll = NaN(6,size(params.timeBins,1));
movementClass        = {};

figs = [];

figs_subj = struct('Spect_Freq',[],'Spect_Time',[],...
    'Spect_TxF',[],'Region_Bar',[],'nElecs',[]);

% load the weights for this subject
out = [];

nFeatures = length(res.A{1});
nFreqs    = size(params.freqBins,1);
nTimes    = size(params.timeBins,1);
nElecs = nFeatures/nFreqs;

% % pval of different time bins
% p = 1-perc;
% thresh = p < .05;
% sigTimes = p < thresh;
% out.sigTimes = sigTimes;


% average across all folds and reshape into freq x elec x time
meanTmp = mean(horzcat(res.A{:}),2);
meanTmp = reshape(meanTmp,nFreqs,nElecs);
meanWeightsPerTime = meanTmp;

% filter by regions
% A-A
indAA = intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2)));
% A-B
indAB = union(intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2))), ...
    intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2))));
% A-C
indAC = union(intersect(strmatch('A',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2))), ...
    intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('A',statL.labelcmb(:,2))));
% B-B
indBB = intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2)));
% B-C
indBC = union(intersect(strmatch('B',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2))), ...
    intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('B',statL.labelcmb(:,2))));
% C-C
indCC = intersect(strmatch('C',statL.labelcmb(:,1)),strmatch('C',statL.labelcmb(:,2)));
out.meanWeightsPerTime   = meanWeightsPerTime;
out.meanWeightsPerTimeAA = meanWeightsPerTime(:,indAA);
out.meanWeightsPerTimeAB = meanWeightsPerTime(:,indAB);
out.meanWeightsPerTimeAC = meanWeightsPerTime(:,indAC);
out.meanWeightsPerTimeBB = meanWeightsPerTime(:,indBB);
out.meanWeightsPerTimeBC = meanWeightsPerTime(:,indBC);
out.meanWeightsPerTimeCC = meanWeightsPerTime(:,indCC);


% elec_order = [find(hipp_elecs) find(ec_elecs) find(mtl_elecs) ...
%               find(temp_elecs) find(frontal_elecs) find(occ_elecs) ...
%               find(par_elecs) find(other_elecs)];
% out.meanWeightsPerTimeSort = meanWeightsPerTime(:,elec_order,:);

regions     = {'H','EC','MTL','TC','FC','OC','PC','X'};
regionCount = sum([hipp_elecs' ec_elecs' ...
    mtl_elecs' temp_elecs' frontal_elecs'...
    occ_elecs' par_elecs' other_elecs']);
noElecs     = regionCount == 0;


out.regionCutoffs = cumsum(regionCount(~noElecs));
out.regions       = regions(~noElecs);

if ~isempty(out)
    figs_subj.subj = subj;
    figs_subj.Region_Bar = {};
    %         movementClass{s} = out.movementClass;
    
    % for each subject, we will plot:
    %    elec x freq spectrogram
    %    elec x time spectrogram
    %    time x freq spectrogram
    
    %% FIGURE 1 - elec x freq (averaged across time)
    fname = fullfile(figDir,[subj 'Spect_Freq.png']);
    figs_subj.Spect_Freq = fname;
    
    % create y labels based on freq bins
    yBinsStrF = {};
    for x = 1:size(params.freqBins,1)
        yBinsStrF{x} = [num2str(params.freqBins(x,1)), '-', num2str(params.freqBins(x,2))];
    end
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(1)
        clf
        imagesc(squeeze(mean(abs(out.meanWeightsPerTimeSort),3)))
        colormap('jet')
        colorbar
        set(gca,'ytick',1:nFreqs)
        set(gca,'yticklabel',yBinsStrF)
        set(gca,'xtick',[])
        set(gca,'fontsize',16)
        set(gca,'xtick',out.regionCutoffs);
        set(gca,'xticklabel','|');
        set(gca,'xticklabel',out.regions);
        set(gca,'XTickLabelRotation',270)
        xlabel('Electrodes','fontsize',20)
        ylabel('Frequency','fontsize',20)
        print('-dpng','-loose',fname);
    end
    
    
    %% FIGURE 2 - elec x time (averaged across freq)
    fname = fullfile(figDir,[subj 'Spect_Time.png']);
    figs_subj.Spect_Time = fname;
    
    % create y labels based on freq bins
    yBins    = round((params.timeBins) / 10) / 100;
    yBinsStrT = {};
    for x = 1:size(yBins,1)
        yBinsStrT{x} = [num2str(yBins(x,1)), '-', num2str(yBins(x,2))];
    end
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        imagesc(squeeze(mean(abs(out.meanWeightsPerTimeSort),1))');
        colormap('jet')
        colorbar
        set(gca,'ytick',1:nTimes)
        set(gca,'yticklabel',yBinsStrT)
        set(gca,'xtick',[])
        set(gca,'fontsize',16)
        set(gca,'xtick',out.regionCutoffs);
        set(gca,'xticklabel','|');
        set(gca,'xticklabel',out.regions);
        set(gca,'XTickLabelRotation',270)
        xlabel('Electrodes','fontsize',20)
        ylabel('Time Bin','fontsize',20)
        print('-dpng','-loose',fname);
        
    end
    
    
    %% FIGURE 3 - time x freq (averaged across elecs)
    fname = fullfile(figDir,[subj 'Spect_TxF.png']);
    figs_subj.Spect_TxF = fname;
    spectTimeFreq(:,:,s) = squeeze(mean(abs(out.meanWeightsPerTimeSort),2));
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(3)
        clf
        imagesc(squeeze(mean(abs(out.meanWeightsPerTimeSort),2)));
        colormap('jet')
        colorbar
        
        set(gca,'xtick',1:nTimes)
        set(gca,'xticklabel',yBinsStrT)
        set(gca,'ytick',1:nFreqs)
        set(gca,'yticklabel',yBinsStrF)
        xlabel('Time Bin','fontsize',20)
        ylabel('Frequency','fontsize',20)
        set(gca,'fontsize',20)
        print('-dpng','-loose',fname);
    end
    
    %% FIGURE 4 - mean absolute weights by region for each time
    labels = params.timeBinLabels;
    if isempty(labels);labels=repmat({''},1,size(params.timeBins,1));end
    regions     = {'H','EC','MTL','TC','FC','OC','PC','X'};
    fieldsToUse = {'meanWeightsPerTimeHipp','meanWeightsPerTimeEC',...
        'meanWeightsPerTimeMTL','meanWeightsPerTimeTC',...
        'meanWeightsPerTimeFC','meanWeightsPerTimeOC',...
        'meanWeightsPerTimePC','meanWeightsPerTimeOth'};
    
    % this is stupid. I'm basically doing this twice so that I can find
    % out what the axis ranges will be and set them all the same
    ylims = [0 0];
    for t = 1:size(params.timeBins,1)
        regionData        = NaN(1,length(fieldsToUse));
        regionDataNonZero = NaN(1,length(fieldsToUse));
        for r = 1:length(fieldsToUse)
            weights       = out.(fieldsToUse{r})(:,:,t);
            isZero        = weights == 0;
            regionData(r)        = mean(nanmean(abs(weights),2));
            regionDataAll(r,s,t) = regionData(r);
            
            regionDataNonZero(r)        = nanmean(abs(weights(~isZero)));
            regionDataNonZeroAll(r,s,t) = regionDataNonZero(r);
        end
        ylims(1) = max([ylims(1) max(regionData)]);
        ylims(2) = max([ylims(2) max(regionDataNonZero)]);
    end
    ylims = ceil(ylims * 100)/100;
    
    % loop over each time bin
    for t = 1:size(params.timeBins,1)
        
        % average weights within region, both including and excluding
        % zero weights
        regionData        = NaN(1,length(fieldsToUse));
        regionDataNonZero = NaN(1,length(fieldsToUse));
        for r = 1:length(fieldsToUse)
            weights       = out.(fieldsToUse{r})(:,:,t);
            isZero        = weights == 0;
            regionData(r)        = mean(nanmean(abs(weights),2));
            regionDataNonZero(r) = nanmean(abs(weights(~isZero)));
        end
        
        % plot bar
        plotData = {regionData,regionDataNonZero};
        ylabels  = {'Mean Abs Weights','Mean Abs NonZero Weights'};
        for panel = 1:2
            
            fname = fullfile(figDir,[subj '_bar_region_' labels{t} '_' strrep(ylabels{panel},' ','_') '.eps']);
            figs_subj.Region_Bar{t,panel} = fname;
            
            if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
                figure(4)
                clf
                
                %                     subplot(2,1,panel)
                h=bar(plotData{panel},'linewidth',2,'facecolor',[.6 .6 .6]);
                set(gca,'xtick',1:length(regions));
                xNames = regions;
                [xNames{isnan(plotData{panel})}] = deal('');
                set(gca,'xticklabel',xNames)
                ylabel(ylabels{panel});
                set(gca,'fontsize',16)
                set(gca,'xlim',[0 length(regions)+1])
                set(gca,'ylim',[0 ylims(panel)])
                grid on
                set(gca,'gridlinestyle',':');
                h=title([labels{t} ' Period'],'fontsize',20);
                set(h,'fontweight','normal');
                print('-depsc2','-loose',fname);
            end
            
        end
    end
    
end
figs = [figs;figs_subj];


