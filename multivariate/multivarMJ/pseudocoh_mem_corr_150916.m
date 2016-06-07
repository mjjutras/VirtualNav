% datdir = 'C:\Data\VR';
datdir = 'R:\Buffalo Lab\Mike\VirtualNav\MAT files\workspace';

load(fullfile(datdir,'cohData150908.mat'))
load(fullfile(datdir,'memsel150908.mat'))
load(fullfile(datdir,'statL_150916.mat'))


%% do univariate correlation to memory

% response data
trialInds = round((1:length(memsel)*2)/2);
Y = memsel(trialInds);

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


for frqlop = 1:size(cohData,2)
    
    [corAA, pAA] = corr(Y,squeeze(cohData(:,frqlop,indAA)));
    [corAB, pAB] = corr(Y,squeeze(cohData(:,frqlop,indAB)));
    [corAC, pAC] = corr(Y,squeeze(cohData(:,frqlop,indAC)));
    [corBB, pBB] = corr(Y,squeeze(cohData(:,frqlop,indBB)));
    [corBC, pBC] = corr(Y,squeeze(cohData(:,frqlop,indBC)));
    [corCC, pCC] = corr(Y,squeeze(cohData(:,frqlop,indCC)));
    
    edg = -0.2:0.01:0.2;
    edx = edg + 0.005;
    [n1,b1] = histc(corAA',edg);
    [n2,b2] = histc(corAB',edg);
    [n3,b3] = histc(corAC',edg);
    [n4,b4] = histc(corBB',edg);
    [n5,b5] = histc(corBC',edg);
    [n6,b6] = histc(corCC',edg);
    
    figure;hold on
%     bar(edx,([n1 n2 n3 n4 n5 n6]),'stacked')
%     legend({'A x A','A x B','A x C','B x B','B x C','C x C'},'Location','northeast')
    h=bar(edx,([n2 n5 n3]),'grouped'); % pay attention to order!
    set(h,'BarWidth',2.5)
    legend({'A x B','B x C','A x C'},'Location','northeast')
    
    h=get(gca);
%     line([median(corAA) median(corAA)],ylim,'Color',([53 42 134])/255,'LineWidth',1.5)
%     line([median(corAB) median(corAB)],ylim,'Color',([12 116 220])/255,'LineWidth',1.5)
%     line([median(corAC) median(corAC)],ylim,'Color',([6 169 193])/255,'LineWidth',1.5)
%     line([median(corBB) median(corBB)],ylim,'Color',([124 191 123])/255,'LineWidth',1.5)
%     line([median(corBC) median(corBC)],ylim,'Color',([240 185 73])/255,'LineWidth',1.5)
%     line([median(corCC) median(corCC)],ylim,'Color',([248 250 13])/255,'LineWidth',1.5)

    line([median(corAB) median(corAB)],ylim,'Color',([12 116 220])/255,'LineWidth',1.5)
    line([median(corBC) median(corBC)],ylim,'Color',([6 169 193])/255,'LineWidth',1.5)
    line([median(corAC) median(corAC)],ylim,'Color',([240 185 73])/255,'LineWidth',1.5)

    line([0 0],ylim,'Color','k','LineStyle','--')

    if frqlop==1
        title({'Theta (3-12) Corr: pseudo-coherence x perf. factor' ...
            ['signtest: AB = ' num2str(signtest(corAB)) '; BC = ' num2str(signtest(corCC)) ...
            '; AC = ' num2str(signtest(corAC))]})
    elseif frqlop==2
        title({'Low Gamma (30-60) Corr: pseudo-coherence x perf. factor' ...
            ['signtest: AB = ' num2str(signtest(corAB)) '; BC = ' num2str(signtest(corCC)) ...
            '; AC = ' num2str(signtest(corAC))]})
    elseif frqlop==3
        title({'High Gamma (60-100) Corr: pseudo-coherence x perf. factor' ...
            ['signtest: AB = ' num2str(signtest(corAB)) '; BC = ' num2str(signtest(corCC)) ...
            '; AC = ' num2str(signtest(corAC))]})
    end
    xlim([-0.2 0.2])
    box on
    grid on
    set(gca,'GridLineStyle','--','FontSize',14)
    set(gcf,'Position',[142   423   630   435])
    xlabel('Correlation (rho)');ylabel('Count')

end

%% do yes/no memory/coherence relationship

% response data
trialInds = round((1:length(memsel)*2)/2);
Y = memsel(trialInds);
Y  = Y > median(Y);

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


for frqlop = 1:size(cohData,2)
    
    figure
    subplot(3,2,1)
    scatter(squeeze(mean(cohData(Y==1,frqlop,indAA),1)),squeeze(mean(cohData(Y==0,frqlop,indAA),1)))
    title('AA')
    subplot(3,2,2)
    scatter(squeeze(mean(cohData(Y==1,frqlop,indAB),1)),squeeze(mean(cohData(Y==0,frqlop,indAB),1)))
    title('AB')
    subplot(3,2,3)
    scatter(squeeze(mean(cohData(Y==1,frqlop,indAC),1)),squeeze(mean(cohData(Y==0,frqlop,indAC),1)))
    title('AC')
    subplot(3,2,4)
    scatter(squeeze(mean(cohData(Y==1,frqlop,indBB),1)),squeeze(mean(cohData(Y==0,frqlop,indBB),1)))
    title('BB')
    subplot(3,2,5)
    scatter(squeeze(mean(cohData(Y==1,frqlop,indBC),1)),squeeze(mean(cohData(Y==0,frqlop,indBC),1)))
    title('BC')
    subplot(3,2,6)
    scatter(squeeze(mean(cohData(Y==1,frqlop,indCC),1)),squeeze(mean(cohData(Y==0,frqlop,indCC),1)))
    title('CC')
    
end

