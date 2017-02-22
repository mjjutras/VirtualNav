ft_defaults

% BRnam = 'JN140812001';
% BRnam = 'JN140813002';
% BRnam = 'JN140814002';
% BRnam = 'JN140815002';
% BRnam = 'JN140818002';
% BRnam = 'JN140911002';
% BRnam = 'JN141013002';
BRnam = 'JN141212003';

%%
% decimate the recording (downsample from 30 kS/s to 1 kS/s
cd('C:\Users\michael.jutras\Documents\MATLAB\VirtualNav\BlackrockRelatedCode')
SF = 30; % skipfactor
[~,hostname]= system('hostname');
if strncmp(hostname,'RBU-MikeJ',9) && ~strncmp(hostname,'RBU-MikeJ2',10)
    decDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated';
elseif strncmp(hostname,'RBU-MikeJ2',10)
    decDir = 'C:\Data\MAT\NS6 - decimated';
end
netDir = 'R:\Mike\VirtualNavigationProject\MATFiles\NS6decimated';

if exist(fullfile(decDir,[BRnam '_NS6_SF30.mat']),'file')
    load(fullfile(decDir,[BRnam '_NS6_SF30.mat']))
elseif exist(fullfile(netDir,[BRnam '_NS6_SF30.mat']),'file')
    copyfile(fullfile(netDir,[BRnam '_NS6_SF30.mat']),fullfile(decDir,[BRnam '_NS6_SF30.mat']))
    load(fullfile(decDir,[BRnam '_NS6_SF30.mat']))
else
    NS6 = decNS6(BRnam,decDir,SF);
end

%%
% xB = [544000:574000]; % JN140812001
% xB = [303211:333211]; % JN140813002
% xB = [1:500000]; % JN140815002
% xB = [172949:202949]; % JN140818002
% xB = [152350:182350]; % JN140911002
% xB = [135207:165207]; % JN141013002
xB = [1:size(NS6.Data,2)]; % JN141212003

channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];

for chnlop = 1:36
    
    if chnlop>=1 && chnlop<=12
        figure(1); set(gcf,'Position',[92 510 560 420])
        yaxplc = chnlop;
    elseif chnlop>=13 && chnlop<=24
        figure(2); set(gcf,'Position',[690 510 560 420])
        yaxplc = chnlop-12;
    elseif chnlop>=25 && chnlop<=36
        figure(3); set(gcf,'Position',[1209 510 560 420])
        yaxplc = chnlop-24;
    end
    
    hold on
    dat = ft_preproc_detrend(NS6.Data(chnlop,xB), 1, length(xB));
    plot(NS6.Data(chnlop,xB)+2000*yaxplc,'Color',channelcolormap(yaxplc,:))
    
%     xlim([0 30000])
    ylim([-5000 3e4])
    
end

% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN140812001_544000-574000_probeA.fig
% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN140818002_172949-202949_probeA.fig
% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN140911002_152350-182350_probeA.fig
% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN141013002_135207-165207_probeB_lines10410-11410.fig

%%
% xB = [554000:555000]; % JN140812001
% xB = [184128:185128]; % JN140818001
% xB = [176301:177301]; % JN140911002 % not as good
% xB = [171584:172584]; % JN140911002 % use this one
% xB = [166225:167225]; % JN140911002 % not as good
xB = [145617:146617]; % JN141013002

channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];

for chnlop = 1:36
    
    if chnlop>=1 && chnlop<=12
        figure(1); set(gcf,'Position',[92 510 560 420])
        yaxplc = chnlop;
    elseif chnlop>=13 && chnlop<=24
        figure(2); set(gcf,'Position',[690 510 560 420])
        yaxplc = chnlop-12;
    elseif chnlop>=25 && chnlop<=36
        figure(3); set(gcf,'Position',[1209 510 560 420])
        yaxplc = chnlop-24;
    end
    
    hold on
%     dat = ft_preproc_detrend(NS6.Data(chnlop,xB), 1, length(xB));
    dat = ft_preproc_highpassfilter(ft_preproc_detrend(NS6.Data(chnlop,xB), 1, length(xB)),1000,0.8);
    plot(dat+500*yaxplc,'Color',channelcolormap(yaxplc,:))
    
    xlim([0 1000])
    ylim([-1000 7000])
    
end

% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN140812001_554000-555000_probeA.fig
% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN140818002_184128-185128_probeA.fig
% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN140911002_171584-172584_probeA.fig
% R:\Buffalo Lab\Flexshaft\Figs_for_ChrisLewis_manuscript\JN141013002_145617-146617_probeB.fig
