function NS6 = decNS6(BRnam,decdir,SF)

% 161026 - revised to require entering decdir, the directory to save the
% resulting downsampled NS6 structure

% load from network by defauly
ns6dir = 'R:\Virtual Navigation\Recording Data\Blackrock Data\';
% load from local drive if available (faster but data must be copied from network first
[~,hostname]= system('hostname');
if strncmp(hostname,'RBU-MikeJ2',10)
    locDir = 'C:\Data\Blackrock_VR';
    if exist(fullfile(locDir,[BRnam '.ns6']),'file')
        ns6dir = locDir;
    end
end

NS6 = openNSx(fullfile(ns6dir,[BRnam '.ns6']),'read','channels',1:36,'skipfactor',SF);

try
    save(fullfile(decdir, [BRnam '_NS6_SF' num2str(SF) '.mat']),'NS6')
catch
    save(fullfile(decdir, [BRnam '_NS6_SF' num2str(SF) '.mat']),'NS6','-v7.3')
end
disp(['Created ' fullfile(decdir, [BRnam '_NS6_SF' num2str(SF) '.mat'])])
