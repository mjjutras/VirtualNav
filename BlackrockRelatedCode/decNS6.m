function NS6 = decNS6(BRnam,decDir,SF)

% 161026 - revised to require entering decdir, the directory to save the
% resulting downsampled NS6 structure

fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';
% fildir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\Data\';
% fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data - Peepers\';

NS6 = openNSx(fullfile(fildir,[BRnam '.ns6']),'read','channels',1:36,'skipfactor',SF);

try
    save(fullfile(decDir, [BRnam '_NS6_SF' num2str(SF) '.mat']),'NS6')
catch
    save(fullfile(decDir, [BRnam '_NS6_SF' num2str(SF) '.mat']),'NS6','-v7.3')
end
disp(['Created ' fullfile(decDir, [BRnam '_NS6_SF' num2str(SF) '.mat'])])
