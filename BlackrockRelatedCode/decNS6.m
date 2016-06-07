function NS6 = decNS6(BRnam, SF)

fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data\';
% fildir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\Data\';
% fildir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data - Peepers\';

NS6 = openNSx(fullfile(fildir,[BRnam '.ns6']),'read','channels',1:36,'skipfactor',SF);

decdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated';

try
    save(fullfile(decdir, [BRnam '_NS6_SF' num2str(SF) '.mat']),'NS6')
catch
    save(fullfile(decdir, [BRnam '_NS6_SF' num2str(SF) '.mat']),'NS6','-v7.3')
end
disp(['Created ' fullfile(decdir, [BRnam '_NS6_SF' num2str(SF) '.mat'])])
