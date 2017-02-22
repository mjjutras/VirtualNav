% processNav: process a navigation recording from beginning to end
% first need to run the appropriat version of MakeParFiles in the directory
% containing the python behavioral data log

BRnam = 'JN140812001';
sesDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\2014\JN_14_08_12_12_26';

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
if exist(fullfile(decDir,[BRnam '_NS6_SF30.mat']),'file')
    load(fullfile(decDir,[BRnam '_NS6_SF30.mat']))
else
    NS6 = decNS6(BRnam,decDir,SF);
end

sesNam = sesDir(find(sesDir=='\',1,'last')+1:end);
trldatFil = fullfile('R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\trldat',[sesNam '_trldat.mat']);
if exist(trldatFil,'file')
    load(trldatFil)
else
    cd('C:\Users\michael.jutras\Documents\MATLAB\VirtualNav\PARtoTRLDAT')
    trldat = PARtoTRLDAT(sesDir,1);
    movefile(fullfile(sesDir,[sesNam '_trldat.mat']),trldatFil)
end

cd('C:\Users\michael.jutras\Documents\MATLAB\VirtualNav\BRtoFT')
data = BRtoFT(sesDir,decDir,BRnam);

