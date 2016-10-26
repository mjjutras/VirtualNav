BRnam = 'JN140812001';

%%
% decimate the recording (downsample from 30 kS/s to 1 kS/s
cd('C:\Users\michael.jutras\Documents\MATLAB\VirtualNav\BlackrockRelatedCode')
SF = 30; % skipfactor
NS6 = decNS6(BRnam, SF);

data = BRtoFT(sesDir,BRnam)