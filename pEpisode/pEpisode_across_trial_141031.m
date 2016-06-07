BRnam = 'JN140825010'; sessionDir = 'JN_14_08_25_13_12';
% BRnam = 'JN140825011'; sessionDir = 'JN_14_08_25_13_57';

BRDir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';
pepdir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\pEpisode';
trlDir = 'C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\trial data';

%%

NS2 = openNSx(fullfile(BRDir,[BRnam '.NS2']),'read');
NEV = openNEV(fullfile(BRDir,[BRnam '.nev']));
% % NS6 = openNSx(fullfile(datdir,[BRnam '.ns6']));
% % load decimated MAT version of NS6 file instead of directly from NS6 file
load(['C:\Users\michael.jutras\Documents\Virtual Navigation Study\MATLAB\MAT files\NS6 - decimated\' BRnam '_NS6_SF30.mat'])

ns2DTR = NS2.MetaTags.DateTimeRaw;
% ns6DTR = NS6.MetaTags.DateTimeRaw; % (NS6 same as NS2, but Timestamp indicates slightly delayed start time)
nevDTR = NEV.MetaTags.DateTimeRaw;

% get date vectors from BR MetaTags
% (bug in NS2 and NS6 files: the hour value is off (value #5), although it
% is correct in the NEV file)
ns2datevec = [ns2DTR([1 2 4]) nevDTR(5) ns2DTR(6) ns2DTR(7)+ns2DTR(8)/1000];
% ns6datevec = [ns6DTR([1 2 4]) nevDTR(5) ns6DTR(6) ns6DTR(7)+ns6DTR(8)/1000+NS6.MetaTags.Timestamp/NS6.MetaTags.SamplingFreq];
nevdatevec = [nevDTR([1 2 4 5 6]) nevDTR(7)+nevDTR(8)/1000];

fs = 1/NS2.MetaTags.SamplingFreq;

%%
load(fullfile(trlDir,[BRnam '_trldat.mat']))
load(fullfile(pepdir,'JN140825010_pEpisode_141031'));
freqs = (2^(1/8)).^(8:42);

trlsiz = nan(size(data.trl,1),1);
for trllop = 1:size(data.trl,1)
    trlsiz(trllop) = length(data.time{trllop});
end

trl = data.trl;
clear data

% chan_freq_time
UV = nan(length(unionvector_all),length(freqs),max(trlsiz));
parfor chnlop = 1:length(unionvector_all)
    UVdum=nan(size(trl,1),length(freqs),max(trlsiz));
    for trllop = 1:size(trl,1)
        UVdum(trllop,:,1:trlsiz(trllop)) = unionvector_all{chnlop}(:,trl(trllop,1):trl(trllop,2));
    end
    UV(chnlop,:,:) = nanmean(UVdum,1);
end
clear UVdum
