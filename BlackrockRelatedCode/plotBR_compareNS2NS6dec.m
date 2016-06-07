BRnam = 'JN150209001';

BRdir = 'R:\Buffalo Lab\Virtual Navigation\Recording Data\Blackrock Data';

NS2 = openNSx(fullfile(BRdir,[BRnam '.NS2']),'read','c:1');
NS6 = openNSx(fullfile(BRdir,[BRnam '.ns6']),'read','c:1');

% load decimated MAT version of NS6 file
dec = load(fullfile(decDir,[BRnam '_NS6_SF30.mat']));

chn = 1;
figure;hold on
plot((1:20)/1000,NS2.Data(chn,1:20))
plot((1:600)/30000+(NS6.MetaTags.Timestamp/30000),NS6.Data(chn,1:600),'r')
plot((1/30000:0.001:0.02)+(NS6.MetaTags.Timestamp/30000),dec.NS6.Data(1,1:20),'g')
