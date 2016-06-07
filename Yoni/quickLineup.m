

%Preparsed data from vr computer log
data = open('C:\Users\Yoni Browning\Documents\SomeTestData_082\JN_14_08_21_14_00\navData.mat')
data = data.data;
%1000 Hz black rock files
openNEV('C:\Users\Yoni Browning\Documents\SomeTestData_082\JN140821002.nev')
openNSx('read','C:\Users\Yoni Browning\Documents\SomeTestData_082\JN140821002.ns2')

%Empty structure
LFPData =[];

%Because position data is highly periodic, we must align files using eye
%data. 

%Number of trials to be aligned. Not sure why, but there seem to be fewer
%trials on bk compomputer than vr computer.
num_trls = 125;

%Build LFPData struct;
% for trl = 1:num_trls;
%     LFPData.time{trl} = data.time{trl};
%     LFPData.eyedat{trl} = data.eyedat{trl};
%     LFPData.posdat{trl} = data.posdat{trl};
%     LFPData.dirdat{trl} = data.dirdat{trl};
%     LFPData.banpos{trl} = data.banpos{trl};
%     LFPData.baneat{trl} = data.baneat{trl};
%     LFPData.neural{trl} = [];
%     LFPData.bk_pos{trl} = [];
%     LFPData.bk_eye{trl} = [];
% end


for trl =1:num_trls;
%Find start of trial on bk computer
Xs = xcorr(data.eyedat{trl}(1,:),double(NS2.Data(37,:)));
[~,Is] = max(Xs);
trl_start = round(length(Xs)/2-Is);

%Find end of trial on bk computer
Xe = xcorr(data.eyedat{trl+1}(1,:),double(NS2.Data(37,:)));
[~,Ie] = max(Xe);
trl_end = round(length(Xe)/2-Ie);

ttl = strcat('trial ',num2str(trl));
disp(ttl)
subplot(2,1,1);
plot(NS2.Data(39,trl_start:trl_end),NS2.Data(40,trl_start:trl_end)); drawnow;
axis square
title(strcat(ttl,' BK'))
subplot(2,1,2);
plot(data.posdat{trl}(1,:),data.posdat{trl}(2,:)); drawnow;
axis square
title(strcat(ttl,' VR'))


% Save as a new file structure containing all the necessary info.
LFPData.time{trl} = data.time{trl};
LFPData.eyedat{trl} = data.eyedat{trl};
LFPData.posdat{trl} = data.posdat{trl};
LFPData.dirdat{trl} = data.dirdat{trl};
LFPData.banpos{trl} = data.banpos{trl};
LFPData.baneat{trl} = data.baneat{trl};
LFPData.neural{trl} = NS2.Data(1:36,trl_start:trl_end);
LFPData.bk_pos{trl} = NS2.Data(39:40,trl_start:trl_end);
LFPData.bk_eye{trl} = NS2.Data(37:38,trl_start:trl_end);
end
save('LFPData_JN140821.mat','LFPData')
