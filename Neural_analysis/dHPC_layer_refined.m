%function dHPC_layer_refined(in)
% here is a refined version of dHPC layer specification based on a series
% of work done by March 2024. 
% Advancement: 1. use the mean of 3 peak DG channels to detect dorsal dentate spikes; 
% 2. classify ds into 2 types to more robustly estimate lpp & mpp sublayers;
% 3. specify dca1 slm layer
% Mengni Wang, 04/2024

clear all;set(0,'DefaultFigureColormap',feval('turbo'));close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet');
in=11;
savedir=SessionSet{in};
cd(savedir);
load('Dentate_Spike_CSD.mat');
load('LFP_CSD.mat');

% 1. specify dca1 pyr, rad, & channel id for dentate spike classification;
dca1_pyr=floor(median(Ripple_Channels));
lfprange=max(lfpcsd_ref,[],1)-min(lfpcsd_ref,[],1);
[~,I]=max(lfprange(SPW_Channels));
dca1_rad=SPW_Channels(I);
[~,I]=max(lfprange(Dorsal_DG_CSD_Channels));
dg=Dorsal_DG_CSD_Channels(I);

figure(in);
subplot(2,3,1);
plot(lfprange,'ko-');xlabel('Channel ID');ylabel('Raw LFP Range (mV)');
xline(SPW_Channels([1,end]),'r--');hold on;plot(dca1_rad,lfprange(dca1_rad),'r*');
xline(Ripple_Channels([1,end]),'m--');hold on;plot(dca1_pyr,lfprange(dca1_pyr),'m*');
xline(Dorsal_DG_CSD_Channels([1,end]),'b--');hold on;plot(dg,lfprange(dg),'b*');
title('Please check if ranges are reasonable to continue :')
%pause;
%-------------------------------------------------------------------------------------
% 2. detect dorsal dentate spikes & calculate csd;
AmpCut=0.5;
dslfp_d=mean(lfpcsd_ref(:,dg-1:dg+1),2);
[dspeaks_d3,ifplot]=Detect_Dentate_Spikes(dslfp_d,AmpCut,LFP_Frequency,lfptime(:,1));
cd(savedir);figure_title='dentate_spikes_dorsal_detection_3channels';save_current_figure(figure_title);
disp(['Previous ds N = ',num2str(size(dspeaks_d,1)),' ; Current ds N = ',num2str(size(dspeaks_d3,1))]);
%pause;
%-------------------------------------------------------------------------------------
v=[-0.23,-0.08,0.62,-0.08,-0.23]; csdwin=0.03;spacing=4*10^(-5);
dnum=round(csdwin*LFP_Frequency);
if size(Filtered_LFP,1)<size(Filtered_LFP,2)
    Filtered_LFP=Filtered_LFP';
end
dDS_CSD=zeros(size(channelcsd,1)-4,2*dnum+1,size(dspeaks_d3,1));
for p=1:size(dspeaks_d3,1)
    ind=[dspeaks_d3(p,3)-dnum:dspeaks_d3(p,3)+dnum];
    dslfp=Filtered_LFP(ind,:);
    for t=1:size(dslfp,1)
        u=dslfp(t,:);
        dDS_CSD(:,t,p) = conv(u,v,'valid');
    end
end
dDS_CSD_Profile=zeros(size(channelcsd,1)-4,size(dspeaks_d3,1));
for p=1:size(dDS_CSD,3)
    dDS_CSD(:,:,p)=imgaussfilt(dDS_CSD(:,:,p),2);
    ind=[dnum+1-(dspeaks_d3(p,3)-dspeaks_d3(p,1)):dnum+1+dspeaks_d3(p,2)-dspeaks_d3(p,3)];
    a=dDS_CSD(:,ind,p);
    dDS_CSD_Profile(:,p)=min(a,[],2);
end

% 3. classify for 2 types of ds;
dds_csd_chrange=[Not_HPC(end)+1:dca1_rad];
ddscsd=dDS_CSD_Profile(dds_csd_chrange-2,:)';
[idx0,C0] = kmeans(ddscsd,2,'Distance','correlation');C0=C0';
maxsink=[];
for i=1:2
    [pk,loc]=findpeaks(-C0(:,i),'MinPeakHeight',0,'MinPeakDistance',5);
    maxsink(i,1)=loc(end);
end
if maxsink(1)>maxsink(2)
    idx=3*ones(size(idx0))-idx0;
    C=C0(:,[2,1]);
else
    idx=idx0;C=C0;
end
dspeaks_d3(:,6)=idx;

% 4. estimate lpp, mpp layer and edge;
dif12=C(:,1)-C(:,2);
[pk1,loc1]=findpeaks(dif12,'MinPeakHeight',0,'MinPeakDistance',5);
[pk2,loc2]=findpeaks(-dif12,'MinPeakHeight',0,'MinPeakDistance',5);
[~,I1]=sort(pk1,'descend');I1=I1(1:2);
[~,I2]=sort(pk2,'descend');I2=I2(1:2);
locpks=[[loc1(I1);loc2(I2)],[pk1(I1);-pk2(I2)]];
[~,I]=sort(locpks(:,1));
locpks=locpks(I,:);
if locpks(1,2)<0
    disp('Error in locpks!');
end
LPP_edge=[];MPP_edge=[];
if isempty(find( C(1:locpks(1,1),2)>0,1,'last' ))
    LPP_edge(1)=1;
else
    LPP_edge(1)=find( C(1:locpks(1,1),2)>0,1,'last' )+1;
end
t1=find(locpks(:,2)<0,1);
LPP_edge(2)=find( dif12(locpks(1,1):locpks(t1,1))>0,1,'last' )+locpks(1,1)-1;
MPP_edge(1)=LPP_edge(2)+1;
MPP_edge(2)=find( C(locpks(t1,1):end,1)>=0,1)+locpks(t1,1)-1;

t2=find(locpks(t1+1:end,2)<0,1)+t1;
MPP_edge(3)=find( C(locpks(t1,1):locpks(t2,1),1)>0,1,'last' )+locpks(t1,1);
p2=find(locpks(t2:end,2)>=0,1)+t2-1;
MPP_edge(4)=find( dif12(locpks(t2,1):locpks(p2,1))<0,1,'last' )+locpks(t2,1)-1;
LPP_edge(3)=MPP_edge(4)+1;
if isempty(find( C( locpks(p2,1):end ,2 )>0,1 ))
    LPP_edge(4)=size(C,1);
else
    LPP_edge(4)=find( C( locpks(p2,1):end ,2 )>0,1 )+locpks(p2,1)-2;
end
%LPP_edge=[33,38];MPP_edge=[26,32];
[~,I]=min(C(LPP_edge(1):LPP_edge(2),2));LPP_center(1)=I+LPP_edge(1)-1;
[~,I]=min(C(LPP_edge(3):LPP_edge(4),2));LPP_center(2)=I+LPP_edge(3)-1;
[~,I]=min(C(MPP_edge(1):MPP_edge(2),1));MPP_center(1)=I+MPP_edge(1)-1;
[~,I]=min(C(MPP_edge(3):MPP_edge(4),1));MPP_center(2)=I+MPP_edge(3)-1;
[~,I]=max(C(MPP_edge(2):MPP_edge(3),1));CA3DG=I+MPP_edge(2)-1;

figure(in);
subplot(2,3,2);imagesc(C0);set(gca,'YDir','normal');xlabel('raw type');title(num2str([sum(idx0==1),sum(idx0==2)]));
subplot(2,3,4);imagesc(C);set(gca,'YDir','normal');xlabel('curated type');title(num2str([sum(idx==1),sum(idx==2)]));
yline(LPP_edge,'b-');yline(MPP_edge,'r-'); 
yline(LPP_center,'b--');yline(MPP_center,'r--');yline(CA3DG,'k--');
subplot(2,3,5);
for i=1:2
    hold on;plot(C(:,i),'o-');
end
hold on;plot(dif12,'ko-');xlabel('DG Channels');ylabel('Min CSD value');
xline(LPP_edge,'b-');xline(MPP_edge,'r-'); yline(0);
xline(LPP_center,'b--');xline(MPP_center,'r--');xline(CA3DG,'k--');
title(savedir);legend('MPP DS','LPP DS','MPP-LPP');

LPP_edge=LPP_edge+dds_csd_chrange(1)-1;
MPP_edge=MPP_edge+dds_csd_chrange(1)-1;
LPP_center=LPP_center+dds_csd_chrange(1)-1;
MPP_center=MPP_center+dds_csd_chrange(1)-1;
CA3DG=CA3DG+dds_csd_chrange(1)-1;
SLM_Channels=[LPP_edge(end)+1:dca1_rad-1];
if ~isempty(SLM_Channels)
    dca1_lm=floor(median(SLM_Channels));
else
    dca1_lm=dca1_rad-1;
end
dHPC_layerchannel=[dca1_pyr,dca1_rad,dca1_lm,max(LPP_center),max(MPP_center),CA3DG,min(MPP_center),min(LPP_center)]';% 8 LAYERS
dHPC_layerchannel(:,2:3)=channelcsd(dHPC_layerchannel(:,1),3:4);
dHPC_layers=struct;
dHPC_layers.dca1_pyr=dca1_pyr;
dHPC_layers.dca1_rad=dca1_rad;
dHPC_layers.dca1_lm=dca1_lm;
dHPC_layers.lpp_center=LPP_center;
dHPC_layers.mpp_center=MPP_center;
dHPC_layers.lpp_edge=LPP_edge;
dHPC_layers.mpp_edge=MPP_edge;
dHPC_layers.ca3dg=CA3DG;
dHPC_layers.ripple_channel=Ripple_Channels;
dHPC_layers.spw_channel=SPW_Channels;
dHPC_layers.slm_channel=SLM_Channels;
dHPC_layers.not_hpc=Not_HPC;
cd(savedir);save('dHPC_layer_refined','channelcsd','dDS_CSD','dDS_CSD_Profile','dspeaks_d3',...
    'C','dds_csd_chrange','Filtered_LFP','dHPC_layers','dHPC_layerchannel','Not_HPC','-v7.3');

[M1,I1]=min(dDS_CSD_Profile(dds_csd_chrange,:),[],1);[b,I]=sort(I1);
subplot(2,3,[3,6]);
imagesc([1 size(dDS_CSD_Profile,2)],[3 size(dDS_CSD_Profile,1)+2],dDS_CSD_Profile(:,I));set(gca,'YDir','normal');
yline(dHPC_layerchannel(:,1),'k--');yline(Not_HPC([1,end]),'k');
yline(MPP_edge,'r-');yline(LPP_edge,'b-');yline(Ripple_Channels([1,end]),'m-');yline(SPW_Channels([1,end]),'r-');
xlabel('Dentate Spike Events');ylabel('Channel ID');title(['Dorsal dentate spikes sorted by sink']);
cd(savedir);figure_title='dHPC_layer_refined';save_current_figure(figure_title);

load('Spike_Data.mat','unit_metric_goodokay');
depth=ceil(max(unit_metric_goodokay(:,5))/40+2)*40;
figure(in+20);histogram(unit_metric_goodokay(:,5),[0:40:depth]);
xline(dHPC_layerchannel(:,2),'r');xline(dHPC_layerchannel(1:3,2),'k');
xline(dHPC_layerchannel(6,2),'m');xline(channelcsd(Not_HPC(end)+1,3),'b');
xline(channelcsd(Not_HPC(1)-1,3),'b');
title(savedir);xlabel('Depth on Probe (um)');ylabel('Unit N');

disp(['dHPC_layer_refined : Session ',num2str(in),' Completed!']);