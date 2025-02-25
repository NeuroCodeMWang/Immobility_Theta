function analyze_thetacycle(in)
% the goal here is to study theta cyles during run vs pause; assymetry of
% theta cycle shape; the relation with sink & source at OML & MML of DG
%--------------------------------------------------------------------------------
close all;set(0,'DefaultFigureColormap',feval('turbo'));%clear all;in=1
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
savedir=SessionSet16{in};
cd(savedir);
load('Position_Data_Maze.mat');
load('InOutBound_Behavior_Analysis.mat','pos_behave_marker');
load('theta_delta_LFP_HPC.mat', 'thetaphase', 'lfptime', 'thetapower');
load('HPC_Layer_Channel_LFP.mat', 'LFP_HPC','LFP_Frequency');
load('CSD_LFP_Visualization.mat', 'CSDall', 'dHPC_layer7channel3', 'channels');
dHPC_layer7channel3(:,4)=0;
for i=1:7
    dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2));
end

% 0. define reftheta as from OML DG & identify individual theta cycles
filted=zscore(LFPfilt(LFP_HPC(:,4),LFP_Frequency));
phase=thetaphase(:,4)-180;phase=phase+(phase<0)*360; % make thetapeak as 180; trough as 0
reftheta=[lfptime(:,1),filted,zscore(thetapower(:,4)),phase];
indtrial=lfptime(:,2)>0;
reftheta=reftheta(indtrial,:);
lfpz_trial=zscore(LFP_HPC(indtrial,:),0,1);
lfptime_trial=lfptime(indtrial,:);
CSDall_trial=CSDall(:,indtrial);
reftheta(:,5)=CSDall_trial(dHPC_layer7channel3(5,4),:)';
reftheta(:,6)=CSDall_trial(dHPC_layer7channel3(5,4)+1,:)';
thetacycle=identify_thetacycles(reftheta);
thetacycle(:,17:20)=thetacycle(:,10:13);
ruler=Position_Data(:,1);template=thetacycle(:,6);[outputindex,error]=match(template,ruler,0);
thetacycle(:,9:12)=[outputindex,Position_Data(outputindex,5),pos_behave_marker(outputindex,:)];
ind=min(thetacycle(:,11:12),[],2)~=-100;
thetacycle=thetacycle(ind,:);
thetacycle(:,13)=[1:size(thetacycle,1)];
clear LFP_HPC thetapower thetaphase 
lfptime_trial(:,5)=0;
for i=1:size(thetacycle,1)
    ind=thetacycle(i,1):thetacycle(i,2);
    lfptime_trial(ind,5)=i;
end
thetacycle(:,14)=(thetacycle(:,5)-thetacycle(:,7)).*thetacycle(:,3); % descending half dur
ind=thetacycle(:,7)>thetacycle(:,5);
thetacycle(ind,14)=(1-thetacycle(ind,7)+thetacycle(ind,5)).*thetacycle(ind,3); % descending half dur;
thetacycle(:,15)=thetacycle(:,3)-thetacycle(:,14);% ascending half dur;
thetacycle(:,16)=thetacycle(:,14)./thetacycle(:,15); % descending/ascending

propset=[3,4,7,5,8,15,14,16:20,10];
propname={'Duration','Theta Power','Peak Percentile','Trough Percentile','Theta Prom','Ascending Half Dur',...
    'Descending Half Dur','Descend/Ascend','MML Sink','MML Sink Percentile','MML Source','MML Source Percentile','Velocity'};
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG IML/GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};

cd(savedir);save('Analyze_Theta_Cycle_Broad','thetacycle',"propname",'propset',"dHPC_layer7channel3",'CSDall_trial',...
    'lfpz_trial','HPC_layer_name','reftheta','lfptime_trial','-v7.3');
disp(['Finish Session : ',num2str(in)]);

if 0
figure;
for p=1:13
    pid=propset(p);
    subplot(3,5,p);histogram(thetacycle(:,pid));
    m1=max(thetacycle(:,pid));m2=min(thetacycle(:,pid));
    xline([m1,m2],'r');title(propname{p});
end

%----------------------------------------------------------------------------------
% 1.below plot lfp of each theta cycle grouped by velocity
thetacycle_lfp=zeros(16,36,11);dv=5;
for i=1:36
    indi=ceil(reftheta(:,4)/10)==i;
    for v=1:80/dv
        indv=find(ceil(thetacycle(:,10)/dv)==v);
        lia=ismember(lfptime_trial(:,5),indv);
        ind= lia & indi;
        thetacycle_lfp(v,i,:)=nanmean(lfpz_trial(ind ,:),1);
    end
end
typenum=80/dv;typeindex=[1:typenum]*floor(256/typenum);
plotcolors=colormap(jet);
plotcolors=plotcolors(typeindex,:);
figure;
for i=1:11
    subplot(3,4,i);
    for v=1:80/dv
        a=squeeze(thetacycle_lfp(v,:,i));
        a=[a,a];
        hold on;plot([5:10:715],a,'Color',plotcolors(v,:));
    end
    xlabel('Theta Phase');ylabel('Mean LFP');title(HPC_layer_name{i});xlim([0 720]);
end
cd(savedir);figure_title='ThetaCycle_LFP_Velocity';save_current_figure(figure_title);
%----------------------------------------------------------------------------------
% 2. below plot theta cycle csd
thetacycle_csd=[];dv=10;
for i=1:36
    indi=ceil(reftheta(:,4)/10)==i;
    for v=1:80/dv
        indv=find(ceil(thetacycle(:,10)/dv)==v);
        lia=ismember(lfptime_trial(:,5),indv);
        ind= lia & indi;
        thetacycle_csd(v,i,:)=nanmean(CSDall_trial(:,ind),2);
    end
end
as=[];sinkloc=[];
for v=1:size(thetacycle_csd,1)
    as=[as,squeeze(thetacycle_csd(v,:,:))'];
    sink1=min(squeeze(thetacycle_csd(v,:,:)),[],2);
    [M,I]=min(sink1);
    sinkloc(v,:)=[M,I];
end
figure;imagesc([0 size(thetacycle_csd,1)],[],as);set(gca,'YDir','normal');
xline([1:size(thetacycle_csd,1)],'k');xlabel('Theta Phase : sorted by velocity low to high');ylabel('Channel #');
title('Mean CSD during Theta Cycles');yline(dHPC_layer7channel3(1:3,4),'b');
yline(dHPC_layer7channel3(4:6,4),'r');yline(dHPC_layer7channel3(7,4),'m');
cd(savedir);figure_title='ThetaCycle_CSD_Velocity';save_current_figure(figure_title);
%----------------------------------------------------------------------------------
% 3. below plot thetacycle prop - velocity 2d dist
figure;
for p=1:12
    pid=propset(p);
    if p~=8
        m2=max(thetacycle(:,pid));
    else
        m2=5;
    end
    m1=min(thetacycle(:,pid));dm=(m2-m1)/50;
    bins=[m1:dm:m2];
    a=histcounts2(thetacycle(:,pid),thetacycle(:,10),bins,[0:5:80]);a=a/sum(a(:));
    subplot(4,4,p);
    imagesc([0 80],[m1 m2],a);set(gca,'YDir','normal');
    xlabel('Velocity');ylabel(propname{p});
end
bins=[0:0.02:1];
ind=thetacycle(:,10)>=10;
subplot(4,4,13);
a=histcounts2(thetacycle(ind,7),thetacycle(ind,18),bins,bins);a=a/sum(a(:));
imagesc([0 1],[0 1],a);set(gca,'YDir','normal');
hold on;plot([0 1],[0 1],'r--');
xlabel('Sink Percentile');ylabel('Theta Peak Percentile');title('Run');
subplot(4,4,14);
a=histcounts2(thetacycle(ind,5),thetacycle(ind,20),bins,bins);a=a/sum(a(:));
imagesc([0 1],[0 1],a);set(gca,'YDir','normal');
hold on;plot([0 1],[0 1],'r--');
xlabel('Source Percentile');ylabel('Theta Trough Percentile');title('Run');
subplot(4,4,15);
a=histcounts2(thetacycle(ind,20),thetacycle(ind,18),bins,bins);a=a/sum(a(:));
imagesc([0 1],[0 1],a);set(gca,'YDir','normal');
hold on;plot([0 1],[0 1],'r--');
ylabel('MML Source Percentile');xlabel('MML Sink Percentile');title('Run');
subplot(4,4,16);
a=histcounts2(thetacycle(ind,19),thetacycle(ind,17),[0:0.005:0.5],[-0.5:0.005:0]);a=a/sum(a(:));
imagesc([-0.5 0],[0 0.5],a);set(gca,'YDir','normal');
hold on;plot([0 1],[0 1],'r--');
ylabel('MML Source');xlabel('MML Sink');title('Run');
cd(savedir);figure_title='ThetaCycleGroup_Dist_Velocity';save_current_figure(figure_title);
%----------------------------------------------------------------------------------
% 4. below quantify correlation btw velocity & thetacycle property by
% thetacycle group
groupN=100;
thetacyclegroup=zeros(groupN,13);
[B,I]=sort(thetacycle(:,10));
I(:,2)=[1:length(I)]/length(I);
for g=1:groupN
    indv=I(ceil(I(:,2)*groupN)==g,1);
    thetacyclegroup(g,:)=nanmean(thetacycle(indv,[10,3,4,8,15,14,16,17,19,18,20,7,5]),1);
end
propnamegroup={'Velocity','Duration','Theta Power','Theta Prom','Ascending Half Dur','Descending Half Dur',...
    'Descend/Ascend','MML Sink','MML Source','Sink Percentile','Source Percentile','Theta Peak Percentile','Theta Trough Percentile'};
figure;
for p=2:11
    subplot(3,4,p-1);
    if p<=9
        plot(thetacyclegroup(:,1),thetacyclegroup(:,p),'.');
        [r,pval]=corrcoef(thetacyclegroup(:,1),thetacyclegroup(:,p));
        xlabel('Velocity');ylabel(propnamegroup{p});
    else
        plot(thetacyclegroup(:,p),thetacyclegroup(:,p+2),'.');
        [r,pval]=corrcoef(thetacyclegroup(:,p),thetacyclegroup(:,p+2));
        xlabel(propnamegroup{p});ylabel(propnamegroup{p+2});
    end
    title(['R, P = ',num2str([r(2),pval(2)],2)]);
end
subplot(3,4,11);p=5;
plot(log2(thetacyclegroup(:,1)),thetacyclegroup(:,p),'.');
[r,pval]=corrcoef(thetacyclegroup(:,1),thetacyclegroup(:,p));
xlabel('log2(Velocity)');ylabel(propnamegroup{p});
title(['R, P = ',num2str([r(2),pval(2)],2)]);
cd(savedir);figure_title='ThetaCycleGroup_Correlation';save_current_figure(figure_title);
%--------------------------------------------------------------------------------------
% to summarize across sessions
% s1. csd across sessions
typenum=16;typeindex=[1:typenum]*floor(256/typenum);
plotcolors=colormap(jet);
plotcolors=plotcolors(typeindex,:);
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Analyze_Theta_Cycle');
    figure(1);subplot(4,4,in);
    a=[squeeze(thetacycle_csd(1,:,:))',squeeze(thetacycle_csd(4,:,:))',squeeze(thetacycle_csd(8,:,:))'];
    imagesc([0 3],[ ],a);set(gca,'YDir','normal');xline(1:2,'k');
    yline(dHPC_layer7channel3(1:3,4),'b');yline(dHPC_layer7channel3(4:6,4),'r');yline(dHPC_layer7channel3(7,4),'m');
    xticks([0.5,1.5,2.5]);xticklabels({'Pause','Low Velocity','High Velocity'});ylabel('Channel #');

    figure(2);subplot(4,4,in);js=[1,4:7];
    for i=1:5:16
        for j=1:4
            a=squeeze(thetacycle_lfp(i,:,js(j)));
            hold on;plot([5:10:715],[a,a]+(7-j)*2,'color',plotcolors(i,:));
        end
    end
    xlabel('Theta Phase');ylabel(' Mean LFP');xlim([0 720]);

    figure(3);
    for p=2:11
        subplot(3,4,p-1);hold on;
        if p<=9
            plot(thetacyclegroup(:,1),thetacyclegroup(:,p),'-','color',plotcolors(in,:));
            %[r,pval]=corrcoef(thetacyclegroup(:,1),thetacyclegroup(:,p));
            xlabel('Velocity');ylabel(propnamegroup{p});
        else
            plot(thetacyclegroup(:,p),thetacyclegroup(:,p+2),'-','color',plotcolors(in,:));
            %[r,pval]=corrcoef(thetacyclegroup(:,p),thetacyclegroup(:,p+2));
            xlabel(propnamegroup{p});ylabel(propnamegroup{p+2});
        end
        %title(['R, P = ',num2str([r(2),pval(2)],2)]);
    end
    subplot(3,4,11);p=5;
    hold on;plot(log2(thetacyclegroup(:,1)),thetacyclegroup(:,p),'-','color',plotcolors(in,:));
    %[r,pval]=corrcoef(thetacyclegroup(:,1),thetacyclegroup(:,p));
    xlabel('log2(Velocity)');ylabel(propnamegroup{p});
    %title(['R, P = ',num2str([r(2),pval(2)],2)]);
end

figure(1);cd(homedir);figure_title='ThetaCycle_CSD_Velocity_Sessions';save_current_figure(figure_title);
figure(2);cd(homedir);figure_title='ThetaCycle_LFP_Velocity_Sessions';save_current_figure(figure_title);
figure(3);cd(homedir);figure_title='ThetaCycleGroup_Correlation_Sessions';save_current_figure(figure_title);
end