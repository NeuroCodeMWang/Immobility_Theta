%function visualize_trial_decode_info_v6(in)
if 0
clear all;close all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
in=10
savedir=SessionSet16{in};
cd(savedir);
factor=50;nn=6;timewin=0.01;
binary=slanCM('binary');
bone=slanCM('bone');bone=bone(256:-1:1,:);
seismic=slanCM('seismic');

load('ThetaCycle_Decode_Info_dHPC2','deseqall');deseqallt=deseqall;
load('CSD_LFP_Visualization.mat', 'CSDall', 'dHPC_layer7channel3', 'channels');
load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','LFP_HPC','lfptime','dHPC_layer7channel3');
load('Position_Data_Maze.mat');
load('Analyze_Theta_Cycle_Broad','reftheta');
load('sorted_spike_decoding_8arm_dHPC.mat','deseqall','dedataall','maze1d','WinAdvance');
load('ThetaCycle_Decode_Info_dHPC2','thetacycle','pksets','inter_run','linear_pos');
load('InOutBound_Behavior_Analysis.mat','pos_behave_marker');
lfp7=zscore(LFP_HPC(:,1:7),0,1);
pksets(:,7)=round(pksets(:,7));
pksets(:,10)=Position_Data(pksets(:,22),1); % reward start time
pksets(:,12)=Position_Data(pksets(:,23),1); % reward end time

ruler=reftheta(:,1);template=lfptime(:,1);[outputindex,error]=match(template,ruler,0);
lfptime(:,5)=reftheta(outputindex,3);
    load('Ripple_Events_PAPER5.mat','ripples','Ripple_Amplitude_CA1');
    ind=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035;
    ripples=ripples(ind,:);
    ruler=thetacycle(:,6);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ruler=thetacycle(:,22);template=ripples(:,3);[outputindex,error1]=match(template,ruler,0);
    ruler=thetacycle(:,23);template=ripples(:,3);[outputindex,error2]=match(template,ruler,0);
    ripples(:,15)=min([abs(error),abs(error1),abs(error2)],[],2); 
    ind=abs(ripples(:,11))==0.5;
    ripples(ind,12)=inter_run(ripples(ind,12),2);
    ripples(:,18)=pksets(ripples(:,12),7);
    indrt=abs(ripples(:,15))<=0.05 ;  
    ripples(:,19)=indrt;
    ind=(ripples(:,18)<=-4 | (ripples(:,18)>0 & ripples(:,18)<=7) );
    ripples=ripples(ind,:);
    ruler=reftheta(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,22:23)=reftheta(outputindex,3:4);
    ruler=reftheta(:,1);template=deseqallt(:,1);[outputindex,error]=match(template,ruler,0);
    deseqallt(:,45)=reftheta(outputindex,3);
    thetacycle(:,79)=0; deseqallt(:,43:44)=0;% associated ripple id
    ripple_decode=nan*ones(size(ripples,1),35);ripples(:,16:17)=0;
    for i=1:size(ripples,1)
        ind= deseqallt(:,1)>=ripples(i,1)-timewin & deseqallt(:,1)<=ripples(i,2)+timewin ;
        deseqallt(ind,42+ripples(i,6))=i;
        ind= thetacycle(:,6)>=ripples(i,1)-0.05 & thetacycle(:,6)<=ripples(i,2)+0.05 ;
        thetacycle(ind,79)=i;
    end
    indtheta= thetacycle(:,79)==0 ;
    liatheta=ismember(deseqallt(:,40),thetacycle(indtheta,21));
    lia_theta_i= liatheta & deseqallt(:,40)>0 & deseqallt(:,43)==0 & deseqallt(:,44)>=0;
    lia_none_i=deseqallt(:,43)==0 & deseqallt(:,44)>=0 & deseqallt(:,40)==0 & deseqallt(:,45)<-1 ;
    deseqallt(:,46)=nan;
    deseqallt(lia_theta_i,46)=1;
    deseqallt(lia_none_i,46)=0;
    ind=deseqallt(:,43)>0;deseqallt(ind,46)=2;
    ind=deseqallt(:,44)>0;deseqallt(ind,46)=3;
    if 0
    load('ThetaCycle_Decode_Info_dHPC2','pksets','inter_run','ripple1_trial','ripple2');
    ripples=[ripple1_trial(:,1:6);ripple2(:,1:6)];
    [~,I]=sort(ripples(:,3));
    ripples=ripples(I,:);
    ruler=Position_Data(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,7:12)=[Position_Data(outputindex,2:5),pos_behave_marker(outputindex,:)];
    ripples(:,13)=vecnorm(ripples(:,7:8),2,2);
    ripples(:,14)=in;
    ind=ripples(:,10)<=10 & ripples(:,11)~=-100 & ripples(:,12)~=-100 & ripples(:,2)-ripples(:,1)>0.035;ripples=ripples(ind,:);
    end
end
if 0
cd(homedir);load('next_arm_deseq_theta_ripple_4','pksetsall','pksetsall3','rippleall'); 
load('Ripple_Events_PAPER3.mat', 'Ripple_Amplitude_CA1')
end
[6477.6 6484]
start_time=6477.5;end_time=6484;
lfpn=1;plotwin=5;
figure;linewidth=1.5;
trialset=unique(pksets(:,1));
for trial=6:length(trialset)
    trialid=trialset(trial);
    trialstarttime=min(Position_Data(Position_Data(:,4)==trialid,1));
    trialendtime=max(Position_Data(Position_Data(:,4)==trialid,1));
    trialdur=trialendtime-trialstarttime;
    ind=deseqall(:,14)==trialid;
    deseq=deseqall(ind,:);
    dedata=dedataall{trial};
for t=1:ceil(trialdur/plotwin)
    start_time=trialstarttime+(t-1)*plotwin;
    end_time=trialstarttime+t*plotwin;
    indpos= Position_Data(:,1)>=start_time & Position_Data(:,1)<=end_time ;
    indlfp= lfptime(:,1)>=start_time & lfptime(:,1)<=end_time ;
    inddecode= deseq(:,1)>=start_time & deseq(:,1)<=end_time ;
    indtheta=thetacycle(:,6)>=start_time & thetacycle(:,6)<=end_time;
    indreftheta=indlfp & lfptime(:,5)>=0;
    inddecode2= deseqallt(:,1)>=start_time & deseqallt(:,1)<=end_time & deseqallt(:,14)==trialid;
    deseqt=deseqallt(inddecode2,:);
    if sum(indpos)>5
    posid=find(Position_Data(:,1)<start_time,1,"last");
    pkstart=find(pksets(:,5)>=start_time & pksets(:,5)<=end_time);
    pkend=find(pksets(:,6)>=start_time & pksets(:,6)<=end_time);
    pkstartreward=find(pksets(:,10)>=start_time & pksets(:,10)<=end_time);
    pkendreward=find(pksets(:,12)>=start_time & pksets(:,12)<=end_time);
    pkid=find(pksets(:,5)<=end_time,1,'last');
    pk=pksets(pksets(:,1)==trialid,:);
    armseq=pk(:,2);

    % below for lfp data
    subplot(4,11,[1:10,12:21]);
    imagesc(lfptime(indlfp,1),[channels(1,3) channels(end,3)],imgaussfilt(CSDall(:,indlfp),2));set(gca,'YDir','normal');
    for j=1:7
        if j<=3
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
        elseif j==4
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
            hold on;plot(lfptime(indreftheta,1),lfp7(indreftheta,j)*factor+dHPC_layer7channel3(j,2),'r.','MarkerSize',10);
        else
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
        end
    end
    hold on;plot(lfptime(indlfp,1),Ripple_Amplitude_CA1(indlfp,2)*factor+dHPC_layer7channel3(1,2)-50,'m');
    yline(3*factor+dHPC_layer7channel3(1,2)-50,'m--');
    yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
    if sum(indtheta)
        xline(thetacycle(indtheta,6),'k--');
    end
    if 1
    rip=ripples(ripples(:,3)>=min(lfptime(indlfp,1)) & ripples(:,3)<=max(lfptime(indlfp,1)), :);
    for i=1:size(rip,1)
        indt2= lfptime(:,1)>=rip(i,1) & lfptime(:,1)<=rip(i,2);
        if rip(i,6)==1
            j=1;hold on;plot(lfptime(indt2,1),lfp7(indt2,j)*factor+dHPC_layer7channel3(j,2),'r');
        elseif rip(i,6)==2
            j=7;hold on;plot(lfptime(indt2,1),lfp7(indt2,j)*factor+dHPC_layer7channel3(j,2),'b');
        end
    end
    ind=rip(:,6)==1;
    if sum(ind)>0
        xline(rip(ind,3),'r','LineWidth',linewidth);
    end
    ind=rip(:,6)==2;
    if sum(ind)>0
        xline(rip(ind,3),'b','LineWidth',linewidth);
    end
    if 0
    ds=dspeaks_d4( dspeaks_d4(:,9)>=min(lfptime(indlfp,1)) & dspeaks_d4(:,9)<=max(lfptime(indlfp,1)), : );
    ind=ds(:,6)==1;
    if sum(ind)>0
        xline(ds(ind,9),'b','LineWidth',linewidth);
    end
    ind=ds(:,6)==2;
    if sum(ind)>0
        xline(ds(ind,9),'r','LineWidth',linewidth);
    end
    ind=ds(:,6)==3;
    if sum(ind)>0
        xline(ds(ind,9),'g','LineWidth',linewidth);
    end
    
    if ~isempty(pkstart)
        xline(pksets(pkstart,5),'m');
    end
    if ~isempty(pkend)
        xline(pksets(pkend,6),'m--');
    end
    end
    if ~isempty(pkstartreward)
        xline(pksets(pkstartreward,10),'y');
    end
    if ~isempty(pkendreward)
        xline(pksets(pkendreward,12),'y--');
    end
    end
    ylabel('Depth on Probe');%xlabel('Time (s)');
    xlim([start_time end_time]);title('CSD + LFP');

    if 1
    % below for decode data
    deseq1=deseq(inddecode,[1,4,5]);
    dedata1=dedata(:,inddecode)';
    [linearized_maze,linear_pos]=linearize_8arm(dedata1,maze1d,deseq1);
    linear_maze=nansum(linearized_maze(:,1:35,:),3)';linear_maze(16:end,:)=0;
    for j=1:8
        linear_maze=[linear_maze;linearized_maze(:,16:end,j)'];
    end
    a=-deseq1(1,1):WinAdvance:-start_time;a=-a(end:-1:1);
    deseqfull=[a(1):WinAdvance:end_time]';
    ruler=deseq1(:,1);template=deseqfull(:,1);[outputindex,error]=match(template,ruler,0);
    lia=abs(error)<WinAdvance/10;
    if sum(lia)~=size(deseq1,1)
        disp('Error deseqfull size!');
    else
        deseqfull(:,2)=lia;
    end
    linear_maze_full=nan*ones(315,size(deseqfull,1));
    linear_maze_full(:,lia)=linear_maze;
    linear_pos_full=nan*ones(size(deseqfull,1),6);
    linear_pos_full(lia,:)=linear_pos;
    ratpos_full=linear_pos_full(:,5)+(linear_pos_full(:,2)-15)/35-0.5;
    ind=linear_pos_full(:,5)==0;
    ratpos_full(ind)=linear_pos_full(ind)/35-0.5;
    end

    subplot(4,11,[23:32,34:43]);
    imagesc([deseqfull(1,1) deseqfull(end,1)],[-0.5 8.5],linear_maze_full);set(gca,'YDir','normal'); 
    hold on;plot(deseqfull(:,1),ratpos_full,'c.');ylabel('Arm #');yline(0.5:7.5,'c');%colormap("hot");
    hold on;plot(deseqt(:,1),deseqt(:,46),'mo');
    yyaxis right
    hold on;plot(Position_Data(indpos,1),Position_Data(indpos,5));yline(10,'r--');ylabel('Velocity (cm/s)');ylim([0 100]);
    %hold on;plot(Position_Data(indpos,1),Position_Data(indpos,5),'m-');yline(10,'m--');ylim([0 100]);
    title(['Trial # , Visit #, Arm # : ',num2str([trialid,pksets(pkid,[7,2])])]);
    xlim([start_time end_time]);xlabel('Time (s)');
    if sum(indtheta)
        xline(thetacycle(indtheta,6),'c--');
    end
    if ~isempty(pkstartreward)
        xline(pksets(pkstartreward,10),'y');
    end
    if ~isempty(pkendreward)
        xline(pksets(pkendreward,12),'y--');
    end
    ind=rip(:,6)==1;
    if sum(ind)>0
        xline(rip(ind,3),'r','LineWidth',linewidth);
    end
    ind=rip(:,6)==2;
    if sum(ind)>0
        xline(rip(ind,3),'b','LineWidth',linewidth);
    end
    if 0
    if ~isempty(pkstart)
        xline(pksets(pkstart,5),'m');
    end
    if ~isempty(pkend)
        xline(pksets(pkend,6),'m--');
    end    
    ind=ds(:,6)==1;
    if sum(ind)>0
        xline(ds(ind,9),'b','LineWidth',linewidth);
    end
    ind=ds(:,6)==2;
    if sum(ind)>0
        xline(ds(ind,9),'r','LineWidth',linewidth);
    end
    ind=ds(:,6)==3;
    if sum(ind)>0
        xline(ds(ind,9),'g','LineWidth',linewidth);
    end
    end

    if 1
    subplot(4,11,[33,44]);
    plot(Position_Data(indpos,2),Position_Data(indpos,3),'.');
    xlim([-110 110]);ylim([-110 110]);title(num2str(armseq'));
    end

    set(findall(gcf, 'Type', 'Line'),'LineWidth',1.5,'MarkerSize', 6);set(gcf,'color','w');set(findobj(gcf,'type','axes'),'FontSize',10);
    end
    pause;clf;
end
end
if 0
    clear all;close all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
in=10
savedir=SessionSet16{in};
cd(savedir);
factor=50;nn=6;timewin=0.01;
binary=slanCM('binary');
bone=slanCM('bone');bone=bone(256:-1:1,:);
seismic=slanCM('seismic');
ifsavefig=0;
load('ThetaCycle_Decode_Info_dHPC2','deseqall');deseqallt=deseqall;
load('CSD_LFP_Visualization.mat', 'CSDall', 'dHPC_layer7channel3', 'channels');
load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','LFP_HPC','lfptime','dHPC_layer7channel3');
load('Position_Data_Maze.mat');
load('Analyze_Theta_Cycle_Broad','reftheta');
load('sorted_spike_decoding_8arm_dHPC.mat','deseqall','dedataall','maze1d','WinAdvance');
load('ThetaCycle_Decode_Info_dHPC2','thetacycle','pksets','inter_run','linear_pos');
load('InOutBound_Behavior_Analysis.mat','pos_behave_marker');
lfp7=zscore(LFP_HPC(:,1:7),0,1);
pksets(:,7)=round(pksets(:,7));
pksets(:,10)=Position_Data(pksets(:,22),1); % reward start time
pksets(:,12)=Position_Data(pksets(:,23),1); % reward end time

ruler=reftheta(:,1);template=lfptime(:,1);[outputindex,error]=match(template,ruler,0);
lfptime(:,5)=reftheta(outputindex,3);
  load('Ripple_Events_PAPER5.mat','ripples','Ripple_Amplitude_CA1');
    ind=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035;
    ripples=ripples(ind,:);
    ruler=thetacycle(:,6);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ruler=thetacycle(:,22);template=ripples(:,3);[outputindex,error1]=match(template,ruler,0);
    ruler=thetacycle(:,23);template=ripples(:,3);[outputindex,error2]=match(template,ruler,0);
    ripples(:,15)=min([abs(error),abs(error1),abs(error2)],[],2); 
    ind=abs(ripples(:,11))==0.5;
    ripples(ind,12)=inter_run(ripples(ind,12),2);
    ripples(:,18)=pksets(ripples(:,12),7);
    indrt=abs(ripples(:,15))<=0.05 ;  
    ripples(:,19)=indrt;
    ind=(ripples(:,18)<=-4 | (ripples(:,18)>0 & ripples(:,18)<=7) );
    ripples=ripples(ind,:);
    ruler=reftheta(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,22:23)=reftheta(outputindex,3:4);
    ruler=reftheta(:,1);template=deseqallt(:,1);[outputindex,error]=match(template,ruler,0);
    deseqallt(:,45)=reftheta(outputindex,3);
    thetacycle(:,79)=0; deseqallt(:,43:44)=0;% associated ripple id
    ripple_decode=nan*ones(size(ripples,1),35);ripples(:,16:17)=0;
    for i=1:size(ripples,1)
        ind= deseqallt(:,1)>=ripples(i,1)-timewin & deseqallt(:,1)<=ripples(i,2)+timewin ;
        deseqallt(ind,42+ripples(i,6))=i;
        ind= thetacycle(:,6)>=ripples(i,1)-0.05 & thetacycle(:,6)<=ripples(i,2)+0.05 ;
        thetacycle(ind,79)=i;
    end
    indtheta= thetacycle(:,79)==0 ;
    liatheta=ismember(deseqallt(:,40),thetacycle(indtheta,21));
    lia_theta_i= liatheta & deseqallt(:,40)>0 & deseqallt(:,43)==0 & deseqallt(:,44)>=0;
    lia_none_i=deseqallt(:,43)==0 & deseqallt(:,44)>=0 & deseqallt(:,40)==0 & deseqallt(:,45)<-1 ;
    deseqallt(:,46)=nan;
    deseqallt(lia_theta_i,46)=1;
    deseqallt(lia_none_i,46)=0;
    ind=deseqallt(:,43)>0;deseqallt(ind,46)=2;
    ind=deseqallt(:,44)>0;deseqallt(ind,46)=3;

    nn=6;factor=50;
    puor=slanCM('PuOr');puor=puor(256:-1:1,:); 
    start_time=6477.7;end_time=6484;trialid=20;trial=20;
    trialstarttime=min(Position_Data(Position_Data(:,4)==trialid,1));
    trialendtime=max(Position_Data(Position_Data(:,4)==trialid,1));
    trialdur=trialendtime-trialstarttime;
    ind=deseqall(:,14)==trialid;
    deseq=deseqall(ind,:);
    dedata=dedataall{trial};
        
    indpos= Position_Data(:,1)>=start_time & Position_Data(:,1)<=end_time ;
    indlfp= lfptime(:,1)>=start_time & lfptime(:,1)<=end_time ;
    inddecode= deseq(:,1)>=start_time & deseq(:,1)<=end_time ;
    posid=find(Position_Data(:,1)<start_time,1,"last");
    pkid=find(pksets(:,5)<=end_time,1,'last');
    pk=pksets(pksets(:,1)==trialid,:);
    armseq=pk(:,2);
    inddecode2= deseqallt(:,1)>=start_time & deseqallt(:,1)<=end_time & deseqallt(:,14)==trialid;
    deseqt=deseqallt(inddecode2,:);

    % below for lfp data
    figure(201);a1=subplot(5,11,[1:10,12:21]);
    imagesc(lfptime(indlfp,1),[channels(1,3) channels(end,3)],imgaussfilt(CSDall(:,indlfp),2));set(gca,'YDir','normal');
    ap=imgaussfilt(CSDall(:,indlfp),2);
    caxis(max(abs(ap(:)))*[-1 1]);colormap(a1,puor);%ylim([4600 7100]);
    for j=1:7
        if j<=3
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
        elseif j==4
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'r');
            %hold on;plot(lfptime(indreftheta,1),lfp7(indreftheta,j)*factor+dHPC_layer7channel3(j,2),'r.','MarkerSize',10);
        else
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
        end
    end
    %hold on;plot(lfptime(indlfp,1),Ripple_Amplitude_CA1(indlfp,2)*factor+dHPC_layer7channel3(1,2)-50,'m');
    %yline(3*factor+dHPC_layer7channel3(1,2)-50,'m--');
    yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
    ylabel('Depth on Probe');%xlabel('Time (s)');
    xlim([start_time end_time]);title('CSD + LFP');
    FIG_INDEX=['sup_denext_method_example_1'];save_fig(FIG_INDEX,ifsavefig,puor,1);

    if 1
    % below for decode data
    deseq1=deseq(inddecode,[1,4,5]);
    dedata1=dedata(:,inddecode)';
    [linearized_maze,linear_pos]=linearize_8arm(dedata1,maze1d,deseq1);
    linear_maze=nansum(linearized_maze(:,1:35,:),3)';linear_maze(16:end,:)=0;
    for j=1:8
        linear_maze=[linear_maze;linearized_maze(:,16:end,j)'];
    end
    a=-deseq1(1,1):WinAdvance:-start_time;a=-a(end:-1:1);
    deseqfull=[a(1):WinAdvance:end_time]';
    ruler=deseq1(:,1);template=deseqfull(:,1);[outputindex,error]=match(template,ruler,0);
    lia=abs(error)<WinAdvance/10;
    if sum(lia)~=size(deseq1,1)
        disp('Error deseqfull size!');
    else
        deseqfull(:,2)=lia;
    end
    linear_maze_full=nan*ones(315,size(deseqfull,1));
    linear_maze_full(:,lia)=linear_maze;
    linear_pos_full=nan*ones(size(deseqfull,1),6);
    linear_pos_full(lia,:)=linear_pos;
    ratpos_full=linear_pos_full(:,5)+(linear_pos_full(:,2)-15)/35-0.5;
    ind=linear_pos_full(:,5)==0;
    ratpos_full(ind)=linear_pos_full(ind)/35-0.5;
    end

    figure(201);a1=subplot(5,11,[23:32,34:43]);
    imagesc([deseqfull(1,1) deseqfull(end,1)],[-0.5 8.5],linear_maze_full);set(gca,'YDir','normal');caxis([0 0.018]);
    hold on;plot(deseqfull(:,1),ratpos_full,'c.');ylabel('Arm #');yline(0.5:7.5,'Color',[1 1 1],'linewidth',2);colormap(a1,'hot');%colormap("hot");
    hold on;plot(Position_Data(indpos,1),Position_Data(indpos,5)*9/100-0.5,'m');%yline(10,'m--');%ylabel('Velocity (cm/s)');ylim([0 100]);
    title(['Trial # , Visit #, Arm # : ',num2str([trialid,pksets(pkid,[7,2])])]);
    xlim([start_time end_time]);xlabel('Time (s)');
    FIG_INDEX=['sup_denext_method_example_2'];save_fig(FIG_INDEX,ifsavefig,'hot');
    figure(201);subplot(5,11,[45:54]);
    hold on;plot(deseqt(:,1),deseqt(:,46),'k.'); xlim([start_time end_time]);
    FIG_INDEX=['sup_denext_method_example_3'];save_fig(FIG_INDEX,ifsavefig,'hot');
    figure_title=['Sup_Figure12_decode_example'];save_current_figure(figure_title);
    set(findall(gcf, 'Type', 'Line'),'LineWidth',1.5,'MarkerSize', 6);set(gcf,'color','w');set(findobj(gcf,'type','axes'),'FontSize',10);
end
cd(homedir);load('next_arm_deseq_theta_ripple_5.mat');
%id=find(pksetsall3(:,45)==in&pksetsall3(:,1)==trialid & pksetsall3(:,2)==pksets(pkid,2));
id=2996;
a=[squeeze(nanmean(deseq_thetaall(id,1,1:3,1:9),3))';squeeze(deseq_rippleall(id,1,1,1:9))';squeeze(nanmean(deseq_noneall(id,1,1:3,1:9),3))';];
cs=[1,0,0;0.1,0.1,0.8;.5 .5 .5];
figure;
for i=1:3
    hold on;plot(a(i,:),'Color',cs(i,:));
end
xlim([0.5 9.5]);
FIG_INDEX=['sup_denext_method_example_4'];save_fig(FIG_INDEX,ifsavefig,[],3);
   
FIG_INDEX=['sup_denext_method_example_4'];save_fig(FIG_INDEX,ifsavefig,[],3);