%function generate_swr_theta_lfp_example(ifsavefig)
%if nargin<1
    ifsavefig=0;
%end
close all;
in=1;timeset=[3666.3 3670.3;3656 3661];
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
figdir='X:\Mengni\Data_Analysis\Paper_Figures\figures';
cd(figdir); 
if exist(['swr_theta_lfp_example_',num2str(in),',mat'])
    load(['swr_theta_lfp_example_',num2str(in)]);
else
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_delta_LFP_HPC.mat', 'lfptime');
    load('HPC_Layer_Channel_LFP.mat', 'LFP_HPC');
    load('CSD_LFP_Visualization.mat', 'CSDall', 'dHPC_layer7channel3', 'channels');
    load('Position_Data_Maze.mat');
    load('Analyze_Theta_Cycle_Broad','reftheta');
    dHPC_layer7channel3(:,4)=0;
    for i=1:7
        dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2));
    end
    lfp7=zscore(LFP_HPC(:,1:7),0,1);
    load('sorted_spike_decoding_8arm_dHPC.mat','deseqall','dedataall','maze1d','WinAdvance');
    load('ThetaCycle_Decode_Info_dHPC2','thetacycle','pksets','inter_run','linear_pos');
    load('InOutBound_Behavior_Analysis.mat','pos_behave_marker');
    pksets(:,7)=round(pksets(:,7));
    pksets(:,10)=Position_Data(pksets(:,22),1); % reward start time
    pksets(:,12)=Position_Data(pksets(:,23),1); % reward end time
    
    load('Ripple_Events_PAPER3.mat','ripples');
    ind=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035;
    ripples=ripples(ind,:);
    thetacycle(:,1)=reftheta(thetacycle(:,1),1);
    thetacycle(:,2)=reftheta(thetacycle(:,2),1);
    ruler=thetacycle(:,6);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ruler=thetacycle(:,1);template=ripples(:,3);[outputindex,error1]=match(template,ruler,0);
    ruler=thetacycle(:,2);template=ripples(:,3);[outputindex,error2]=match(template,ruler,0);
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

    thetacycle(:,79)=0; deseqall(:,43:44)=0;% associated ripple id
    ripples(:,16:17)=0;lfptime(:,5:6)=0;
    for i=1:size(ripples,1)
        ind= deseqall(:,1)>=ripples(i,1) & deseqall(:,1)<=ripples(i,2) ;
        deseqall(ind,42+ripples(i,6))=i;
        ind= lfptime(:,1)>=ripples(i,1) & lfptime(:,1)<=ripples(i,2) ;
        lfptime(ind,4+ripples(i,6))=i;
        ripples(i,16)=sum(ind);
        ind= thetacycle(:,6)>=ripples(i,1)-0.05 & thetacycle(:,6)<=ripples(i,2)+0.05 ;
        thetacycle(ind,79)=i;
    end

    nn=6;factor=60;
    puor=slanCM('PuOr');puor=puor(256:-1:1,:); 
    trialset=unique(deseqall(:,14));

figure(in);
[~,I2]=sort(ripples(:,22),'descend');rtype=1;
for tt=1:size(ripples,1)
    t=I2(tt);
    if ripples(t,6)==rtype
    start_time=ripples(t,1)-1;
    end_time=ripples(t,2)+1;
    indpos= Position_Data(:,1)>=start_time & Position_Data(:,1)<=end_time ;
    indlfp= lfptime(:,1)>=start_time & lfptime(:,1)<=end_time ;
    indlfp1=indlfp & lfptime(:,5)>0;
    indlfp2=indlfp & lfptime(:,6)>0;
    indtheta= thetacycle(:,22)>=start_time & thetacycle(:,23)<=end_time ;
    indrip1= find(ripples(:,3)>=start_time & ripples(:,3)<=end_time & ripples(:,6)==1);
    indrip2= find(ripples(:,3)>=start_time & ripples(:,3)<=end_time & ripples(:,6)==2);
    trialid=unique(Position_Data(indpos,4));
    ind=deseqall(:,14)==trialid;
    deseq=deseqall(ind,:);
    trial=find(trialset==trialid);
    dedata=dedataall{trial};
    inddecode= deseq(:,1)>=start_time & deseq(:,1)<=end_time ;
    pk=pksets(pksets(:,1)==trialid,:);
    armseq=pk(:,2);
    pkid=find(pksets(:,5)<=end_time,1,'last');

    % below for lfp data
    figure(in);subplot(8,11,[1:10,12:21]);
    for j=1:7
        if j==4
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'r');
        else
            hold on;plot(lfptime(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'Color',[1 1 1]*0);
        end
    end
    hold on;plot(lfptime(indlfp1,1),lfp7(indlfp1,1)*factor+dHPC_layer7channel3(1,2),'r');
    hold on;plot(lfptime(indlfp2,1),lfp7(indlfp2,7)*factor+dHPC_layer7channel3(7,2),'b');
    yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
    ylabel('Depth on Probe');xlim([start_time end_time]);title(num2str([t,ripples(t,[6,15,22])]));
    if ~isempty(indrip1)
        xline(ripples(indrip1,3),'r');
    end
    if ~isempty(indrip2)
        xline(ripples(indrip2,3),'b');
    end
    xline(thetacycle(indtheta,22),'c--');xline(thetacycle(indtheta,23),'c--');
    %FIG_INDEX=['sup_decode_example1_',num2str(in)];save_fig(FIG_INDEX,ifsavefig,[],1.5);

    % below for lfp csd data
    figure(in);a1=subplot(8,11,[1:10,12:21]+22);
    a=imgaussfilt(CSDall(:,indlfp),1);
    imagesc(lfptime(indlfp,1),[channels(1,3) channels(end,3)],a);set(gca,'YDir','normal');colormap(a1,puor);
    yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
    ylabel('Depth on Probe');xlim([start_time end_time]);xticks([]);%colorbar;
    %caxis(0.9*max(abs(a(:)))*[-1 1]);
    if ~isempty(indrip1)
        xline(ripples(indrip1,3),'r');
    end
    if ~isempty(indrip2)
        xline(ripples(indrip2,3),'b');
    end
    xline(thetacycle(indtheta,22),'c--');xline(thetacycle(indtheta,23),'c--');
    %FIG_INDEX=['sup_decode_example2_',num2str(in)];save_fig(FIG_INDEX,ifsavefig,puor);

    % below for decode data
    deseq1=deseq(inddecode,[1,4,5]);
    dedata1=dedata(:,inddecode)';
    [linearized_maze,linear_pos]=linearize_8arm(dedata1,maze1d,deseq1);
    ind=linear_pos(:,2)>50;linear_pos(ind,2)=50;
    linear_maze=nanmean(linearized_maze(:,1:35,:),3)';linear_maze(16:end,:)=0;
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

    figure(in);a1=subplot(8,11,[23:32,34:43,45:54,56:65]+22);pn=6;
    imagesc([deseqfull(1,1) deseqfull(end,1)],[-0.5 8.5],linear_maze_full);set(gca,'YDir','normal'); 
    hold on;plot(deseqfull(1:pn:end,1),ratpos_full(1:pn:end),'c.');ylabel('Arm #');yline(0.5:7.5,'w','LineWidth',2);colormap(a1,"hot");
    hold on;plot(Position_Data(indpos,1),Position_Data(indpos,5)*9/100-0.5,'m');%colorbar;%yline(4*9/100-0.5,'m--');
    title(['Trial # , Visit #, Arm # : ',num2str([trialid,pksets(pkid,[7,2])])]);
    xlim([start_time end_time]);xlabel('Time (s)');
    if ~isempty(indrip1)
        xline(ripples(indrip1,3),'r');
    end
    if ~isempty(indrip2)
        xline(ripples(indrip2,3),'b');
    end
    xline(thetacycle(indtheta,22),'c--');xline(thetacycle(indtheta,23),'c--');
    %FIG_INDEX=['sup_decode_example3_',num2str(in)];save_fig(FIG_INDEX,ifsavefig,'hot');

    if 1
    figure(in);subplot(8,11,[2:4]*11+44);
    plot(Position_Data(indpos,2),Position_Data(indpos,3),'k.');
    xlim([-110 110]);ylim([-110 110]);title(num2str(armseq'));
    %FIG_INDEX=['sup_decode_example4_',num2str(in)];save_fig(FIG_INDEX,ifsavefig);
    end
    cd(figdir);figure_title=['fig4_decode_example',num2str(in)];save_current_figure(figure_title);
    set(findall(gcf, 'Type', 'Line'),'LineWidth',0.6,'MarkerSize', 6);set(gcf,'color','w');set(findobj(gcf,'type','axes'),'FontSize',15);
    pause;clf;
    end
end
lfpexample=[lfptime(indlfp,1),lfp7(indlfp,1:7)];
csdexample=CSDall(:,indlfp);
posexample=Position_Data(indpos,:);
thetacycle1=thetacycle(indtheta,:);
cd(figdir);save(['decode_example_',num2str(in)],'lfpexample','csdexample','posexample','deseqfull','linear_maze_full','ratpos_full', ...
    'dHPC_layer7channel3','channels','armseq','trialid','pkid','pksets','timeset','thetacycle1');

    for i=1:size(timeset,1)
        indlfp= lfptime(:,1)>=timeset(i,1) & lfptime(:,1)<=timeset(i,2);
        lfpexample=[lfptime(indlfp,1),lfp7(indlfp,:)];
        csdexample=CSDall(:,indlfp);
        indpos= Position_Data(:,1)>=timeset(i,1) & Position_Data(:,1)<=timeset(i,2) ;
        posexample=Position_Data(indpos,:);
        lfpexampleset{i}=lfpexample;
        csdexampleset{i}=csdexample;
        posexampleset{i}=posexample;
    end
    cd(figdir);save(['immobile_theta_lfp_example_',num2str(in)],'timeset',"lfpexampleset","csdexampleset",'posexampleset','channels','dHPC_layer7channel3','-v7.3');
end
puor=slanCM('PuOr');puor=puor(256:-1:1,:);factor=40;
    for i=1%:size(timeset,1)
        lfpexample=lfpexampleset{i};
        csdexample=csdexampleset{i};
        posexample=posexampleset{i};
        % remove one point of outlier in dCA1
        indout=find(zscore(lfpexample(:,2))>5);
        lfpexample(indout,2)=nanmean(lfpexample([indout-1,indout+1],2),1);
        figure(i);subplot(5,1,[1,2]);        
        for j=1:7
            if j==4
                hold on;plot(lfpexample(:,1),lfpexample(:,1+j)*factor+dHPC_layer7channel3(j,2),'r');
            elseif j==1
                hold on;plot(lfpexample(:,1),lfpexample(:,1+j)*factor+dHPC_layer7channel3(j,2),'k');
            else
                hold on;plot(lfpexample(:,1),lfpexample(:,1+j)*factor+dHPC_layer7channel3(j,2),'k');
            end
        end
        xlim([min(lfpexample(:,1)),max(lfpexample(:,1))]);title(num2str(in));
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        FIG_INDEX=['immobile_theta_lfp_example10_',num2str(i)];save_fig(FIG_INDEX,ifsavefig);
        
        figure(i);subplot(5,1,[3 4]);
        a=imgaussfilt(csdexample,1);
        imagesc(lfpexample(:,1),[channels(1,3) channels(end,3)],a);set(gca,'YDir','normal');
        caxis(0.10*[-1 1]);colorbar;caxis(0.8*max(abs(a(:)))*[-1 1]);
        xlim([min(lfpexample(:,1)),max(lfpexample(:,1))]);colormap(puor);
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        FIG_INDEX=['immobile_theta_csd_example10_',num2str(i)];save_fig(FIG_INDEX,ifsavefig,puor);

        figure(i);subplot(5,1,5);
        hold on;plot(posexample(:,1),posexample(:,5),'m-');ylabel('Velocity (cm/s)');
        hold on;plot(posexample(:,1),vecnorm(posexample(:,2:3),2,2),'k-');ylim([0 100]);%yline(4,'r--');
        xlim([min(lfpexample(:,1)),max(lfpexample(:,1))]);xlabel('Time (s)');legend('Velocity','Distance to Center');
        FIG_INDEX=['immobile_theta_vel_example10_',num2str(i)];save_fig(FIG_INDEX,ifsavefig,[],3);
        cd(figdir);figure_title='fig1_theta_example';save_current_figure(figure_title);
        set(findall(gcf, 'Type', 'Line'),'LineWidth',2,'MarkerSize', 6);set(gcf,'color','w');set(findobj(gcf,'type','axes'),'FontSize',15);
    end

