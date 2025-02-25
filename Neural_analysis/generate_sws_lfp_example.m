function generate_sws_lfp_example(ifsavefig,in2)
if nargin<1
    ifsavefig=0;
end
close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);
figdir='X:\Mengni\Data_Analysis\Paper_Figures\figures';
load('SessionSet16');
load('SlowWaveSleep_Example',"timeset");
puor=slanCM('PuOr');puor=puor(256:-1:1,:);
for in=1:16
    if timeset(in,1)>0 & ( nargin<2 | (nargin==2 & in==in2) )
    savedir=SessionSet16{in};
    cd(savedir);
    %load('theta_delta_LFP_HPC.mat', 'lfptime');
    %load('HPC_Layer_Channel_LFP.mat', 'LFP_HPC');
    load('CSD_LFP_Visualization.mat', 'CSDall', 'channels');
    load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','LFP_HPC','lfptime','dHPC_layer7channel3');
    load('Position_Data_Raw.mat');
    Position_Data_Raw(:,5)=vecnorm(Position_Data_Raw(:,2:3),2,2);
    Position_Data_Raw(1:end-1,6)=diff(Position_Data_Raw(:,5))/0.0335;

    dHPC_layer7channel3(:,4)=0;
    for i=1:7
        dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2));
    end
    LFP_HPC=zscore(LFP_HPC,0,1);
    indsleep=lfptime(:,2)<0;
    lfptime_sleep=lfptime(indsleep,:);
    CSDall_sleep=CSDall(:,indsleep);
    lfp7=LFP_HPC(indsleep,1:7);

    factor=20;
    if in==10
        timeset(in,:)=[949.7 956.2];
    end
    for i=1%:size(timeset,1)
        indlfp= lfptime_sleep(:,1)>=timeset(in,1) & lfptime_sleep(:,1)<=timeset(in,2);
        indpos= Position_Data_Raw(:,1)>=timeset(in,1) & Position_Data_Raw(:,1)<=timeset(in,2);
        figure(in);subplot(5,1,[1,2]);        
        for j=1:7
            if j==1
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
            elseif j==4
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'r');
            else
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
            end
        end
        xlim([min(lfptime_sleep(indlfp,1)),max(lfptime_sleep(indlfp,1))]);title(num2str(in));
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        FIG_INDEX=['SupFig4_sws_lfp_example1_',num2str(in)];save_fig_epsc(FIG_INDEX,ifsavefig,[],1.5);
        
        figure(in);subplot(5,1,[3,4]);
        a=imgaussfilt(CSDall_sleep(:,indlfp),1);
        imagesc(lfptime_sleep(indlfp,1),[channels(1,3) channels(end,3)],a);set(gca,'YDir','normal');
        caxis(0.15*[-1 1]);colorbar;caxis(0.7*max(abs(a(:)))*[-1 1]);
        xlabel('Time in Slow Wave Sleep (s)');xlim([min(lfptime_sleep(indlfp,1)),max(lfptime_sleep(indlfp,1))]);colormap(puor);
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        FIG_INDEX=['SupFig4_sws_lfp_example2_',num2str(in)];save_fig_epsc(FIG_INDEX,ifsavefig,puor);
        
        figure(in);subplot(5,1,5);
        hold on;plot(Position_Data_Raw(indpos,1),Position_Data_Raw(indpos,6),'m-');ylabel('Velocity (cm/s)');
        %hold on;plot(posexample(:,1),vecnorm(posexample(:,2:3),2,2),'k-');ylim([0 100]);%yline(4,'r--');
        xlim([min(Position_Data_Raw(indpos,1)),max(Position_Data_Raw(indpos,1))]);yticks([]);
        xlabel('Time (s)');%legend('Velocity','Distance to Center');
        FIG_INDEX=['SupFig4_sws_lfp_example3_',num2str(in)];save_fig_epsc(FIG_INDEX,ifsavefig,[],5);
        cd(figdir);figure_title=['SupFigure4_sws_example',num2str(in)];save_current_figure(figure_title);
    end
    set(findall(gcf, 'Type', 'Line'),'LineWidth',2,'MarkerSize', 6);set(gcf,'color','w');set(findobj(gcf,'type','axes'),'FontSize',15);
    end
end
