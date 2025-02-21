function generate_SupFigure_5_Jan25(ifsavedata)
if nargin<1
    ifsavedata=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
puor=slanCM('PuOr');puor=puor(256:-1:1,:);
plotcolors=slanCM('dense',16);
cd(datadir);load('Figure2_theta_quantification','thetacycle_csdset','thetacycle_lfp','HPC_layer_name');
% below sup figure 5
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Analyze_Theta_Cycle_Broad','dHPC_layer7channel3');
    thetacycle_csd=thetacycle_csdset{in};
    figure(1);subplot(4,4,in);
    a=[squeeze(thetacycle_csd(1,:,:))',squeeze(thetacycle_csd(2,:,:))',squeeze(thetacycle_csd(3,:,:))',squeeze(thetacycle_csd(4,:,:))'];
    b=[squeeze(thetacycle_lfp(1,:,1:7,in))',squeeze(thetacycle_lfp(2,:,1:7,in))',squeeze(thetacycle_lfp(3,:,1:7,in))',squeeze(thetacycle_lfp(4,:,1:7,in))'];
    imagesc([0 4*36],[ ],a);set(gca,'YDir','normal');xline([1:4]*36,'k');colormap(puor);
    for i=1:7
        hold on;plot([1:36*4],b(i,:)*2+dHPC_layer7channel3(i,4),'k');
    end
    hold on;plot([0 0],[0 3]*2+dHPC_layer7channel3(2,4),'r');
    yticks(dHPC_layer7channel3(end:-1:1,4));yticklabels(HPC_layer_name(7:-1:1));ylabel('Channel on Probe');%set(gca,'TickLength',[0 0]);
    xticks([0.5,1.5,2.5,3.5]*36);xticklabels({'0-1','1-5','5-10','> 10'});xlabel('Theta cycles grouped by rat velocity (cm/s)');
    caxis(0.8*max(abs(a(:)))*[-1 1]);colorbar;
    
    if ifsavedata
        rowname={'dCA1 pyr_LFP','dCA1 st rad_LFP','dCA1 slm_LFP','DG OML_LFP','DG MML_LFP','DG GCL_LFP','dCA3_LFP'};
        for i=1:size(thetacycle_csd,3)
            rowname{7+i}=['CSD_',num2str(i)];
        end
        figdata=[b;a];
        T = array2table(figdata,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Sup Figure 5_',num2str(in)],'WriteRowNames',true);
    end
end
figure(1);figure_title='Sup_Figure_5_Jan2025';save_current_figure(figure_title);