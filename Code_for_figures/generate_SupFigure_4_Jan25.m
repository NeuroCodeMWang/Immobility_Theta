function generate_SupFigure_4_Jan25(ifsavedata)
if nargin<1
    ifsavedata=0;
end
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
cd(datadir);load('Figure2_theta_quantification','thetacycle_lfp','HPC_layer_name');
plotcolors4=slanCM('matter',5);lettername={'a','b','c','d','e','f','g','h','i'};
rowname=[];vname={'0-1 cm/s','1-5 cm/s','5-10 cm/s','> 10 cm/s'};
for i=1:4
    for in=1:16
        rowname{in+(i-1)*16}=[vname{i},'_',num2str(in)];
    end
end
for h=1:7
    figure(4);subplot(4,2,h);
    for v=1:4
        a=squeeze(thetacycle_lfp(v,:,h,:));
        hold on;shaded_errbar([5:10:715],[a;a],plotcolors4(1+v,:));
        %a=squeeze(nanmean(thetacycle_lfp(v,:,h,:),4));
        %hold on;plot([5:10:715],[a(:);a(:)]','Color',plotcolors4(1+v,:));
    end
    xticks([0:360:720]);
    xlabel('Theta Phase');ylabel('Mean LFP Voltage (zscore)');title(HPC_layer_name(h));xlim([0 720]);
    if ifsavedata
        figdata=[squeeze(thetacycle_lfp(1,:,h,:))';squeeze(thetacycle_lfp(2,:,h,:))';squeeze(thetacycle_lfp(3,:,h,:))';squeeze(thetacycle_lfp(4,:,h,:))'];
        T = array2table(figdata,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Sup Figure 6',lettername{h}],'WriteRowNames',true);
    end
end
legend('0-1 cm/s','','1-5 cm/s','','5-10 cm/s','','> 10 cm/s','');
figure(4);figure_title='SupFigure_4_Jan2025';save_current_figure(figure_title);    