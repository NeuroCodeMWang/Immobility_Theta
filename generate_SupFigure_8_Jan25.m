function generate_SupFigure_8_Jan25(ifsavedata)
if nargin<1
    ifsavedata=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
if 0
ifpeak=0;pos1=1;pos2=101;h=1;thetaseq_align=[];ns=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    if ifpeak
        load('theta_peak_peak');
    else
        load('theta_trough_trough');
    end
    as=[];
    for k=1:2 % outbound/ inbound
        indk=thetacycle2(:,62)==5-k & thetacycle2(:,10)>=10;ns(in,k)=sum(indk);
        c=thetacycle2(indk,:);
        a=thetaseq2(indk,:,:,h);
        [~,I]=sort(c(:,10));
        I(:,2)=ceil(5*[1:length(I)]/length(I));
        a1=[];
        for i=1:5
            ind=I(:,2)==i;
            a1=[a1,squeeze(nanmean(a(ind,:,:),1))];
        end
        as=[as;a1(31:71,:)];%
    end
    thetaseq_align(:,:,in)=as;
end
ns(:,3)=sum(ns,2);% outbound/ inbound / sum
cd(datadir);save('SupFigure_8_sessions','thetaseq_align','ns','-v7.3');
else
    cd(datadir);load('SupFigure_8_sessions');
end
figure(4);
for in=1:16
    as=thetaseq_align(:,:,in);
    subplot(4,4,in);imagesc([ ],[-1 1]*40,as);set(gca,'YDir','normal');xticks(0:36:36*5);yticks(-40:20:40);%colorbar;: out|in
    xline([36.5:36:36*5],'w');yline(0,'w');xlabel('Theta phase');ylabel('Distance to rat (cm)');title(['N = ',num2str(ns(in,:))]);
    if ifsavedata
        T = array2table(as);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Sup Figure 8_',num2str(in)]);
    end
end
colormap('hot');
figure(4);figure_title='SupFigure_8_Jan2025';save_current_figure(figure_title);