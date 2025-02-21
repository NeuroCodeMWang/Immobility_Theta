function generate_Figure_2_theta_Jan25(ifsavedata)
if nargin<1
    ifsavedata=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
puor=slanCM('PuOr');puor=puor(256:-1:1,:);
plotcolors=slanCM('dense',16);
if 0
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
propnamegroup={'Velocity','Duration','Theta Peak Power','Theta Prominence','Ascending Half Dur','Descending Half Dur',...
    'Descend/Ascend','MML Sink','MML Source','Sink Percentile','Source Percentile','Theta Peak Percentile','Theta Trough Percentile'};
thetacycle_csdset=[];thetacycle_N=[];groupN=100;
thetacyclegroup=zeros(groupN,13,16,2);corr_csd=[];thetacycle_lfp=[];corr_csd_shuffle=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_trough_trough','thetacycle');
    load('Analyze_Theta_Cycle_Broad','CSDall_trial','reftheta','lfptime_trial', 'dHPC_layer7channel3','lfpz_trial');
    ind=thetacycle(:,79)==0 ;
    thetacycle=thetacycle(ind,:);
    indstill0=thetacycle(:,10)<1;
    indstill1=thetacycle(:,10)<5 & thetacycle(:,10)>=1;
    indstill2= thetacycle(:,10)<10 & thetacycle(:,10)>=5;
    indlowv= thetacycle(:,10)>=10 ;
    indhighv= thetacycle(:,10)>=30 ;
    thetacycle_N(in,:)=[sum(indstill0),sum(indstill1),sum(indstill2),sum(indlowv),sum(indhighv)];
    thetacycle_csd=[];
    for i=1:36
        indi=ceil(reftheta(:,4)/10)==i;
        for v=1:4
            if v==1
                indv=indstill0;
            elseif v==2
                indv=indstill1;
            elseif v==3
                indv=indstill2;
            elseif v==4
                indv=indlowv;
            elseif v==5
                indv=indhighv;
            end
            lia=ismember(lfptime_trial(:,5),thetacycle(indv,13));
            ind= lia & indi;
            thetacycle_csd(v,i,:)=nanmean(CSDall_trial(:,ind),2);
            thetacycle_lfp(v,i,:,in)=nanmean(lfpz_trial(ind ,:),1);
        end
    end
    thetacycle_csdset{in}=thetacycle_csd;
    a1=squeeze(thetacycle_csd(1,:,:));a1=a1(:);
    for v=2:4
        a=squeeze(thetacycle_csd(v,:,:));a=a(:);
        [r,p]=corrcoef(a,a1);
        corr_csd(in,v-1,:)=[r(2),p(2)];
    end
    for sh=1:100
        psh=randperm(length(a1));
        a1_sh=a1(psh);
        for v=2:4
            a=squeeze(thetacycle_csd(v,:,:));a=a(:);
            [r,p]=corrcoef(a,a1_sh);
            corr_csd_shuffle(sh,v-1,:,in)=[r(2),p(2)];
        end
    end
    [B,I]=sort(thetacycle(:,10));
    I(:,2)=[1:length(I)]/length(I);
    for g=1:groupN
        indv=I(ceil(I(:,2)*groupN)==g,1);
        thetacyclegroup(g,:,in,1)=nanmean(thetacycle(indv,[10,3,4,8,15,14,16,17,19,18,20,7,5]),1);
        thetacyclegroup(g,:,in,2)=nanmedian(thetacycle(indv,[10,3,4,8,15,14,16,17,19,18,20,7,5]),1);
    end
end
cd(datadir);save('Figure2_theta_quantification','corr_csd_shuffle','corr_csd','thetacyclegroup','thetacycle_csdset','thetacycle_N',...
    'thetacycle_lfp','HPC_layer_name','propnamegroup','-v7.3');
else
    cd(datadir);load('Figure2_theta_quantification','corr_csd_shuffle','corr_csd','thetacyclegroup','thetacycle_csdset','thetacycle_N',...
    'thetacycle_lfp','HPC_layer_name','propnamegroup');
end
%---------------------------------------------------------------------------------
% 1. figure2 theta quantification
pn1=2;pn2=5;  
figure(3);subplot(pn1,pn2,1);
for h=1:7
    as=[];
    for v=1:4
        a=squeeze(thetacycle_lfp(v,:,h,:));
        as=[as;a];
    end
    hold on;plot([5:10:4*360-5],8-h+nanmean(as,2),'k');%shaded_errbar([5:10:4*360-5],8-h+as,'k');
    %hold on;plot([5:10:715],[a(:);a(:)]','Color',plotcolors4(1+v,:));
end
xticks([0:360:360*4]);xlabel('Theta Phase');ylabel('Mean LFP Voltage (z)');xlim([0 360*4]);
yticks([1:7]);yticklabels(HPC_layer_name(7:-1:1));xline([0:360:4*360],'k');

a1=subplot(pn1,pn2,2);in=1;  
savedir=SessionSet16{in};
cd(savedir);load('Analyze_Theta_Cycle_Broad','dHPC_layer7channel3');
thetacycle_csd=thetacycle_csdset{in};
a=[squeeze(thetacycle_csd(1,:,:))',squeeze(thetacycle_csd(2,:,:))',squeeze(thetacycle_csd(3,:,:))',squeeze(thetacycle_csd(4,:,:))'];
b=[squeeze(thetacycle_lfp(1,:,1:7,in))',squeeze(thetacycle_lfp(2,:,1:7,in))',squeeze(thetacycle_lfp(3,:,1:7,in))',squeeze(thetacycle_lfp(4,:,1:7,in))'];
imagesc([0 36*4],[ ],a);set(gca,'YDir','normal');xline([1:4]*36,'k');colormap(a1,puor);
for i=1:7
    hold on;plot([1:36*4],b(i,:)*3+dHPC_layer7channel3(i,4),'k');
end
hold on;plot([0 0],[0 1]*3+dHPC_layer7channel3(2,4),'r');
yticks(dHPC_layer7channel3(end:-1:1,4));yticklabels(HPC_layer_name(7:-1:1));ylabel('Channel on Probe');%set(gca,'TickLength',[0 0]);
xticks([0.5,1.5,2.5,3.5]*36);xticklabels({'0-1','1-5','5-10','> 10'});xlabel('Theta cycles grouped by rat velocity (cm/s)');
title('Mean Theta CSD for an Example Session');colorbar;caxis(max(abs(a(:)))*[-1 1]);

figure(3);subplot(pn1,pn2,3);
x = [ones(1,1600) 2*ones(1,1600) 3*ones(1,1600)];
y1 = corr_csd_shuffle(:,1,1,:);
y2 = corr_csd_shuffle(:,2,1,:);
y3 = corr_csd_shuffle(:,3,1,:);
y = [y1(:)' y2(:)' y3(:)'];
swarmchart(x,y,5,[1 1 1]*0.4,'.');
for in=1:16
    hold on;plot([1:3],corr_csd(in,:,1),'ko-');
end
xlim([0.5 3.5]);xticks([1:3]);xticklabels({'0-1 vs 1-5','0-1 vs 5-10','0-1 vs > 10'});ylabel('Correlation coefficient');ylim([-0.2 1.1])
xlabel('Theta CSD Pair: 0-1 vs ');title('Correlation with theta CSD at 0-1');legend('Shuffle','Session Data');
if ifsavedata
    figdata=[squeeze(corr_csd(:,:,1))',[y1(:), y2(:), y3(:)]']';
    T = array2table(figdata,'VariableNames',{'1-5','5-10','> 10'});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2b']);
end

propset=[2,5:7,3,8:10];mid=2;letternum={'c','d','e','f','g','h','i'};
thetacycleprop_vel_corr=[];
for pid=1:7%8
    p=propset(pid);
    figure(3);subplot(pn1,pn2,pid+3);hold on;
    for in=1:16
        if p<=9
            plot(log10(thetacyclegroup(:,1,in,mid)),thetacyclegroup(:,p,in,mid),'-','color',[plotcolors(in,:),0.3]);%[0.6431 0.8078 0.9176]);%
            xlabel('Velocity (cm/s)');ylabel(propnamegroup{p});xlim([-1 2.3]);
            xticks(log10([0.1 1 5 10 30 80]));xticklabels({'0','1','5','10','30','80'});
            hold on;plot(log10(nanmean(thetacyclegroup(:,1,:,mid),3)),nanmean(thetacyclegroup(:,p,:,mid),3),'k-');
            if in==16
                [r,pval]=corrcoef(log10(nanmean(thetacyclegroup(:,1,:,mid),3)),nanmean(thetacyclegroup(:,p,:,mid),3));
                thetacycleprop_vel_corr(pid,:)=[r(2),pval(2)];
                title(num2str([r(2),pval(2)],2));
            end
        else
            plot(thetacyclegroup(:,p,in,mid),thetacyclegroup(:,p+2,in,mid),'-','color',[plotcolors(in,:),0.3]);%[0.6431 0.8078 0.9176]);%,
            hold on;plot(nanmean(thetacyclegroup(:,p,:,mid),3),nanmean(thetacyclegroup(:,p+2,:,mid),3),'k-');
            xlabel(propnamegroup{p});ylabel(propnamegroup{p+2});
            if in==16
                [r,pval]=corrcoef(nanmean(thetacyclegroup(:,p,:,mid),3),nanmean(thetacyclegroup(:,p+2,:,mid),3));
                thetacycleprop_vel_corr(pid,:)=[r(2),pval(2)];
            end
        end
    end
    if pid==6
        ylim([-inf 0]);
    elseif pid<8
        ylim([0, inf]);
    end
    if pid<=7 & ifsavedata
        figdata=[squeeze(thetacyclegroup(:,1,:,mid))';squeeze(thetacyclegroup(:,p,:,mid))'];
        rowname={'X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9','X_10','X_11','X_12','X_13','X_14','X_15','X_16'...
            'Y_1','Y_2','Y_3','Y_4','Y_5','Y_6','Y_7','Y_8','Y_9','Y_10','Y_11','Y_12','Y_13','Y_14','Y_15','Y_16'};
        T = array2table(figdata,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2',letternum{pid}],'WriteRowNames',true);
    end
end
figure(3);figure_title='Figure_2_theta_quantification_Jan2025';save_current_figure(figure_title);
%save('Figure_2_theta_statistics','thetacycleprop_vel_corr');