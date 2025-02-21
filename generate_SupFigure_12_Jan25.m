%function SupFigure_13
%clear all;close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
velocity_dif=1;thetapower_dif=0.1;thetaphase_dif=5;halfwinN=150;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('fig_ripple_theta_csd_lfp_v2','ripple1');
    load('theta_delta_LFP_HPC.mat','thetapower', 'thetaphase','lfptime');
    load('CSD_LFP_Visualization.mat', 'CSDall');
    load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','LFP_HPC');
    phase=thetaphase(:,4)-180;phase=phase+(phase<0)*360; % make thetapeak as 180; trough as 0
    ripple1=ripple1(ripple1(:,6)==1,:);
    LFP_HPC_z=zscore(LFP_HPC(:,1:7),0,1);
    lfptime(:,5:6)=[zscore(thetapower(:,4)),phase];
    % below randomly match a non-SWR period to each ripple with same velocity,
    % theta phase and theta power
    lfptime(:,7:8)=0;
    for i=1:size(ripple1,1)
        ind= lfptime(:,1)>=ripple1(i,1) & lfptime(:,1)<=ripple1(i,2) ;
        lfptime(ind,7)=lfptime(ind,7)+1;
    end
    lfptime(halfwinN+1:end,9)=lfptime(halfwinN+1: end,1)-lfptime(1:end-halfwinN,1);
    lfptime(1:end-halfwinN,10)=lfptime(halfwinN+1:end,1)-lfptime(1:end-halfwinN,1);
    lfptime(:,11)= (lfptime(:,9)*625/(halfwinN)>=0.9 & lfptime(:,9)*625/(halfwinN)<=1.1) & ...
        (lfptime(:,10)*625/(halfwinN)>=0.9 & lfptime(:,10)*625/(halfwinN)<=1.1);   
    ripple1match_csd=nan*ones(size(ripple1,1),size(CSDall,1),2*halfwinN+1);
    ripple1match_lfp=nan*ones(size(ripple1,1),7,2*halfwinN+1);
    ripple1match_lfp_z=nan*ones(size(ripple1,1),7,2*halfwinN+1);
    ripple1(:,26:27)=0;
    for i=1:size(ripple1,1)
        indmatch=lfptime(:,8)==0 & lfptime(:,11)==1 & lfptime(:,7)==0 & lfptime(:,3)>=0 & ...
            abs(ripple1(i,10)-lfptime(:,3))<velocity_dif & ...
            abs(ripple1(i,16)-lfptime(:,5))<thetapower_dif & ( abs(ripple1(i,17)-lfptime(:,6))<thetaphase_dif | abs(ripple1(i,17)-lfptime(:,6))>360-thetaphase_dif );
        ripple1(i,26)=sum(indmatch);
        if ripple1(i,26)>0
            lfpid=find(indmatch>0);
            shid=randperm(sum(indmatch),1);
            ripple1(i,27:28)=[lfpid(shid),lfptime(lfpid(shid),6)];
            ind= [ripple1(i,27)-halfwinN:ripple1(i,27)+halfwinN];
            ind2= [ripple1(i,27)-50:ripple1(i,27)+50];
            lfptime(ind2,8)=lfptime(ind2,8)+1;
            ripple1match_csd(i,:,:)=CSDall(:,ind);
            ripple1match_lfp(i,:,:)=LFP_HPC(ind,1:7)';
            ripple1match_lfp_z(i,:,:)=LFP_HPC_z(ind,1:7)';
        end
    end  
    cd(savedir);save('CA1_ripple1_match','ripple1match_lfp_z','ripple1match_lfp','ripple1match_csd','ripple1',...
        'halfwinN','velocity_dif','thetaphase_dif','thetapower_dif','lfptime','-v7.3');
    in
end

homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
HPC_layer_name={'CA1 pyr','CA1 st rad','CA1 slm','DG OML','DG MML','DG GCL','CA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
puor=slanCM('PuOr');puor=puor(256:-1:1,:); 
theta_time_cutoff=0.05;nn=[];mean_lfpmatch=[];mean_lfp_zmatch=[];mean_csdmatch=[];ripple1all=[];
figure(1);
for in=1%:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('CA1_ripple1_match');
    load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','dHPC_layer7channel3');
    ripple1(:,22)=0;
    ind= ripple1(:,16)>0;%ripple1(:,15)<=theta_time_cutoff &
    ripple1(ind,22)=1;
    ripple1(:,29)=in;
    ripple1all=[ripple1all;ripple1];
    a1=[];a2=[];a3=[];a4=[];
    for r=1%:2
        for t=1:2
            if t==1
                ind=ripple1(:,22)==1 & ripple1(:,6)==r & ripple1(:,26)>0;
            elseif t==2
                ind=ripple1(:,22)==0 & ripple1(:,6)==r & ripple1(:,26)>0;
            end
            nn(in,r,t)=[sum(ind)];
            a=squeeze(nanmean(ripple1match_csd(ind,:,:),1));
            a1=[a1,a];mean_csdmatch{in,r,t}=a;
            a=squeeze(nanmean(ripple1match_lfp(ind,:,:),1));
            a2=[a2,a];mean_lfpmatch(in,r,t,:,:)=a;
            a=squeeze(nanmean(ripple1match_lfp_z(ind,:,:),1));
            a4=[a4,a];mean_lfp_zmatch(in,r,t,:,:)=a;
        end
    end
    a=squeeze(nn(in,:,:));
    subplot(4,4,in);
    imagesc(a1);set(gca,'YDir','normal');title(num2str(a(:)'));
    xline([1:6]*(2*halfwinN+1));colormap(puor);caxis(0.9*max(abs(a1(:)))*[-1 1]);colorbar;
    for k=1:7
        hold on;plot(a4(k,:)*6+dHPC_layer7channel3(k,4),'k');
    end
    yticks(dHPC_layer7channel3(end:-1:1,4));yticklabels(HPC_layer_name(7:-1:1));ylabel('Channel on Probe');%set(gca,'TickLength',[0 0]);
    xticks([0.5:1.5]*(2*halfwinN+1));xticklabels({'High theta-match','Low theta-match'});
    xlabel('Time around match peak (s)');
end
cd(datadir);save('SupFigure_SWR_match','mean_lfp_zmatch',"mean_csdmatch",'mean_lfpmatch','ripple1all','nn','-v7.3');

cd(homedir);load('fig5_ripple_theta_csd_lfp_mean_v2');nn1=nn;
cd(datadir);load('SupFigure_SWR_match');
% Sup Figure 10
m=2;
for r=1%:2
    for in=1:16
        savedir=SessionSet16{in};
        cd(savedir);
        load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','dHPC_layer7channel3');
        n1=squeeze(nn1(in,r,:));n2=squeeze(nn(in,r,:));n=[n1(:);n2(:)];
        a1=[];a2=[];
        for t=1:3
            a=mean_csd{in,r,t};
            a=squeeze(a(:,:,m));a1=[a1,a];
            a=squeeze(mean_lfp_z(in,r,t,:,:,m));
            a2=[a2,a];
        end
        for t=1:2
            a=mean_csdmatch{in,r,t};a1=[a1,a];
            a=squeeze(mean_lfp_zmatch(in,r,t,:,:));
            a2=[a2,a];
        end
        figure(r);subplot(4,4,in);
        imagesc(a1);set(gca,'YDir','normal');title(num2str(n(:)'));
        xline([1:5]*(2*halfwinN+1));colormap(puor);caxis(0.7*max(abs(a1(:)))*[-1 1]);colorbar;
        for k=1:7
            if ceil(in/2)==2 | ceil(in/2)==5
                hold on;plot(a2(k,:)*3+dHPC_layer7channel3(k,4),'k');
            else
                hold on;plot(a2(k,:)*3+dHPC_layer7channel3(k,4),'k');
            end
        end
        hold on;plot([0 0]+3,[0 3]*3+dHPC_layer7channel3(4,4),'r');
        ylim([0.5 ceil(max(a2(1,:))*3)+dHPC_layer7channel3(1,4)+0.5]);
        yticks(dHPC_layer7channel3(end:-1:1,4));yticklabels(HPC_layer_name(7:-1:1));ylabel('Channel on Probe');%set(gca,'TickLength',[0 0]);
        xticks([0.5:4.5]*(2*halfwinN+1));
        if r==1
            xticklabels({'theta-CA1','none-CA1','SWS-CA1','High theta-match','Low theta-match'});
        else
            xticklabels({'theta-CA3','none-CA3','SWS-CA3'});
        end
        xlabel('Time around ripple peak (s)');
        if 0
        % below save figdata
        figdata=[a2(1:7,:);a1];
        rowname={'CA1 pyr_LFP','CA1 st rad_LFP','CA1 slm_LFP','DG OML_LFP','DG MML_LFP','DG GCL_LFP','CA3_LFP'};
        for i=1:size(a1,1)
            rowname{7+i}=['CSD_',num2str(i)];
        end
        T = array2table(figdata,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['SupFigure 10_',num2str(in)],'WriteRowNames',true);
        end
    end
    figure(r);figure_title=['SupFigure_10_csd_lfp_match_16sessions_SWRs'];save_current_figure(figure_title);
end