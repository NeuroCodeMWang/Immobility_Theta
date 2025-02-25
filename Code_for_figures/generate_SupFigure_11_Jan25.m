%function generate_SupFigure_11_Jan25(ifsavedata)
%if nargin<1
    ifsavedata=0;
%end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
vgroup={'0-1','1-5','5-10','> 10'};
vname={'Pause','Run'};
bname={'Reward','Center'};
behave={'Reward','Center','In','Out'};
localremote={'local','remote'};
regionrepname={'Center','Stem','Reward'};
rname={'CA1 ripples','CA3 ripples'};
HPC_layer_name={'CA1 pyr','CA1 st rad','CA1 slm','DG OML','DG MML','DG GCL','CA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
plotcolors3=slanCM('bone',8);plotcolors3=plotcolors3(6:-2:2,:);
plotcolors4=slanCM('Paired',12);
defaultcolor4=slanCM('YlGnBu',8);defaultcolor4=defaultcolor4([8:-2:1],:);defaultcolor4(4,:)=[0.8828 0.7305 0.2656];%[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560]];
plotcolors5=slanCM('bone',6);plotcolors5=plotcolors5(6:-1:1,:);
reward_dis_cut=70;
cd(homedir);load('next_arm_deseq_theta_ripple_6');
cd(datadir);load('Figure4_next_shuffle_theta_ripple_none_data_v6','next_shuffle','session_rep_n','ripple_percent','ripple_percent16',...
    'theta_dist','theta_dist16','rippleall','reward_center_correlation');
pksetsall(:,7)=round(pksetsall(:,7));
pksetsall3(:,7)=round(pksetsall3(:,7));
mu_theta=next_shuffle.theta.mu;
err_theta=next_shuffle.theta.err;
ns_theta=next_shuffle.theta.ns;
ps_theta=next_shuffle.theta.ps;
ps_theta_r=next_shuffle.theta.ps_r;
ts_theta_r=next_shuffle.theta.ts_r;
ps2_theta=next_shuffle.theta.ps_forcefree;
mu_none=next_shuffle.none.mu;
err_none=next_shuffle.none.err;
ns_none=next_shuffle.none.ns;
ps_none=next_shuffle.none.ps;
ps_none_r=next_shuffle.none.ps_r;
ts_none_r=next_shuffle.none.ts_r;
ps2_none=next_shuffle.none.ps_forcefree;
mu_ripple=next_shuffle.ripple.mu;
err_ripple=next_shuffle.ripple.err;
ns_ripple=next_shuffle.ripple.ns;
ps_ripple=next_shuffle.ripple.ps;
ps_ripple_r=next_shuffle.ripple.ps_r;
ts_ripple_r=next_shuffle.ripple.ts_r;
ps2_ripple=next_shuffle.ripple.ps_forcefree;
mu_ripple_pk=next_shuffle.ripple_pk.mu;
err_ripple_pk=next_shuffle.ripple_pk.err;
ns_ripple_pk=next_shuffle.ripple_pk.ns;
deseq_theta_mean=reward_center_correlation.deseq_theta_mean;
deseq_none_mean=reward_center_correlation.deseq_none_mean;
deseq_ripple_mean=reward_center_correlation.deseq_ripple_mean;
% sup figure 11
% below for last arm prediction
fign=10;pn1=2;pn2=3;h=3;hname={'next','current','last','center'}; 
for b=1:2
    figure(fign);subplot(pn1,pn2,b); % for theta
    hold on;ba=bar([1:3],squeeze(mu_theta(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_theta(b,1,1:3,h)),squeeze(err_theta(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_theta(b,2,1:3,h)),squeeze(err_theta(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_theta(b,1:2,1:3,h);ar=squeeze(ps_theta_r(b,1:2,1:3,h))';tr=squeeze(ts_theta_r(b,1:2,1:3,h))';
   % title({ ['p:',num2str(a(:)',2)];['p2:',num2str(ps2_theta(b,1:3,h),2)];['p:',num2str(ar(:)',2)]});
    title({ ['t:',num2str(tr(:)',4)];['p:',num2str(ar(:)',2)]});
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel([bname{b},' theta']);xlim([0.5 3.5]);
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');

    figure(fign);subplot(pn1,pn2,b+2);  % for none
    hold on;ba=bar([1:3],squeeze(mu_none(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_none(b,1,1:3,h)),squeeze(err_none(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_none(b,2,1:3,h)),squeeze(err_none(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_none(b,1:2,1:3,h);n=ns_none(b,:,1:3,h);ar=squeeze(ps_none_r(b,1:2,1:3,h))';tr=squeeze(ts_none_r(b,1:2,1:3,h))';
    %title({['p:',num2str(a(:)',2)];['p2:',num2str(ps2_none(b,:,h),2)];['p:',num2str(ar(:)',2)]});
    title({ ['t:',num2str(tr(:)',4)];['p:',num2str(ar(:)',2)]});
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel([bname{b},' none']);xlim([0 3.5]);
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');
    
    figure(fign);subplot(pn1,pn2,b+4);  % for99 ripple
    hold on;ba=bar(1,squeeze(mu_ripple(b,1:2,1,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar(1-bar_pos,squeeze(mu_ripple(b,1,1,h)),squeeze(err_ripple(b,1,1,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar(1+bar_pos,squeeze(mu_ripple(b,2,1,h)),squeeze(err_ripple(b,2,1,h)),'k.','LineWidth',errbar_width);
    ps_r1=squeeze(ps_ripple(b,1:2,1,h))';ar=squeeze(ps_ripple_r(b,1:2,1,h));tr=squeeze(ts_ripple_r(b,1:2,1,h));
    %title({['p:',num2str(ps_r1(:)',2)];['p2:',num2str(ps2_ripple(b,1,h),2)];['p:',num2str(ar(:)',2)]});%['n = ',num2str(ns_r1(:)')];
    title({ ['t:',num2str(tr(:)',4)];['p:',num2str(ar(:)',2)]});
    xticks(1);xticklabels({'CA1 ripples'});xlabel([bname{b},' ripple']);
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');
end
figure(fign);figure_title='SupFigure11_Jan2025';save_current_figure(figure_title);