%function generate_Figure_4_Jan25(ifsavedata)
%if nargin<1
    ifsavedata=0;
%end

% 0. variables for plotting
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
puor=slanCM('PuOr');puor=puor(256:-1:1,:); 
reward_dis_cut=70;
%---------------------------------------------------------------------------------------
% fig 4-1. below quantify ripple rate during theta vs non-theta at reward pause
cd(datadir);load('SupFigure_SWR_match');nnmatch=nn;
cd(homedir);load('pk_behaveall','pk_behaveall','indpk','rippleall','pksetsall','ca1ca3_overlap');
load('fig5_ripple_theta_csd_lfp_mean_v2');
pksetsall(:,7)=round(pksetsall(:,7));
if 0
indi= (pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60) | ((pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3) ;
behave_time=[];
for b=1:4
    for in=1:16
        indin= pksetsall(:,42)==in & indi & indpk; 
        behave_time(in,b,1:4,:)=nanmean(pk_behaveall(indin,b,:,[1:2,13]),1);
    end
end
pn1=3;pn2=4;fign=1;
figure(fign);subplot(pn1,pn2,1);
for b=1:4
    y=squeeze(behave_time(:,b,:,1))';
    hold on;shaded_errbar([1:4],y,defaultcolor4(b,:));
end
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Velocity Ranges');
ylabel('Total time (s)');legend({behave{1},'',behave{2},'',behave{3},'',behave{4},''});%ylim([0 50]);

figure(fign);subplot(pn1,pn2,2);
for b=1:4
    y=squeeze(behave_time(:,b,:,2))';
    hold on;shaded_errbar([1:4],y,defaultcolor4(b,:));
end
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Velocity Ranges');
ylabel('Total theta time (s)');legend({behave{1},'',behave{2},'',behave{3},'',behave{4},''});%ylim([0 50]);

figure(fign);subplot(pn1,pn2,3);
for b=1:4
    y=squeeze(behave_time(:,b,:,3))';
    hold on;shaded_errbar([1:4],y,defaultcolor4(b,:));
end
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Velocity Ranges');ylim([0 1]);
ylabel('Theta time proportion');legend({behave{1},'',behave{2},'',behave{3},'',behave{4},''});%ylim([0 50]);

for r=1:2
    figure(fign);subplot(pn1,pn2,3+r);
    for b=1:4
        y=squeeze(pk_behaveall(indpk,b,:,13+r))';
        hold on;shaded_errbar([1:4],y,defaultcolor4(b,:));
    end
    xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Velocity Ranges');title(rname{r});
    ylabel('Ripple rate (/s)');legend({behave{1},'',behave{2},'',behave{3},'',behave{4},''});%ylim([0 50]);
end
figure(fign);subplot(pn1,pn2,6);b=1;
for r=1:2
    y=squeeze(pk_behaveall(indpk,b,:,15+r))';
    hold on;shaded_errbar([1:4],y,plotcolors5(r*2+1,:));
end
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Velocity Ranges');
ylabel('Ripple rate (/s)');legend({'theta-ca3 ripple','','Nontheta-ca3 ripple',''});%ylim([0 50]);

edges=[-2.4:0.2:5.2];edges_plot=[-2.3:0.2:5.1];
for r=1:2
    indr=rippleall(:,6)==r;
    figure(fign);subplot(pn1,pn2,6+r);
    ind0=rippleall(:,11)==1;
    ind1=(rippleall(:,11))==2;
    ind2=rippleall(:,11)==3;
    ind3=rippleall(:,11)==4;
    forpie=[sum(ind0& indr),sum(ind1& indr),sum(ind2& indr),sum(ind3& indr)];
    h=pie(forpie);patchHand = findobj(h, 'Type', 'Patch'); 
    set(patchHand, {'FaceColor'}, mat2cell(defaultcolor4, ones(size(defaultcolor4,1),1), 3));
    set(h,'LineWidth',2,'EdgeColor','k');
    legend('Reward','Inter','In','Out');
    title([rname{r},' : N = ',num2str(sum(forpie))]);
end
figure(fign);subplot(pn1,pn2,9);
for r=1%:2
    indr=rippleall(:,6)==r;
    [a]=histcounts(rippleall(indr,16),edges);a=a/sum(a);
    if r==2
        hold on;bar(edges_plot,a,'FaceAlpha',0.3,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1);
    else
        hold on;bar(edges_plot,a,'FaceAlpha',0.7,'FaceColor',[1 1 1]*0.5,'EdgeColor',[0 0 0],'LineWidth',1);
    end
end
xline(0,'k');xlabel('Theta power at ripple peak (z-score)');xlim([-2.4 4]);
ylabel('Ripple Proportion');%legend(rname);
figdata=rippleall(indr,16);
T = array2table(figdata);
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4f']);

figure(fign);subplot(pn1,pn2,10);
for r=1%:2
    indr=rippleall(:,6)==r;
    [a]=histcounts(rippleall(indr,17),[0:10:360]);a=a(:)/sum(a);
    if r==2
        hold on;bar([5:10:715],[a;a],'FaceAlpha',0.3,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1);
    else
        hold on;bar([5:10:715],[a;a],'FaceAlpha',0.7,'FaceColor',[1 1 1]*0.5,'EdgeColor',[0 0 0],'LineWidth',1);
    end
end
xlabel('Theta phase at ripple peak');xlim([0 720]);ylabel('Ripple Proportion');%

m=2;r=1;a1=[];
for t=1:3
    a=squeeze(nanmean(mean_lfp_z(:,r,t,:,:,m),1));
    a1=[a1,a];
end
for t=1:2
    a=squeeze(nanmean(mean_lfp_zmatch(:,r,t,:,:),1));
    a1=[a1,a];
end
figure(fign);subplot(pn1,pn2,11);
for k=1:7
    hold on;plot(a1(k,:)+8-k,'k');
end
yticks([1:7]);yticklabels(HPC_layer_name(7:-1:1));ylabel('Channel on Probe');%set(gca,'TickLength',[0 0]);
xticks([0.5:4.5]*(2*halfwinN+1));xline([1:5]*(2*halfwinN+1),'k');xlim([0 5*(2*halfwinN+1)]);
xticklabels({'theta-CA1','none-CA1','SWS-CA1','High theta-match','Low theta-match'});
xlabel('Time around ripple peak (s)');

for r=2
    a1=[];
    for t=1:3
        a=squeeze(nanmean(mean_lfp_z(:,r,t,:,:,m),1));
        a1=[a1,a];
    end
    figure(fign);subplot(pn1,pn2,10+r);
    for k=1:7
        hold on;plot(a1(k,:)+8-k,'k');
    end
    yticks([1:7]);yticklabels(HPC_layer_name(7:-1:1));ylabel('Channel on Probe');%set(gca,'TickLength',[0 0]);
    xticks([0.5:2.5]*(2*halfwinN+1));xline([1:3]*(2*halfwinN+1),'k');xlim([0 3*(2*halfwinN+1)]);
    if r==1
        xticklabels({'theta-CA1','none-CA1','SWS-CA1'});
    else
        xticklabels({'theta-CA3','none-CA3','SWS-CA3'});
    end
    xlabel('Time around ripple peak (s)');
    %FIG_INDEX=['fig4_1_',num2str(9+r)];save_fig(FIG_INDEX,ifsavefig,[],1);
end
figure(fign);figure_title='Figure4_ripple_basic';save_current_figure(figure_title);
end
%-------------------------------------------------------------------------------------
% fig4 major plots
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

fign=1;pn1=4;pn2=4;errbar_width=2;bar_pos=0.14;
% 'fig4_2_next_shuffle'
for b=1:2
    figure(fign);subplot(pn1,pn2,b);
    ba=bar([1:4],squeeze(theta_dist(1:4,b,1:3)),'stacked','EdgeColor',[0 0 0],'LineWidth',2);
    for i=1:3
        ba(i).FaceColor = plotcolors5(i*2-1,:);
    end
    xticks([1:4]);xticklabels(vgroup(1:4));xlim([.5 4.5]);xlabel('Thetaseqs of Velocity Ranges');%xlim([0 4]); 
    ylim([0 1]);ylabel('Percent');title([bname{b},' thetacycle N = ',num2str(squeeze(theta_dist(1:4,b,4))')]);
    if b==1
        legend('Center','Local arm','Remote arm');
    else
        legend('Center','','Arm');
    end
end

% ripple percent
figure(fign);subplot(pn1,pn2,3);
ba=bar([1:2],squeeze(ripple_percent(1,1:2,1:3)),'stacked','EdgeColor',[0 0 0],'LineWidth',2);
for i=1:3
    ba(i).FaceColor = plotcolors5(i*2-1,:);
end
xticks([1:2]);xticklabels({'ca1 ripple at reward','ca1 ripple at center'});xlim([0.5 2.5]); 
ylabel('Ripple proportion');title(['ripple N = ',num2str(squeeze(ripple_percent(1,:,4)))]);
legend('Center','Local arm','Remote arm');

% below compare theta & ripple duration
pksetsall(:,7)=round(pksetsall(:,7));
thetarippledur_session=[];
for i=1:3
    if i==1 % correct forced visit
        indi= pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60 ;
    elseif i==2 % free visit
        indi= (pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3;
    else
        indi= (pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60) | ((pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3) ;
    end
    for b=1:2
        % below for theta
        for in=1:16
            indin= pksetsall3(:,45)==in & indi ; 
            thetarippledur_session(in,b,i,1:4)=mean(pk_behaveall(indin,b,:,2),1);
            thetarippledur_session(in,b,i,5:8)=mean(sum(pk_behaveall(indin,b,1:3,18:21),3),1);
        end
    end
end
thetaripple_mu=[];thetaripple_err=[];
for b=1:2
    for i=1:3
        thetaripple_mu(b,i,:)=mean(thetarippledur_session(:,b,i,:),1);
        thetaripple_err(b,i,:)=std(thetarippledur_session(:,b,i,:),0,1)/4;
    end
end

xdif=0.14;n1=0.05;d=1;n=16;errbar_width=2;bar_pos=0.14;
% below for next coding proportion
figure(fign);subplot(pn1,pn2,4);i=3;m=1;h=2;comp=[];indpp=[1:3,4.5];
mu=[];err=[];
for b=1:2
    a=[squeeze(session_rep_n(b,i,1,1:3,h,:,m));squeeze(session_rep_n(b,i,2,1,h,:,m))'];
    mu(:,b)=nanmean(a,2);
    err(:,b)=nanstd(a,0,2)/4;
end
hold on;ba=bar(indpp,mu,'FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
ba(1).FaceColor =[1 1 1];ba(2).FaceColor =defaultcolor4(1,:);ba(2).FaceAlpha=0.3;
hold on;errorbar(indpp-bar_pos,mu(:,1),err(:,1),'k.','LineWidth',errbar_width);
hold on;errorbar(indpp+bar_pos,mu(:,2),err(:,2),'k.','LineWidth',errbar_width);
dat1=[];stats=[];
for b=1:2
    f0=squeeze(session_rep_n(b,i,2,1,h,:,m));
    for v=1:3
        f1=squeeze(session_rep_n(b,i,1,v,h,:,m));
        [p,~,stat]=ranksum(f1,f0,"Tail","right",'method','approximate');%p=signrank(f1-f0,0,"Tail","right");%p=signrank(f1-f0);%[~,p]=ttest(f1-f0,0,"Tail","right");%
        comp(b,v)=p;dat1(:,v,b)=f1-f0;stats(b,v)=stat.zval;
        hold on;
        if b==1
            plot(d*v*ones(n,1)+n1*randn(n,1)-xdif,f1,'k.');
            %hold on;plot(d*v*ones(n,1)-xdif+[-0.1 0.1],mean(f1)*[1 1],'-');
        else
            plot(d*v*ones(n,1)+n1*randn(n,1)+xdif,f1,'.','color',defaultcolor4(1,:));
            %hold on;plot(d*v*ones(n,1)+xdif+[-0.1 0.1],mean(f1)*[1 1],'-','color',defaultcolor4(2,:));
        end
    end
    hold on;
    if b==1
        plot(d*4.5*ones(n,1)+n1*randn(n,1)-xdif,f0,'k.');
        %hold on;plot(d*5*ones(n,1)-xdif+[-0.1 0.1],mean(f0)*[1 1],'k-');
    else
        plot(d*4.5*ones(n,1)+n1*randn(n,1)+xdif,f0,'.','color',defaultcolor4(1,:));
        %hold on;plot(d*5*ones(n,1)+xdif+[-0.1 0.1],mean(f0)*[1 1],'-','color',defaultcolor4(2,:));
    end
end
xlim([0 5.5]);xticks(d*[1:3,4.5]);ylabel('Next arm coding proportion of visits');
xticklabels({'Theta: 0-1 cm/s','Theta: 1-5 cm/s','Theta: 5-10 cm/s','CA1 SWRs'});
title(num2str(comp,2));xlabel(num2str(stats,2));

% below for duration
figure(fign);subplot(pn1,pn2,5);i=3;comp=[];indp=[1:3,5];
ba=bar(indpp,squeeze(thetaripple_mu(:,3,indp))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
ba(1).FaceColor =[1 1 1];ba(2).FaceColor =defaultcolor4(1,:);ba(2).FaceAlpha=0.3;
hold on;errorbar(indpp-bar_pos,squeeze(thetaripple_mu(1,3,indp)),squeeze(thetaripple_err(1,3,indp)),'k.','LineWidth',errbar_width);
hold on;errorbar(indpp+bar_pos,squeeze(thetaripple_mu(2,3,indp)),squeeze(thetaripple_err(2,3,indp)),'k.','LineWidth',errbar_width);
dat=[];stats=[];
for b=1:2
    f0=squeeze(thetarippledur_session(:,b,i,5));
    for v=1:3
        f1=squeeze(thetarippledur_session(:,b,i,v));
        [p,~,stat]=ranksum(f1,f0,"Tail","right",'method','approximate');%p=signrank(f1-f0,0,"Tail","right");%[~,p]=ttest(f1-f0,0,"Tail","right");%
        comp(b,v)=p;dat(:,v,b)=f1-f0;stats(b,v)=stat.zval;
        hold on;
        if b==1
            plot(d*v*ones(n,1)+n1*randn(n,1)-xdif,f1,'k.');
            %hold on;plot(d*v*ones(n,1)-xdif+[-0.1 0.1],mean(f1)*[1 1],'-');
        else
            plot(d*v*ones(n,1)+n1*randn(n,1)+xdif,f1,'.','color',defaultcolor4(1,:));
            %hold on;plot(d*v*ones(n,1)+xdif+[-0.1 0.1],mean(f1)*[1 1],'-','color',defaultcolor4(2,:));
        end
    end    
    hold on;
    if b==1
        plot(d*4.5*ones(n,1)+n1*randn(n,1)-xdif,f0,'k.');
        %hold on;plot(d*5*ones(n,1)-xdif+[-0.1 0.1],mean(f0)*[1 1],'k-');
    else
        plot(d*4.5*ones(n,1)+n1*randn(n,1)+xdif,f0,'.','color',defaultcolor4(1,:));
        %hold on;plot(d*5*ones(n,1)+xdif+[-0.1 0.1],mean(f0)*[1 1],'-','color',defaultcolor4(2,:));
    end
end
ylabel('Duration per visit (s)');title(num2str(comp,2));xlim([0 6]);xlabel(num2str(stats,2));
xticks([1:3,4.5]);xticklabels({'Theta: 0-1 cm/s','Theta: 1-5 cm/s','Theta: 5-10 cm/s','CA1 SWRs'});
legend('Rat at reward','Rat at center');


% below for next arm prediction
h=1;hname={'next','current','last','center'}; % next arm
for b=1:2
    figure(fign);subplot(pn1,pn2,b+5); % for theta
    hold on;ba=bar([1:3],squeeze(mu_theta(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_theta(b,1,1:3,h)),squeeze(err_theta(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_theta(b,2,1:3,h)),squeeze(err_theta(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_theta(b,1:2,1:3,h);ar=squeeze(ps_theta_r(b,1:2,1:3,h))';tr=squeeze(ts_theta_r(b,1:2,1:3,h))';
    title({ ['t:',num2str(tr(:)',4)];['p:',num2str(ar(:)',2)]});%['p2:',num2str(ps2_theta(b,1:3,h),2)];
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel([bname{b},' theta']);xlim([0.5 3.5]);
    ylabel('P(next > shuffle)-P(next < shuffle)');legend('Forced','Free');

    figure(fign);subplot(pn1,pn2,b+7);  % for none
    hold on;ba=bar([1:3],squeeze(mu_none(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_none(b,1,1:3,h)),squeeze(err_none(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_none(b,2,1:3,h)),squeeze(err_none(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_none(b,1:2,1:3,h);n=ns_none(b,:,1:3,h);ar=squeeze(ps_none_r(b,1:2,1:3,h))';tr=squeeze(ts_none_r(b,1:2,1:3,h))';
    title({['t:',num2str(tr(:)',4)];['p:',num2str(ar(:)',2)]});%['p2:',num2str(ps2_none(b,:,h),2)];
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel([bname{b},' none']);xlim([0 3.5]);
    ylabel('P(next > shuffle)-P(next < shuffle)');legend('Forced','Free');
    
    figure(fign);subplot(pn1,pn2,b+9);  % for ripple
    hold on;ba=bar(1,squeeze(mu_ripple(b,1:2,1,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar(1-bar_pos,squeeze(mu_ripple(b,1,1,h)),squeeze(err_ripple(b,1,1,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar(1+bar_pos,squeeze(mu_ripple(b,2,1,h)),squeeze(err_ripple(b,2,1,h)),'k.','LineWidth',errbar_width);
    ps_r1=squeeze(ps_ripple(b,1:2,1,h))';ar=ps_ripple_r(b,1:2,1,h);tr=squeeze(ts_ripple_r(b,1:2,1,h));
    title({['t:',num2str(tr(:)',4)];['p:',num2str(ar(:)',2)]});%['p2:',num2str(ps2_ripple(b,1,h),2)];['n = ',num2str(ns_r1(:)')];
    xticks(1);xticklabels({'CA1 ripples'});xlabel([bname{b},' ripple']);
    ylabel('P(next > shuffle)-P(next < shuffle)');legend('Forced','Free');
end
% below for correlation btw center & reward
nn=[];binum=20;ifbin=1;n_cutoff=10;markersize=15;
for ii=1:2
    i=3-ii;
    if i==1 % correct forced visit
        indi= pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60 ;
    elseif i==2 % free visit
        indi= (pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3;
    end
    %indi=indi & pksetsall(:,3)>=reward_dis_cut ;%& pksetsall(:,19)<30;
    for h=1%:4
        if h==1
            ind0=indi & pksetsall3(:,44)>0; % next
        elseif h==3
            ind0=indi & pksetsall3(:,43)>0; % last
        else
            ind0=indi;
        end
        % below for theta
        ind=ind0 & squeeze(min(deseq_theta_mean(:,:,5),[],2))>=n_cutoff;
        nn(1,1,i,h)=sum(ind0&squeeze(deseq_theta_mean(:,1,5))>=n_cutoff);
        nn(1,2,i,h)=sum(ind0&squeeze(deseq_theta_mean(:,2,5))>=n_cutoff);
        nn(1,3,i,h)=sum(ind);
        nn(1,6,i,h)=sum(ind0);
        a=squeeze(deseq_theta_mean(ind,1:2,h));
        if sum(isnan(a(:))) 
            'error theta!'
        end
        if ifbin
        [~,I]=sort(a(:,1));
        a=a(I,:);
        a(:,3)=ceil(binum*[1:size(a,1)]/size(a,1));
        ab=[];
        for j=1:binum
            indj=a(:,3)==j;
            ab(j,:)=mean(a(indj,1:2),1);
        end
        else
            ab=a;
        end
        figure(fign);subplot(pn1,pn2,12);
        if i==2
            hold on;plot(ab(:,1),ab(:,2),'o','MarkerFaceColor',plotcolors5(3,:),'MarkerEdgeColor','k');
        else
            hold on;plot(ab(:,1),ab(:,2),'ko');
        end
        [rh,pv]=corr(a(:,1),a(:,2),'type','Spearman');%[rh,pv]=corrcoef(ab(:,1),ab(:,2));
        nn(1,4:5,i,h)=[rh,pv];
        rp=squeeze(nn(1,4:5,:,h));
        title(['Theta: R,P= ',num2str([rp(:)'],2)]);
        %a=histcounts2(ab(:,1),ab(:,2),histn);imagesc(a');set(gca,'YDir','normal');
        xlabel({['prob at reward ',num2str(squeeze(nn(1,[6,1:3],1,h)))];[num2str(squeeze(nn(1,[6,1:3],2,h)))]});
        ylabel([hname{h},' prob at center']);
        if 0%ii==2
        Y_Limits=(ylim);ylim([0 ceil(Y_Limits(2)*10)/10]);
        X_Limits=(xlim);xlim([-0.1 ceil(X_Limits(2)*10)/10]);
        %set(findall(gcf, 'Type', 'Line'),'MarkerSize', markersize);
        %FIG_INDEX=['fig4_2_reward_center_corr_theta'];save_fig(FIG_INDEX,ifsavefig);
        end

        % below for none
        ind=ind0 & squeeze(min(deseq_none_mean(:,:,5),[],2))>=n_cutoff;
        nn(2,1,i,h)=sum(ind0&squeeze(deseq_none_mean(:,1,5))>=n_cutoff);
        nn(2,2,i,h)=sum(ind0&squeeze(deseq_none_mean(:,2,5))>=n_cutoff);
        nn(2,3,i,h)=sum(ind);
        nn(2,6,i,h)=sum(ind0);
        a=squeeze(deseq_none_mean(ind,1:2,h));
        if sum(isnan(a(:))) 
            'error none!'
        end
        if ifbin
        [~,I]=sort(a(:,1));
        a=a(I,:);
        a(:,3)=ceil(binum*[1:size(a,1)]/size(a,1));
        ab=[];
        for j=1:binum
            indj=a(:,3)==j;
            ab(j,:)=mean(a(indj,1:2),1);
        end
        else
            ab=a;
        end
        figure(fign);subplot(pn1,pn2,13);
        if i==2
            hold on;plot(ab(:,1),ab(:,2),'o','MarkerFaceColor',plotcolors5(3,:),'MarkerEdgeColor','k');
        else
            hold on;plot(ab(:,1),ab(:,2),'ko');
        end
        [rh,pv]=corr(a(:,1),a(:,2),'type','Spearman');%[rh,pv]=corrcoef(ab(:,1),ab(:,2));
        nn(2,4:5,i,h)=[rh,pv];
        rp=squeeze(nn(2,4:5,:,h));
        title(['None: R,P= ',num2str([rp(:)'],2)]);
        xlabel({['prob at reward ',num2str(squeeze(nn(2,[6,1:3],1,h)))];[num2str(squeeze(nn(2,[6,1:3],2,h)))]});
        ylabel([hname{h},' prob at center']);
        if 0%ii==2
        Y_Limits=(ylim);ylim([0 ceil(Y_Limits(2)*10)/10]);
        X_Limits=(xlim);xlim([-0.1 ceil(X_Limits(2)*10)/10]);
        %set(findall(gcf, 'Type', 'Line'),'MarkerSize', markersize);
        %FIG_INDEX=['fig4_2_reward_center_corr_none'];save_fig(FIG_INDEX,ifsavefig);
        end

        % below for ripple
        for r=1%:2
        ind=ind0 & squeeze(min(deseq_ripple_mean(:,:,r,5),[],2))>=n_cutoff;
        nn(2+r,1,i,h)=sum(ind0&squeeze(deseq_ripple_mean(:,1,r,5))>=n_cutoff);
        nn(2+r,2,i,h)=sum(ind0&squeeze(deseq_ripple_mean(:,2,r,5))>=n_cutoff);
        nn(2+r,3,i,h)=sum(ind);
        nn(2+r,6,i,h)=sum(ind0);
        a=squeeze(deseq_ripple_mean(ind,1:2,r,h));
        if sum(isnan(a(:))) 
            'error none!'
        end
        if ifbin
        [~,I]=sort(a(:,1));
        a=a(I,:);
        a(:,3)=ceil(binum*[1:size(a,1)]/size(a,1));
        ab=[];
        for j=1:binum
            indj=a(:,3)==j;
            ab(j,:)=mean(a(indj,1:2),1);
        end
        else
            ab=a;
        end
        figure(fign);subplot(pn1,pn2,13+r);
        if i==2
            hold on;plot(ab(:,1),ab(:,2),'o','MarkerFaceColor',plotcolors5(3,:),'MarkerEdgeColor','k');
        else
            hold on;plot(ab(:,1),ab(:,2),'ko');
        end
        [rh,pv]=corr(a(:,1),a(:,2),'type','Spearman');%[rh,pv]=corrcoef(ab(:,1),ab(:,2));
        nn(2+r,4:5,i,h)=[rh,pv];
        rp=squeeze(nn(2+r,4:5,:,h));
        title([rname{r},': R,P= ',num2str([rp(:)'],2)]);
        xlabel({['prob at reward ',num2str(squeeze(nn(2+r,[6,1:3],1,h)))];[num2str(squeeze(nn(2+r,[6,1:3],2,h)))]});
        ylabel([hname{h},' prob at center']);
        legend('Free phase','Forced phase');
        if 0%ii==2
        Y_Limits=(ylim);ylim([-0.05 ceil(Y_Limits(2)*10)/10]);
        X_Limits=(xlim);xlim([-0.1 ceil(X_Limits(2)*11)/10]); 
        set(findall(gcf, 'Type', 'Line'),'MarkerSize', markersize);
        %FIG_INDEX=['fig4_2_reward_center_corr_ripple'];save_fig(FIG_INDEX,ifsavefig);
        end
        end
    end
end
% below for dot plot remote arm coding proportion
xdif=0.2;n1=0.08;d=1;n=16;
figure(fign);subplot(pn1,pn2,15);corrp=[];
for b=1:2
    f1=squeeze(theta_dist16(:,:,b,3));
    f2=ones(16,1)*[1 2 3 4];
    [rho,pv]=corrcoef(f1(:),f2(:));
    corrp(b,:)=[rho(2),pv(2)];
end
for b=1:2
    for v=1:4
        hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)-xdif,squeeze(theta_dist16(:,v,b,3)),'k.');
    end
    hold on;plot((d*[1:4]+d*(b-1)*4.5)-xdif,squeeze(nanmean(theta_dist16(:,:,b,3),1)),'-','Color',[1 1 1]/2);
end
v=4;
for b=3:4
    hold on;plot((d*(b-2)*1.5+d*9-.5)*ones(n,1)+n1*randn(n,1)-xdif,squeeze(theta_dist16(:,v,b,3)),'k.');a=squeeze(nanmean(theta_dist16(:,v,b,3)));
    hold on;plot(d*[-0.15 0.15]+d*9+d*(b-2)*1.5-xdif-.5,a*[1 1],'Color',[1 1 1]/2);
end
xlim([0 12]);xticks(d*[1:4,5.5:8.5,10,11.5]);ylabel('Remote arm coding proportion');
xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: >10','out: >10'});
title(num2str(corrp,3));

if ifsavedata
figdata=[squeeze(theta_dist16(:,:,1,3)),squeeze(theta_dist16(:,:,2,3)),squeeze(theta_dist16(:,4,3,3)),squeeze(theta_dist16(:,4,4,3))];
T = array2table(figdata,'VariableNames',{'Reward_theta 0-1cm/s','Reward_theta 1-5cm/s','Reward_theta 5-10cm/s','Reward_theta >10cm/s',...
    'Center_theta 0-1cm/s','Center_theta 1-5cm/s','Center_theta 5-10cm/s','Center_theta >10cm/s','Inbound run','Outbound run'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 3i']);
end

% below for dot plot center coding proportion
xdif=0.2;n1=0.08;d=1;n=16;
figure(fign);subplot(pn1,pn2,16);comp=[];stats=[];
for b=1:2
    f0=squeeze(ripple_percent16(:,1,b,1));
    for v=1:3
        f1=squeeze(theta_dist16(:,v,b,1));
        [p,~,stat]=ranksum(f0,f1,"Tail","right",'method', 'approximate');%signrank(f1-f0);
        comp(b,v)=p;stats(b,v)=stat.zval;
        hold on;
        if b==1
            plot(d*v*ones(n,1)+n1*randn(n,1)-xdif,f1,'k.');
            %hold on;plot(d*v*ones(n,1)-xdif+[-0.1 0.1],mean(f1)*[1 1],'-');
        else
            plot(d*v*ones(n,1)+n1*randn(n,1)+xdif,f1,'.','color',defaultcolor4(1,:));
            %hold on;plot(d*v*ones(n,1)+xdif+[-0.1 0.1],mean(f1)*[1 1],'-','color',defaultcolor4(2,:));
        end
    end
    hold on;
    if b==1
        plot(d*4.5*ones(n,1)+n1*randn(n,1)-xdif,f0,'k.');
        %hold on;plot(d*5*ones(n,1)-xdif+[-0.1 0.1],mean(f0)*[1 1],'k-');
    else
        plot(d*4.5*ones(n,1)+n1*randn(n,1)+xdif,f0,'.','color',defaultcolor4(1,:));
        %hold on;plot(d*5*ones(n,1)+xdif+[-0.1 0.1],mean(f0)*[1 1],'-','color',defaultcolor4(2,:));
    end
end
b=1;hold on;plot(d*[1:3,4.5]-xdif,[squeeze(mean(theta_dist16(:,1:3,b,1),1)),nanmean(ripple_percent16(:,1,b,1))],'-','Color',[1 1 1]/2);
b=2;hold on;plot(d*[1:3,4.5]+xdif,[squeeze(mean(theta_dist16(:,1:3,b,1),1)),nanmean(ripple_percent16(:,1,b,1))],'-','Color',defaultcolor4(2,:));
xlim([0 5.5]);xticks(d*[1:3,4.5]);ylabel('Center coding proportion');
xticklabels({'Theta: 0-1 cm/s','Theta: 1-5 cm/s','Theta: 5-10 cm/s','CA1 SWRs'});
xlabel(num2str(stats,2));title(num2str(comp,2));

if ifsavedata
figdata=[squeeze(theta_dist16(:,1:3,1,1)),squeeze(ripple_percent16(:,1,1,1)),squeeze(theta_dist16(:,1:3,2,1)),squeeze(ripple_percent16(:,1,2,1))];
T = array2table(figdata,'VariableNames',{'Reward_theta 0-1cm/s','Reward_theta 1-5cm/s','Reward_theta 5-10cm/s','Reward_SWR',...
    'Center_theta 0-1cm/s','Center_theta 1-5cm/s','Center_theta 5-10cm/s','Center_SWR'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4l']);
end
figure(fign);figure_title=['Figure4_Jan2025'];save_current_figure(figure_title);
% -----------------------------------------------------------------------------
if 0
if 0
% below fig4_sup for last arm & current arm
fign=6;pn1=2;pn2=3;errbar_width=2;bar_pos=0.14;
% 'fig5sup_last_current_center'
h=3;hname={'next','current','last','center'}; % next arm
for b=1:2
    figure(fign);subplot(pn1,pn2,b); % for theta
    hold on;ba=bar([1:3],squeeze(mu_theta(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_theta(b,1,1:3,h)),squeeze(err_theta(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_theta(b,2,1:3,h)),squeeze(err_theta(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_theta(b,1:2,1:3,h);
    title({[bname{b},' theta : p = ',num2str(a(:)',2)];['p2 = ',num2str(ps2_theta(b,1:3,h),2)];});% ['n = ',num2str(n(:)')];
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});legend('Forced','Free');xlim([0.5 3.5]);
    ylabel('P(last > shuffle)-P(last < shuffle)');
    FIG_INDEX=['fig4sup_last_current_center_theta_',bname{b}];save_fig(FIG_INDEX,ifsavefig);

    figure(fign);subplot(pn1,pn2,b+2);  % for none
    hold on;ba=bar([1:3],squeeze(mu_none(b,1:2,1:3,h)),'FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_none(b,1,1:3,h)),squeeze(err_none(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_none(b,2,1:3,h)),squeeze(err_none(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_none(b,1:2,1:3,h);
    title({[bname{b},' none : p = ',num2str(a(:)',2)];['p2 = ',num2str(ps2_none(b,:,h),2)];});% ['n = ',num2str(n(:)')];
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});legend('Forced','Free');
    ylabel('P(last > shuffle)-P(last < shuffle)');xlim([0 3.5]);
    FIG_INDEX=['fig4sup_last_current_center_none_',bname{b}];save_fig(FIG_INDEX,ifsavefig);
    
    figure(fign);subplot(pn1,pn2,b+4);  % for ripple
    hold on;ba=bar(1,squeeze(mu_ripple(b,1:2,1,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar(1-bar_pos,squeeze(mu_ripple(b,1,1,h)),squeeze(err_ripple(b,1,1,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar(1+bar_pos,squeeze(mu_ripple(b,2,1,h)),squeeze(err_ripple(b,2,1,h)),'k.','LineWidth',errbar_width);
    ps_r1=squeeze(ps_ripple(b,1:2,1,h))';
    title({[bname{b},' ripple p = ',num2str(ps_r1(:)',2)];['p2 = ',num2str(ps2_ripple(b,1,h),2)];});
    xticks(1);xticklabels({'CA1 ripples'});
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');
    FIG_INDEX=['fig4sup_last_current_center_ripple_',bname{b}];save_fig(FIG_INDEX,ifsavefig);
end
if 0
nn=[];binum=20;ifbin=1;n_cutoff=7;
for ii=1:2
    i=3-ii;
    if i==1 % correct forced visit
        indi= pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60 ;
    elseif i==2 % free visit
        indi= (pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3;
    end
    %indi=indi & pksetsall(:,3)>=reward_dis_cut ;%& pksetsall(:,19)<30;
    for h=2:4
        if h==1
            ind0=indi & pksetsall3(:,44)>0; % next
        elseif h==3
            ind0=indi & pksetsall3(:,43)>0; % last
        else
            ind0=indi;
        end
        % below for theta
        ind=ind0 & squeeze(min(deseq_theta_mean(:,:,5),[],2))>=n_cutoff;
        nn(1,1,i,h)=sum(ind0&squeeze(deseq_theta_mean(:,1,5))>=n_cutoff);
        nn(1,2,i,h)=sum(ind0&squeeze(deseq_theta_mean(:,2,5))>=n_cutoff);
        nn(1,3,i,h)=sum(ind);
        a=squeeze(deseq_theta_mean(ind,1:2,h));
        if sum(isnan(a(:))) 
            'error theta!'
        end
        if ifbin
        [~,I]=sort(a(:,1));
        a=a(I,:);
        a(:,3)=ceil(binum*[1:size(a,1)]/size(a,1));
        ab=[];
        for j=1:binum
            indj=a(:,3)==j;
            ab(j,:)=mean(a(indj,1:2),1);
        end
        else
            ab=a;
        end
        figure(fign);subplot(pn1,pn2,7+(h-2)*4);
        if i==2
            hold on;plot(ab(:,1),ab(:,2),'o','MarkerFaceColor',plotcolors5(3,:),'MarkerEdgeColor','k');
        else
            hold on;plot(ab(:,1),ab(:,2),'ko');
        end
        [rh,pv]=corr(a(:,1),a(:,2),'type','Spearman');%[rh,pv]=corrcoef(ab(:,1),ab(:,2));
        nn(1,4:5,i,h)=[rh,pv];
        rp=squeeze(nn(1,4:5,:,h));
        title(['Theta: R,P= ',num2str([rp(:)'],2)]);
        %a=histcounts2(ab(:,1),ab(:,2),histn);imagesc(a');set(gca,'YDir','normal');
        xlabel(['prob at reward ',num2str([squeeze(nn(1,1:3,1,h)),squeeze(nn(1,1:3,2,h))])]);ylabel([hname{h},' prob at center']);
        if ii==2
            Y_Limits=(ylim);ylim([-0.05 ceil(Y_Limits(2)*10)/10]);
            X_Limits=(xlim);xlim([-0.1 ceil(X_Limits(2)*10)/10]);
            FIG_INDEX=['fig4sup_last_current_center_',num2str(7+(h-2)*4)];save_fig(FIG_INDEX,ifsavefig);
        end

        % below for none
        ind=ind0 & squeeze(min(deseq_none_mean(:,:,5),[],2))>=n_cutoff;
        nn(2,1,i,h)=sum(ind0&squeeze(deseq_none_mean(:,1,5))>=n_cutoff);
        nn(2,2,i,h)=sum(ind0&squeeze(deseq_none_mean(:,2,5))>=n_cutoff);
        nn(2,3,i,h)=sum(ind);
        a=squeeze(deseq_none_mean(ind,1:2,h));
        if sum(isnan(a(:))) 
            'error none!'
        end
        if ifbin
        [~,I]=sort(a(:,1));
        a=a(I,:);
        a(:,3)=ceil(binum*[1:size(a,1)]/size(a,1));
        ab=[];
        for j=1:binum
            indj=a(:,3)==j;
            ab(j,:)=mean(a(indj,1:2),1);
        end
        else
            ab=a;
        end
        figure(fign);subplot(pn1,pn2,8+(h-2)*4);
        if i==2
            hold on;plot(ab(:,1),ab(:,2),'o','MarkerFaceColor',plotcolors5(3,:),'MarkerEdgeColor','k');
        else
            hold on;plot(ab(:,1),ab(:,2),'ko');
        end
        [rh,pv]=corr(a(:,1),a(:,2),'type','Spearman');%[rh,pv]=corrcoef(ab(:,1),ab(:,2));
        nn(2,4:5,i,h)=[rh,pv];
        rp=squeeze(nn(2,4:5,:,h));
        title(['None: R,P= ',num2str([rp(:)'],2)]);
        xlabel(['prob at reward ',num2str([squeeze(nn(2,1:3,1,h)),squeeze(nn(2,1:3,2,h))])]);ylabel([hname{h},' prob at center']);
        if ii==2
            Y_Limits=(ylim);ylim([-0.05 ceil(Y_Limits(2)*10)/10]);
            X_Limits=(xlim);xlim([-0.1 ceil(X_Limits(2)*10)/10]);
            FIG_INDEX=['fig4sup_last_current_center_',num2str(8+(h-2)*4)];save_fig(FIG_INDEX,ifsavefig);
        end

        % below for ripple
        for r=1%:2
        ind=ind0 & squeeze(min(deseq_ripple_mean(:,:,r,5),[],2))>=n_cutoff;
        nn(2+r,1,i,h)=sum(ind0&squeeze(deseq_ripple_mean(:,1,r,5))>=n_cutoff);
        nn(2+r,2,i,h)=sum(ind0&squeeze(deseq_ripple_mean(:,2,r,5))>=n_cutoff);
        nn(2+r,3,i,h)=sum(ind);
        a=squeeze(deseq_ripple_mean(ind,1:2,r,h));
        if sum(isnan(a(:))) 
            'error none!'
        end
        if ifbin
        [~,I]=sort(a(:,1));
        a=a(I,:);
        a(:,3)=ceil(binum*[1:size(a,1)]/size(a,1));
        ab=[];
        for j=1:binum
            indj=a(:,3)==j;
            ab(j,:)=mean(a(indj,1:2),1);
        end
        else
            ab=a;
        end
        figure(fign);subplot(pn1,pn2,8+r+(h-2)*4);
        if i==2
            hold on;plot(ab(:,1),ab(:,2),'o','MarkerFaceColor',plotcolors5(3,:),'MarkerEdgeColor','k');
        else
            hold on;plot(ab(:,1),ab(:,2),'ko');
        end
        [rh,pv]=corr(a(:,1),a(:,2),'type','Spearman');%[rh,pv]=corrcoef(ab(:,1),ab(:,2));
        nn(2+r,4:5,i,h)=[rh,pv];
        rp=squeeze(nn(2+r,4:5,:,h));
        title([rname{r},': R,P= ',num2str([rp(:)'],2)]);
        xlabel(['prob at reward ',num2str([squeeze(nn(2+r,1:3,1,h)),squeeze(nn(2+r,1:3,2,h))])]);ylabel([hname{h},' prob at center']);
        legend('Free phase','Forced phase');
        if ii==2
            Y_Limits=(ylim);ylim([-0.05 ceil(Y_Limits(2)*10)/10]);
            X_Limits=(xlim);xlim([-0.1 ceil(X_Limits(2)*10)/10]);
            FIG_INDEX=['fig4sup_last_current_center_',num2str(8+r+(h-2)*4)];save_fig(FIG_INDEX,ifsavefig);
        end
        end
    end
end
end
figure_title=['fig4sup_last_current_center'];save_current_figure(figure_title);
end
%------------------------------------------------------------------------------------
% fig 4-3, theta duration & others 
% below look at overlap in visits for theta and ripples
overlap_session=[];
for r=1%:4
for h=1:5
for i=1:3
    if i==1 % correct forced visit
        indi= pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60 ;
    elseif i==2 % free visit
        indi= (pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3;
    else
        indi= (pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60) | ((pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3);
    end
    for b=1:2
        if b==1
            indb=pksetsall(:,3)>=reward_dis_cut;
        else
            indb=1;
        end
        % below for theta
        for in=1:16
            indin= pksetsall3(:,45)==in & indi & indb ; 
            overlap_session(in,b,i,1:4,h,r)=[sum(indin),sum(deseq_theta_mean(indin,b,h)>0),sum(deseq_ripple_mean(indin,b,r,h)>0),sum(deseq_theta_mean(indin,b,h)>0 & deseq_ripple_mean(indin,b,r,h)>0)];
            %overlap_session(in,b,i,1:4,h)=[sum(indin),sum(deseq_theta_mean(indin,b,h)>0),sum(deseq_none_mean(indin,b,h)>0),sum(deseq_theta_mean(indin,b,h)>0 & deseq_none_mean(indin,b,h)>0)];
        end
    end
end
end
end
overlap_session(:,:,:,5,:,:)=overlap_session(:,:,:,2,:,:).*overlap_session(:,:,:,3,:,:)./overlap_session(:,:,:,1,:,:);
overlap_session(:,:,:,4,:,:)=overlap_session(:,:,:,4,:,:)./overlap_session(:,:,:,1,:,:);
overlap_session(:,:,:,5,:,:)=overlap_session(:,:,:,5,:,:)./overlap_session(:,:,:,1,:,:);

overlap_session=overlap_session(:,:,:,:,[1,5],:);
hname2={'center','next','last','current','unvisited','visited'};
fign=7;pn1=3;pn2=5;errbar_width=2;bar_pos=0.14;
for b=1:2
    for h=1:2
        for m=1:2
            mu=[];err=[];
            for i=1:2
                a=[];
                for k=1:3
                    a=[a;squeeze(session_rep_n(b,i,k,:,h,:,m))];
                end
                mu(:,i)=nanmean(a,2);
                err(:,i)=nanstd(a,0,2)/4;
            end
            figure(fign);subplot(pn1,pn2,b+(m-1)*2+(h-1)*4);indp=[1:3,5];indpp=[1:3,4.5];
            hold on;ba=bar(indpp,mu(indp,:),'FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
            ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
            hold on;errorbar(indpp-bar_pos,mu(indp,1),err(indp,1),'k.','LineWidth',errbar_width);
            hold on;errorbar(indpp+bar_pos,mu(indp,2),err(indp,2),'k.','LineWidth',errbar_width);
            title([bname{b},' ',hname2{h}]);legend('Forced','Free');xlim([0 5.5]);
            xticks(indpp);xticklabels({'theta: 0-1','1-5','5-10','CA1 ripples'});
            if m==1
                ylabel(['Proportion of all ',bname{b}]);
            else
                ylabel(['Proportion of coding ',bname{b}]);
            end
            %FIG_INDEX=['fig4_3_',bname{b},'_',hname2{h},'_',num2str(m)];save_fig(FIG_INDEX,ifsavefig);
        end
    end
end
% below save data
if 0
m=1;h=1;b=1;      
a=[];
for i=1:2
    a=[a,squeeze(session_rep_n(b,i,1,1:3,h,:,m))',squeeze(session_rep_n(b,i,2,1,h,:,m))'];
end
T = array2table(a,'VariableNames',{'Forced_theta 0-1cm/s','Forced_theta 1-5cm/s','Forced_theta 5-10cm/s','Forced_CA1 SWRs','Choice_theta 0-1cm/s','Choice_theta 1-5cm/s','Choice_theta 5-10cm/s','Choice_CA1 SWRs'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4l']);
end
m=1;h=2;b=1;      
a=[];
for i=1:2
    a=[a,squeeze(session_rep_n(b,i,1,1:3,h,:,m))',squeeze(session_rep_n(b,i,2,1,h,:,m))];
end
T = array2table(a,'VariableNames',{'Reward_Forced_theta 0-1cm/s','Reward_Forced_theta 1-5cm/s','Reward_Forced_theta 5-10cm/s',...
    'Reward_Forced_CA1 SWRs','Reward_Choice_theta 0-1cm/s','Reward_Choice_theta 1-5cm/s','Reward_Choice_theta 5-10cm/s','Reward_Choice_CA1 SWRs'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4m'],'Range','A1');
m=1;h=2;b=2;      
a=[];
for i=1:2
    a=[a,squeeze(session_rep_n(b,i,1,1:3,h,:,m))',squeeze(session_rep_n(b,i,2,1,h,:,m))];
end
T = array2table(a,'VariableNames',{'Center_Forced_theta 0-1cm/s','Center_Forced_theta 1-5cm/s','Center_Forced_theta 5-10cm/s',...
    'Center_Forced_CA1 SWRs','Center_Choice_theta 0-1cm/s','Center_Choice_theta 1-5cm/s','Center_Choice_theta 5-10cm/s','Center_Choice_CA1 SWRs'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4m'],'Range','K1');
                    
hname3={'next','maze'};%'current','last','center',
ps=[];r=1;
for b=1:2
    for i=3
        figure(fign);subplot(pn1,pn2,b+8);
        for h=1:2
            a=squeeze(overlap_session(:,b,i,4,h,r)-overlap_session(:,b,i,5,h,r));
            hold on;plot(h*ones(16,1),a,'ko');
            ps(b,i,h)=signrank(a);
        end
        title([bname{b}, ' P = ',num2str(squeeze(ps(b,i,:))',2)]);xlim([0.5 2.5]);xticks([1:2]);xticklabels(hname3);
        Y_Limits=(ylim);ylim([-0.01 ceil(Y_Limits(2)*10)/10]);
        ylabel('Proportion - chance');
        %FIG_INDEX=['fig4_3_overlap_',bname{b}];save_fig(FIG_INDEX,ifsavefig);
    end
end
% below compare theta & ripple duration
pksetsall(:,7)=round(pksetsall(:,7));
thetarippledur_session=[];
for i=1:3
    if i==1 % correct forced visit
        indi= pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60 ;
    elseif i==2 % free visit
        indi= (pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3;
    else
        indi= (pksetsall(:,7)>0 & pksetsall(:,7)<=3 & pksetsall(:,16)<=60) | ((pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3) ;
    end
    for b=1:2
        % below for theta
        for in=1:16
            indin= pksetsall3(:,45)==in & indi ; 
            thetarippledur_session(in,b,i,1:4)=mean(pk_behaveall(indin,b,:,2),1);
            thetarippledur_session(in,b,i,5:8)=mean(sum(pk_behaveall(indin,b,1:3,18:21),3),1);
        end
    end
end
thetaripple_mu=[];thetaripple_err=[];
for b=1:2
    for i=1:3
        thetaripple_mu(b,i,:)=mean(thetarippledur_session(:,b,i,:),1);
        thetaripple_err(b,i,:)=std(thetarippledur_session(:,b,i,:),0,1)/4;
    end
end
for b=1:2
    figure(fign);subplot(pn1,pn2,10+b);indp=[1:3,5];
    hold on;ba=bar(indpp,squeeze(thetaripple_mu(b,:,indp))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar(indpp-bar_pos,squeeze(thetaripple_mu(b,1,indp)),squeeze(thetaripple_err(b,1,indp)),'k.','LineWidth',errbar_width);
    hold on;errorbar(indpp+bar_pos,squeeze(thetaripple_mu(b,2,indp)),squeeze(thetaripple_err(b,2,indp)),'k.','LineWidth',errbar_width);
    title([bname{b}]);legend('Forced','Free');
    xticks(indpp);xticklabels({'theta: 0-1','1-5','5-10','CA1 ripples'});
    ylabel(['Duration per visit (s)']);xlim([0 5.5]);
    %FIG_INDEX=['fig4_3_codingduration_',bname{b}];save_fig(FIG_INDEX,ifsavefig);
end
% below save data
b=1;a=[];
for i=1:2
    a=[a,squeeze(thetarippledur_session(:,b,i,indp))];
end
T = array2table(a,'VariableNames',{'Reward_Forced_theta 0-1cm/s','Reward_Forced_theta 1-5cm/s','Reward_Forced_theta 5-10cm/s',...
    'Reward_Forced_CA1 SWRs','Reward_Choice_theta 0-1cm/s','Reward_Choice_theta 1-5cm/s','Reward_Choice_theta 5-10cm/s','Reward_Choice_CA1 SWRs'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4n'],'Range','A1');
b=2;a=[];
for i=1:2
    a=[a,squeeze(thetarippledur_session(:,b,i,indp))];
end
T = array2table(a,'VariableNames',{'Center_Forced_theta 0-1cm/s','Center_Forced_theta 1-5cm/s','Center_Forced_theta 5-10cm/s',...
    'Center_Forced_CA1 SWRs','Center_Choice_theta 0-1cm/s','Center_Choice_theta 1-5cm/s','Center_Choice_theta 5-10cm/s','Center_Choice_CA1 SWRs'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4n'],'Range','K1'); 

xdif=0.14;n1=0.05;d=1;n=16;errbar_width=2;bar_pos=0.14;
% below for next coding proportion
figure(fign);subplot(pn1,pn2,13);i=3;m=1;h=2;comp=[];indpp=[1:3,4.5];
mu=[];err=[];
for b=1:2
    a=[squeeze(session_rep_n(b,i,1,1:3,h,:,m));squeeze(session_rep_n(b,i,2,1,h,:,m))'];
    mu(:,b)=nanmean(a,2);
    err(:,b)=nanstd(a,0,2)/4;
end
hold on;ba=bar(indpp,mu,'FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
ba(1).FaceColor =[1 1 1];ba(2).FaceColor =defaultcolor4(1,:);ba(2).FaceAlpha=0.3;
hold on;errorbar(indpp-bar_pos,mu(:,1),err(:,1),'k.','LineWidth',errbar_width);
hold on;errorbar(indpp+bar_pos,mu(:,2),err(:,2),'k.','LineWidth',errbar_width);
dat1=[];
for b=1:2
    f0=squeeze(session_rep_n(b,i,2,1,h,:,m));
    for v=1:3
        f1=squeeze(session_rep_n(b,i,1,v,h,:,m));
        p=ranksum(f1,f0,"Tail","right");%p=signrank(f1-f0,0,"Tail","right");%p=signrank(f1-f0);%[~,p]=ttest(f1-f0,0,"Tail","right");%
        comp(b,v)=p;dat1(:,v,b)=f1-f0;
        hold on;
        if b==1
            plot(d*v*ones(n,1)+n1*randn(n,1)-xdif,f1,'k.');
            %hold on;plot(d*v*ones(n,1)-xdif+[-0.1 0.1],mean(f1)*[1 1],'-');
        else
            plot(d*v*ones(n,1)+n1*randn(n,1)+xdif,f1,'.','color',defaultcolor4(1,:));
            %hold on;plot(d*v*ones(n,1)+xdif+[-0.1 0.1],mean(f1)*[1 1],'-','color',defaultcolor4(2,:));
        end
    end
    hold on;
    if b==1
        plot(d*4.5*ones(n,1)+n1*randn(n,1)-xdif,f0,'k.');
        %hold on;plot(d*5*ones(n,1)-xdif+[-0.1 0.1],mean(f0)*[1 1],'k-');
    else
        plot(d*4.5*ones(n,1)+n1*randn(n,1)+xdif,f0,'.','color',defaultcolor4(1,:));
        %hold on;plot(d*5*ones(n,1)+xdif+[-0.1 0.1],mean(f0)*[1 1],'-','color',defaultcolor4(2,:));
    end
end
xlim([0 5.5]);xticks(d*[1:3,4.5]);ylabel('Next arm coding proportion of visits');
xticklabels({'Theta: 0-1 cm/s','Theta: 1-5 cm/s','Theta: 5-10 cm/s','CA1 SWRs'});
title(num2str(comp,2));
figdata=[squeeze(session_rep_n(1,i,1,1:3,h,:,m))',squeeze(session_rep_n(1,i,2,1,h,:,m)),squeeze(session_rep_n(2,i,1,1:3,h,:,m))',squeeze(session_rep_n(2,i,2,1,h,:,m))];
T = array2table(figdata,'VariableNames',{'Reward_theta 0-1cm/s','Reward_theta 1-5cm/s','Reward_theta 5-10cm/s','Reward_SWR',...
    'Center_theta 0-1cm/s','Center_theta 1-5cm/s','Center_theta 5-10cm/s','Center_SWR'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4m']);

% below for duration
figure(fign);subplot(pn1,pn2,14);i=3;comp=[];indp=[1:3,5];
ba=bar(indpp,squeeze(thetaripple_mu(:,3,indp))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
ba(1).FaceColor =[1 1 1];ba(2).FaceColor =defaultcolor4(1,:);ba(2).FaceAlpha=0.3;
hold on;errorbar(indpp-bar_pos,squeeze(thetaripple_mu(1,3,indp)),squeeze(thetaripple_err(1,3,indp)),'k.','LineWidth',errbar_width);
hold on;errorbar(indpp+bar_pos,squeeze(thetaripple_mu(2,3,indp)),squeeze(thetaripple_err(2,3,indp)),'k.','LineWidth',errbar_width);
dat=[];
for b=1:2
    f0=squeeze(thetarippledur_session(:,b,i,5));
    for v=1:3
        f1=squeeze(thetarippledur_session(:,b,i,v));
        p=ranksum(f1,f0,"Tail","right");%p=signrank(f1-f0,0,"Tail","right");%[~,p]=ttest(f1-f0,0,"Tail","right");%
        comp(b,v)=p;dat(:,v,b)=f1-f0;
        hold on;
        if b==1
            plot(d*v*ones(n,1)+n1*randn(n,1)-xdif,f1,'k.');
            %hold on;plot(d*v*ones(n,1)-xdif+[-0.1 0.1],mean(f1)*[1 1],'-');
        else
            plot(d*v*ones(n,1)+n1*randn(n,1)+xdif,f1,'.','color',defaultcolor4(1,:));
            %hold on;plot(d*v*ones(n,1)+xdif+[-0.1 0.1],mean(f1)*[1 1],'-','color',defaultcolor4(2,:));
        end
    end    
    hold on;
    if b==1
        plot(d*4.5*ones(n,1)+n1*randn(n,1)-xdif,f0,'k.');
        %hold on;plot(d*5*ones(n,1)-xdif+[-0.1 0.1],mean(f0)*[1 1],'k-');
    else
        plot(d*4.5*ones(n,1)+n1*randn(n,1)+xdif,f0,'.','color',defaultcolor4(1,:));
        %hold on;plot(d*5*ones(n,1)+xdif+[-0.1 0.1],mean(f0)*[1 1],'-','color',defaultcolor4(2,:));
    end
end
ylabel('Duration per visit (s)');title(num2str(comp,2));xlim([0 6]);
xticks([1:3,4.5]);xticklabels({'Theta: 0-1 cm/s','Theta: 1-5 cm/s','Theta: 5-10 cm/s','CA1 SWRs'});
legend('Rat at reward','Rat at center');

figdata=[squeeze(thetarippledur_session(:,1,i,[1:3,5])),squeeze(thetarippledur_session(:,2,i,[1:3,5]))];
T = array2table(figdata,'VariableNames',{'Reward_theta 0-1cm/s','Reward_theta 1-5cm/s','Reward_theta 5-10cm/s','Reward_SWR',...
    'Center_theta 0-1cm/s','Center_theta 1-5cm/s','Center_theta 5-10cm/s','Center_SWR'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 4n']);
figure(fign);figure_title='Figure4_pkN_theta_ripple';save_current_figure(figure_title);
%--------------------------------------------------------------------------------------------------
% sup figure 11
% below for last arm prediction
fign=10;pn1=2;pn2=3;h=3;hname={'next','current','last','center'}; 
for b=1:2
    figure(fign);subplot(pn1,pn2,b); % for theta
    hold on;ba=bar([1:3],squeeze(mu_theta(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_theta(b,1,1:3,h)),squeeze(err_theta(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_theta(b,2,1:3,h)),squeeze(err_theta(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_theta(b,1:2,1:3,h);ar=squeeze(ps_theta_r(b,1:2,1:3,h))';tr=squeeze(tstat_theta(b,1:2,1:3,h))';
   % title({ ['p:',num2str(a(:)',2)];['p2:',num2str(ps2_theta(b,1:3,h),2)];['p:',num2str(ar(:)',2)]});
   title({ ['t:',num2str(tr(:)',2)];['p:',num2str(ar(:)',2)]});
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel([bname{b},' theta']);xlim([0.5 3.5]);
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');
    %FIG_INDEX=['fig4_2_next_theta_',bname{b}];save_fig(FIG_INDEX,ifsavefig,[],[],1);

    figure(fign);subplot(pn1,pn2,b+2);  % for none
    hold on;ba=bar([1:3],squeeze(mu_none(b,1:2,1:3,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-bar_pos,squeeze(mu_none(b,1,1:3,h)),squeeze(err_none(b,1,1:3,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:3]+bar_pos,squeeze(mu_none(b,2,1:3,h)),squeeze(err_none(b,2,1:3,h)),'k.','LineWidth',errbar_width);
    a=ps_none(b,1:2,1:3,h);n=ns_none(b,:,1:3,h);ar=squeeze(ps_none_r(b,1:2,1:3,h))';tr=squeeze(tstat_none(b,1:2,1:3,h))';
    %title({['p:',num2str(a(:)',2)];['p2:',num2str(ps2_none(b,:,h),2)];['p:',num2str(ar(:)',2)]});
    title({ ['t:',num2str(tr(:)',2)];['p:',num2str(ar(:)',2)]});
    xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel([bname{b},' none']);xlim([0 3.5]);
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');
    %FIG_INDEX=['fig4_2_next_none_',bname{b}];save_fig(FIG_INDEX,ifsavefig,[],[],1);
    
    figure(fign);subplot(pn1,pn2,b+4);  % for ripple
    hold on;ba=bar(1,squeeze(mu_ripple(b,1:2,1,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar(1-bar_pos,squeeze(mu_ripple(b,1,1,h)),squeeze(err_ripple(b,1,1,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar(1+bar_pos,squeeze(mu_ripple(b,2,1,h)),squeeze(err_ripple(b,2,1,h)),'k.','LineWidth',errbar_width);
    ps_r1=squeeze(ps_ripple(b,1:2,1,h))';ar=squeeze(ps_ripple_r(b,1:2,1,h));tr=squeeze(tstat_ripple(b,1:2,1,h));
    %title({['p:',num2str(ps_r1(:)',2)];['p2:',num2str(ps2_ripple(b,1,h),2)];['p:',num2str(ar(:)',2)]});%['n = ',num2str(ns_r1(:)')];
    title({ ['t:',num2str(tr(:)',2)];['p:',num2str(ar(:)',2)]});
    xticks(1);xticklabels({'CA1 ripples'});xlabel([bname{b},' ripple']);
    ylabel('P(last > shuffle)-P(last < shuffle)');legend('Forced','Free');
    %FIG_INDEX=['fig4_2_next_ripple_',bname{b}];save_fig(FIG_INDEX,ifsavefig,[],[],1);
end
figure(fign);figure_title='SupFigure11_last_arm';save_current_figure(figure_title);
%------------------------------------------------------------------------------------
% 'fig4sup_deseqN_'
fign=8;pn1=2;pn2=3;errbar_width=2;bar_pos=0.14;
deseq_num_mu_thetanone=next_shuffle.deseq_num_thetanone.mu;
deseq_num_err_thetanone=next_shuffle.deseq_num_thetanone.err;
deseq_num_ns_thetanone=next_shuffle.deseq_num_thetanone.ns;
deseq_num_mu_ripple=next_shuffle.deseq_num_ripple.mu;
deseq_num_err_ripple=next_shuffle.deseq_num_ripple.err;
deseq_num_ns_ripple=next_shuffle.deseq_num_ripple.ns;
mu_ripple_pk=next_shuffle.ripple_pk.mu;
err_ripple_pk=next_shuffle.ripple_pk.err;
ns_ripple_pk=next_shuffle.ripple_pk.ns;
% 'deseq_n_next_shuffle'
for h=1%:3
for b=1:2
    figure(fign);subplot(pn1,pn2,b); % for theta
    hold on;ba=bar([1:4],squeeze(deseq_num_mu_thetanone(1,b,1:2,:,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:4]-bar_pos,squeeze(deseq_num_mu_thetanone(1,b,1,:,h)),squeeze(deseq_num_err_thetanone(1,b,1,:,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:4]+bar_pos,squeeze(deseq_num_mu_thetanone(1,b,2,:,h)),squeeze(deseq_num_err_thetanone(1,b,2,:,h)),'k.','LineWidth',errbar_width);
    n=deseq_num_ns_thetanone(1,b,1:2,:,h);
    title({[bname{b},' theta : n = ',num2str(n(:)')]});
    xticks([1:4]);xticklabels({'0-1','1-5','5-10','>10'});legend('Forced','Free');
    ylabel([hname{h},': deseq N']);
    FIG_INDEX=['fig4sup_deseqN_',hname{h},'_theta_',bname{b}];save_fig(FIG_INDEX,ifsavefig);

    figure(fign);subplot(pn1,pn2,b+2); % for none
    hold on;ba=bar([1:4],squeeze(deseq_num_mu_thetanone(2,b,1:2,:,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:4]-bar_pos,squeeze(deseq_num_mu_thetanone(2,b,1,:,h)),squeeze(deseq_num_err_thetanone(2,b,1,:,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:4]+bar_pos,squeeze(deseq_num_mu_thetanone(2,b,2,:,h)),squeeze(deseq_num_err_thetanone(2,b,2,:,h)),'k.','LineWidth',errbar_width);
    n=deseq_num_ns_thetanone(2,b,1:2,:,h);
    title({[bname{b},' none : n = ',num2str(n(:)')]});
    xticks([1:4]);xticklabels({'0-1','1-5','5-10','>10'});legend('Forced','Free');
    ylabel([hname{h},': deseq N']);
    FIG_INDEX=['fig4sup_deseqN_',hname{h},'_none_',bname{b}];save_fig(FIG_INDEX,ifsavefig);

    figure(fign);subplot(pn1,pn2,b+4); % for ripple
    hold on;ba=bar(1,squeeze(deseq_num_mu_ripple(b,1:2,1,h))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar(1-bar_pos,squeeze(deseq_num_mu_ripple(b,1,1,h)),squeeze(deseq_num_err_ripple(b,1,1,h)),'k.','LineWidth',errbar_width);
    hold on;errorbar(1+bar_pos,squeeze(deseq_num_mu_ripple(b,2,1,h)),squeeze(deseq_num_err_ripple(b,2,1,h)),'k.','LineWidth',errbar_width);
    n=squeeze(deseq_num_ns_ripple(b,1:2,1,h));
    title({[bname{b},' ripple : n = ',num2str(n(:)')]});
    xticks(1);xticklabels({'CA1 ripples'});
    legend('Forced','Free');ylabel([hname{h},': deseq N']);
    FIG_INDEX=['fig4sup_deseqN_',hname{h},'_ripple_',bname{b}];save_fig(FIG_INDEX,ifsavefig);
end
figure(fign);figure_title=['fig4sup_deseqN_',hname{h}];save_current_figure(figure_title);
end
%------------------------------------------------------------------------------------
end