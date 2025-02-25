function generate_SupFigure_11_Jan25(ifsavedata)
%if nargin<1
    ifsavedata=0;
%end

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