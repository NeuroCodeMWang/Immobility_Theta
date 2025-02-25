%function generate_SupFigure_9_Jan25(ifsavedata)
%if nargin<1
    ifsavedata=0;
%end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
cd(datadir);load('Figure3_thetaseq_figure_v4');
spikepropname={'phase onset','max phase','sig seq percent','seq length (cm)','dir, out-in','theta freq','theta power',...
    'major','minor','minor/major','ca3 p'};
bname={'reward','center','in','out'};mname={'local','remote'};

% below statistic test
fign=20;pn1=2;pn2=4;mu=[];xdif=0.2;n1=0.08;c1=0.4;d=1;thetacycleN_cut=50;stds=[];
lettername={'g','f','a','b','c','h','i','d','e','f','g','j','k'};inter=1;fn={'velocity','local/remote'};
bname={'reward','center','in','out'};mname={'local','remote'};
if 0
    
    % sequence length at reward
    b=1;v=1:4;m=1:2;k=4;reps=16;
    anovadata=[];
    for v=1:4
        a=squeeze(spikeprop(1:16,b,v,m,k));
        anovadata=[anovadata;a];
    end
    [p,table] = anova2(anovadata,reps);

    % sig sequence fraction
    b=1:2;v=1:4;m=1:2;ks=[3,4,5,7,6];anovap=[];
    for kk=5%:5
        k=ks(kk);ps=[];
        for b=1:2
            anovadata=[];mlabels=[];vlabels=[];
            for v=1:4
                a=squeeze(spikeprop(1:16,b,v,m,k));
                anovadata=[anovadata;a];
                mlabels=[mlabels;ones(16,1)*[1,2]];
                vlabels=[vlabels;ones(16,2)*v];
            end
            anovadata=anovadata(:);mlabels=mlabels(:);vlabels=vlabels(:);
            ind=~isnan(anovadata);anovadata=anovadata(ind);mlabels=mlabels(ind);vlabels=vlabels(ind);
            p = anovan(anovadata,{vlabels,mlabels},"Model","interaction");
            ps(:,b)=p;
            anovap(b,:,k)=p;
        end
    end

    pss=[];
for k=1%:2
    for b=2%:2%:4
        alpha=[];idp=[];idq=[];
        for m=1:2
            for v=2:4
                ind1=find((~isnan(spikeprop(1:16,b,v,m,k))) & nn_session(1:16,b,v,m)>=thetacycleN_cut);
                n=length(ind1);
                if n>0
                    alpha=[alpha;squeeze(spikeprop(ind1,b,v,m,k))];
                    idp=[idp;ones(n,1)*v];
                    idq=[idq;ones(n,1)*(m)];
                end
            end
        end
        [pval table] = circ_hktest(alpha*pi/180, idp, idq, inter, fn);
        pss(k,b,:)=pval
    end
end

circ_linear_correlation=[];
for k=1:2
    for b=1:2
        for m=1:2 
            alpha=[];idp=[];idq=[];
            for v=1:4
                ind1=find((~isnan(spikeprop(1:16,b,v,m,k))) & nn_session(1:16,b,v,m)>thetacycleN_cut);
                n=length(ind1);
                if n>0
                    alpha=[alpha;squeeze(spikeprop(ind1,b,v,m,k))];
                    idp=[idp;ones(n,1)*v];
                    idq=[idq;ones(n,1)*m];
                end
            end
            [rho pval] = circ_corrcl(alpha*pi/180, idp);
            rho
            pval
            circ_linear_correlation(k,b,m,:)=[rho pval];
        end
    end
end
end

%------------------------------------------------------------------------------
% below to plot figures
thetacycleN_cut=50;ns=[];
for b=1:4
    for m=1:2
        for v=1:4
            for k=1:7%11
                ind1=find((~isnan(spikeprop(1:16,b,v,m,k))) & nn_session(1:16,b,v,m)>=thetacycleN_cut);
                n=length(ind1);
                if n>0
                    if k<3
                        mu(b,v,m,k)=rem(circ_meand(squeeze(spikeprop(ind1,b,v,m,k))),360);
                        [~,~,err]=circ_meand(squeeze(spikeprop(ind1,b,v,m,k)));
                        stds(b,v,m,k)=err;%*n^(-0.5);
                    else
                        mu(b,v,m,k)=squeeze(nanmean(spikeprop(ind1,b,v,m,k),1));
                        stds(b,v,m,k)=std(squeeze(spikeprop(ind1,b,v,m,k)));%*n^(-0.5);
                    end
                    figure(fign);subplot(pn1,pn2,k);
                    if m==1 & b<3
                        figure(fign);subplot(pn1,pn2,k);
                        hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)-xdif,squeeze(spikeprop(ind1,b,v,m,k)),'.','Color',[1 1 1]*c1);
                        if k==2
                            hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)-xdif,squeeze(spikeprop(ind1,b,v,m,k))+360,'.','Color',[1 1 1]*c1);
                        end
                    elseif m==2 & b<3
                        figure(fign);subplot(pn1,pn2,k);
                        hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)+xdif,squeeze(spikeprop(ind1,b,v,m,k)),'.','Color',[0.8500 0.3250 0.0980]);
                        if k==2
                            hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)+xdif,squeeze(spikeprop(ind1,b,v,m,k))+360,'.','color',[0.8500 0.3250 0.0980]);
                        end
                    elseif m==1 & b>=3 & v==4
                        figure(fign);subplot(pn1,pn2,k);
                        hold on;plot((d*10+d*(b-3)*1.5)*ones(n,1)+n1*randn(n,1),squeeze(spikeprop(ind1,b,v,m,k)),'.','Color',[1 1 1]*c1);
                        if k==2
                            hold on;plot((d*10+d*(b-3)*1.5)*ones(n,1)+n1*randn(n,1),squeeze(spikeprop(ind1,b,v,m,k))+360,'.','Color',[1 1 1]*c1);
                        end
                    end
                else
                    mu(b,v,m,k)=nan;stds(b,v,m,k)=nan;
                end
            end
        end
    end
end
k2phase=240;k1phase=20;
for k=1:7%11
    figure(fign);subplot(pn1,pn2,k);
    for b=1:2
        if k==2 
            a=squeeze(mu(b,:,1,k))';a(a<k2phase)=a(a<k2phase)+360;hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'Color',[1 1 1]*c1*0);
            a=squeeze(mu(b,:,2,k))';a(a<k2phase)=a(a<k2phase)+360;hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'color',[0.8500 0.3250 0.0980]);
        elseif k==1 %& b==2
            a=squeeze(mu(b,:,1,k))';a(a<k1phase)=a(a<k1phase)+360;hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'Color',[1 1 1]*c1*0);
            a=squeeze(mu(b,:,2,k))';a(a<k1phase)=a(a<k1phase)+360;hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'color',[0.8500 0.3250 0.0980]);  
            %a=squeeze(mu(b,1:3,2,k))';hold on;plot(d*[1:3]+d*(b-1)*4.5,a,'color',[0.8500 0.3250 0.0980]);
        else
            a=squeeze(mu(b,:,1,k))';hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'Color',[1 1 1]*c1*0);
            a=squeeze(mu(b,:,2,k))';hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'color',[0.8500 0.3250 0.0980]);
        end
        if k==2
            a=squeeze(mu(2+b,4,1,k))';a(a<k2phase)=a(a<k2phase)+360;hold on;plot(d*[-0.15 0.15]+d*10+d*(b-1)*1.5,a*[1 1],'Color',[1 1 1]*c1*0);
        elseif k==1
            a=squeeze(mu(2+b,4,1,k))';a(a<k1phase)=a(a<k1phase)+360;hold on;plot(d*[-0.15 0.15]+d*10+d*(b-1)*1.5,a*[1 1],'Color',[1 1 1]*c1*0);
        else
            a=squeeze(mu(2+b,4,1,k))';hold on;plot(d*[-0.15 0.15]+d*10+d*(b-1)*1.5,a*[1 1],'Color',[1 1 1]*c1*0);
        end
    end
    if k==1
        ylim([0 360]+k1phase);
    elseif k==2
        ylim([0 360]+k2phase);%ylim([180 540]);
    elseif k==3
        ylim([0 1.05]);
    elseif k==4
        ylim([0 45]);
    elseif k==6
        ylim([6 10]);
    elseif k==7
        ylim([0 2.5]);
    elseif k==9
        ylim([-0.55 1.1]);
    end
    xticks(d*[1:4,5.5:8.5,10,11.5]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: >10','out: >10'});
    ylabel(spikepropname{k});xlim(d*[0 12.5]);
    if k>=9
        set(findall(gcf, 'Type', 'Line'),'MarkerSize', 20);
    end
end
figure(fign);figure_title='SupFigure9_Jan2025';save_current_figure(figure_title);
% below to generate heatmap for sup 9
fign=4;pn1=4;pn2=9;ifdist=1;thetacycleN_cut=50;mnum=2;lettername={'c','d','e'};gaustd2=1;
% below mean seq part
for b=1%:2
    for m=1:2
        as=[];bs=[];lfp=[];cens=[];
        for v=1:4
           ind=1:16;%
           if ifdist
               a=imgaussfilt(squeeze(nanmean(deposphase_dist(ind,b,v,m,:,:),1)),gaustd2);
           else
               a=squeeze(nanmean(seqmean_session(ind,b,v,m,:,:)./nansum(seqmean_session(ind,b,v,m,:,:),5),1));
           end
           as=[as,a];
           b1=squeeze(nanmean(firemean_session(ind,b,v,m,:,:,1),1));
           for k=1:6
                b1(k,:)=gaussian_smooth(b1(k,:));
                b1(k,:)=(b1(k,:)-min(b1(k,:)))/(max(b1(k,:))-min(b1(k,:)));
           end
           bs=[bs,b1];
           lfp1=squeeze(nanmean(lfpmean_session(ind,b,v,m,:,:),1));
           lfp=[lfp,lfp1]; 
           cen=ceil(mu(b,v,m,1)/10);%cen=rem(circ_meand(seqcut(1:16,b,v,m,5)),360)/10;
           cens=[cens,[cen,cen+36]+(v-1)*72];
        end
        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2)]+(m-1)*floor(pn2/2)+(b-1)*pn2*4);
        hold on;plot(lfp(1,:)+2,'k');hold on;plot(lfp(4,:),'k');xline([1:4]*72,'k');
        xticks([0:36:4*72]);xlim([0 4*72]);ylabel('zscore lfp');xticks([]);hold on;xline(cens,'c--');%plot(cens,4,'k+');%

        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2)]+(m-1)*floor(pn2/2)+(b-1)*pn2*4+pn2);
        hold on;plot(bs(3,:)+1,'k');%'color',[0.8500 0.3250 0.0980]);
        hold on;plot(bs(4,:),'k');ylim([-0.5 2.5]);
        %hold on;plot(bs(3,:)-bs(4,:));
        xticks([0:36:4*72]);xline([1:4]*72,'k');xlim([0 4*72]);xline(cens,'c--');%xline([36:36:36*4],'k');xline([52:72:72*4],'m');
        ylabel('zscore spike N');yticks([0:3]*5);yticklabels({'CA1 P','CA3 P','CA1 I','CA3 I'});

        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2),1+pn2:floor(pn2/2)+pn2]+(m-1)*floor(pn2/2)+(b-1)*pn2*4+pn2*2);
        imagesc([],[-100 100],as);set(gca,'YDir','normal');
        if b==1
            ylim([-30 30]);
        elseif b==2
            ylim([-40 90]);
        elseif b==3
            ylim([-60 60]);
        elseif b==4
            ylim([-60 60]);
        end
        xticks([0:36:4*72]);xline([1:4]*72,'k');colormap('hot');colorbar;xline(cens,'c--');
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
        if ifsavedata
            figdata=[lfp([1,4],:);bs(3:4,:);as];
            T = array2table(figdata,'RowNames',rowname);
            writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 3',lettername{b},'_',num2str(m)],'WriteRowNames',true);
        end
    end
end
figure(fign);figure_title='SupFigure9_heatmap_Jan2025';save_current_figure(figure_title);


if ifsavedata
% below to save data
bname={'reward','center','in','out'};mname={'local','remote'};
for k=1:7
    figdata=[[squeeze(spikeprop(1:16,1,:,1,k));squeeze(spikeprop(1:16,1,:,2,k))],...
        [squeeze(spikeprop(1:16,2,:,1,k));squeeze(spikeprop(1:16,2,:,2,k))]];
    T = array2table(figdata,'VariableNames',{'reward: 0-1 cm/s','reward: 1-5 cm/s','reward: 5-10 cm/s','reward: > 10 cm/s',...
        'center: 0-1 cm/s','center: 1-5 cm/s','center: 5-10 cm/s','center: > 10 cm/s'},'RowNames',...
        {'local_1','local_2','local_3','local_4','local_5','local_6','local_7','local_8','local_9','local_10',...
        'local_11','local_12','local_13','local_14','local_15','local_16','remote_1','remote_2','remote_3','remote_4',...
        'remote_5','remote_6','remote_7','remote_8','remote_9','remote_10','remote_11','remote_12','remote_13','remote_14','remote_15','remote_16'});
    if k<=2
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 3',lettername{k}],'WriteRowNames',true);
    else
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['SupFigure 9',lettername{k}],'WriteRowNames',true);
    end
    T = array2table(squeeze(spikeprop(1:16,3:4,4,1,k)),'VariableNames',{'Inbound run','Outbound run'});
    if k<=2
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 3',lettername{k}],'Range','K1');
    else
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['SupFigure 9',lettername{k}],'Range','K1');
    end
end
colname=[];vname={'0-1 cm/s','1-5 cm/s','5-10 cm/s','> 10 cm/s'};
for v=1:4
    for i=1:61
        colname{i+(v-1)*61}=[vname{v},'_',num2str(i)];
    end
end
for b=1:2
    figdata=[[squeeze(velmean_session(:,b,1,1,:));squeeze(velmean_session(:,b,2,1,:))],[squeeze(velmean_session(:,b,1,2,:));squeeze(velmean_session(:,b,2,2,:))],...
        [squeeze(velmean_session(:,b,1,3,:));squeeze(velmean_session(:,b,2,3,:))],[squeeze(velmean_session(:,b,1,4,:));squeeze(velmean_session(:,b,2,4,:))]];
    T = array2table(figdata,'VariableNames',colname,'RowNames',...
        {'local_1','local_2','local_3','local_4','local_5','local_6','local_7','local_8','local_9','local_10',...
        'local_11','local_12','local_13','local_14','local_15','local_16','remote_1','remote_2','remote_3','remote_4',...
        'remote_5','remote_6','remote_7','remote_8','remote_9','remote_10','remote_11','remote_12','remote_13','remote_14','remote_15','remote_16'});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['SupFigure 9',lettername{b+11}],'WriteRowNames',true);
end
end