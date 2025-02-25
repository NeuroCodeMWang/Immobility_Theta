function theta_peak_peak
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
% 1.below plot lfp of each theta cycle grouped by velocity to check peak
% peak theta cycle
thetacycle_lfp=zeros(4,36,11,16);
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Analyze_Theta_Cycle_Broad_peak_peak','thetacycle','reftheta','lfptime_trial','lfpz_trial');
    for i=1:36
        indi=ceil(reftheta(:,4)/10)==i;
        for v=1:4
            if v==1
                indv=thetacycle(:,10)>=0 & thetacycle(:,10)<1;
            elseif v==2
                indv=thetacycle(:,10)>=1 & thetacycle(:,10)<5;
            elseif v==3
                indv=thetacycle(:,10)>=5 & thetacycle(:,10)<10;
            elseif v==4
                indv=thetacycle(:,10)>=10;
            end
            lia=ismember(lfptime_trial(:,5),find(indv));
            ind= lia & indi;
            thetacycle_lfp(v,i,:,in)=nanmean(lfpz_trial(ind ,:),1);
        end
    end
end
typenum=4;typeindex=[1:typenum]*floor(256/typenum);
plotcolors=colormap(jet);
plotcolors=plotcolors(typeindex,:);
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG IML/GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
figure;
for in=1:16
for i=1:7
    subplot(2,4,i);
    for v=1:4
        a=squeeze(thetacycle_lfp(v,:,i,in));
        a=[a,a];
        hold on;plot([5:10:715],a,'Color',plotcolors(v,:));
    end
    xlabel('Theta Phase');ylabel('Mean LFP');title(HPC_layer_name{i});xlim([0 720]);
end
pause;clf;
end
figure;
for i=1:7
    subplot(2,4,i);
    for v=1:4
        a=squeeze(nanmean(thetacycle_lfp(v,:,i,:),4));
        a=[a,a];
        hold on;plot([5:10:715],a,'Color',plotcolors(v,:));
    end
    xlabel('Theta Phase');ylabel('Mean LFP');title(HPC_layer_name{i});xlim([0 720]);
end
% below first calcualte data for individual sessions for fig4
for in=1:16
    analyze_theta_sequence_at_pause(in,ifdataall);
end

% below quantify thetaseq1
figdir='X:\Mengni\Data_Analysis\Paper_Figures';
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
thetacycleall=[];depos_distall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_dHPC_peak_peak.mat','thetacycle','thetacycle_depos_dist');
    thetacycleall=[thetacycleall;thetacycle];
    depos_distall=cat(1,depos_distall,thetacycle_depos_dist);
end

clear all;close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
figdir='X:\Mengni\Data_Analysis\Paper_Figures';
timewin=0.01;reward_dis_cut=60;%thetaseqshift=zeros(101,36,4,4,3,16);thetaseq=zeros(50,36,4,4,3,16);
min_dewinN=14;ifpeak=0;phasebin=10;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    if ifpeak
    load('ThetaCycle_Decode_Info_dHPC_peak_peak.mat','linear_maze','linear_pos','deseqall','thetacycle','thetacycle_depos_dist');
    load('Analyze_Theta_Cycle_Broad_peak_peak','reftheta');
    else
    load('ThetaCycle_Decode_Info_dHPC2.mat','linear_maze','linear_pos','deseqall','thetacycle','thetacycle_depos_dist');
    load('Analyze_Theta_Cycle_Broad','reftheta');
    end
    load('Ripple_Events_PAPER5.mat','ripples');
    load('Position_Data_Maze.mat');

    ind=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035;
    ripples=ripples(ind,:);
    thetacycle(:,79)=0; deseqall(:,43:44)=0;% associated ripple id
    for i=1:size(ripples,1)
        ind= deseqall(:,1)>=ripples(i,1)-timewin & deseqall(:,1)<=ripples(i,2)+timewin ;
        deseqall(ind,42+ripples(i,6))=i;
        ind= thetacycle(:,6)>=ripples(i,1)-0.05 & thetacycle(:,6)<=ripples(i,2)+0.05 ;
        thetacycle(ind,79)=i;
    end
    ruler=reftheta(:,1);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
    deseqall(:,45:46)=reftheta(outputindex,3:4);
    deseqall(:,47)=linear_pos(:,3);
    deseqall(:,48)=linear_pos(:,3)==linear_pos(:,1); % if local
    [M,I]=max(squeeze(thetacycle_depos_dist(:,1,:))./thetacycle(:,24),[],2); 
    thetacycle(:,59:60)=[I-1,M];
    [M,I]=max(squeeze(thetacycle_depos_dist(:,1,2:9))./thetacycle(:,24),[],2); 
    thetacycle(:,63)=I;
    thetacycle(:,64)=M+squeeze(thetacycle_depos_dist(:,1,1))./thetacycle(:,24);
    thetacycle(:,61)=thetacycle(:,26)==thetacycle(:,59); % if local
    thetacycle(:,62)=nan;
    for b=1:4
        if b==1
            indb= thetacycle(:,27)>reward_dis_cut & thetacycle(:,11)==0;
        elseif b==2
            indb=thetacycle(:,27)<30 & abs(thetacycle(:,11))==0.5;
        elseif b==3
            indb=thetacycle(:,11)==-1; % inbound
        elseif b==4
            indb=thetacycle(:,11)==1; % outbound
        end
        thetacycle(indb,62)=b;
    end
    deseqall(:,27)=ceil(deseqall(:,27)/phasebin);
    despikez=zscore(deseqall(:,[6:7,19:20,23:24]),0,1);
    thetacycle2=thetacycle(thetacycle(:,79)==0 & thetacycle(:,24)>=min_dewinN & thetacycle(:,64)>=0.8 & thetacycle(:,62)>0,:);
    thetaid=thetacycle2(:,21);
    thetaseq=zeros(size(thetaid,1),50,36,2);
    thetaseq2=zeros(size(thetaid,1),101,36);
    thetafire=zeros(size(thetaid,1),6,36);
    lia=ismember(deseqall(:,40),thetaid(:,1));
    deseqall0=deseqall(lia,:);
    despikez0=despikez(lia,:);
    linear_maze0=linear_maze(lia,:,:);
    linear_pos0=linear_pos(lia,:);
    linear_maze_rat0=linearize_8arm_rat_centered(linear_pos0,linear_maze0);
    seqn=zeros(size(thetaid,1),36,2);
    for j=1:size(deseqall0,1)
        id=deseqall0(j,40);
        i=find(thetaid(:,1)==id);
        phase=deseqall0(j,27);           
        a=squeeze(linear_maze0(j,:,thetacycle2(i,63)));
        thetaseq(i,:,phase,1)=thetaseq(i,:,phase,1)+a;
        thetafire(i,:,phase)=thetafire(i,:,phase)+despikez0(j,:);
        a=squeeze(nanmean(linear_maze0(j,:,:),3));
        thetaseq(i,:,phase,2)=thetaseq(i,:,phase,2)+a;
        seqn(i,phase)=seqn(i,phase)+1;
        thetaseq2(i,:,phase)=thetaseq2(i,:,phase)+linear_maze_rat0(j,:);
    end
    for i=1:size(thetaid,1)
        for phase=1:36
            thetafire(i,:,phase)=thetafire(i,:,phase)/seqn(i,phase);
            thetaseq2(i,:,phase)=thetaseq2(i,:,phase)/seqn(i,phase);
            for h=1:2
                thetaseq(i,:,phase,h)=thetaseq(i,:,phase,h)/seqn(i,phase);
            end
        end
    end
    thetaseqslope=[];
    for h=1:3
        if h<3
            thetaseq1=squeeze(thetaseq(:,:,:,h));
        else
            thetaseq1=squeeze(thetaseq2(:,:,:));
        end
        seqslope=quantify_shuffle_thetaseq(thetaseq1);
        seqslope(:,5)=seqslope(:,1); % start location
        seqslope(:,6)=seqslope(:,1)+seqslope(:,2)*360; % end location
        seqslope(:,7)=seqslope(:,6)-seqslope(:,5); % path length
        thetaseqslope(:,:,h)=seqslope;
    end
    if ifpeak
        cd(savedir);save('theta_peak_peak2','thetaseqslope','thetafire','thetaseq','thetaseq2','thetacycle2','thetacycle','deseqall','-v7.3');
    else
        cd(savedir);save('theta_trough_trough','thetaseqslope','thetafire','thetaseq','thetaseq2','thetacycle2','thetacycle','deseqall','-v7.3');
    end
    in
end
ifpeak=0;
thetacycle2all=[];thetafireall=[];thetaseqslopeall=[];thetaseqall=[];thetacyclelfpall=[];%thetaseq2all=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    if ifpeak
        load('theta_peak_peak2');
    else
        load('theta_trough_trough');
    end
    thetacycle2(:,65)=in;
    thetacycle2all=[thetacycle2all;thetacycle2];
    thetafireall=cat(1,thetafireall,thetafire);
    thetaseqslopeall=cat(1,thetaseqslopeall,thetaseqslope);
    thetaseqall=cat(1,thetaseqall,thetaseq);
    %thetaseq2all=cat(1,thetaseq2all,thetaseq2);
    load('Analyze_Theta_Cycle_Broad','reftheta','lfptime_trial','lfpz_trial');
    lia=ismember(lfptime_trial(:,5),thetacycle2(:,13));
    lfp=[lfptime_trial(lia,5),ceil(reftheta(lia,4)/10),lfpz_trial(lia,1:7)];
    thetacycle_lfp=zeros(size(thetacycle2,1),36,7);nn=zeros(size(thetacycle2,1),36);
    for i=1:size(lfp,1)
        ind=find(thetacycle2(:,13)==lfp(i,1));
        thetacycle_lfp(ind,lfp(i,2),:)=squeeze(thetacycle_lfp(ind,lfp(i,2),:))'+lfp(i,3:9);
        nn(ind,lfp(i,2))=nn(ind,lfp(i,2))+1;
    end
    thetacycle_lfp=thetacycle_lfp./nn;
    thetacyclelfpall=cat(1,thetacyclelfpall,thetacycle_lfp);
    in
end
thetaseq3all=nan*ones(size(thetaseqall,1),101,36);
for i=1:50
    ind=ceil(thetacycle2all(:,27)/2)==i;
    startind=52-i;
    endind=101-i;
    thetaseq3all(ind,startind:endind,:)=thetaseqall(ind,:,:,1);
end

ifconnect=thetacycle2all(2:end,1)-thetacycle2all(1:end-1,2);
indconnect=find(ifconnect==1);
thetacycleconnect=thetacycle2all(indconnect,:);
thetacycleconnect2=thetacycle2all(indconnect+1,:);
thetaseqconnect=cat(3,thetaseqall(indconnect,:,:,1),thetaseqall(indconnect+1,:,:,1));
thetaseq3connect=cat(3,thetaseq3all(indconnect,:,:),thetaseq3all(indconnect+1,:,:));
thetafireconnect=cat(3,thetafireall(indconnect,:,:),thetafireall(indconnect+1,:,:));
thetalfpconnect=cat(2,thetacyclelfpall(indconnect,:,:),thetacyclelfpall(indconnect+1,:,:));
thetacycleconnectv=(thetacycleconnect(:,10)+thetacycleconnect2(:,10))/2;
indsame=thetacycleconnect(:,61)==thetacycleconnect2(:,61) & min(thetacycleconnect(:,64),thetacycleconnect2(:,64))>0.8 & thetacycleconnect(:,3)<0.15;
figure;depath=50;h=1;nn=[];seqmean=[];firemean=[];
for b=1:4
    for m=1:2
        as=[];bs=[];lfp=[];
        for v=1:4
           if v==1
               indv=thetacycleconnectv>=0 & thetacycleconnectv<1;
           elseif v==2
               indv=thetacycleconnectv>=1 & thetacycleconnectv<5;
           elseif v==3
               indv=thetacycleconnectv>=5 & thetacycleconnectv<10;
           else
               indv=thetacycleconnectv>=10 ;
           end
           ind00= indv & indsame & thetacycleconnect(:,61)==2-m & thetacycleconnect(:,62)==b ;
           nn(b,v,m)=sum(ind00);
           if b>2
               a=squeeze(nanmean(thetaseq3connect(ind00,:,:),1));
           else
               a=squeeze(nanmean(thetaseqconnect(ind00,:,:),1));
           end
           as=[as,[a]];
           b1=squeeze(nanmean(thetafireconnect(ind00,:,:),1));
           bs=[bs,[b1]];
           f1=squeeze(nanmean(thetafireconnect(ind00,:,:),1));
           bs=[bs,[b1]];
           %seqmean(b,m,v,:,:)=a;
           firemean(b,m,v,:,:)=b1;
           lfp1=squeeze(nanmean(thetalfpconnect(ind00,:,:),1))';
           lfp=[lfp,lfp1];
        end
        subplot(4,2,b+2*(b>2)+(m-1)*2);
        imagesc([as]);set(gca,'YDir','normal');xline([1:4]*72,'k');
        for f=3:4
            hold on;plot(gaussian_smooth(bs(f,:))*40+10*f-15,'k');
        end        
        for f=1:3:6
            hold on;plot(lfp(f,:)*10+10*f,'r');
        end 
        title([bname{b},mname{m},num2str([sum(ind00)])]);
    end
end

figure;
for b=1:4
    for m=1:2
        subplot(2,2,b);%subplot(1,2,m);
        ind=thetacycle2all(:,62)==b & thetacycle2all(:,61)==2-m;
        hold on;histogram(thetacycle2all(ind,3));
    end
end

cyclelength=[];mu=[];err=[];n=[];
for b=1:4
    for m=1:2
        for v=1:5
            if v==1
                indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
            elseif v==2
                indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
            elseif v==3
                indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
            elseif v==4
                indv=thetacycle2all(:,10)>=10 ;
            else
                indv=1;
            end
            ind=thetacycle2all(:,62)==b & thetacycle2all(:,61)==2-m & indv;
            mu(b,v,m)=nanmean(thetacycle2all(ind,3));
            n(b,v,m)=sum(ind);
            err(b,v,m)=nanstd(thetacycle2all(ind,3))*(n(b,v,m)^(-0.5));
        end
    end
end
plotcolors5=slanCM('bone',6);plotcolors5=plotcolors5(6:-1:1,:);bar_pos=0.14;errbar_width=1;
figure;
for v=1:4
    subplot(2,2,v);
    bar([1:4],squeeze(mu(:,v,:)),'FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:4]-bar_pos,squeeze(mu(:,v,1)),squeeze(err(:,v,1)),'k.','LineWidth',errbar_width);
    hold on;errorbar([1:4]+bar_pos,squeeze(mu(:,v,2)),squeeze(err(:,v,2)),'k.','LineWidth',errbar_width);
    xticks([1:4]);xticklabels(bname);legend('Local','Remote');title(vgroup{v});
end

x=thetacycle2all(:,62)*2-thetacycle2all(:,61);
y=thetacycle2all(:,3);
figure;
for v=1:4
    if v==1
        indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
    elseif v==2
        indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
    elseif v==3
        indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
    else
        indv=thetacycle2all(:,10)>=10 ;
    end
    subplot(2,2,v);swarmchart(x(indv),y(indv),[],[1 1 1]*0.5,'.');
end

fname={'ca1 p','ca3 p','ca1 i','ca3 i'};
figure;
for b=1:4
    for m=1:2
        for f=3:4
            subplot(4,2,b+(f-3)*4);
            bs=squeeze(firemean(b,m,:,f,:))';
            hold on;plot(gaussian_smooth(bs(:)),'-');xline([1:4]*72,'k');
            title([bname{b},mname{m},fname{f}]);
        end
    end
end

bname={'reward ','center ','in ','out '};mname={'local ','remote '};
figure;depath=50;h=1;nn=[];
for b=1:4
    for m=1:2
        as=[];bs=[];
        for v=1:4
           if v==1
               indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
           elseif v==2
               indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
           elseif v==3
               indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
           else
               indv=thetacycle2all(:,10)>=10 ;
           end
           ind= indv & thetacycle2all(:,64)>0.9 & thetacycle2all(:,61)==2-m & thetacycle2all(:,62)==b ;%& thetaseqslopeall(:,4,h)<0.05;
           ind00=ind & abs(thetaseqslopeall(:,7,h))<=depath ;%& thetacycle2all(:,65)==in;
           nn(b,v,m)=sum(ind00);
           if b>2
               a=squeeze(nanmean(thetaseq3all(ind00,:,:),1));
           else
               a=squeeze(nanmean(thetaseqall(ind00,:,:,h),1));
           end
           as=[as,imgaussfilt([a,a])];
           b1=squeeze(nanmean(thetafireall(ind00,:,:),1));
           bs=[bs,[b1,b1]];
        end
        subplot(2,4,b+(m-1)*4);
        imagesc([as]);set(gca,'YDir','normal');xline([1:4]*72,'k');
        for f=3:4
            hold on;plot(gaussian_smooth(bs(f,:))*40+10*f-20,'k');
        end
        title([bname{b},mname{m},num2str([sum(ind),in])]);
    end
end





bname={'reward ','center ','in ','out '};mname={'local ','remote '};
figure;depath=50;h=1;
for b=1:4
    for m=1:2
        as=[];bs=[];
        for v=1:4
           if v==1
               indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
           elseif v==2
               indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
           elseif v==3
               indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
           else
               indv=thetacycle2all(:,10)>=10 ;
           end
           ind= indv & thetacycle2all(:,64)>0.5 & thetacycle2all(:,61)==2-m & thetacycle2all(:,62)==b & thetaseqslopeall(:,4,h)<0.05;
           ind00=ind & abs(thetaseqslopeall(:,7,h))<=depath ;
           ind1=ind & thetaseqslopeall(:,7,h)>0 & thetaseqslopeall(:,7,h)<depath;
           ind2=ind & thetaseqslopeall(:,7,h)>-depath & thetaseqslopeall(:,7,h)<0;
           if b>2
               a=[squeeze(nanmean(thetaseq3all(ind00,:,:),1))];%squeeze(nanmean(thetaseq2all(ind1,:,:),1));squeeze(nanmean(thetaseq2all(ind2,:,:),1))];
           else
               a=[squeeze(nanmean(thetaseqall(ind00,:,:,h),1));squeeze(nanmean(thetaseqall(ind1,:,:,h),1));squeeze(nanmean(thetaseqall(ind2,:,:,h),1))];
           end
           as=[as,[a,a]];
           b1=[squeeze(nanmean(thetafireall(ind00,:,:),1));squeeze(nanmean(thetafireall(ind1,:,:),1));squeeze(nanmean(thetafireall(ind2,:,:),1))];
           bs=[bs,[b1,b1]];
        end
        subplot(2,4,b+(m-1)*4);
        imagesc([as]);set(gca,'YDir','normal');xline([1:4]*72,'k');
        for f=3:4
            hold on;plot(bs(f,:)*40+10*f-20,'k');
            hold on;plot(bs(f+6,:)*40+50+10*f-20,'r');
            hold on;plot(bs(f+12,:)*40+100+10*f-20,'k');
        end
        title([bname{b},mname{m},num2str([sum(ind),in])]);
    end
end


figure;depath=50;h=1;
for b=1:4
    for m=1:2
        for v=1:4
           if v==1
               indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
           elseif v==2
               indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
           elseif v==3
               indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
           else
               indv=thetacycle2all(:,10)>=10 ;
           end
           ind= indv & thetacycle2all(:,64)>0.5 & thetacycle2all(:,61)==2-m & thetacycle2all(:,62)==b & thetaseqslopeall(:,4,h)<0.05;
           as=[];
           for f=3:6
               a=squeeze(nanmean(thetafireall(ind,f,:),1))';
               as=[as,[a,a]];
           end
           subplot(2,4,b+(m-1)*4);
           hold on;plot(as);xline([1:4]*72,'k');
        end
        ylabel('zscore fire rate');
        title([bname{b},mname{m},num2str([sum(ind),in])]);
    end
end

figure;depath=50;h=1;
for b=1:4
    for m=1:2
        for f=3:4
           as=[];
           for v=1:4
           if v==1
               indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
           elseif v==2
               indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
           elseif v==3
               indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
           else
               indv=thetacycle2all(:,10)>=10 ;
           end
           ind= indv & thetacycle2all(:,64)>0.5 & thetacycle2all(:,61)==2-m & thetacycle2all(:,62)==b & thetaseqslopeall(:,4,h)<0.05;
           a=squeeze(nanmean(thetafireall(ind,f,:),1))';
           as=[as,[a,a]];
           end
           subplot(2,4,b+(m-1)*4);
           hold on;plot(as);xline([1:4]*72,'k');
        end
        ylabel('zscore fire rate');
        title([bname{b},mname{m},num2str([sum(ind),in])]);
    end
end

vgroup={'0-1','1-5','5-10','> 10'};
fname={'ca1 p','ca3 p','ca1 i','ca3 i'};
figure;depath=50;h=1;
for v=1:4
    if v==1
        indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
    elseif v==2
        indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
    elseif v==3
        indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
    else
        indv=thetacycle2all(:,10)>=10 ;
    end   
    for m=1:2
        for b=1:4
            as=[];
            for f=3:4
                ind= indv & thetacycle2all(:,64)>0.5 & thetacycle2all(:,61)==2-m & thetacycle2all(:,62)==b & thetaseqslopeall(:,4,h)<0.05;
                a=gaussian_smooth(squeeze(nanmean(thetafireall(ind,f,:),1)))';
                as=[as,[a,a]];
            end
            subplot(2,4,v+(m-1)*4);
            hold on;plot(as);
        end
        ylabel('zscore fire rate');
        title([vgroup{v},mname{m}]);xline([1:4]*72,'k');
        xticks([1:4]*72-36);xticklabels(fname)
    end
end
legend(bname);

figure;
for b=1:2
    for m=1:2
        as=[];
        for v=1:4
           if v==1
               indv=thetacycle2all(:,10)>=0 & thetacycle2all(:,10)<1;
           elseif v==2
               indv=thetacycle2all(:,10)>=1 & thetacycle2all(:,10)<5;
           elseif v==3
               indv=thetacycle2all(:,10)>=5 & thetacycle2all(:,10)<10;
           else
               indv=thetacycle2all(:,10)>=10 ;
           end
           ind=ind0 & indv & thetacycle2all(:,60)>0.9 & thetacycle2all(:,61)==2-m & thetacycle2all(:,62)==b & thetaseqslopeall(:,4,h)<0.05;
           ind00=ind & abs(thetaseqslopeall(:,7,1))<=10 ;
           ind1=ind & thetaseqslopeall(:,7,1)>5 ;
           ind2=ind & thetaseqslopeall(:,7,1)<-5 ;
           subplot(2,2,b+(m-1)*2+(h-1)*4);
           hold on;histogram(thetacycle2all(ind,3));
           title([bname{b},mname{m},num2str(sum(ind))]);
        end
    end
end
      
cd(homedir);save('thetaseq_2way',"thetaseq",'thetaseqshift')

cd(homedir);load('thetaseq_2way',"thetaseq",'thetaseqshift')
figure;b=1;indpos=36:65;
for in=1:16
    for b=1:4
        subplot(1,4,b);
        as=[];
        for m=1:3
        a=[];
        for v=1:4
            a=[a;[squeeze(thetaseqshift(indpos,:,b,v,m,in)),squeeze(thetaseqshift(indpos,:,b,v,m,in))]];
        end
        as=[as,a];
        end
        imagesc(as);set(gca,'YDir','normal');xline([1:6]*36);
    end
    pause;clf;
end


figure;b=1;indpos=36:65;
for b=1:2
    subplot(1,2,b);
    as=[];
    for m=1:3
        a=[];
        for v=1:4
            a1=squeeze(nanmean(thetaseqshift(indpos,:,b,v,m,:),6));
            a=[a;[a1,a1]];
        end
        as=[as,a];
    end
    imagesc(as);set(gca,'YDir','normal');xline([1:6]*36);
end

figure;
for b=1:2
    subplot(1,2,b);
    as=[];
    for m=1:3
        a=[];
        for v=1:4
            a1=squeeze(nanmean(thetaseq(:,:,b,v,m,:),6));
            a1=a1./(ones(50,1)*sum(a1,1));
            %a1=a1./(sum(a1,2)*ones(1,36));
            a=[a;[a1,a1]];
        end
        as=[as,a];
    end
    imagesc(as);set(gca,'YDir','normal');xline([1:6]*36);
end