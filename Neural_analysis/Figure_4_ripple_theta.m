% 1. code for theta, nontheta ripple csd, lfp, spectrogram, generate
% 'fig_ripple_theta_csd_lfp_v2' for each session
clear all;close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
cd(homedir);load('Rem_SW_sleep_analysis_2s','nonremall');
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_dHPC2','pksets','inter_run','thetacycle');
    load('CSD_LFP_Visualization.mat', 'CSDall');
    load('dentate_spike_spectrum_analysis_2.mat', 'spec_wavelet','fre');  
    load('Ripple_Events_PAPER5.mat','ripples','ripples_sleep_maze');
    load('dHPC_layer7channel3_ca1_adjusted_LFP_v2','LFP_HPC','lfptime','dHPC_layer7channel3');
    LFP_HPC_z=zscore(LFP_HPC(:,1:7),0,1);
   
    ruler=thetacycle(:,6);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ruler=thetacycle(:,22);template=ripples(:,3);[outputindex,error1]=match(template,ruler,0);
    ruler=thetacycle(:,23);template=ripples(:,3);[outputindex,error2]=match(template,ruler,0);
    ripples(:,15)=min([abs(error),abs(error1),abs(error2)],[],2);  
    ind=abs(ripples(:,11))==0.5;
    ripples(ind,12)=inter_run(ripples(ind,12),2);
    ripples(:,18)=round(pksets(ripples(:,12),7));
    ripples(:,19)=in;

    ind1=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035;
    ind2=(ripples(:,18)<=-4 | (ripples(:,18)>0 & ripples(:,18)<=8) );
    ripples=ripples(ind1 & ind2,:);
    ruler=lfptime(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,20)=outputindex;
    [~,I]=sort(ripples(:,7));ripples=ripples(I,:);
    ruler=lfptime(:,1);template=ripples(:,7);[outputindex,error]=match(template,ruler,0);
    ripples(:,21)=outputindex;
    [~,I]=sort(ripples(:,8));ripples=ripples(I,:);
    ruler=lfptime(:,1);template=ripples(:,8);[outputindex,error]=match(template,ruler,0);
    ripples(:,22)=outputindex;
    [~,I]=sort(ripples(:,3));ripples=ripples(I,:); % sort ripples according to peak time

    % below to calculate for ripples during maze 
    ripple1=ripples(ripples(:,5)==0,:);halfwinN=150;
    ripple1_csd=zeros(size(ripple1,1),size(CSDall,1),2*halfwinN+1,3);
    ripple1_lfp=zeros(size(ripple1,1),7,2*halfwinN+1,3);
    ripple1_lfp_z=zeros(size(ripple1,1),7,2*halfwinN+1,3);
    ripple1_spectrogram=zeros(size(ripple1,1),length(fre),2*halfwinN+1,7,3);
    ripple1(:,23:25)=0;
    for i=1:size(ripple1,1)
        for m=1:3
        ind= [ripple1(i,19+m)-halfwinN:ripple1(i,19+m)+halfwinN];
        if ind(1)>0 & ind(end)<=size(lfptime,1)
            ripple1(i,22+m)=max(lfptime(ind,1))-min(lfptime(ind,1));
            ripple1_csd(i,:,:,m)=CSDall(:,ind);
            ripple1_lfp(i,:,:,m)=LFP_HPC(ind,1:7)';
            ripple1_lfp_z(i,:,:,m)=LFP_HPC_z(ind,1:7)';
            for k=1:7
                spec=spec_wavelet{k};
                ripple1_spectrogram(i,:,:,k,m)=spec(:,ind);
            end
        end
        end
    end  

    % below for sws
    ripplesleep=ripples_sleep_maze(ripples_sleep_maze(:,9)<0,:);
    sws=nonremall(nonremall(:,8)==in,:);
    ripplesleep(:,10)=0;% if in sws
    for i=1:size(sws,1)
        ind= ripplesleep(:,3)>=sws(i,4) & ripplesleep(:,3)<=sws(i,5) ;
        ripplesleep(ind,10)=i;
    end
    ripplesleep(:,14)=in;
    ind1= ripplesleep(:,2)-ripplesleep(:,1)>0.035 & ripplesleep(:,10)>0;
    ripplesleep=ripplesleep(ind1,:);
    if ~isempty(ripplesleep)
    ruler=lfptime(:,1);template=ripplesleep(:,3);[outputindex,error]=match(template,ruler,0);
    ripplesleep(:,20)=outputindex;
    [~,I]=sort(ripplesleep(:,7));ripplesleep=ripplesleep(I,:);
    ruler=lfptime(:,1);template=ripplesleep(:,7);[outputindex,error]=match(template,ruler,0);
    ripplesleep(:,21)=outputindex;
    [~,I]=sort(ripplesleep(:,8));ripplesleep=ripplesleep(I,:);
    ruler=lfptime(:,1);template=ripplesleep(:,8);[outputindex,error]=match(template,ruler,0);
    ripplesleep(:,22)=outputindex;
    [~,I]=sort(ripplesleep(:,3));ripplesleep=ripplesleep(I,:);
    ripplesleep(:,23:25)=0;
    end
    ripple1sleep=ripplesleep(ripplesleep(:,5)==0,:);
    ripple1sleep_csd=nan*ones(size(ripple1sleep,1),size(CSDall,1),2*halfwinN+1,3);
    ripple1sleep_lfp=nan*ones(size(ripple1sleep,1),7,2*halfwinN+1,3);
    ripple1sleep_lfp_z=nan*ones(size(ripple1sleep,1),7,2*halfwinN+1,3);
    ripple1sleep_spectrogram=nan*ones(size(ripple1sleep,1),length(fre),2*halfwinN+1,7,3);
    
    for i=1:size(ripple1sleep,1)
        for m=1:3
        ind= [ripple1sleep(i,19+m)-halfwinN:ripple1sleep(i,19+m)+halfwinN];
        if ind(1)>0 & ind(end)<=size(lfptime,1)
            ripple1sleep(i,22+m)=max(lfptime(ind,1))-min(lfptime(ind,1));
            ripple1sleep_csd(i,:,:,m)=CSDall(:,ind);
            ripple1sleep_lfp(i,:,:,m)=LFP_HPC(ind,1:7)';
            ripple1sleep_lfp_z(i,:,:,m)=LFP_HPC_z(ind,1:7)';
            for k=1:7
                spec=spec_wavelet{k};
                ripple1sleep_spectrogram(i,:,:,k,m)=spec(:,ind);
            end
        end
        end
    end  
    cd(savedir);save('fig_ripple_theta_csd_lfp_v2','ripples','ripple1','ripple1_spectrogram',"ripple1_lfp","ripple1_lfp_z","ripple1_csd",'halfwinN',...
        'ripplesleep','ripple1sleep','ripple1sleep_spectrogram',"ripple1sleep_lfp","ripple1sleep_lfp_z","ripple1sleep_csd",'fre','-v7.3');
    in
end
 
% 2. code for pk_behaveall, quantify ripple rate during theta vs non-theta at reward pause
clear all;close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
posdt=0.0335;rippleall=[];pk_behaveall=[];pksetsall=[];ca1ca3_overlap=[];reward_dis_cut=70;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('fig_ripple_theta_csd_lfp','ripples');
    load('ThetaCycle_Decode_Info_dHPC2','pksets','inter_run','thetacycle');
    load('Position_Data_Maze.mat','Position_Data');
    load('new_reward_center_in_out','pos_behave_marker');
    % below use new version of behavior epochs
    ruler=Position_Data(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,11:12)=pos_behave_marker(outputindex,3:4);
    ruler=Position_Data(:,1);template=thetacycle(:,6);[outputindex,error]=match(template,ruler,0);
    thetacycle(:,11:12)=pos_behave_marker(outputindex,3:4);

    ind1=ripples(:,6)==1;ind2=ripples(:,6)==2;ripples(:,21)=0;
    ruler=ripples(ind2,3);template=ripples(ind1,3);[outputindex,error]=match(template,ruler,0);
    ripples(ind1,21)=abs(error);
    ruler=ripples(ind1,3);template=ripples(ind2,3);[outputindex,error]=match(template,ruler,0);
    ripples(ind2,21)=abs(error);
    ca1ca3_overlap(in,:)=[sum(ind1),sum(ind2),nanmean(abs(ripples(ind1,21))<0.01),nanmean(abs(ripples(ind2,21))<0.01)];
    ripples(:,25)=ripples(:,2)-ripples(:,1);

    thetacycle(:,82)=0; % associated ripple id
    for i=1:size(ripples,1)
        if ripples(i,6)==1
        ind= thetacycle(:,6)>=ripples(i,1)-0.05 & thetacycle(:,6)<=ripples(i,2)+0.05 ;
        thetacycle(ind,82)=i;
        end
    end
    indrt=ripples(:,16)>0 & abs(ripples(:,15))<=0.05;
    pk_behave=[];
    for i=1:size(pksets,1)
        for b=1:4
            indb= pos_behave_marker(:,3)==b & pos_behave_marker(:,4)==i ;
            indt= thetacycle(:,11)==b & thetacycle(:,12)==i & thetacycle(:,82)==0;
            indr= ripples(:,11)==b & ripples(:,12)==i ;
            for v=1:4
                if v==1
                    indv= Position_Data(:,5)>=0 & Position_Data(:,5)<1;
                    indtv= thetacycle(:,10)>=0 & thetacycle(:,10)<1;
                    indrv= ripples(:,10)>=0 & ripples(:,10)<1;
                elseif v==2
                    indv= Position_Data(:,5)>=1 & Position_Data(:,5)<5;
                    indtv= thetacycle(:,10)>=1 & thetacycle(:,10)<5;
                    indrv= ripples(:,10)>=1 & ripples(:,10)<5;
                elseif v==3
                    indv= Position_Data(:,5)>=5 & Position_Data(:,5)<10;
                    indtv= thetacycle(:,10)>=5 & thetacycle(:,10)<10;
                    indrv= ripples(:,10)>=5 & ripples(:,10)<10;
                elseif v==4
                    indv= Position_Data(:,5)>=10;
                    indtv= thetacycle(:,10)>=10;
                    indrv= ripples(:,10)>=10;
                end
                pk_behave(i,b,v,1)=sum(indv & indb)*posdt; % total pos time
                pk_behave(i,b,v,2)=sum(thetacycle( indtv & indt, 3)); % total thetacycle time
                pk_behave(i,b,v,3)=pk_behave(i,b,v,1)-pk_behave(i,b,v,2);
                pk_behave(i,b,v,4)=sum(indrv & indr); % total ripple N
                pk_behave(i,b,v,5)=sum(indrv & indr & indrt); % ripple N associated with thetacycles
                pk_behave(i,b,v,6)=pk_behave(i,b,v,4)-pk_behave(i,b,v,5);
                pk_behave(i,b,v,7)=sum(indrv & indr & ripples(:,6)==1); % total ca1 ripple N
                pk_behave(i,b,v,8)=sum(indrv & indr & ripples(:,6)==1 & indrt); % ca1 ripple N associated with thetacycles
                pk_behave(i,b,v,9)=pk_behave(i,b,v,7)-pk_behave(i,b,v,8);
                pk_behave(i,b,v,10)=sum(indrv & indr & ripples(:,6)==2); % total ca3 ripple N
                pk_behave(i,b,v,11)=sum(indrv & indr & ripples(:,6)==2 & indrt); % ca3 ripple N associated with thetacycles
                pk_behave(i,b,v,12)=pk_behave(i,b,v,10)-pk_behave(i,b,v,11);
                pk_behave(i,b,v,18)=sum(ripples(indrv & indr & ripples(:,6)==1,25)); % total ca1 ripple duration
                pk_behave(i,b,v,19)=sum(ripples(indrv & indr & ripples(:,6)==2,25)); % ca3 ripple duration 
                pk_behave(i,b,v,20)=sum(ripples(indrv & indr & ripples(:,6)==2 & (indrt),25)); % ca3 ripple duration with theta
                pk_behave(i,b,v,21)=sum(ripples(indrv & indr & ripples(:,6)==2 & (~indrt),25)); % ca3 ripple duration without theta
            end
        end
    end   
    pk_behaveall=cat(1,pk_behaveall,pk_behave);
    rippleall=[rippleall;ripples];
    pksets(:,42)=in;
    pksetsall=[pksetsall;pksets];
    in
end
pk_behaveall(:,:,:,13)=pk_behaveall(:,:,:,2)./pk_behaveall(:,:,:,1); % thetacycle time percent
pk_behaveall(:,:,:,14)=pk_behaveall(:,:,:,7)./pk_behaveall(:,:,:,1); % ca1 ripple rate
pk_behaveall(:,:,:,15)=pk_behaveall(:,:,:,10)./pk_behaveall(:,:,:,1); % ca3 ripple rate
pk_behaveall(:,:,:,16)=pk_behaveall(:,:,:,11)./pk_behaveall(:,:,:,1); % theta ca3 ripple rate
pk_behaveall(:,:,:,17)=pk_behaveall(:,:,:,12)./pk_behaveall(:,:,:,1); % not theta ca3 ripple rate
indexcludepk= ( round(pksetsall(:,7))>-4 & round(pksetsall(:,7))<=0 ) | ( round(pksetsall(:,7))~=pksetsall(:,7) );
indpk=~indexcludepk;
cd(homedir);save('pk_behaveall','pk_behaveall','indpk','rippleall','pksetsall','ca1ca3_overlap');


if 0
puor=slanCM('PuOr');puor=puor(256:-1:1,:); colormap(puor);
theta_time_cutoff=0.05;nn=[];mean_lfp=[];mean_csd=[];mean_spectrogram=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('fig4_ripple_4');
    ripple1(:,22)=0;
    ind=ripple1(:,15)<=theta_time_cutoff & ripple1(:,16)>0;
    ripple1(ind,22)=1;
    a1=[];a2=[];a3=[];
    for r=1:2
        for t=1:3
            if t==1
                ind=ripple1(:,22)==1 & ripple1(:,6)==r;
            elseif t==2
                ind=ripple1(:,22)==0 & ripple1(:,6)==r;
            else
                ind=ripple1sleep(:,6)==r;
            end
            nn(in,r,t)=sum(ind);
            
            if t<3
                a=squeeze(nanmean(ripple1_csd(ind,:,:),1));
                a1=[a1,a];mean_csd{in,r,t}=a;
                a=squeeze(nanmean(ripple1_lfp(ind,:,:),1));
                a2=[a2,a];mean_lfp(in,r,t,:,:)=a;
                a=squeeze(nanmean(ripple1_spectrogram(ind,:,:,:),1));
                a3=[a3,a];mean_spectrogram(in,r,t,:,:,:)=a;
            else
                a=squeeze(nanmean(ripple1sleep_csd(ind,:,:),1));
                a1=[a1,a];mean_csd{in,r,t}=a;
                a=squeeze(nanmean(ripple1sleep_lfp(ind,:,:),1));
                a2=[a2,a];mean_lfp(in,r,t,:,:)=a;
                a=squeeze(nanmean(ripple1sleep_spectrogram(ind,:,:,:),1));
                a3=[a3,a];mean_spectrogram(in,r,t,:,:,:)=a;
            end
        end
    end
    a=squeeze(nn(in,:,:));
    figure(1);subplot(4,4,in);imagesc(a1);set(gca,'YDir','normal');title(num2str(a(:)'));
    xline([1:6]*(2*halfwinN+1));colormap(puor);caxis(max(abs(a1(:)))*[-1 1]);
    figure(2);subplot(4,4,in);title(num2str(a(:)'));xline([1:6]*(2*halfwinN+1));
    for k=1:7
        hold on;plot(a2(k,:)+8-k,'k');
    end
    %figure(3);subplot(4,4,in);imagesc(a3);set(gca,'YDir','normal');title(num2str(a(:)'));
end
figure;
a1=[];
for r=1:2
    for t=1:3
        a=squeeze(nanmean(mean_lfp(:,r,t,:,:),1));
        a1=[a1,a];
    end
end
subplot(4,2,1);
xline([1:6]*(2*halfwinN+1));xlim([0 6*(2*halfwinN+1)]);title('CA1: theta | non-theta | sws | CA3: theta | non-theta | sws ');
for k=1:7
    hold on;plot(a1(k,:)+8-k,'k');
end
for k=1:7
    a1=[];
    for r=1:2
        for t=1:3
            a=squeeze(nanmean(mean_spectrogram(:,r,t,:,:,k),1));
            a1=[a1,a];
        end
    end
    subplot(4,2,1+k);imagesc([],fre,a1);set(gca,'YDir','normal');xline([1:6]*(2*halfwinN+1),'k');
end

cd(homedir);load('next_arm_deseq_theta_ripple_4','pksetsall','pksetsall3','rippleall');  
ind=~( rippleall(:,18)>-4 & rippleall(:,18)<=0 );
rippleall=rippleall(ind,:);
ind=abs(rippleall(:,11))~=1;
rippleall=rippleall(ind,:);
pksetsall(:,7)=round(pksetsall(:,7));
runduration_cutoff=10;
indb1= pksetsall(:,17)<30 & pksetsall(:,19)<30 & pksetsall(:,11)>=60 & pksetsall(:,13)>=60 & max(pksetsall(:,14:15),[],2)<=runduration_cutoff & pksetsall(:,16)<=60 & pksetsall(:,43)==1 ;
indb2= pksetsall(:,17)<30 & pksetsall(:,19)<30 & pksetsall(:,11)>=60 & pksetsall(:,13)>=60 & max(pksetsall(:,14:15),[],2)<=runduration_cutoff & pksetsall(:,44)<=15 & pksetsall(:,43)==1 ;
%indb1= pksetsall(:,16)<=60*2 & pksetsall(:,43)==1 ;
%indb2= pksetsall(:,44)<=15*2 & pksetsall(:,43)==1 ;
%indb1= pksetsall(:,43)==1 ;
%indb2= pksetsall(:,43)==1 ;
%indb1=1;
%indb2=1;
rippleall(:,20:21)=0;
for in=1:16
    indin=pksetsall3(:,45)==in;
    indin2=rippleall(:,19)==in;
    for b=1:2
        if b==1
            indb=indb1;indbr=rippleall(:,11)==0;
        else
            indb=indb2;indbr=abs(rippleall(:,11))==0.5;
        end
        rippleall(indin2 & indbr,21)=b;
        lia=ismember(rippleall(:,12),pksetsall3(indin & indb,31));
        rippleall(lia & indin2 & indbr,20)=b;
    end
end
rname={'ca1 ripple','ca3 ripple'};bname={'reward','center'};
figure;
for r=1:2
    for b=1:2
        ind0=rippleall(:,6)==r & rippleall(:,21)==b;
        ind=rippleall(:,6)==r & rippleall(:,20)==b;
        ind3=rippleall(:,6)==r & rippleall(:,20)==b & rippleall(:,10)<=5;
        ind4=rippleall(:,6)==r & rippleall(:,20)==b & rippleall(:,4)>4;
        subplot(2,4,b+(r-1)*2);
        histogram(rippleall(ind0,22));hold on;histogram(rippleall(ind,22));
        hold on;histogram(rippleall(ind3,22));hold on;histogram(rippleall(ind4,22));
        title([rname{r},' + ',bname{b}]);
        subplot(2,4,b+(r-1)*2+4);
        histogram(rippleall(ind0,23));hold on;histogram(rippleall(ind,23));
        hold on;histogram(rippleall(ind3,23));hold on;histogram(rippleall(ind4,23));
        title(sum(ind)/sum(ind0));
    end
end



cd(homedir);load('next_arm_deseq_theta_ripple_3');
% below shuffle next arm id
ripple_shuffle=[];
for h=1:2
    for i=1:size(rippleall,1)
        sh=randi(6,500,1);
        a=squeeze(ripple_decode6all(i,:,h));
        ripple_shuffle(i,:,h)=a(sh);
    end
end
cd(figdir);save('ripple_theta','ripple_decodeall','rippleall','timewin','ripple_decode6all','ripple_shuffle');

indrt=abs(rippleall(:,15))<=0.05;
indreward=rippleall(:,11)==0 ;
indforce= rippleall(:,24)<=3 & rippleall(:,24)>0 ;
indfree= abs(rippleall(:,24))>=4 & abs(rippleall(:,24))<=7;

figure;histogram(rippleall(indrt & indreward,6));
hold on;histogram(rippleall((~indrt) & indreward,6));

ripple_comp=[];num=[];dif=rippleall(:,22)-ripple_shuffle(:,:,1);
for r=1:2
    indr=rippleall(:,6)==r;
for t=1:3
    if t==1
        indt=indrt;
    elseif t==2
        indt=~indrt;
    else
        indt=1;
    end
    for h=1:2
        if h==1
            indh=indforce;
        else
            indh=indfree;
        end
        for v=1:3
            if v==1
                indv= rippleall(:,10)>=0 & rippleall(:,10)<1;
            elseif v==2
                indv= rippleall(:,10)>=1 & rippleall(:,10)<5;
            elseif v==3
                indv= rippleall(:,10)>=5 & rippleall(:,10)<10;
            end
            ind=indr & indt & indh & indv & indreward;
            a=nanmean(dif(ind,:),1);
            num(r,t,h,v)=sum(~isnan(a));
            ripple_comp(r,t,h,v,1)=nanmean(a);
            ripple_comp(r,t,h,v,2)=nanstd(a)*(num(r,t,h,v)^(-0.5));
            %p=signrank(dif(ind));%
            [~,p]=ttest(a);
            ripple_comp(r,t,h,v,3)=p;
        end
    end
end
end
%cd(figdir);save('ripple_theta','ripple_decodeall','rippleall','timewin','ripple_comp','num');
plotcolors5=slanCM('bone',6);plotcolors5=plotcolors5(6:-1:1,:);
tname={'Theta state ','Non-theta ','All '};rname={'CA1 ripple','CA3 ripple'};
figure;
for r=1:2
for t=1:3
    subplot(2,3,t+(r-1)*3);
    hold on;ba=bar([1:3],squeeze(ripple_comp(r,t,:,:,1))','FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(3,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:3]-0.15,squeeze(ripple_comp(r,t,1,:,1)),squeeze(ripple_comp(r,t,1,:,2)),'k.','LineWidth',2);
    hold on;errorbar([1:3]+0.15,squeeze(ripple_comp(r,t,2,:,1)),squeeze(ripple_comp(r,t,2,:,2)),'k.','LineWidth',2);
    a=squeeze(ripple_comp(r,t,:,:,3));n=squeeze(num(r,t,:,:));
    title({[tname{t},rname{r}];[': p = ',num2str(a(:)',2)];['n = ',num2str(n(:)')];});
    xticks([1:4]);xticklabels({'0-1','1-5','5-10'});legend('Forced','Free');
    ylabel('Next - Other Unvisited Mean Prob');
end
end

% 5. below quantify ripple rate during theta vs non-theta at reward pause
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
figdir='X:\Mengni\Data_Analysis\Paper_Figures';
posdt=0.0335;rippleall=[];pk_rewardall=[];pksetsall=[];trialsetall=[];trial_rewardall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_2','ripple1_trial','ripple2');
    load('InOutBound_Behavior_Analysis.mat','pos_behave_marker');
    load('Position_Data_Maze.mat');
    ripples=[ripple1_trial(:,1:6);ripple2(:,1:6)];
    ripples=ripples(ripples(:,6)==1,:);
    [~,I]=sort(ripples(:,3));
    ripples=ripples(I,:);
    ruler=Position_Data(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,7:12)=[Position_Data(outputindex,2:5),pos_behave_marker(outputindex,:)];
    ripples(:,13)=vecnorm(ripples(:,7:8),2,2);
    ripples(:,14)=in;

    load('ThetaCycle_Decode_Info_2','pksets','trialset','thetacycle');
    thetacycle(:,82)=0; % associated ripple id
    for i=1:size(ripples,1)
        ind= thetacycle(:,6)>=ripples(i,1)-0.05 & thetacycle(:,6)<=ripples(i,2)+0.05 ;
        thetacycle(ind,82)=i;
    end
    ruler=thetacycle(:,6);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    ripples(:,15)=error;
    indrt=abs(ripples(:,15))<=0.05;
    pk_reward=[];
    for i=1:size(pksets,1)
        indp= pos_behave_marker(:,1)==0 & pos_behave_marker(:,2)==i ;
        indt= thetacycle(:,11)==0 & thetacycle(:,12)==i ;
        indr= ripples(:,11)==0 & ripples(:,12)==i ;
        for v=1:3
            if v==1
                indv= Position_Data(:,5)>=0 & Position_Data(:,5)<1;
                indtv= thetacycle(:,10)>=0 & thetacycle(:,10)<1;
                indrv= ripples(:,10)>=0 & ripples(:,10)<1;
            elseif v==2
                indv= Position_Data(:,5)>=1 & Position_Data(:,5)<5;
                indtv= thetacycle(:,10)>=1 & thetacycle(:,10)<5;
                indrv= ripples(:,10)>=1 & ripples(:,10)<5;
            elseif v==3
                indv= Position_Data(:,5)>=5 & Position_Data(:,5)<10;
                indtv= thetacycle(:,10)>=5 & thetacycle(:,10)<10;
                indrv= ripples(:,10)>=5 & ripples(:,10)<10;
            end
            pk_reward(i,v,1)=sum(indv & indp)*posdt; % total pos time
            pk_reward(i,v,2)=sum(thetacycle( indtv & indt, 3)); % total thetacycle time
            pk_reward(i,v,3)=pk_reward(i,v,1)-pk_reward(i,v,2);
            pk_reward(i,v,4)=sum(indrv & indr); % total ripple N
            pk_reward(i,v,5)=sum(indrv & indr & indrt); % ripple N associated with thetacycles
            pk_reward(i,v,6)=pk_reward(i,v,4)-pk_reward(i,v,5);
        end
    end
    pk_rewardall=cat(1,pk_rewardall,pk_reward);
    rippleall=[rippleall;ripples];
    pksets(:,42)=in;
    pksets(:,7)=round(pksets(:,7));
    pksetsall=[pksetsall;pksets];
    trial_reward=[];
    for i=1:size(trialset,1)
        ind=pksets(:,1)==trialset(i,1);
        trial_reward(i,:,:)=nansum(pk_reward(ind,:,:),1);
    end
    trialset(:,7)=in;
    trialsetall=[trialsetall;trialset];
    trial_rewardall=cat(1,trial_rewardall,trial_reward);
end
trial_rewardall(:,:,7)=trial_rewardall(:,:,2)./trial_rewardall(:,:,1); % thetacycle time percent
trial_rewardall(:,:,8)=trial_rewardall(:,:,5)./trial_rewardall(:,:,2); % thetacycle associated ripple rate
trial_rewardall(:,:,9)=trial_rewardall(:,:,6)./trial_rewardall(:,:,3); % non thetacycle associated ripple rate
session_reward=[];session_reward2=[];
for in=1:16
    ind=trialsetall(:,7)==in;
    session_reward(in,:,:)=nanmean(trial_rewardall(ind,:,7:9),1);
    
    ind=pksetsall(:,42)==in;
    pksets=pksetsall(ind,:);
    pksets(:,7)=round(pksets(:,7));
    pk_reward=pk_rewardall(ind,:,:);
    ind1= pksets(:,7)>0 & pksets(:,7)<=4;
    ind2= pksets(:,7)>4 & pksets(:,7)<=8;
    ind3= pksets(:,7)<=-4 ;
    session_reward2(in,:,1)=squeeze(nanmean(nansum(pk_reward(ind1,:,4:6),2),1));
    session_reward2(in,:,2)=squeeze(nanmean(nansum(pk_reward(ind2,:,4:6),2),1));
    session_reward2(in,:,3)=squeeze(nanmean(nansum(pk_reward(ind3,:,4:6),2),1));
end
session_stas=[];
for i=1:3
    for j=1:3
        session_stas(i,j,1)=nanmean(session_reward2(:,i,j),1);
        session_stas(i,j,2)=nanstd(session_reward2(:,i,j),0,1)/4;
    end
end
cd(homedir);save('fig4_thetaseq_ripple',"session_reward",'rippleall','trial_rewardall','pksetsall','pk_rewardall','trialsetall');

if 0
figure;iname={'Forced','Correct Free','Error Free'};
for i=1:3
    subplot(1,3,i);
    ba=bar([1:3],squeeze(session_stas(:,i,1)),'FaceAlpha',0.6,'EdgeColor',[0 0 0],'LineWidth',2);
    hold on;errorbar([1:3],squeeze(session_stas(:,i,1)),squeeze(session_stas(:,i,2)),'k.','LineWidth',2);
    xticks([1:3]);xticklabels({'All ','Theta state ','Non-theta '});title(iname{i});ylabel('Ripple N per reward');
end
end
cd(homedir);
load('fig4_thetaseq_ripple',"session_reward",'rippleall','trial_rewardall','pksetsall','pk_rewardall','trialsetall');
figure;
subplot(2,3,1);
ind0=rippleall(:,11)==0;
ind1=abs(rippleall(:,11))==0.5;
ind2=rippleall(:,11)==-1;
ind3=rippleall(:,11)==1;
forpie=[sum(ind0),sum(ind1),sum(ind2),sum(ind3)];
pie(forpie);legend('Reward','Inter','In','Out');title(['Ripple: N = ',num2str(sum(forpie))]);
subplot(2,3,2);
y=squeeze(session_reward(:,:,1));
shaded_errbar([1:3],y');
xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel('Velocity ranges (cm/s)');
ylabel('Theta State Time Percent');xlim([0.5 3.5]);title('Reward Pause Time Percent with Theta, N = 16 sessions');
subplot(2,3,3);
ps=[];mu=[];
for v=1:3
    a=session_reward(:,v,2)-session_reward(:,v,3);
    ps(v)=signrank(a);
    mu(v)=nanmean(a);
end
x = [1 2 3];vals = squeeze(nanmean(session_reward(:,:,2:3),1))';
bar(x,vals,'FaceAlpha',0.7);
for v=1:3
    for in=1:16
        hold on;plot(v+[-0.15 0.15],squeeze(session_reward(in,v,2:3)),'k');
    end
end
xticks([1:3]);xticklabels({'0-1','1-5','5-10'});xlabel('Velocity ranges (cm/s)');
ylabel('Ripple Rate(/s)');legend('Theta','Non-theta');title({['Mu = ',num2str(mu,2)];['P = ',num2str(ps,2)]});

iname={'All ','Theta state ','Non-theta '};
for i=1:3
    subplot(2,3,i+3);
    ba=bar([1:3],squeeze(session_stas(i,:,1)),'FaceAlpha',0.6,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',2);
    hold on;errorbar([1:3],squeeze(session_stas(i,:,1)),squeeze(session_stas(i,:,2)),'r.','LineWidth',2);
    for in=1:16
        %hold on;plot([1:3],squeeze(session_reward2(in,i,:)),'o','Color',[1 1 1]*0.3);
    end
    xticks([1:3]);xticklabels({'Forced','Correct Free','Error Free'});title(iname{i});ylabel('Ripple N per reward');
end
cd(figdir);figure_title='fig4sup_thetapause_ripple';save_current_figure(figure_title);


homedir='X:\Mengni\Data_Analysis\Session_combined_0324';
cd(homedir);load('SessionSet16');load('fig4_thetaseq_behave_correlate_dHPC_v','pksetsall');
figdir='X:\Mengni\Data_Analysis\Paper_Figures';cd(figdir);load('ripple_theta')
pk_rippleall=[];pksetsall(:,7)=round(pksetsall(:,7));
for in=1:16
    ind=pksetsall(:,42)==in;pksets=pksetsall(ind,:);
    ripples=rippleall(rippleall(:,14)==in,:);
    pk_ripple=[];
    for i=1:size(pksets,1)
        indr=ripples(:,11)==0 & ripples(:,12)==i ;
        for v=1:3
            if v==1
                indv= ripples(:,10)>=0 & ripples(:,10)<1;
            elseif v==2
                indv= ripples(:,10)>=1 & ripples(:,10)<5;
            elseif v==3
                indv= ripples(:,10)>=5 & ripples(:,10)<10;
            end
            indrv=indv & indr ;
            pk_ripple(i,1,v)=sum(indrv);
            for k=1:6
                indk=indrv & ripples(:,17+k)>0;
                pk_ripple(i,1+k,v)=sum(indk)/sum(indrv);
                if k==2
                    pk_ripple(i,1+k,v)=pk_ripple(i,1+k,v)/(abs(pksets(i,7))-2);
                elseif k==6
                    pk_ripple(i,1+k,v)=pk_ripple(i,1+k,v)/(8-abs(pksets(i,7)));
                    end
                end
            end
        end
    end
    pk_rippleall=cat(1,pk_rippleall,pk_ripple);

pksetsall3(:,7)=round(pksetsall3(:,7));
ind=pksetsall3(:,7)==1;pk_rippleall(ind,2:3,:,:)=nan;
ind=pksetsall3(:,7)==2;pk_rippleall(ind,2,:,:)=nan;
ind=abs(pksetsall3(:,7))==7;pk_rippleall(ind,6,:,:)=nan;
ind=abs(pksetsall3(:,7))==8;pk_rippleall(ind,5:6,:,:)=nan;
cd(homedir);
load('pksets_interrun.mat');
pksetsall2=[];interrunall(:,23)=interrunall(:,7)-interrunall(:,5);
for in=1:16
    ind=interrunall(:,21)==in;inter_run=interrunall(ind,:);
    ind=pksetsall(:,42)==in;pksets=pksetsall(ind,:);
    pksets(:,44:45)=0;
    for i=1:size(pksets,1)
        id=find(inter_run(:,2)==i);
        if length(id)==1
            pksets(i,44:45)=[interrunall(id,23),length(id)];
        else
            pksets(i,44:45)=[0,length(id)];
        end
    end
    pksetsall2=[pksetsall2;pksets];
end
pksetsall=pksetsall2;
cd(homedir);save('fig4_thetaseq_behave_correlate_dHPC_v','thetacycleall3','pksetsall3','pksetsall',"pk_rippleall",'-v7.3');


pksetsall(:,7)=round(pksetsall(:,7));
runduration_cutoff=10;mu=[];err=[];ns=[];ps=[];ps2=[];
for v=1:4
for b=1:2
    if b==1
        indb= pksetsall(:,17)<30 & pksetsall(:,19)<30 & pksetsall(:,11)>=60 & pksetsall(:,13)>=60 & max(pksetsall(:,14:15),[],2)<=runduration_cutoff & pksetsall(:,16)<=60 & pksetsall(:,43)==1 ;
    else
        indb= pksetsall(:,17)<30 & pksetsall(:,19)<30 & pksetsall(:,11)>=60 & pksetsall(:,13)>=60 & max(pksetsall(:,14:15),[],2)<=runduration_cutoff & pksetsall(:,44)<=15 & pksetsall(:,43)==1 ;
    end
    for i=1:2
        if i==1 % correct forced visit
            indi= pksetsall(:,7)>0 & pksetsall(:,7)<=3  ;
        elseif i==2 % free visit
            indi= (pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3;
        end
        ind=indi & indb;
        a=squeeze(pk_rippleall(ind,5,b,v)-pk_rippleall(ind,6,b,v));
        mu(b,i,v)=nanmean(a);
        ns(i,b,v)=sum(~isnan(a));
        err(b,i,v)=nanstd(a).*(ns(i,b,v).^(-0.5));
        [h,p]=ttest(a);
        ps(i,b,v)=p;
        %subplot(2,3,b+(i-1)*3);hold on;histogram(a);p1=signrank(a);title(num2str([p1,p],2));
    end
    ind1= pksetsall(:,7)>0 & pksetsall(:,7)<=3 & indb ;
    ind2= ((pksetsall(:,7)>=4 & pksetsall(:,7)<=7) | pksetsall(:,7)<-3) & indb ;
    a=squeeze(pk_rippleall(:,5,b)-pk_rippleall(:,6,b));
    ps2(b,v)=ranksum(a(ind1),a(ind2));
end
end
figure;bname={'Reward','Center'};
for b=1:2
    subplot(1,2,b);
    hold on;ba=bar([1:4],squeeze(mu(b,:,:)),'FaceAlpha',0.6);
    ba(1).FaceColor =[1 1 1];ba(2).FaceColor =plotcolors5(5,:);ba(2).FaceAlpha=1;
    hold on;errorbar([1:4]-0.15,squeeze(mu(b,1,:)),squeeze(err(b,1,:)),'k.');
    hold on;errorbar([1:4]+0.15,squeeze(mu(b,2,:)),squeeze(err(b,2,:)),'k.');
    a=ps(:,b,:);n=ns(:,b,:);
    title({[bname{b},': p = ',num2str(a(:)',2)];['n = ',num2str(n(:)')];['p2 = ',num2str(ps2(b,:),2)];});
    xticks([1:4]);xticklabels({'0-1','1-5','5-10','>10'});legend('Forced','Free');
    ylabel('Next - Other Unvisited Mean Prob');%FIG_INDEX=['fig4_13'];save_fig(FIG_INDEX,ifsavefig);
end
end