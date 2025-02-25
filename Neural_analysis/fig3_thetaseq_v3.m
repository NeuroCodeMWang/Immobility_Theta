function fig3_thetaseq_v3(ifsavefig)

%clear all;close all;
if nargin<1
    ifsavefig=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\figure_v3_data';
figdir='X:\Mengni\Data_Analysis\Paper_Figures\figures_v3';
timewin=0.01;reward_dis_cut=70;min_dewinN=1;ifpeak=0;phasebin=10;
if 0
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    if ifpeak
        load('ThetaCycle_Decode_Info_dHPC_peak_peak.mat','linear_maze','linear_pos','deseqall','thetacycle','thetacycle_depos_dist');
        load('Analyze_Theta_Cycle_Broad_peak_peak','reftheta','lfptime_trial','lfpz_trial');
    else
        load('ThetaCycle_Decode_Info_dHPC2.mat','linear_maze','linear_pos','deseqall','thetacycle','thetacycle_depos_dist');
        load('Analyze_Theta_Cycle_Broad','reftheta','lfptime_trial','lfpz_trial');
    end
    load('Ripple_Events_PAPER5.mat','ripples');
    load('Position_Data_Maze.mat');
    load('new_reward_center_in_out','pos_behave_marker');
   
    % below use new version of behavior epochs 
    ruler=Position_Data(:,1);template=thetacycle(:,6);[outputindex,error]=match(template,ruler,0);
    thetacycle(:,[62,68])=pos_behave_marker(outputindex,3:4);
    %ruler=Position_Data(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
    %ripples(:,11:12)=pos_behave_marker(outputindex,3:4);
    %ruler=Position_Data(:,1);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
    %deseqall(:,28:29)=pos_behave_marker(outputindex,3:4);

    ind=thetacycle(:,27)<=30;thetacycle(ind,26)=0;
    ind=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035 & ripples(:,6)==1;
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
    %deseqall(:,47)=linear_pos(:,3); % decoded arm no center
    %deseqall(:,48)=linear_pos(:,5)==linear_pos(:,6); % if local
    thetacycle_depos_dist(:,1,:)=squeeze(thetacycle_depos_dist(:,1,:))./thetacycle(:,24);
    [M,I]=max(squeeze(thetacycle_depos_dist(:,1,:)),[],2); 
    thetacycle(:,59:60)=[I-1,M];
    [M,I]=max(squeeze(thetacycle_depos_dist(:,1,2:9)),[],2); 
    thetacycle(:,63)=I;
    thetacycle(:,64)=M+squeeze(thetacycle_depos_dist(:,1,1));
    thetacycle(:,61)=thetacycle(:,26)==thetacycle(:,59) ;%| thetacycle(:,26)==thetacycle(:,63); % if local
    for i=1:size(thetacycle,1)     
        thetacycle(i,67)=squeeze(thetacycle_depos_dist(i,1,thetacycle(i,26)+1));
    end
    deseqall(:,27)=ceil(deseqall(:,27)/phasebin);
    despikez=zscore(deseqall(:,[6:7,19:20,23:24]),0,1);
    despike=deseqall(:,[6:7,19:20,23:24]);
    thetacycle2=thetacycle(thetacycle(:,79)==0 & thetacycle(:,24)>=min_dewinN & thetacycle(:,64)>=0.5 & thetacycle(:,62)>0,:);
    thetaid=thetacycle2(:,21);
    thetaseq=zeros(size(thetaid,1),50,36,3);
    thetaseq2=zeros(size(thetaid,1),101,36,3);
    thetafire=zeros(size(thetaid,1),6,36,2);
    lia=ismember(deseqall(:,40),thetaid(:,1));
    deseqall0=deseqall(lia,:);
    despikez0=despikez(lia,:);
    despike0=despike(lia,:);
    linear_maze0=linear_maze(lia,:,:);
    linear_pos0=linear_pos(lia,:);
    linear_pos0(:,2)=ceil(linear_pos0(:,2)/2);
    %ind=linear_pos0(:,2)>50;
    %linear_pos0(ind,2)=50;
    seqn=zeros(size(thetaid,1),36);
    for j=1:size(deseqall0,1)
        ratbin=linear_pos0(j,2); % shift from the rat's current position 
        if ratbin<=50
        startind=52-ratbin;
        endind=101-ratbin;  
        id=deseqall0(j,40);
        i=find(thetaid(:,1)==id);
        phase=deseqall0(j,27);  
        seqn(i,phase)=seqn(i,phase)+1;
        thetafire(i,:,phase,1)=thetafire(i,:,phase,1)+despikez0(j,:);    
        thetafire(i,:,phase,2)=thetafire(i,:,phase,2)+despike0(j,:);
        for h=1:3
            if h==1
                a=squeeze(linear_maze0(j,:,thetacycle2(i,63)));
            elseif h==2
                a=squeeze(linear_maze0(j,:,thetacycle2(i,63)));
                a(1:15)=squeeze(nanmean(linear_maze0(j,1:15,:),3));
            else
                a=squeeze(nanmean(linear_maze0(j,:,:),3));
            end
            thetaseq(i,:,phase,h)=thetaseq(i,:,phase,h)+a;
            thetaseq2(i,startind:endind,phase,h)=thetaseq2(i,startind:endind,phase,h)+a;
        end
        end
    end
    for i=1:size(thetaid,1)
        for phase=1:36
            thetafire(i,:,phase,1)=thetafire(i,:,phase,1)/seqn(i,phase);
            thetafire(i,:,phase,2)=thetafire(i,:,phase,2)/seqn(i,phase);
            for h=1:3
                thetaseq(i,:,phase,h)=thetaseq(i,:,phase,h)/seqn(i,phase);
                thetaseq2(i,:,phase,h)=thetaseq2(i,:,phase,h)/seqn(i,phase);
            end
        end
    end
    lia=ismember(lfptime_trial(:,5),thetacycle2(:,13));
    lfp=[lfptime_trial(lia,5),ceil(reftheta(lia,4)/10),lfpz_trial(lia,1:7)];
    thetacycle_lfp=zeros(size(thetacycle2,1),7,36);nn=zeros(size(thetacycle2,1),36);
    for i=1:size(lfp,1)
        ind=find(thetacycle2(:,13)==lfp(i,1));
        thetacycle_lfp(ind,:,lfp(i,2))=squeeze(thetacycle_lfp(ind,:,lfp(i,2)))+lfp(i,3:9);
        nn(ind,lfp(i,2))=nn(ind,lfp(i,2))+1;
    end
    for h=1:7
        thetacycle_lfp(:,h,:)=squeeze(thetacycle_lfp(:,h,:))./nn;
    end
    if ifpeak
        cd(savedir);save('theta_peak_peak','thetacycle_lfp','thetafire','thetaseq','thetaseq2','thetacycle2','thetacycle','deseqall','seqn','-v7.3');%'thetaseqslope',
    else
        cd(savedir);save('theta_trough_trough','thetacycle_lfp','thetafire','thetaseq','thetaseq2','thetacycle2','thetacycle','deseqall','seqn','-v7.3');%'thetaseqslope',
    end
    in
end

homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\figure_v2_data';
figdir='X:\Mengni\Data_Analysis\Paper_Figures\figures_v2';
timewin=0.01;reward_dis_cut=70;min_dewinN=1;ifpeak=0;phasebin=10;
thetacycle2all=[];thetacycleconnect=[];thetacycleconnect2=[];thetafireconnect=[];thetalfpconnect=[];
thetaseqconnect=[];thetaseq2connect=[];thetacycleconnectv=[];
pos1=1;pos2=101;h=1;duon=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    if ifpeak
        load('theta_peak_peak');
    else
        load('theta_trough_trough');
    end
    thetacycle2(:,65)=in;
    %thetacycle2(:,61)= thetacycle2(:,26)==thetacycle2(:,59) ;
    %thetacycle2(:,61)=thetacycle2(:,61) | thetacycle2(:,26)==thetacycle2(:,63) ;%| (thetacycle2(:,26)==0 & thetacycle2(:,67)>0.1);
    %thetacycle2(:,61)=thetacycle2(:,61) | (thetacycle2(:,59)==0 & thetacycle2(:,67)>0.2);| thetacycle2(:,59)==0;%
    %ind=thetacycle2(:,67)>0.1 & thetacycle2(:,61)==0;thetacycle2(ind,61)=2;
    ifconnect=thetacycle2(2:end,1)-thetacycle2(1:end-1,2);
    indconnect=find(ifconnect==1);
    cycleconnect=thetacycle2(indconnect,:);
    cycleconnect2=thetacycle2(indconnect+1,:);
    %indsame= cycleconnect(:,61)==cycleconnect2(:,61) & cycleconnect(:,68)==cycleconnect2(:,68) ;
    indsame= ( cycleconnect(:,63)==cycleconnect2(:,63) | cycleconnect(:,59)==0 | cycleconnect2(:,59)==0 ) ...
     & cycleconnect(:,62)==cycleconnect2(:,62) ;
    cycleconnect=cycleconnect(indsame,:);
    cycleconnect2=cycleconnect2(indsame,:);
    cycleconnectv=(cycleconnect(:,10)+cycleconnect2(:,10))/2;
    seqconnect=cat(3,thetaseq(indconnect(indsame),:,:,h),thetaseq(indconnect(indsame)+1,:,:,h));
    seq2connect=cat(3,thetaseq2(indconnect(indsame),pos1:pos2,:,h),thetaseq2(indconnect(indsame)+1,pos1:pos2,:,h));
    fireconnect=cat(3,thetafire(indconnect(indsame),:,:,:),thetafire(indconnect(indsame)+1,:,:,:));
    lfpconnect=cat(3,thetacycle_lfp(indconnect(indsame),:,:),thetacycle_lfp(indconnect(indsame)+1,:,:));
    indboth=unique([indconnect(indsame);indconnect(indsame)+1]);
    thetacycle2=thetacycle2(indboth,:);
    thetacycle2all=[thetacycle2all;thetacycle2];
    thetafireconnect=cat(1,thetafireconnect,fireconnect);
    thetalfpconnect=cat(1,thetalfpconnect,lfpconnect);
    thetaseqconnect=cat(1,thetaseqconnect,seqconnect);
    thetaseq2connect=cat(1,thetaseq2connect,seq2connect);
    thetacycleconnect=[thetacycleconnect;cycleconnect];
    thetacycleconnect2=[thetacycleconnect2;cycleconnect2];
    thetacycleconnectv=[thetacycleconnectv;cycleconnectv];
    duon(in,:)=[mean(indsame),mean(indboth)];
    in
end
thetaseq2connect_pos=nan*ones(size(thetaseq2connect,1),72);
for phase=1:72
    [M,I]=max(thetaseq2connect(:,:,phase),[],2);
    ind=M>0;
    thetaseq2connect_pos(ind,phase)=I(ind)-51+pos1-1;
end
if ifpeak
    cd(datadir);save(['fig3_thetaseq_peak_peak'],"thetaseq2connect","thetaseqconnect","thetalfpconnect","thetafireconnect","thetacycle2all",...
        'thetaseq2connect_pos','thetacycleconnect','thetacycleconnect2',"thetacycleconnectv",'duon','pos1','pos2','h','-v7.3');
else
    cd(datadir);save(['fig3_thetaseq'],"thetaseq2connect","thetaseqconnect","thetalfpconnect","thetafireconnect","thetacycle2all",...
        'thetaseq2connect_pos','thetacycleconnect','thetacycleconnect2',"thetacycleconnectv",'duon','pos1','pos2','h','-v7.3');
end
%else
    %cd(datadir);load('fig3_thetaseq');

%--------------------------------------------------------------------------------
% below quantify decode location of thetaseq
thecyclecurrent=(thetacycleconnect(:,67)+thetacycleconnect2(:,67))/2;
deposphase_dist=nan*ones(16,4,4,3,101,72);deposphase=nan*ones(16,4,4,3,72,5);
nn_session=[];seqmean_session=[];firemean_session=[];lfpmean_session=[];thetadur_session=[];
for in=1:17
    for b=1:4
        for m=1:3
            indm= ceil(3*thecyclecurrent)==m;
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
                if in<17
                    ind0= indv & indm & thetacycleconnect(:,62)==b & thetacycleconnect(:,42)==in;
                else
                    ind0= indv & indm & thetacycleconnect(:,62)==b ;
                end
                nn_session(in,b,v,m)=sum(ind0);
                thetadur_session(in,b,v,m)=nanmean(thetacycleconnect(ind0,3)+thetacycleconnect2(ind0,3),1)/2;
                seqmean_session(in,b,v,m,:,:)=squeeze(nanmean(thetaseq2connect(ind0,:,:),1));
                firemean_session(in,b,v,m,:,:,:)=squeeze(nanmean(thetafireconnect(ind0,:,:,:),1));
                lfpmean_session(in,b,v,m,:,:)=squeeze(nanmean(thetalfpconnect(ind0,:,:),1));
                for phase=1:72
                    a=histcounts(thetaseq2connect_pos(ind0,phase),[-50.5:1:50.5]);a=a/sum(a);
                    deposphase_dist(in,b,v,m,:,phase)=a;
                    deposphase(in,b,v,m,phase,1:2)=[mean(~isnan(thetaseq2connect_pos(ind0,phase))),nanstd(thetaseq2connect_pos(ind0,phase))];
                end
                seq1=imgaussfilt(squeeze(deposphase_dist(in,b,v,m,:,:)),3);
                for phase=1:72
                    a=seq1(41:61,phase);a=a/sum(a);
                    deposphase(in,b,v,m,phase,5)=[41:61]*a-51;
                    a=seq1(1:101,phase);[M,I]=max(a);
                    deposphase(in,b,v,m,phase,3:4)=[I-51,M];
                end
            end
        end
    end
end
bname={'reward ','center ','in ','out '};mname={'low','mid','high'};%'local ','remote '};
seqcut=[];
figure;
for in=1:17
    for b=1:4
        for m=1:3
            as=[];bs=[];cs=[];cut=[];lfp=[];as2=[];
            for v=1:4
                lfp1=squeeze(lfpmean_session(in,b,v,m,:,:));
                lfp=[lfp,lfp1]; 
                bs=[bs,squeeze(deposphase(in,b,v,m,:,2))'];
                cs=[cs,squeeze(deposphase(in,b,v,m,:,4))'];
                cuts=[];
                for i=1:36
                    cuts(i)=deposphase(in,b,v,m,i,4)+deposphase(in,b,v,m,i+36,4);
                end
                %cuts=gaussian_smooth(cuts);
                [M,I]=min(cuts);
                if M==0
                    I=nan;
                end
                seqcut(in,b,v,m,1:3)=[I,M,(max(cuts)-min(cuts))/max(cuts)];
                cut(v,:)=[I+72*(v-1),I+36+72*(v-1)];

                a=imgaussfilt(squeeze(deposphase_dist(in,b,v,m,:,:)),2);as=[as,a];
                a=a(:,I+1:I+36);
                [r0,c0]=find(a==max(a(:)));
                [M,I2]=max(a(51,:));
                seqcut(in,b,v,m,4:6)=[rem((mean(I2)+I)*10-5,360),rem((mean(c0)+I)*10-5,360),(mean(r0)-51)*2];
                a=squeeze(seqmean_session(in,b,v,m,:,:));
                as2=[as2,a./nansum(a,1)];
            end
            subplot(4,3,m+(b-1)*3);
            %imagesc([],[ ],[imgaussfilt(as(31:81,:),2);imgaussfilt(as2(31:81,:),0.1)]);
            imagesc([],[-50 50],as);
            set(gca,'YDir','normal');xline([1:4]*72,'k');
            %hold on;plot(lfp(4,:)*3,'r');
            for v=1:4
                cut1=seqcut(in,b,v,m,1);
                %hold on;plot([cut1+1,cut1+36]+(v-1)*72,squeeze(bestfit_session(in,b,v,m,[2,3]))/2,'m');
            end
            hold on;plot(bs,'b');hold on;
            plot(cs*500,'r');xline(cut(:),'r');
            title([bname{b},mname{m},num2str([in,squeeze(nn_session(in,b,:,m))'])]);%,squeeze(seqcut(in,b,:,m,3))'])]);
            if b~=2
                ylim([-10 10]*3);
            else
                ylim([-30 50]);
            end
        end
    end
    pause;clf;
end
seqcut(:,:,:,:,7)=seqcut(:,:,:,:,1)*10-5;

% below quantify seq length & quantify mean fire phase
seqlength_session=[];meanfirephase_session=[];
for in=1:16
    for b=1:4
        for m=1:2
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
                if in<17
                    ind0= indv & thetacycleconnect(:,61)==2-m & thetacycleconnect(:,62)==b & thetacycleconnect(:,42)==in;
                else
                    ind0= indv & thetacycleconnect(:,61)==2-m & thetacycleconnect(:,62)==b ;
                end
                cut1=seqcut(in,b,v,m,1);
                thetaseq1=squeeze(thetaseqconnect(ind0,:,cut1+1:cut1+36));
                thetafire1=squeeze(thetafireconnect(ind0,:,cut1+1:cut1+36,2));
                indq=sum(mean(~isnan(thetaseq1),2)>0.5,3)>20;
                thetaseq1=thetaseq1(indq,:,:);
                thetafire1=thetafire1(indq,:,:);
                if sum(indq)>0
                [seqslope]=quantify_shuffle_thetaseq(thetaseq1);
                %[seqslope,seqscore]=quantify_thetaseq(thetaseq1);seqslope(:,3:4)=seqscore;
                seqslope(:,5)=max(0,seqslope(:,1)); % start location
                seqslope(:,5)=min(50,seqslope(:,5)); % start location
                seqslope(:,6)=max(0,seqslope(:,1)+seqslope(:,2)*360); % end location
                seqslope(:,6)=min(50,seqslope(:,6)); % end location
                seqslope(:,7)=seqslope(:,6)-seqslope(:,5); % path length
                indq2=seqslope(:,4)<=0.05;
                a=histcounts(seqslope(indq2,7),[-50.5:50.5]);a=a(:)'/sum(a);
                seqlength_session(in,b,v,m,:)=[sum(indq),mean(indq2),2*nanmean(abs(seqslope(indq2,7))),sum(a(52:end))-sum(a(1:50)),a];
                alpha=([cut1+1:cut1+36]*10-5)*pi/180;   
                for h=1:6
                    meanfire1=[];
                    for i=1:size(thetafire1,1)
                        w=squeeze(thetafire1(i,h,:));
                        ind=~isnan(w);
                        [mu,r] = circ_mean(alpha(ind), w(ind));
                        meanfire1(i,:)=[mu*180/pi,r];
                    end
                    ind=~isnan(meanfire1(:,2));
                    meanfirephase_session(in,b,v,m,h,:)=[sum(ind),sum(ind&indq2),circ_meand(meanfire1(ind,1)),nanmean(meanfire1(ind,2)),...
                        circ_meand(meanfire1(indq2&ind,1)),nanmean(meanfire1(indq2&ind,2))];
                end
                else
                    seqlength_session(in,b,v,m,:)=nan;
                    meanfirephase_session(in,b,v,m,h,:)=nan;
                end
            end
        end
    end
    in
end

cd(datadir);save(['fig3_thetaseq_figure'],'nn_session','firemean_session','lfpmean_session','seqmean_session','meanfirephase_session',...
        "seqlength_session",'thetadur_session','deposphase','deposphase_dist','seqcut','-v7.3');
end
%--------------------------------------------------------------------------------
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\figure_v2_data';
figdir='X:\Mengni\Data_Analysis\Paper_Figures\figures_v2';
cd(datadir);load('fig3_thetaseq_figure');
% below plot figures: mean thetaseq
pn1=12;pn2=8;ifdist=1;thetacycleN_cut=200;
for b=1:4
    for m=1:2
        as=[];bs=[];lfp=[];
        for v=1:4
           ind=17;
           if ifdist
               a=imgaussfilt(squeeze(nanmean(deposphase_dist(ind,b,v,m,:,:),1)),2);
           else
               a=squeeze(nanmean(seqmean_session(ind,b,v,m,:,:)./nansum(seqmean_session(ind,b,v,m,:,:),5),1));
           end
           as=[as,a];
           b1=squeeze(nanmean(firemean_session(ind,b,v,m,:,:,1),1));
           bs=[bs,b1];
           lfp1=squeeze(nanmean(lfpmean_session(ind,b,v,m,:,:),1));
           lfp=[lfp,lfp1]; 
        end
        figure(1);subplot(pn1,pn2,[1:pn2/2]+(m-1)*pn2/2+(b-1)*pn2*3);
        hold on;plot(lfp(1,:)*5+5,'k');%hold on;plot(lfp(4,:),'k');xline([1:4]*72,'k');
        hold on;plot(bs(3,:)*15+5,'r');%hold on;plot(bs(4,:)*15,'r');
        xlim([0 4*72]);ylabel('zscore lfp');xticks([]);%title('DG oml mean lfp');
        FIG_INDEX=['fig3sup_mean_seq_',bname{b},'_',mname{m},"_1"];save_fig(FIG_INDEX,ifsavefig);

        figure(1);subplot(pn1,pn2,[1:pn2/2,1+pn2:pn2/2+pn2]+(m-1)*pn2/2+(b-1)*pn2*3+pn2);
        imagesc([],[-100 100],as);set(gca,'YDir','normal');
        if b==1
            ylim([-30 30]);
        elseif b==2
            ylim([-40 80]);
        elseif b==3
            ylim([-60 60]);
        elseif b==4
            ylim([-60 60]);
        end
        xline([1:4]*72,'k');colormap('hot');%colorbar;
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
        FIG_INDEX=['fig3sup_mean_seq_',bname{b},'_',mname{m},"_2"];save_fig(FIG_INDEX,ifsavefig,'hot');
    end
end
figure(1);figure_title='fig3sup_mean_seq';save_current_figure(figure_title);
figure;
fireindex_session=[];pkall=[];firecut=[];bimodality=[];peakfire=[];mnum=3;
for in=1:17
    for b=1:4
        for m=1:mnum
            for v=1:4
                for h=1:6
                    if nn_session(in,b,v,m)>50
                        for i=1:2
                            a=squeeze(firemean_session(in,b,v,m,h,:,3-i));a=gaussian_smooth(a);
                            fireindex_session(in,b,v,m,h,:,i)=a;
                        end
                        ai=(a-min(a))/(max(a)-min(a));
                        fireindex_session(in,b,v,m,h,:,3)=ai;
                        [M,I]=min(a);
                        peakfire(in,b,v,m,h,:)=[M,rem(I*10-5,360)];
                        if h==3
                        a=squeeze(fireindex_session(in,b,v,m,h,:,:));
                        if I>=36 
                            cutf=I-35;
                        else
                            cutf=I;
                        end
                        a1=a(cutf:cutf+35,:);
                        [pks,locs]=findpeaks(a1(:,2),'MinPeakDistance',7);%,'MinPeakWidth',2
                        pk1=[locs+cutf-1,pks,locs,a1(locs,1)];
                        lc=pk1(:,1)-36*(pk1(:,1)>36);
                        pk1(:,5)=lc;
                        pk1(:,6)=0;
                        ind= lc>=7 & lc<=28 ;pk1(ind,6)=1; % major peak
                        pk1(~ind,6)=2; % minor peak
                        ind=pk1(:,6)>0;pk1=pk1(ind,:);
                        for k=1:2
                            ind=pk1(:,6)==k;
                            if sum(ind)==1
                                bimodality(in,b,v,m,k,1:5)=pk1(ind,1:5);
                            elseif sum(ind)==0
                                bimodality(in,b,v,m,k,1:5)=nan;
                            elseif sum(ind)>1
                                pk2=pk1(ind,:);
                                [M,I]=max(pk2(:,2));
                                bimodality(in,b,v,m,k,1:5)=pk2(I,1:5);
                            end
                        end
                        bimodality(in,b,v,m,:,7)=[nanmean(a(46:59,2)),nanmean(a(33:39,2))];
                        bimodality(in,b,v,m,:,6)=bimodality(in,b,v,m,:,4)./bimodality(in,b,v,m,1,4);
                        %nanmean(a(33:39,1))./nanmean(a(46:59,1));%bimodality(in,b,v,m,:,4)./bimodality(in,b,v,m,1,4);
                        pkall=[pkall;pk1];
                        firecut(in,b,v,m,:)=[nn_session(in,b,v,m),cutf,length(pks)];
                        if 0
                        subplot(4,2,b+(m-1)*4);
                        hold on;plot([1:72]+(v-1)*72,a,'ko-');xline(I+(v-1)*72,'m');
                        hold on;plot(pk1(:,1)+(v-1)*72,pk1(:,2),'r*');xline([1:4]*72);xlim([0 72*4]);
                        title([bname{b},mname{m},num2str([in,squeeze(firecut(in,b,:,m,3))'])]);%xline([1:4]*72,'m');xlim([0 72*4]);
                        %yticks([0 1]);yticklabels({'CA3 P','CA1 P'});
                        end
                        end
                    end
                end
            end
        end
    end
    %pause;clf;
end

spikeprop=cat(5,seqcut(:,:,:,:,[7,5]),1./thetadur_session,squeeze(bimodality(:,:,:,:,:,2)),peakfire(:,:,:,:,4:6,1),squeeze(bimodality(:,:,:,:,2,6)));
figure(4);mu=[];
for b=1:4
    for m=1:mnum
        for v=1:4
            for k=1:9
                ind1=find((~isnan(spikeprop(1:16,b,v,m,k))) & nn_session(1:16,b,v,m)>50);
                n=length(ind1);
                if n>0
                    if k<3
                        mu(b,v,m,k)=rem(circ_meand(squeeze(spikeprop(ind1,b,v,m,k))),360);
                    else
                        mu(b,v,m,k)=squeeze(nanmean(spikeprop(ind1,b,v,m,k),1));
                    end
                    subplot(5,2,k);
                    if m==1
                        hold on;plot((v+(b-1)*4)*ones(n,1)-0.2,squeeze(spikeprop(ind1,b,v,m,k)),'bo');
                    elseif m==2
                        hold on;plot((v+(b-1)*4)*ones(n,1),squeeze(spikeprop(ind1,b,v,m,k)),'ro');
                    elseif m==3
                        hold on;plot((v+(b-1)*4)*ones(n,1)+0.2,squeeze(spikeprop(ind1,b,v,m,k)),'ko');
                    end
                else
                    mu(b,v,m,k)=nan;
                end
            end
        end
    end
end
kname={'phase onset ','max phase ','theta freq ','major ','minor ','ca3 p','ca1 i','ca3 i','minor/major'};
for k=1:9
    figure(4);subplot(5,2,k);
    xline([1:4]*4+0.5,'k');
    a=squeeze(mu(:,:,1,k))';hold on;plot([1:16],a(:),'b-');
    a=squeeze(mu(:,:,2,k))';hold on;plot([1:16],a(:),'r-');
    a=squeeze(mu(:,:,3,k))';hold on;plot([1:16],a(:),'k-');
    xticks([1:16]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: 0-1','1-5','5-10','>10',...
        'out: 0-1','1-5','5-10','>10'});
    title([kname{k}]);xlim([0 17]);
end

figure;mu=[];
for b=1:4
    for m=1:2
        for v=1:4
            for h=4:6
                ind1=find((~isnan(peakfire(1:16,b,v,m,h,1))) & nn_session(1:16,b,v,m)>50);
                n=length(ind1);
                if n>0
                    mu(b,v,m,h,:)=[squeeze(nanmean(peakfire(ind1,b,v,m,h,1),1))',...
                    rem(circ_meand(squeeze(peakfire(ind1,b,v,m,h,2))),360)];
                    for k=1:2
                    subplot(3,2,k+(h-4)*2);
                    if m==1
                        hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(peakfire(ind1,b,v,m,h,k)),'bo');
                    else 
                        hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(peakfire(ind1,b,v,m,h,k)),'ro');
                    end
                    end
                else
                    mu(b,v,m,h,:)=nan;
                end
            end
        end
    end
end
hname={'all spike n','all cell n','ca1 p','ca3 p','ca1 i','ca3 i'};
kname={'peak rate ','peak phase '};
for h=4:6
    for k=1:2
        subplot(3,2,k+(h-4)*2);
        xline([1:4]*4+0.5,'k');
        a=squeeze(mu(:,:,1,h,k))';hold on;plot([1:16],a(:),'b-');
        a=squeeze(mu(:,:,2,h,k))';hold on;plot([1:16],a(:),'r-');
        xticks([1:16]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: 0-1','1-5','5-10','>10',...
        'out: 0-1','1-5','5-10','>10'});
        title([kname{k},hname{h}]);xlim([0 17]);
    end
    FIG_INDEX=['fig3sup_seq_spike_only_',hname{h},"_",kname{k}];save_fig(FIG_INDEX,ifsavefig);
end
figure(4);figure_title='fig3sup_seq_spike_only';save_current_figure(figure_title);

% below quanfify major minor firings and correlation with thetaseq
thetafire=squeeze(firemean_session(17,4,4,1,3,:,1));
[M,I]=min(thetafire);thetafire1=thetafire(28:63);
indmajor=[46:58];

    
    figure;ns=[];
for b=1:4
    for m=1:2
        for h=3:4
            bs=[];err=[];cuts=[];cen=[];
            for v=1:4
                ind=find(nn_session(1:16,b,v,m)>thetacycleN_cut & (~isnan(fireindex_session(1:16,b,v,m,h,1,1))));
                b1=squeeze(nanmean(fireindex_session(ind,b,v,m,h,:,2),1))';n=length(ind);ns(b,v,m,h)=n;
                bs=[bs,b1];
                err=[err,squeeze(nanstd(fireindex_session(ind,b,v,m,h,:,2),0,1))'*n^(-0.5)];
                %cut1=seqcut(ind,b,v,m,1);
                %cuts=[cuts,[cut1+1,cut1+36]+(v-1)*72];
                cen1=circ_meand(seqcut(ind,b,v,m,5))/10;
                cen=[cen,[cen1,cen1+36]+(v-1)*72];
            end
            subplot(4,1,b);
            if m==1
                hold on;errorbar([1:72*4],bs+(h),err,'b');xline(cen,'b');%xline(cuts,'m');
            else
                hold on;errorbar([1:72*4],bs+(h),err,'r');xline(cen,'r');%xline(cuts,'m--');
            end
            yticks([3 4]);yticklabels({'CA1 P','CA3 P'});
            %yticks([1:6]*20);yticklabels({'all spike n','all cell n','ca1 p','ca3 p','ca1 i','ca3 i'});
            title(bname{b});xline([1:4]*72,'k');
        end
    end
end

figure;ns=[];hname={'all spike n','all cell n','ca1 p','ca3 p','ca1 i','ca3 i'};
for b=1:4
    for m=1:2
        for h=3:4
            bs=[];err=[];cuts=[];cen=[];
            for v=1:4
                ind=find(nn_session(1:16,b,v,m)>thetacycleN_cut & (~isnan(fireindex_session(1:16,b,v,m,h,1,1))));
                n=length(ind);ns(b,v,m,h)=n;
                if n>0
                b1=squeeze(nanmean(fireindex_session(ind,b,v,m,h,:,2),1))';
                bs=[bs,b1];
                err1=squeeze(nanstd(fireindex_session(ind,b,v,m,h,:,2),0,1))'*n^(-0.5);
                err=[err,err1];
                %cut1=seqcut(ind,b,v,m,1);
                %cuts=[cuts,[cut1+1,cut1+36]+(v-1)*72];
                cen1=ceil(circ_meand(seqcut(ind,b,v,m,5))/10);
                cen=[cen1,cen1+36];
           
                subplot(4,4,(b-1)*4+m+(h-3)*2);
                hold on;errorbar([1:72],b1+v,err1,'k');
                hold on;plot(cen,b1(cen)+v,'r*');
                end
            end
            yticks([1:4]);yticklabels({'0-1','1-5','5-10','>10'});
            title([bname{b},mname{m},hname{h}]);
        end
    end
end

figure;ns=[];hname={'all spike n','all cell n','ca1 p','ca3 p','ca1 i','ca3 i'};
for b=1:4
    for m=1:2
        for h=3%:4
            bs=[];err=[];cuts=[];cen=[];
            for v=1:4
                ind=find(nn_session(1:16,b,v,m)>thetacycleN_cut & (~isnan(fireindex_session(1:16,b,v,m,h,1,1))));
                n=length(ind);ns(b,v,m,h)=n;
                if n>0
                b1=squeeze(nanmean(fireindex_session(ind,b,v,m,h,:,2),1))';
                bs=[bs,b1];
                err1=squeeze(nanstd(fireindex_session(ind,b,v,m,h,:,2),0,1))'*n^(-0.5);
                err=[err,err1];
                %cut1=seqcut(ind,b,v,m,1);
                %cuts=[cuts,[cut1+1,cut1+36]+(v-1)*72];
                cen1=ceil(circ_meand(seqcut(ind,b,v,m,5))/10);
                cen=[cen1,cen1+36];
           
                subplot(4,2,(b-1)*2+m+(h-3)*4);
                hold on;plot([1:72],b1+v/5,'o-');
                hold on;plot(cen,b1(cen),'r*');
                end
            end
            yticks([1:4]);yticklabels({'0-1','1-5','5-10','>10'});
            title([bname{b},mname{m},hname{h}]);xlim([0 72]);
        end
    end
end
legend('0-1','1-5','5-10','>10');


figure;thetacycleN_cut=10;err=[];
for in=1:16
for b=1:4
    for m=1:2
        bs=[];m1=[];m2=[];
        for v=1:4
           b1=squeeze(fireindex_session(in,b,v,m,h,:,1));
           bs=[bs,b1];
           m1=[m1,bimodality(in,b,v,m,3,:)]
        end
        subplot(4,2,b+(m-1)*2);
        hold on;plot(bs(3,:)*20+20,'k');
        hold on;plot(bs(4,:)*20,'k');
        xline([1:4]*72,'k');
    end
end
end

figure;
for b=1:4
    for m=1:2
        bs=[];lfp=[];
        for v=1:4
           ind=17;%find(nn_session(1:16,b,v,m)>thetacycleN_cut );%17;
           b1=squeeze(nanmean(firemean_session(ind,b,v,m,:,:,1),1));%b1=gaussian_smooth(b1(3,:))';%b1=(b1-min(b1))/(max(b1)-min(b1));
           bs=[bs,b1];
           lfp1=squeeze(nanmean(lfpmean_session(ind,b,v,m,:,:),1));%lfp1=gaussian_smooth(lfp1(1,:))';lfp1=(lfp1-min(lfp1))/(max(lfp1)-min(lfp1));
           lfp=[lfp,lfp1]; 
        end
        subplot(4,1,b);
        if m==1
            hold on;plot(bs(3,:)*20+20,'b');
            hold on;plot(bs(4,:)*20,'b');
        else
            hold on;plot(bs(3,:)*20+20,'r');
            hold on;plot(bs(4,:)*20,'r');
        end
        yticks([0 20]);yticklabels({'CA3 P','CA1 P'});title(bname{b});
        xline([1:4]*72,'k');
    end
end


%------

% below phase onset
mu=[];thetadur_session2=1./thetadur_session;
for b=1:4
    for m=1:2
        for v=1:4
            ind=find(nn_session(1:16,b,v,m)>thetacycleN_cut );
            n=length(ind);
            if n>0
                mu(b,v,m,:)=[circ_meand(squeeze(seqcut(ind,b,v,m,7))),squeeze(nanmean(thetadur_session2(ind,b,v,m,1))),...
                    squeeze(nanmean(nn_session(ind,b,v,m))),circ_meand(squeeze(seqcut(ind,b,v,m,4))),...
                    circ_meand(squeeze(seqcut(ind,b,v,m,5))),squeeze(nanmean(seqlength_session(ind,b,v,m,2))),...
                    squeeze(nanmean(seqlength_session(ind,b,v,m,3))),squeeze(nanmean(seqlength_session(ind,b,v,m,4)))];
            if m==1
                figure(2);subplot(4,2,1);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,7)),'bo');
                %hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,1))*10+355,'bo');
                figure(2);subplot(4,2,2);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(thetadur_session2(ind,b,v,m,1)),'bo');
                figure(2);subplot(4,2,3);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(nn_session(ind,b,v,m)),'bo');
                figure(2);subplot(4,2,4);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,4)),'bo');
                figure(2);subplot(4,2,5);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,5)),'bo');
                figure(2);subplot(4,2,6);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqlength_session(ind,b,v,m,2)),'bo');
                figure(2);subplot(4,2,7);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqlength_session(ind,b,v,m,3)),'bo');
                figure(2);subplot(4,2,8);
                hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(seqlength_session(ind,b,v,m,4)),'bo');
            else
                figure(2);subplot(4,2,1);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,7)),'ro');
                %hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,1))*10+355,'ro');
                figure(2);subplot(4,2,2);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(thetadur_session2(ind,b,v,m,1)),'ro');
                figure(2);subplot(4,2,3);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(nn_session(ind,b,v,m)),'ro');
                figure(2);subplot(4,2,4);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,4)),'ro');
                figure(2);subplot(4,2,5);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,5)),'ro');
                figure(2);subplot(4,2,6);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqlength_session(ind,b,v,m,2)),'ro');
                figure(2);subplot(4,2,7);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqlength_session(ind,b,v,m,3)),'ro');
                figure(2);subplot(4,2,8);
                hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(seqlength_session(ind,b,v,m,4)),'ro');
            end
            else
                mu(b,v,m,:)=[nan,nan,n,nan,nan,nan,nan,nan];
            end
        end
    end
end
hname={'thetaseq onset phase','theta freq','thetacycle N','local phase','max phase','sig seq percent','seq length (cm)','dir: out-in'};
for h=1:8
    figure(2);subplot(4,2,h);
    xline([1:4]*4+0.5,'k');
    a=squeeze(mu(:,:,1,h))';hold on;plot([1:16],a(:),'b-');
    a=squeeze(mu(:,:,2,h))';hold on;plot([1:16],a(:),'r-');
    if 0%h==1
    a=squeeze(seqcut(17,:,:,1,7))';hold on;plot([1:16],a(:),'b--');
    a=squeeze(seqcut(17,:,:,2,7))';hold on;plot([1:16],a(:),'r--');
    end
    xticks([1:16]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: 0-1','1-5','5-10','>10',...
        'out: 0-1','1-5','5-10','>10'});
    ylabel(hname{h});xlim([0 17]);
    FIG_INDEX=['fig3sup_seq_quantification_',hname{h}];save_fig(FIG_INDEX,ifsavefig);
end
figure(2);figure_title='fig3sup_seq_quantification';save_current_figure(figure_title);

% below fire phase
figure(3);mu=[];
for b=1:4
    for m=1:2
        for v=1:4
            for h=1:6
                ind=find(meanfirephase_session(1:16,b,v,m,h,2)>100 );
                n=length(ind);
                if n>0
                    mu(b,v,m,h,:)=[rem(circ_meand(squeeze(meanfirephase_session(ind,b,v,m,h,5))),360),...
                    squeeze(nanmean(meanfirephase_session(ind,b,v,m,h,6)))];
                    if m==1
                        subplot(6,2,1+(h-1)*2);
                        hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5)),'bo');
                        %hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5))+360,'bo');
                        subplot(6,2,2+(h-1)*2);
                        hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(meanfirephase_session(ind,b,v,m,h,6)),'bo');
                    else 
                        subplot(6,2,1+(h-1)*2);
                        hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5)),'ro');
                        %hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5))+360,'ro');
                        subplot(6,2,2+(h-1)*2);
                        hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(meanfirephase_session(ind,b,v,m,h,6)),'ro');
                    end
                else
                    mu(b,v,m,h,:)=[nan,nan];
                end
            end
        end
    end
end
hname={'all spike n','all cell n','ca1 p','ca3 p','ca1 i','ca3 i'};
kname={'mean spike phase','resultant length'};
for h=1:6
    for k=1:2
        figure(3);subplot(6,2,k+(h-1)*2);
        xline([1:4]*4+0.5,'k');
        a=squeeze(mu(:,:,1,h,k))';hold on;plot([1:16],a(:),'b-');
        a=squeeze(mu(:,:,2,h,k))';hold on;plot([1:16],a(:),'r-');
        xticks([1:16]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: 0-1','1-5','5-10','>10',...
        'out: 0-1','1-5','5-10','>10'});
        title([hname{h},kname{k}]);xlim([0 17]);
        if 0%k==1
            a=squeeze( meanfirephase_session(17,:,:,1,h,3))';hold on;plot([1:16],a(:),'b--');
            a=squeeze( meanfirephase_session(17,:,:,2,h,3))';hold on;plot([1:16],a(:),'r--');
        end
    end
    FIG_INDEX=['fig3sup_seq_spike_',hname{h},"_",kname{k}];save_fig(FIG_INDEX,ifsavefig);
end
figure(3);figure_title='fig3sup_seq_spike';save_current_figure(figure_title);
figure(4);mu=[];
for b=1:4
    for m=1:2
        for v=1:4
            for h=1:6
                ind=find(meanfirephase_session(1:16,b,v,m,h,2)>100 );
                n=length(ind);
                if n>0
                    mu(b,v,m,h,:)=[rem(circ_meand(squeeze(meanfirephase_session(ind,b,v,m,h,5))),360),...
                    squeeze(nanmean(meanfirephase_session(ind,b,v,m,h,6)))];
                    if m==1
                        subplot(3,2,h);
                        hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5)),'bo');
                        %hold on;plot((v+(b-1)*4)*ones(n,1)-0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5))+360,'bo');
                    else 
                        subplot(3,2,h);
                        hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5)),'ro');
                        %hold on;plot((v+(b-1)*4)*ones(n,1)+0.1,squeeze(meanfirephase_session(ind,b,v,m,h,5))+360,'ro');
                    end
                else
                    mu(b,v,m,h,:)=[nan,nan];
                end
            end
        end
    end
end
hname={'all spike n','all cell n','ca1 p','ca3 p','ca1 i','ca3 i'};
kname={'mean spike phase','resultant length'};
for h=1:6
    for k=1
        figure(4);subplot(3,2,h);
        xline([1:4]*4+0.5,'k');
        a=squeeze(mu(:,:,1,h,k))';hold on;plot([1:16],a(:),'b-');
        a=squeeze(mu(:,:,2,h,k))';hold on;plot([1:16],a(:),'r-');
        xticks([1:16]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: 0-1','1-5','5-10','>10',...
        'out: 0-1','1-5','5-10','>10'});
        title([hname{h},kname{k}]);xlim([0 17]);
        if 0%k==1
            a=squeeze( meanfirephase_session(17,:,:,1,h,3))';hold on;plot([1:16],a(:),'b--');
            a=squeeze( meanfirephase_session(17,:,:,2,h,3))';hold on;plot([1:16],a(:),'r--');
        end
    end
    FIG_INDEX=['fig3sup_seq_spike_only_',hname{h},"_",kname{k}];save_fig(FIG_INDEX,ifsavefig);
end
figure(4);figure_title='fig3sup_seq_spike_only';save_current_figure(figure_title);
%------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------
% below generate paper fig 3
fign=5;pn1=9;pn2=9;ifdist=1;thetacycleN_cut=200;
% below mean seq part
for b=1:2
    for m=1:2
        as=[];bs=[];lfp=[];cens=[];
        for v=1:4
           ind=17;
           if ifdist
               a=imgaussfilt(squeeze(nanmean(deposphase_dist(ind,b,v,m,:,:),1)),2);
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
           cen=rem(circ_meand(seqcut(ind,b,v,m,5)),360)/10;
           cens=[cens,[cen,cen+36]+(v-1)*72];
        end
        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2)]+(m-1)*floor(pn2/2)+(b-1)*pn2*4);
        hold on;plot(lfp(1,:)*5+5,'k');hold on;plot(lfp(4,:),'k');xline([1:4]*72,'k');
        xlim([0 4*72]);ylabel('zscore lfp');xticks([]);xline(cens,'c');
        FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},"_1"];save_fig(FIG_INDEX,ifsavefig);

        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2)]+(m-1)*floor(pn2/2)+(b-1)*pn2*4+pn2);
        for k=3:4
            hold on;plot(bs(k,:));
        end
        %hold on;plot(bs(3,:)-bs(4,:));
        xline([1:4]*72,'k');xlim([0 4*72]);xline(cens,'c');xline([36:36:36*4],'k');xline([52:72:72*4],'m');
        ylabel('zscore spike N');yticks([0:3]*5);yticklabels({'CA1 P','CA3 P','CA1 I','CA3 I'});
        FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},"_2"];save_fig(FIG_INDEX,ifsavefig);

        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2),1+pn2:floor(pn2/2)+pn2]+(m-1)*floor(pn2/2)+(b-1)*pn2*4+pn2*2);
        imagesc([],[-100 100],as);set(gca,'YDir','normal');
        if b==1
            ylim([-30 30]);
        elseif b==2
            ylim([-40 80]);
        elseif b==3
            ylim([-60 60]);
        elseif b==4
            ylim([-60 60]);
        end
        xline([1:4]*72,'k');colormap('hot');xline(cens,'c');%colorbar; 
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
        FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},"_3"];save_fig(FIG_INDEX,ifsavefig,'hot');
    end
end
for b=3:4
    for m=1
        as=[];bs=[];lfp=[];cens=[];
        for v=4
           ind=17;
           if ifdist
               a=imgaussfilt(squeeze(nanmean(deposphase_dist(ind,b,v,m,:,:),1)),2);
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
           cen=rem(circ_meand(seqcut(ind,b,v,m,5)),360)/10;
           cens=[cens,[cen,cen+36]];
        end
        figure(fign);subplot(pn1,pn2,pn2+(b-3)*pn2*4);
        hold on;plot(lfp(1,:)*5+5,'k');hold on;plot(lfp(4,:),'k');
        ylabel('zscore lfp');xticks([]);xline(cens,'c');
        FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},"_1"];save_fig(FIG_INDEX,ifsavefig);

        figure(fign);subplot(pn1,pn2,pn2*2+(b-3)*pn2*4);
        for k=3:4
            hold on;plot(bs(k,:));
        end
        %hold on;plot(bs(3,:)-bs(4,:));
        xlim([0 1*72]);xline(cens,'c');xline([36:36:36*4],'k');xline([52:72:72*4],'m');
        ylabel('zscore spike N');yticks([0:3]*5);yticklabels({'CA1 P','CA3 P','CA1 I','CA3 I'});
        FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},"_2"];save_fig(FIG_INDEX,ifsavefig);

        figure(fign);subplot(pn1,pn2,[pn2*3,pn2*4]+(b-3)*pn2*4);
        imagesc([],[-100 100],as);set(gca,'YDir','normal');
        if b==1
            ylim([-30 30]);
        elseif b==2
            ylim([-40 80]);
        elseif b==3
            ylim([-60 60]);
        elseif b==4
            ylim([-60 60]);
        end
        xline([1:4]*72,'k');colormap('hot');xline(cens,'c');%colorbar;   
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
        FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},"_3"];save_fig(FIG_INDEX,ifsavefig,'hot');
    end
end
% below seq property part
hname={'thetaseq onset phase','max phase','sig seq percent','seq length (cm)','dir: out-in','theta freq'};
mu=[];thetadur_session2=1./thetadur_session;
for b=1:4
    for m=1:2
        for v=1:4
            ind=find(nn_session(1:16,b,v,m)>thetacycleN_cut );
            n=length(ind);
            if n>0
                mu(b,v,m,:)=[circ_meand(squeeze(seqcut(ind,b,v,m,7))),circ_meand(squeeze(seqcut(ind,b,v,m,5))),...
                    squeeze(nanmean(seqlength_session(ind,b,v,m,2:4),1))',squeeze(nanmean(thetadur_session2(ind,b,v,m,1)))];
                if m==1 & b<3
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]);
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,7)),'bo');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+2);
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,5)),'bo');
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(seqcut(ind,b,v,m,5))+360,'bo');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+4);
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(seqlength_session(ind,b,v,m,2)),'bo');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+6);
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(seqlength_session(ind,b,v,m,3)),'bo');
                    figure(fign);subplot(pn1,pn2,7*pn2+[1,2]);
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(seqlength_session(ind,b,v,m,4)),'bo');
                    figure(fign);subplot(pn1,pn2,7*pn2+[1,2]+2);
                    hold on;plot((v+(b-1)*5)*ones(n,1)-0.1,squeeze(thetadur_session2(ind,b,v,m,1)),'bo');
                elseif m==2 & b<3
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]);
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,7)),'ro');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+2);
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,5)),'ro');
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(seqcut(ind,b,v,m,5))+360,'ro');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+4);
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(seqlength_session(ind,b,v,m,2)),'ro');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+6);
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(seqlength_session(ind,b,v,m,3)),'ro');
                    figure(fign);subplot(pn1,pn2,7*pn2+[1,2]);
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(seqlength_session(ind,b,v,m,4)),'ro');
                    figure(fign);subplot(pn1,pn2,7*pn2+[1,2]+2);
                    hold on;plot((v+(b-1)*5)*ones(n,1)+0.1,squeeze(thetadur_session2(ind,b,v,m,1)),'ro');
                elseif m==1 & b>=3 & v==4
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]);
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(seqcut(ind,b,v,m,7)),'bo');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+2);
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(seqcut(ind,b,v,m,5)),'bo');
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(seqcut(ind,b,v,m,5))+360,'bo');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+4);
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(seqlength_session(ind,b,v,m,2)),'bo');
                    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+6);
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(seqlength_session(ind,b,v,m,3)),'bo');
                    figure(fign);subplot(pn1,pn2,7*pn2+[1,2]);
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(seqlength_session(ind,b,v,m,4)),'bo');
                    figure(fign);subplot(pn1,pn2,7*pn2+[1,2]+2);
                    hold on;plot((11+(b-3)*2)*ones(n,1),squeeze(thetadur_session2(ind,b,v,m,1)),'bo');
                end
            else
                mu(b,v,m,:)=nan;
            end
        end
    end
end
for h=1:6
    figure(fign);subplot(pn1,pn2,6*pn2+[1,2]+(h-1)*2+(h>4));
    for b=1:2
        if h==2 & b==1
            a=squeeze(mu(b,:,1,h))';a(a<180)=a(a<180)+360;hold on;plot([1:4]+(b-1)*5,a,'b-');
            a=squeeze(mu(b,:,2,h))';hold on;plot([1:4]+(b-1)*5,a,'r-');
        elseif h==2 & b==2
            a=squeeze(mu(b,:,1,h))'+360;hold on;plot([1:4]+(b-1)*5,a,'b-');
            a=squeeze(mu(b,:,2,h))';a(a<180)=a(a<180)+360;hold on;plot([1:4]+(b-1)*5,a,'r-');
        else
            a=squeeze(mu(b,:,1,h))';hold on;plot([1:4]+(b-1)*5,a,'b-');
            a=squeeze(mu(b,:,2,h))';hold on;plot([1:4]+(b-1)*5,a,'r-');
        end
        if h==2
            a=squeeze(mu(2+b,4,1,h))'+360;hold on;plot([-0.15 0.15]+11+(b-1)*2,a*[1 1],'b-');
        else
            a=squeeze(mu(2+b,4,1,h))';hold on;plot([-0.15 0.15]+11+(b-1)*2,a*[1 1],'b-');
        end
    end
    if h==2
        ylim([180 540]);
    end
    xticks([1:4,6:9,11,13]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: >10','out: >10'});
    ylabel(hname{h});xlim([0 14]);
    FIG_INDEX=['fig3_thetaseq_',hname{h}];save_fig(FIG_INDEX,ifsavefig);
end
figure(fign);figure_title='fig3_thetaseq';save_current_figure(figure_title);