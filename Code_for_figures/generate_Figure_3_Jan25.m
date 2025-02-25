%function generate_Figure_3_Jan25(ifsavedata)

%clear all;close all;
%if nargin<1
    ifsavedata=0;
%end
close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
timewin=0.01;reward_dis_cut=70;min_dewinN=1;ifpeak=0;phasebin=10;
if 0
generate_SupFigure_7_decode_example(10,1);
inset=[1,3,5,8,10,12,13,16];
for i=1:8
    in=inset(i);
    generate_SupFigure_7_decode_example(in);
    i
end
end

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
        thetafire(i,:,phase,1)=thetafire(i,:,phase,1)+despike0(j,:);
        thetafire(i,:,phase,2)=thetafire(i,:,phase,2)+despikez0(j,:);    
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
% below for N
min_dewinN=14;
majorN=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_trough_trough','thetacycle');
    ind1=thetacycle(:,79)==0 & thetacycle(:,62)>0;
    ind=thetacycle(:,79)==0 & thetacycle(:,24)>=min_dewinN & thetacycle(:,62)>0;
    majorN(in,:)=[mean(thetacycle(ind1,24)>=min_dewinN),mean(thetacycle(ind,64)>=0.5)];
end
majorN=[];
for in=1:16
    ind=thetacycleconnect(:,42)==in;
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_trough_trough','thetacycle2');
    ifconnect=thetacycle2(2:end,1)-thetacycle2(1:end-1,2);
    indconnect=find(ifconnect==1);
    majorN(in,:)=[length(indconnect),sum(ind)]/size(thetacycle2,1);
end
majorN(:,3)=majorN(:,2)./majorN(:,1);
majorN=[];
for in=1:16
    ind=thetacycleconnect(:,42)==in;
    majorN(in,:)=squeeze(nansum(nansum(seqlength_session(in,:,:,:,1),3),4))/sum(ind);
end
majorN(:,2)=majorN(:,2)/2;
a=squeeze(nansum(nansum(seqlength_session(:,1:2,1,:,1).*seqlength_session(:,1:2,1,:,2),3),4))./squeeze(nansum(nansum(seqlength_session(:,1:2,1,:,1),3),4));

if 0
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
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
    thetacycle2(:,61)= thetacycle2(:,26)==thetacycle2(:,59) ;
    ifconnect=thetacycle2(2:end,1)-thetacycle2(1:end-1,2);
    indconnect=find(ifconnect==1);
    cycleconnect=thetacycle2(indconnect,:);
    cycleconnect2=thetacycle2(indconnect+1,:);
    indsame= ( cycleconnect(:,63)==cycleconnect2(:,63) | cycleconnect(:,59)==0 | cycleconnect2(:,59)==0 ) ...
     & cycleconnect(:,61)==cycleconnect2(:,61) & cycleconnect(:,62)==cycleconnect2(:,62) ;
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

velmean_session=[];velmean_session2=[];thetacycleconnectv2=[];pn=30;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    ind1=thetacycleconnect(:,42)==in & thetacycleconnect(:,62)<=4 ;
    cycleconnect=thetacycleconnect(ind1,:);
    cycleconnect2=thetacycleconnect2(ind1,:);
    cycleconnectv=thetacycleconnectv(ind1);
    cycle=(cycleconnect(:,6)+cycleconnect2(:,6))/2;
    ruler=Position_Data(:,1);template=cycle;[outputindex,error]=match(template,ruler,0);
    cycle(:,2)=outputindex;
    cycle_vel=nan*ones(size(cycle,1),pn*2+1);
    for i=1:size(cycleconnect,1)
        trialid=Position_Data(cycle(i,2),4);
        p1=max(1,cycle(i,2)-pn):cycle(i,2);
        p2=cycle(i,2)+1:min(size(Position_Data,1),cycle(i,2)+pn);
        p01=find(Position_Data(p1,4)==trialid,1);p1=p1(p01:end);
        p02=find(Position_Data(p2,4)==trialid,1,'last');p2=p2(1:p02);
        if max(Position_Data(p1,1))-min(Position_Data(p1,1))<1.2
            cycle_vel(i,32-length(p1):31)=Position_Data(p1,5);
        end
        if max(Position_Data(p2,1))-min(Position_Data(p2,1))<1.2
            cycle_vel(i,32:31+length(p2))=Position_Data(p2,5);
        end
    end
    a=[nanmean(cycle_vel(:,16:46),2),nanmean(cycle_vel,2)];
    thetacycleconnectv2=[thetacycleconnectv2;a];
    for b=1:4
        for m=1:2
            for v=1:4
                if v==1
                    indv=cycleconnectv>=0 & cycleconnectv<1;
                elseif v==2
                    indv=cycleconnectv>=1 & cycleconnectv<5;
                elseif v==3
                    indv=cycleconnectv>=5 & cycleconnectv<10;
                else
                    indv=cycleconnectv>=10 ;
                end
                ind0= indv & cycleconnect(:,61)==2-m & cycleconnect(:,62)==b ;
                velmean_session(in,b,m,v,:)=nanmean(cycle_vel(ind0,:),1);
                velmean_session2(in,b,m,v)=nanmean(nanmean(cycle_vel(ind0,:),2),1);%16:46
            end
        end
    end
    in
end
end
if ifpeak
    cd(datadir);save(['Figure3_thetaseq_peak_peak'],"thetaseq2connect","thetaseqconnect","thetalfpconnect","thetafireconnect","thetacycle2all",...
        'thetaseq2connect_pos','thetacycleconnect','thetacycleconnect2',"thetacycleconnectv",'duon','pos1','pos2','h',...
        'velmean_session2','velmean_session','thetacycleconnectv2','-v7.3');
else
    cd(datadir);save(['Figure3_thetaseq'],"thetaseq2connect","thetaseqconnect","thetalfpconnect","thetafireconnect","thetacycle2all",...
        'thetaseq2connect_pos','thetacycleconnect','thetacycleconnect2',"thetacycleconnectv",'duon','pos1','pos2','h',...
        'velmean_session2','velmean_session','thetacycleconnectv2','-v7.3');
end

%--------------------------------------------------------------------------------
% below quantify decode location of thetaseq
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
timewin=0.01;reward_dis_cut=70;min_dewinN=1;ifpeak=0;phasebin=10;
cd(datadir);load(['Figure3_thetaseq']);
deposphase_dist=nan*ones(16,4,4,2,101,72);deposphase=nan*ones(16,4,4,2,72,5);deposphase_dist_n=nan*ones(16,4,4,2,101,72);
nn_session=[];seqmean_session=[];firemean_session=[];lfpmean_session=[];thetaprop_session=[];
for in=1:16
    for b=1:4
        for m=1:3
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
                if m<3
                    ind0= indv & thetacycleconnect(:,61)==2-m & thetacycleconnect(:,62)==b & thetacycleconnect(:,42)==in;
                else
                    ind0= indv & thetacycleconnect(:,62)==b & thetacycleconnect(:,42)==in;
                end
                nn_session(in,b,v,m)=sum(ind0);
                thetaprop_session(in,b,v,m,:)=[nanmean(1./thetacycleconnect(ind0,3)+1./thetacycleconnect2(ind0,3),1)/2,nanmean(thetacycleconnect(ind0,4)+thetacycleconnect2(ind0,4),1)/2];
                seqmean_session(in,b,v,m,:,:)=squeeze(nanmean(thetaseq2connect(ind0,:,:),1));
                firemean_session(in,b,v,m,:,:,:)=squeeze(nanmean(thetafireconnect(ind0,:,:,:),1));
                lfpmean_session(in,b,v,m,:,:)=squeeze(nanmean(thetalfpconnect(ind0,:,:),1));
                for phase=1:72
                    a=histcounts(thetaseq2connect_pos(ind0,phase),[-50.5:1:50.5]);
                    deposphase_dist_n(in,b,v,m,:,phase)=a;
                    a=a/sum(a);
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
bname={'reward ','center ','in ','out '};mname={'local ','remote '};
seqcut=[];
figure;
for in=1:16
    for b=1:4
        for m=1:2
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
                seqcut(in,b,v,m,1:3)=[I,M,(max(cuts)-min(cuts))/max(cuts)];
                cut(v,:)=[I+72*(v-1),I+36+72*(v-1)];

                if 1
                    [M,I]=max(cuts);
                    seqcut(in,b,v,m,4:6)=[M,I*10-5,(max(cuts)-min(cuts))/max(cuts)];
                else
                    a=imgaussfilt(squeeze(deposphase_dist(in,b,v,m,:,:)),2);as=[as,a];
                    a=a(:,I+1:I+36);
                    [r0,c0]=find(a==max(a(:)));
                    [M,I2]=max(a(51,:));
                    seqcut(in,b,v,m,4:6)=[rem((mean(I2)+I)*10-5,360),rem((mean(c0)+I)*10-5,360),(mean(r0)-51)*2];
                    a=squeeze(seqmean_session(in,b,v,m,:,:));
                    as2=[as2,a./nansum(a,1)];

                    if ~(M>0)
                        I=nan;
                    end
                end

            end
            subplot(4,2,b+(m-1)*4);
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
    %pause;clf;
end
seqcut(:,:,:,:,7)=seqcut(:,:,:,:,1)*10-5;

% below quantify seq length & quantify mean fire phase
seqlength_session=[];%meanfirephase_session=[];
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
                ind0= indv & thetacycleconnect(:,61)==2-m & thetacycleconnect(:,62)==b & thetacycleconnect(:,42)==in;
                cut1=seqcut(in,b,v,m,1);
                thetaseq1=squeeze(thetaseqconnect(ind0,:,cut1+1:cut1+36));
                %thetafire1=squeeze(thetafireconnect(ind0,:,cut1+1:cut1+36,2));
                indq=sum(mean(~isnan(thetaseq1),2)>0.5,3)>20;
                thetaseq1=thetaseq1(indq,:,:);
                %thetafire1=thetafire1(indq,:,:);
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
                if 0
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
                end
                else
                    seqlength_session(in,b,v,m,:)=nan;
                    %meanfirephase_session(in,b,v,m,h,:)=nan;
                end
            end
        end
    end
    in
end

% below quantify spike prop
fireindex_session=[];pkall=[];firecut=[];bimodality=[];peakfire=[];mnum=2;
for in=1:16
    for b=1:4
        for m=1:mnum
            for v=1:4
                for h=1:6
                    if nn_session(in,b,v,m)>=50
                        for i=1:2
                            a=squeeze(firemean_session(in,b,v,m,h,:,i));a=gaussian_smooth(a);
                            fireindex_session(in,b,v,m,h,:,i)=a;
                        end
                        ai=(a-min(a))/(max(a)-min(a));
                        fireindex_session(in,b,v,m,h,:,3)=ai;
                        a=squeeze(fireindex_session(in,b,v,m,h,:,1));
                        if max(a)>0
                            a=squeeze(fireindex_session(in,b,v,m,h,:,2));
                            [M,I]=max(a);
                            peakfire(in,b,v,m,h,:)=[M,rem(I*10-5,360)];
                        else
                            peakfire(in,b,v,m,h,:)=nan;
                        end
                        if h==3
                        a=squeeze(fireindex_session(in,b,v,m,h,:,:));
                        [M,I]=min(a(:,1));
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
                        bimodality(in,b,v,m,1,7:8)=nanmean(a(46:59,1:2),1);
                        bimodality(in,b,v,m,2,7:8)=nanmean(a(33:39,1:2),1);
                        bimodality(in,b,v,m,:,6)=bimodality(in,b,v,m,:,7)./bimodality(in,b,v,m,1,7);
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
bimodality(:,:,:,:,:,9)=rem(bimodality(:,:,:,:,:,1)*10-5,360);

spikepropname={'phase onset','max phase','sig seq percent','seq length (cm)','dir, out-in','theta freq','theta power',...
    'major','minor','minor,major','ca3 p'};
spikeprop=cat(5,seqcut(1:16,:,:,:,[7,5]),seqlength_session(:,:,:,:,2:4),thetaprop_session(1:16,:,:,1:2,:),squeeze(bimodality(1:16,:,:,:,:,8)),...
    squeeze(bimodality(1:16,:,:,:,2,6)),peakfire(1:16,:,:,:,4,1));
spikeprop2name={'major phase','minor phase','ca3 p peak phase','ca1 i peak phase','ca3 i peak phase','ca1 i','ca3 i'};
spikeprop2=cat(5,squeeze(bimodality(1:16,:,:,:,:,9)),squeeze(peakfire(1:16,:,:,:,4:6,2)),squeeze(peakfire(1:16,:,:,:,5:6,1)));
cd(datadir);save(['Figure3_thetaseq_figure_v4'],'nn_session','firemean_session','lfpmean_session','seqmean_session',...
        "seqlength_session",'thetaprop_session','deposphase','deposphase_dist','deposphase_dist_n','seqcut','fireindex_session',...
        'pkall','firecut','bimodality','peakfire','spikeprop2','spikeprop2name','spikeprop','spikepropname',...
        'velmean_session2','velmean_session','-v7.3');
end
%------------------------------------------------------------------------------------
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
cd(datadir);load('Figure3_thetaseq_figure_v4');
timewin=0.01;reward_dis_cut=70;min_dewinN=1;ifpeak=0;phasebin=10;
spikepropname={'phase onset','max phase','sig seq percent','seq length (cm)','dir, out-in','theta freq','theta power',...
    'major','minor','minor/major','ca3 p'};
bname={'reward','center','in','out'};mname={'local','remote'};
rowname={'CA1 pyr LFP','DG oml LFP','CA1 principal spike','CA3/DG principal spike'};
for i=1:101
    rowname{4+i}=['Dist: ',num2str(2*(i-51))];
end
% below generate paper fig 3
fign=4;pn1=6;pn2=8;ifdist=1;thetacycleN_cut=50;mnum=2;lettername={'c','d','e'};m=1;d1=5;seqcut2=[];
seqdensity=[];gaustd=2;
if 1
% below mean seq part
for b=1:2
    as=[];bs=[];lfp=[];cens=[];
    for v=1:4
        if ifdist
            if b==1
                am=[];
                for in=1:16
                    a=[squeeze(deposphase_dist_n(in,b,v,1,66:-1:d1,:));squeeze(deposphase_dist_n(in,b,v,2,d1:66,:))];
                    a=a./(ones(size(a,1),1)*sum(a,1));
                    am(in,:,:)=a;%imgaussfilt(a,2);
                    a=imgaussfilt(a,gaustd);
                    [M,~]=max(a,[],1);
                    M2=M(1:36)+M(37:72);
                    M2=gaussian_smooth(M2);
                    [M,I]=min(M2);
                    seqcut2(in,b,v,1:2)=[I*10-5,M];
                    [M,I]=max(M2);
                    seqcut2(in,b,v,3:4)=[I*10-5,M];
                end
            else
                am=[];
                for in=1:16
                    a=[squeeze(nansum(deposphase_dist_n(in,b,v,:,:,:),4))];
                    a=a./(ones(size(a,1),1)*sum(a,1));
                    am(in,:,:)=a;%imgaussfilt(a,gaustd);
                    a=imgaussfilt(a,3);
                    [M,~]=max(a,[],1);
                    M2=M(1:36)+M(37:72);
                    M2=gaussian_smooth(M2);
                    [M,I]=min(M2);
                    seqcut2(in,b,v,1:2)=[I*10-5,M];
                    [M,I]=max(M2);
                    seqcut2(in,b,v,3:4)=[I*10-5,M];
                end
            end
            as=[as,imgaussfilt(squeeze(nanmean(am,1)),gaustd)];
            b1=squeeze(firemean_session(:,b,v,3,:,:,1));
            for k=1:6
                for in=1:16 
                    b4=gaussian_smooth(squeeze(b1(in,k,:)));
                    b1(in,k,:)=(b4-min(b4))/(max(b4)-min(b4));
                end
            end
            bs=[bs,squeeze(nanmean(b1,1))];
            lfp1=squeeze(nanmean(lfpmean_session(:,b,v,3,:,:),1));
            lfp=[lfp,lfp1]; 
            ins=squeeze(nn_session(:,b,v,3))>=thetacycleN_cut;
            cen=rem(circ_meand(squeeze(seqcut2(ins,b,v,1))),360)/10;
            cens=[cens,[cen,cen+36]+(v-1)*72];
        end
    end
    figure(fign);subplot(pn1,pn2,[1:floor(pn2/2)]+(b-1)*floor(pn2/2));
    hold on;plot(lfp(1,:)+2,'k');hold on;plot(lfp(4,:),'k');xline([1:4]*72,'k');
    xticks([0:36:4*72]);xlim([0 4*72]);ylabel('zscore lfp');xticks([]);hold on;xline(cens,'c--');%plot(cens,4,'k+');%

    figure(fign);subplot(pn1,pn2,[1:floor(pn2/2)]+(b-1)*floor(pn2/2)+pn2);
    hold on;plot(bs(3,:)+1,'k');%'color',[0.8500 0.3250 0.0980]);
    hold on;plot(bs(4,:),'k');ylim([-0.5 2.5]);
    %hold on;plot(bs(3,:)-bs(4,:));
    xticks([0:36:4*72]);xline([1:4]*72,'k');xlim([0 4*72]);xline(cens,'c--');%xline([36:36:36*4],'k');xline([52:72:72*4],'m');
    ylabel('zscore spike N');yticks([0:3]*5);yticklabels({'CA1 P','CA3 P','CA1 I','CA3 I'});
    if b==1
        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2),1+pn2:floor(pn2/2)+pn2,1+2*pn2:floor(pn2/2)+2*pn2,1+3*pn2:floor(pn2/2)+3*pn2]+(b-1)*floor(pn2/2)+pn2*2);
        imagesc([],[-30 234-4d1],as);set(gca,'YDir','normal');
        xticks([0:36:4*72]);xline([1:4]*72,'k');colormap('hot');colorbar;xline(cens,'c--');
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
    else
        figure(fign);subplot(pn1,pn2,[1:floor(pn2/2),1+pn2:floor(pn2/2)+pn2]+(b-1)*floor(pn2/2)+pn2*2);
        imagesc([],[-100 100],as);set(gca,'YDir','normal');ylim([-40 90]);
        xticks([0:36:4*72]);xline([1:4]*72,'k');colormap('hot');colorbar;xline(cens,'c--');
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
    end
    if 0
        figdata=[lfp([1,4],:);bs(3:4,:);as];
        T = array2table(figdata,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 3',lettername{b},'_',num2str(m)],'WriteRowNames',true);
    end
end
for b=3:4
    for m=1
        as=[];bs=[];lfp=[];cens=[];
        for v=4
           ind=1:16;%17;%17;
           if ifdist
               a=squeeze(nanmean(deposphase_dist(ind,b,v,m,:,:),1));a=imgaussfilt(a,gaustd);%a=imgaussfilt(squeeze(nanmean(deposphase_dist(ind,b,v,m,:,:),1)),1);
           else
               a=squeeze(nanmean(seqmean_session(ind,b,v,m,:,:)./nansum(seqmean_session(ind,b,v,m,:,:),5),1));
           end
           as=[as,a];
           b1=squeeze(firemean_session(ind,b,v,m,:,:,1));
           for k=1:6
               for in=1:16
                    b1(in,k,:)=gaussian_smooth(squeeze(b1(in,k,:)));
                    b1(in,k,:)=(b1(in,k,:)-min(b1(in,k,:)))/(max(b1(in,k,:))-min(b1(in,k,:)));
               end
           end
           b1=squeeze(nanmean(b1,1));
           bs=[bs,b1];
           
           lfp1=squeeze(nanmean(lfpmean_session(ind,b,v,m,:,:),1));
           lfp=[lfp,lfp1]; 
           %cen=rem(circ_meand([circ_meand(squeeze(seqcut(1:16,b,v,1,5))),circ_meand(squeeze(seqcut(1:16,b,v,2,5)))]),360)/10;
           ins=squeeze(nn_session(:,b,v,1))>=thetacycleN_cut;
           cen=rem(circ_meand(squeeze(seqcut(ins,b,v,1,7))),360)/10;
           cens=[cens,[cen,cen+36]];
        end
        figure(fign);subplot(pn1,pn2,5+(b-3)*2+pn2*4);
        hold on;plot(lfp(1,:)+2,'k');hold on;plot(lfp(4,:),'k');xticks([0:36:72]);
        ylabel('zscore lfp');xticks([]);hold on;xline(cens,'c--');%plot(cens,4,'k+');%
        %FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},'_1'];save_fig(FIG_INDEX,ifsavefig,[],5);

        figure(fign);subplot(pn1,pn2,5+(b-3)*2+pn2*5);
        hold on;plot(bs(3,:)+1,'k');%hold on;plot(bs(3,:),'color',[0.8500 0.3250 0.0980]);
        hold on;plot(bs(4,:),'k');ylim([-0.5 2.5]);
        %hold on;plot(bs(3,:)-bs(4,:));
        xticks([0:36:72]);xlim([0 1*72]);xline(cens,'c--');%xline([36:36:36*4],'k');xline([52:72:72*4],'m');
        ylabel('zscore spike N');yticks([0:3]*5);yticklabels({'CA1 P','CA3 P','CA1 I','CA3 I'});
        %FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},'_2'];save_fig(FIG_INDEX,ifsavefig,[],5);

        figure(fign);subplot(pn1,pn2,[6,14]+(b-3)*2+pn2*4);
        imagesc([],[-100 100],as);set(gca,'YDir','normal');
        if b==1
            ylim([-30 30]);
        elseif b==2
            ylim([-40 80]);
        elseif b==3
            ylim([-50 50]);
        elseif b==4
            ylim([-50 50]);
        end
        xline([1:4]*72,'k');colormap('hot');colorbar;xline(cens,'c--');xticks([0:36:72]);   
        title([bname{b},mname{m}]);ylabel('Distance to rat (cm)');
        %FIG_INDEX=['fig3_thetaseq_',bname{b},'_',mname{m},'_3'];save_fig(FIG_INDEX,ifsavefig,'hot');
        if 0
        figdata=[lfp([1,4],:);bs(3:4,:);as];
        T = array2table(figdata,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 3e_',num2str(b-2)],'WriteRowNames',true);
        end
    end
end
figure(fign);figure_title='Figure3_thetaseq_Jan2025';save_current_figure(figure_title);
end

% below for phase onset and core phase and combine local vs remote
spikeprop2=[];nn_session2=[];mu=[];stds=[];
fign=21;pn1=1;pn2=2;xdif=0.2;n1=0.08;c1=0.4;d=1;
for b=1:4
    for v=1:4
        for k=1:2
            for in=1:16
                if b<3
                    spikeprop2(in,b,v,k)=squeeze(seqcut2(in,b,v,k*2-1));
                    nn_session2(in,b,v,k)=nn_session(in,b,v,3);
                else
                    spikeprop2(in,b,v,k)=squeeze(spikeprop(in,b,v,1,k));
                    nn_session2(in,b,v,k)=nn_session(in,b,v,1);
                end
            end
            ind1=find(~isnan(spikeprop2(1:16,b,v,k)) & nn_session2(1:16,b,v,k)>=thetacycleN_cut);
            n=length(ind1);
            if n>0
                mu(b,v,k)=rem(circ_meand(squeeze(spikeprop2(ind1,b,v,k))),360);
                [~,~,err]=circ_meand(squeeze(spikeprop2(ind1,b,v,k)));
                stds(b,v,k)=err;%*n^(-0.5);
                figure(fign);subplot(pn1,pn2,k);
                if b<3
                    figure(fign);subplot(pn1,pn2,k);
                    hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)-xdif,squeeze(spikeprop2(ind1,b,v,k)),'.','Color',[1 1 1]*c1);
                    if k==2
                        hold on;plot((d*v+d*(b-1)*4.5)*ones(n,1)+n1*randn(n,1)-xdif,squeeze(spikeprop2(ind1,b,v,k))+360,'.','Color',[1 1 1]*c1);
                    end
                elseif b>=3 & v==4
                    figure(fign);subplot(pn1,pn2,k);
                    hold on;plot((d*10+d*(b-3)*1.5)*ones(n,1)+n1*randn(n,1),squeeze(spikeprop2(ind1,b,v,k)),'.','Color',[1 1 1]*c1);
                    if k==2
                        hold on;plot((d*10+d*(b-3)*1.5)*ones(n,1)+n1*randn(n,1),squeeze(spikeprop2(ind1,b,v,k))+360,'.','Color',[1 1 1]*c1);
                    end
                end
            else
                mu(b,v,k)=nan;stds(b,v,k)=nan;
            end
        end
    end
end
% below for statistics
circ_linear_correlation=[];
for k=1:2
    for b=1:2
        alpha=[];idp=[];idq=[];
        for v=1:4
            ind1=find((~isnan(spikeprop2(1:16,b,v,k))) & nn_session2(1:16,b,v,k)>=thetacycleN_cut);
            n=length(ind1);
            if n>0
                alpha=[alpha;squeeze(spikeprop2(ind1,b,v,k))];
                idp=[idp;ones(n,1)*v];
            end
        end
        [rho pval] = circ_corrcl(alpha*pi/180, idp);
        circ_linear_correlation(k,b,:)=[rho pval];
    end
end
k2phase=240;
for k=1:2
    figure(fign);subplot(pn1,pn2,k);
    for b=1:2
        if k==2 
            a=squeeze(mu(b,:,k))';a(a<k2phase)=a(a<k2phase)+360;hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'Color',[1 1 1]*c1*0);
            a=squeeze(mu(2+b,4,k));a(a<k2phase)=a(a<k2phase)+360;hold on;plot(d*[-0.15 0.15]+d*10+d*(b-1)*1.5,a*[1 1],'Color',[1 1 1]*c1*0);
        else
            a=squeeze(mu(b,:,k))';hold on;plot(d*[1:4]+d*(b-1)*4.5,a,'Color',[1 1 1]*c1*0);
            a=squeeze(mu(2+b,4,k))';hold on;plot(d*[-0.15 0.15]+d*10+d*(b-1)*1.5,a*[1 1],'Color',[1 1 1]*c1*0);
        end
    end
    if k==1
        ylim([0 360]);
    elseif k==2
        ylim([0 360]+k2phase);%ylim([180 540]);
    end
    xticks(d*[1:4,5.5:8.5,10,11.5]);xticklabels({'reward: 0-1','1-5','5-10','>10','center: 0-1','1-5','5-10','>10','in: >10','out: >10'});
    ylabel(spikepropname{k});xlim(d*[0 12.5]);title(num2str(squeeze(circ_linear_correlation(k,:,:)),2));
end
figure(fign);figure_title='Figure3_thetaseq_prop_Jan2025';save_current_figure(figure_title);

%------------------------------------------------------------------------------
% below Sup figure 9
% below statistic test
if 0
   % cd(datadir);load('Figure3_thetaseq_figure');
    % sequence length at reward
    b=1;v=1:4;m=1:2;k=4;reps=16;
    anovadata=[];
    for v=1:4
        a=squeeze(spikeprop(1:16,b,v,m,k));
        anovadata=[anovadata;a];
    end
    [p,table] = anova2(anovadata,reps);

    % sig sequence fraction
   
    b=1:2;v=1:4;m=1:2;reps=16;ks=[3,4,5,7,6];anovap=[];
    for kk=1%:5
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
            p = anovan(anovadata,{mlabels,vlabels});
            ps(:,b)=p;
            anovap(b,:,k)=p;
        end
    end

    pss=[];
for k=1:2
    for b=1:2%:4
        alpha=[];idp=[];idq=[];
        for m=1:2%(3-b)
            for v=1:4
                ind1=find((~isnan(spikeprop(1:16,b,v,m,k))) & nn_session(1:16,b,v,m)>thetacycleN_cut);
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
            circ_linear_correlation(k,b,m,:)=[rho pval];
        end
    end
end
end
thetacycleN_cut=50;fign=20;pn1=2;pn2=4;mu=[];xdif=0.2;n1=0.08;c1=0.4;d=1;stds=[];
lettername={'g','f','a','b','c','h','i','d','e','f','g','j','k'};inter=1;fn={'velocity','local/remote'};
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
    %FIG_INDEX=['fig3_thetaseq_prop_',spikepropname{k}];save_fig(FIG_INDEX,ifsavefig,[],3);
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
for k=1:11
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
%-------------------------------------------------------------------------------