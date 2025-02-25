function figure_3_analyze_theta_sequence_at_pause(in,ifdataall)
% the goal here is to study theta cyles during pause; 
%--------------------------------------------------------------------------------
%close all;set(0,'DefaultFigureColormap',feval('turbo'));
if nargin<2
    ifdataall=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
savedir=SessionSet16{in};min_dewinN=14;
cd(savedir);
load('Analyze_Theta_Cycle_Broad','thetacycle','reftheta','lfptime_trial');
load('Position_Data_Maze.mat');
load('Arm_Visit_Confirmed.mat', 'includetrials');
if ifdataall==1
    load('sorted_spike_decoding_8arm6.mat','deseqall','dedataall','maze1d');
else
    load('sorted_spike_decoding_8arm_dHPC.mat','deseqall','dedataall','maze1d');
end
load('InOutBound_Behavior_Analysis.mat','pksets','inter_run','pos_behave_marker');
load('Ripple_Events_PAPER5.mat','ripples');
ind=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>0.035 ;
ripples=ripples(ind,:);

%  below combine pksets with conjunctions
for i=1:size(pksets,1)-1
    if pksets(i,1)==pksets(i+1,1) & pksets(i,9)>30 & pksets(i,2)==pksets(i+1,2)
        if pksets(i,7)>=0
            if pksets(i,7)==round(pksets(i,7))
                pksets(i,7)=pksets(i,7)+0.1;
                pksets(i+1,7)=pksets(i,7)+0.1;
            else
                pksets(i+1,7)=pksets(i,7)+0.1;
            end
        else
            if pksets(i,7)==round(pksets(i,7))
                pksets(i,7)=pksets(i,7)-0.1;
                pksets(i+1,7)=pksets(i,7)-0.1;
            else
                pksets(i+1,7)=pksets(i,7)-0.1;
            end
        end
    end
end
trialset=unique(pksets(:,1));trialset=trialset(:);
for trial=1:size(trialset,1)
    trialid=trialset(trial,1);
    id=find(includetrials(:,1)==trialid);
    if ~isempty(id)
        trialset(trial,2:3)=includetrials(id,[3,5]);
    end
    pk=pksets(pksets(:,1)==trialid,:);
    trialset(trial,4:5)=[round(max(pk(:,7))),sum(pk(:,7)<=-4)];
end
trialset(:,6)=0;
ind=trialset(:,4)==8 & trialset(:,5)==0;
trialset(ind,6)=1;

% 1. prepare data
ruler=lfptime_trial(:,1);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
deseqall(:,27)=reftheta(outputindex,4);figure;histogram(error);
deseqall(:,42)=outputindex;
ruler=Position_Data(:,1);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
deseqall(:,28:30)=[pos_behave_marker(outputindex,:),outputindex];
deseqall0=deseqall;
trialfull=unique(deseqall(:,14));
if length(trialfull)~=length(dedataall)
    disp(['Error dedata size : ',num2str(in)]);
end
linear_maze=[];linear_pos=[];deseqall=[];
for i=1:length(dedataall)
    dedata=dedataall{i};
    deseq=deseqall0(deseqall0(:,14)==trialfull(i),:);
    ind= deseq(:,28)~=-100 & deseq(:,29)~=-100 ;
    if sum(ind)>0
        deseq=deseq(ind,:);
        dedata=dedata(:,ind);
        [linear_maze1,linear_pos1]=linearize_8arm(dedata',maze1d,deseq(:,[1,4,5]));
        linear_maze=cat(1,linear_maze,linear_maze1);
        linear_pos=[linear_pos;linear_pos1];
        deseqall=[deseqall;deseq];
    end
end
clear dedataall
linear_pos(:,[2,4])=linear_pos(:,[2,4])*2; % in cm
linear_pos(:,7)=1; % local encoding
indremote=linear_pos(:,5)~=linear_pos(:,6);
linear_pos(indremote,7)=2; % remote encoding
deseqall(:,31)=linear_pos(:,2);% rat distance to center
deseqall(:,32:36)=0;
for k=1:2
    indk=ripples(:,6)==k;
    ruler=ripples(indk,3);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
    deseqall(:,31+k)=error;
end
if 0
for k=1:3
    indk=dspeaks_d_trial(:,6)==k;
    ruler=dspeaks_d_trial(indk,3);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
    deseqall(:,33+k)=error;
end
end
deseqall(:,37)=pksets(deseqall(:,29),7);% current visit ID
ind=abs(deseqall(:,28))==0.5;
deseqall(ind,37)=inter_run(deseqall(ind,29),4);

% below merge pksets
pksets_m=pksets;
pksets_m(:,42:43)=0;
indmerge=find(pksets(:,7)-round(pksets(:,7))~=0); % pk id
indmerge(:,2)=pksets(indmerge(:,1),1); % trial id
indmerge(:,3)=pksets(indmerge(:,1),7); % visit id
indmerge(:,4)=0;
for i=1:size(indmerge,1)
    if indmerge(i,4)==0
        ind=(indmerge(:,2)==indmerge(i,2) & round(indmerge(:,3))==round(indmerge(i,3)));
        indmerge(ind,4)=1;
        pkids=indmerge(ind,1);
        pksets_m(pkids(2:end),:)=0;
        pksets_m(pkids(1),42)=length(pkids);
    end
end
indinclu=pksets_m(:,1)>0 & (~( pksets_m(:,7)>-4 & pksets_m(:,7)<=0 ));
pksets_m=pksets_m(indinclu,:);
pksets_m(1:end-1,43)=pksets_m(2:end,7); % next visit id
label=pksets_m(1:end-1,1)-pksets_m(2:end,1);label(end+1)=1;
ind=label~=0;
pksets_m(ind,43)=0; % no next arm
pksets_m(:,7)=round(pksets_m(:,7));
deseqall(:,38:39)=0;% decoded visit id | trial correct/error
deseqposindex=deseqall(:,30);
for trial=1:size(trialset,1)
    trialid=trialset(trial,1);
    indtrial=Position_Data(deseqposindex,4)==trialid;
    pk=round(pksets_m(pksets_m(:,1)==trialid,[2,7]));
    deseqall(indtrial,39)=trialset(trial,6);
    if trialset(trial,6)==1 % correct trials
        [~,I]=sort(pk(:,1));pk=pk(I,:);
        pkind=[1;diff(pk(:,1))]==1;pk=pk(pkind,:);
        if size(pk,1)~=8
            disp(['Error in pk size!']);
        else
            deseqall(indtrial,38)=pk(linear_pos(indtrial,3),2);% decoded visit id
        end
    else % error trials
        for i=1:8
            ind=linear_pos(:,3)==i & indtrial;
            visit=unique(pk(pk(:,1)==i,2));
            if isempty(visit) % skipped visit
                deseqall(ind,38)=nan;
            elseif length(visit)==1 % correct visit
                deseqall(ind,38)=visit;
            elseif length(visit)>1 % error visit
                deseqall(ind,38)=0;
            end
        end
    end 
end
ind=linear_pos(:,6)==0;
deseqall(ind,38)=0;

despikez=zscore(deseqall(:,[6:7,19:26]),0,1);
lfpdt=1/625;
% below analyze theta seq associated with thetacyles
ind=thetacycle(:,4)>0 & thetacycle(:,3)<=0.2;
thetacycle=thetacycle(ind,:);
thetacycle(:,21:57)=0;
thetacycle(:,21)=[1:size(thetacycle,1)]; % thetacycle id
thetacycle(:,22)=reftheta(thetacycle(:,1),1)-0.5*lfpdt; % start time
thetacycle(:,23)=reftheta(thetacycle(:,2),1)+0.5*lfpdt; % end time
ruler=deseqall(:,1);template=thetacycle(:,6);[outputindex,error]=match(template,ruler,0);
thetacycle(:,26:28)=[linear_pos(outputindex,[5,2]),deseqall(outputindex,37)]; % current arm | distance to center | current visit id
thetacycle_depos_dist=nan*ones(size(thetacycle,1),3,9);thetacycle_depos_mean=nan*ones(size(thetacycle,1),17);
deseqall(:,40:41)=0;
for i=1:size(thetacycle,1)
    ind= deseqall(:,42)>=thetacycle(i,1) & deseqall(:,42)<=thetacycle(i,2) ; 
    %ind= deseqall(:,1)>=thetacycle(i,22) & deseqall(:,1)<thetacycle(i,23) ; % [ )
    thetacycle(i,24)=sum(ind); % N decode win
    if sum(ind)>0
        deseqall(ind,40)=i; % belong to thetacycle id
        deseqall(ind,41)=deseqall(ind,41)+1; % belong to num of thetacycles
        thetacycle_depos_dist(i,1,:)=histcounts(linear_pos(ind,6),[-0.5:8.5]); % decoded arm
        thetacycle_depos_dist(i,2,:)=histcounts(deseqall(ind,38),[-0.5:8.5]); % decoded visit id
        thetacycle_depos_dist(i,3,:)=histcounts(linear_pos(ind,4),[0:10:80,100]); % decoded distance to center
        thetacycle_depos_mean(i,17)=squeeze(nanmean(nanmean(nanmean(linear_maze(ind,1:15,:),1),2),3));
        thetacycle_depos_mean(i,1:8)=squeeze(nanmean(nanmean(linear_maze(ind,16:30,:),1),2));
        thetacycle_depos_mean(i,9:16)=squeeze(nanmean(nanmean(linear_maze(ind,31:50,:),1),2));
        thetacycle(i,25)=mean(linear_pos(ind,7)==2); % remote encoding percent
        thetacycle(i,29:38)=mean(deseqall(ind,[6:7,19:26]),1);
        thetacycle(i,43:52)=mean(despikez(ind,:),1);
    end
    thetacycle(:,55:56)=0;
    thetacycle(:,55)=Position_Data(thetacycle(:,9),4);
    lia=ismember(thetacycle(:,55),trialset(trialset(:,6)==1,1));
    thetacycle(lia,56)=1;

    thetacycle(:,39)=0; % group
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)<10 & thetacycle(:,27)>60;
    thetacycle(ind,39)=1; % reward pause
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)>=10 & thetacycle(:,27)>60;
    thetacycle(ind,39)=2; % reward run
    
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)<10 & thetacycle(:,27)<30;
    thetacycle(ind,39)=3; % inter pause
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)>=10 & thetacycle(:,27)<30;
    thetacycle(ind,39)=4; % inter run
    
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)<10 & thetacycle(:,27)>=30 & thetacycle(:,27)<=60 & thetacycle(:,11)<0;
    thetacycle(ind,39)=5; % in pause
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)>=10 & thetacycle(:,27)>=30 & thetacycle(:,27)<=60 & thetacycle(:,11)<0;
    thetacycle(ind,39)=6; % in run

    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)<10 & thetacycle(:,27)>=30 & thetacycle(:,27)<=60 & thetacycle(:,11)>0;
    thetacycle(ind,39)=7; % out pause
    ind=thetacycle(:,24)>=min_dewinN & thetacycle(:,10)>=10 & thetacycle(:,27)>=30 & thetacycle(:,27)<=60 & thetacycle(:,11)>0;
    thetacycle(ind,39)=8; % out run
    
    thetacycle(:,40)=sum(thetacycle_depos_dist(:,1,:)>0,3); % encoded arm N
    a=squeeze(thetacycle_depos_dist(:,1,:));a=a./(sum(a,2)*ones(1,9));
    thetacycle(:,41)=sum(abs(a-1/9),2)*9/16;% concentration index [0 1]
    thetacycle(:,42)=in;
    thetacycle(:,57)=sum(squeeze(thetacycle_depos_dist(:,1,2:9))>0,2);
    thetacycle(:,58)=max(squeeze(thetacycle_depos_dist(:,1,:)),[],2)./thetacycle(:,24); % dominance index
end
if ifdataall==1
    cd(savedir);save('ThetaCycle_Decode_Info_2','deseqall','linear_maze','linear_pos','pksets','pksets_m','inter_run','trialset',...
    'thetacycle','thetacycle_depos_dist','maze1d','ripple1_trial','dspeaks_d_trial','ripple2','thetacycle_depos_mean','-v7.3');
else
    cd(savedir);save('ThetaCycle_Decode_Info_dHPC2','deseqall','linear_maze','linear_pos','pksets','pksets_m','inter_run','trialset',...
    'thetacycle','thetacycle_depos_dist','maze1d','ripple1_trial','dspeaks_d_trial','ripple2','thetacycle_depos_mean','-v7.3');
end
disp(['Finish Session : ',num2str(in)]);
if 0
close all;set(0,'DefaultFigureColormap',feval('turbo'));%clear all;in=1
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
unittype10 = {'All','All Cell','dCA1 Pyramidal','dCA3DG Pyramidal','vDGCA3 Pyramidal','v subiculum Pyramidal','dCA1 Interneuron','dCA3DG Interneuron','vDGCA3 Interneuron','v subiculum Interneuron'};
velname={'Pause','Run'};
typenum=16;typeindex=[typenum:-1:1]*floor(256/typenum);
plotcolors=colormap(jet);
plotcolors=plotcolors(typeindex,:);
figure;phasebin=10;thetacycle_spikeall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_2','deseqall','thetacycle','thetacycle_depos_dist','trialset');
    deseqall(:,27)=ceil(deseqall(:,27)/phasebin); 
    despikez=zscore(deseqall(:,[6:7,19:26]),0,1);
    thetacycle_spike=[];
    ind=thetacycle(:,10)<10;
    lia=ismember(deseqall(:,40),thetacycle(ind,21));
    ind=thetacycle(:,10)>=10;
    lia2=ismember(deseqall(:,40),thetacycle(ind,21));
    for phase=1:360/phasebin
        indph=deseqall(:,27)==phase;
        thetacycle_spike(phase,:,1)=squeeze(nanmean(despikez( lia & indph, :),1));
        thetacycle_spike(phase,:,2)=squeeze(nanmean(despikez( lia2 & indph, :),1));
    end
    thetacycle_spikeall(:,:,:,in)=thetacycle_spike;
    for p=1:10
        for i=1:2
            subplot(4,5,p+(i-1)*10);
            a=squeeze(thetacycle_spike(:,p,i))';
            hold on;plot([5:10:715],[a,a],'o-','color',plotcolors(in,:));
            xlabel(unittype10{p});title(velname{i});
        end
    end
end
figure;
for p=1:10
    subplot(2,5,p);
    as=[];
    for i=1:2
        a=squeeze(thetacycle_spikeall(:,p,i,:))';
        as=[as,a,a];
    end
    imagesc([0 720]*2,[1 16],as);set(gca,'YDir','normal'); 
    xlabel('Pause | Run');title(unittype10{p});xline(720,'k');
end
cd(homedir);figure_title='thetaseq_spike';save_current_figure(figure_title);
% below look at spike dist as a function of theta power; velocity; cycle
% length

%---------------------------------------------------------------------------------
% below focus on theta sequence presentation
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
thetacycleall=[];depos_distall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_2.mat','deseqall','thetacycle','thetacycle_depos_dist','trialset');
    load('Position_Data_Maze.mat');
    [velocity,vel_angle,locdis,loc_angle,locdis_change]=calculate_location_movement(Position_Data);
    thetacycle(:,41)=locdis_change(thetacycle(:,9)); % movement direction
    thetacycleall=[thetacycleall;thetacycle];
    depos_distall=cat(1,depos_distall,thetacycle_depos_dist);
end
groupname={'Reward Pause','Reward Run','Inter Pause','Inter Run','In Pause','In Run','Out Pause','Out run'};

% below plot basic properties
figure;
for g=1:4
    for v=1:2
        group=v+(g-1)*2;
        ind=thetacycleall(:,39)==group;
        subplot(4,5,1+(g-1)*5);hold on;histogram(thetacycleall(ind,24));xlabel('De Win N');
        subplot(4,5,2+(g-1)*5);hold on;histogram(thetacycleall(ind,25));xlabel('Remote Percent');    
        subplot(4,5,3+(g-1)*5);hold on;histogram(thetacycleall(ind,40));xlabel('De Region N');  
        subplot(4,5,4+(g-1)*5);hold on;histogram(thetacycleall(ind,57));xlabel('De Arm N'); 
        subplot(4,5,5+(g-1)*5);hold on;histogram(thetacycleall(ind,58));xlabel('Dominance Index');
    end
    legend(groupname{1+(g-1)*2},groupname{2+(g-1)*2});
end

[M,I]=max(squeeze(depos_distall(:,1,:)),[],2);
thetacycleall(:,59)=(squeeze(depos_distall(:,1,1))+M)./thetacycleall(:,24); % (center+max arm) decode percent
ind=I==1;
thetacycleall(ind,59)=thetacycleall(ind,58);
a=squeeze(depos_distall(:,1,:))./(thetacycleall(:,24)*ones(1,9));
thetacycleall(:,60)=sum(a(:,2:9)>=0.2,2); % N of arm with >=0.1 decode percent
[M,I]=max(a,[],2);
thetacycleall(:,61:62)=[I-1,M];
for i=1:9
    ind=I==i;
    a(ind,i)=-a(ind,i);
end
[M2,I2]=max(a,[],2);
thetacycleall(:,63:64)=[I2-1,M2];
for i=1:9
    ind=I2==i;
    a(ind,i)=-a(ind,i);
end
[M3,I3]=max(a,[],2);
thetacycleall(:,65:66)=[I3-1,M3];
for i=1:9
    ind=I3==i;
    a(ind,i)=-a(ind,i);
end
[M4,I4]=max(a,[],2);
thetacycleall(:,67:68)=[I4-1,M4];

a=squeeze(depos_distall(:,1,2:9))./(thetacycleall(:,24)*ones(1,8));
[M,I]=max(a,[],2);
thetacycleall(:,69:70)=[I,M];
for i=1:8
    ind=I==i;
    a(ind,i)=-a(ind,i);
end
[M2,I2]=max(a,[],2);
thetacycleall(:,71:72)=[I2,M2];
for i=1:8
    ind=I2==i;
    a(ind,i)=-a(ind,i);
end
[M3,I3]=max(a,[],2);
thetacycleall(:,73:74)=[I3,M3];

figure;
for k=1:4
    subplot(1,2,1);hold on;histogram(thetacycleall(:,61+(k-1)*2));
    subplot(1,2,2);hold on;histogram(thetacycleall(:,62+(k-1)*2));
end
figure;kname={'Max decode percent','Arm or Center with max decode percent','N of arm with >=0.1 decode percent'};
for k=1:3
    subplot(1,3,k);histogram(thetacycleall(:,58+k));xlabel(kname{k});ylabel('Thetacycle N');
end

ind0=thetacycleall(:,60)==0;depos_dist0=depos_distall(ind0,:,:);thetacycle0=thetacycleall(ind0,:);
ind1=thetacycleall(:,60)==1;depos_dist1=depos_distall(ind1,:,:);thetacycle1=thetacycleall(ind1,:);
ind2=thetacycleall(:,60)==2;depos_dist2=depos_distall(ind2,:,:);thetacycle2=thetacycleall(ind2,:);

% below separately study thetaseq when col 60 = 0,1,2
thetacycle0=[];thetacycle1=[];thetacycle2=[];
thetaseq0=[];thetaseq1=[];thetaseq2=[];phasebin=10;thetaspike0=[];thetaspike1=[];thetaspike2=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_dHPC','linear_maze','deseqall','thetacycle_depos_dist');
    deseqall(:,27)=ceil(deseqall(:,27)/phasebin);
    despikez=zscore(deseqall(:,[6:7,19:26]),0,1);
    indin=thetacycleall(:,42)==in;
    for k=0:2
        indk=thetacycleall(:,60)==k & indin;
        cycle=thetacycleall(indk,:);
        seq=[];spike=[];
        for phase=1:360/phasebin
            indph=deseqall(:,27)==phase;
            for i=1:size(cycle,1)
                lia=ismember(deseqall(:,40),cycle(i,21));
                if k==0
                    if cycle(i,26)>0
                        seq(i,phase,1:15)=squeeze(nanmean(nanmean(linear_maze( lia & indph, 1:15, :),1),3));
                        seq(i,phase,16:50)=squeeze(nanmean(linear_maze( lia & indph, 16:50, cycle(i,26)),1));
                    else
                        seq(i,phase,:)=squeeze(nanmean(nanmean(linear_maze( lia & indph, :, :),1),3));
                    end
                elseif k==1
                    seq(i,phase,1:15)=squeeze(nanmean(nanmean(linear_maze( lia & indph, 1:15, :),1),3));
                    seq(i,phase,16:50)=squeeze(nanmean(linear_maze( lia & indph, 16:50, cycle(i,69)),1));
                elseif k==2
                    seq(i,phase,1:15)=squeeze(nanmean(nanmean(linear_maze( lia & indph, 1:15, :),1),3));
                    seq(i,phase,16:50)=squeeze(nanmean(linear_maze( lia & indph, 16:50, cycle(i,69)),1));
                    seq(i,phase,51:85)=squeeze(nanmean(linear_maze( lia & indph, 16:50, cycle(i,71)),1));
                end
                spike(i,phase,:,1)=squeeze(nanmean(deseqall(lia & indph,[6:7,19:26]),1));
                spike(i,phase,:,2)=squeeze(nanmean(despikez(lia & indph,:),1));
            end
        end
        if k==0
            thetacycle0=[thetacycle0;cycle];
            thetaseq0=cat(1,thetaseq0,seq);
            thetaspike0=cat(1,thetaspike0,spike);
        elseif k==1
            thetacycle1=[thetacycle1;cycle];
            thetaseq1=cat(1,thetaseq1,seq);
            thetaspike1=cat(1,thetaspike1,spike);
        elseif k==2
            thetacycle2=[thetacycle2;cycle];
            thetaseq2=cat(1,thetaseq2,seq);
            thetaspike2=cat(1,thetaspike2,spike);
        end
    end
    in
end
cd(homedir);save('thetacycle_thetaseq_dHPC','thetacycleall','depos_distall','thetacycle0','thetacycle1','thetacycle2','thetaseq0',...
'thetaseq1','thetaseq2','depos_dist0','depos_dist1','depos_dist2','phasebin','thetaspike0','thetaspike1','thetaspike2','-v7.3');
thetaseq1(:,:,1:15)=thetaseq1(:,:,1:15)*8;
% below visualize individual thetaseq1
figure;
for i=1:size(thetaseq1,1)
    a=squeeze(thetaseq1(i,:,:))';
    subplot(1,2,1);imagesc(a);set(gca,'YDir','normal');title(num2str(thetacycle1(i,26)));
    subplot(1,2,2);a(1:15,:)=8*a(1:15,:);imagesc(a);set(gca,'YDir','normal');
    pause;clf;
end
%-----------------------------------------------------------------------------------
% below further look into local vs remote theta sequences
close all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
load('thetacycle_thetaseq2.mat')
figure;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_2','deseqall','thetacycle');
    a=diff(deseqall(:,40));
    a=[0;a];
    ind=deseqall(:,41)==1 & a>0;
    subplot(4,4,in);histogram(deseqall(ind,27));
    ind=deseqall(:,27)>=300 & deseqall(:,41)==1 & a>0;title(num2str(size(unique(deseqall(ind,40)))));
    if 0
    deseq=deseqall(ind,:);
    id=unique(deseqall(ind,40));
    the=thetacycle(id,:);
    deseq(:,42)=deseq(:,1)-the(:,22);
    subplot(4,4,in);histogram(deseq(:,42));title(num2str(mean(deseq(:,42)<0)));
    end
end

% 1. thetaseq in/out, run, local/remote
indv=thetacycle1(:,10)>=10;
indin=thetacycle1(:,11)==-1;
indout=thetacycle1(:,11)==1;
indlocal=thetacycle1(:,69)==thetacycle1(:,26);
indremote=thetacycle1(:,69)~=thetacycle1(:,26);
thetaseq_align1=[];thetaseq_noalign1=[];thetaspike_align1=[];
for dir=1:2
    if dir==1
        indd=indout;
    else
        indd=indin;
    end
    for k=1:2
        if k==1
            ind=indv & indd & indlocal;% 
        else
            ind=indv & indd & indremote;% 
        end
        thetacycle=thetacycle1(ind,:);
        thetaseq=thetaseq1(ind,:,:);
        thetaseq_align=align_thetaseq(thetacycle,thetaseq);
        thetaseq_align1{dir,k}=thetaseq_align;
        thetacycle_align1{dir,k}=thetacycle;
        thetaseq_noalign1{dir,k}=thetaseq;
        thetaspike_align1{dir,k}=thetaspike1(ind,:,:,:);
    end
end
figure;
for in=1:16
    ass=[];
    for k=1:2
        as=[];
        for dir=1:2
            a=thetaseq_align1{dir,k};
            thetacycle=thetacycle_align1{dir,k};
            ind=thetacycle(:,42)==in;
            a1=squeeze(nanmean(a(ind,:,:),1))';
            as=[as,a1,a1];
        end
        ass=[ass,as];
    end
    subplot(4,4,in);imagesc([0 8]*360,[-100 100],ass);set(gca,'YDir','normal');
    xline(360*2*[1:4],'k');ylim([-1 1]*40);
end
unittype10 = {'All','All Cell','dCA1 Pyramidal','dCA3DG Pyramidal','vDGCA3 Pyramidal','v subiculum Pyramidal','dCA1 Interneuron','dCA3DG Interneuron','vDGCA3 Interneuron','v subiculum Interneuron'};

% below plot theta seq & spike as a function of velocity
figure;
for in=1:16
    ass=[];
    for k=1
        as=[];b3=[];
        for dir=1:2
            a=thetaseq_align1{dir,k};
            b=thetaspike_align1{dir,k};
            thetacycle=thetacycle_align1{dir,k};
            ind=thetacycle(:,42)==in;
            c=thetacycle(ind,:);
            a=a(ind,:,:);
            b=b(ind,:,:,:);
            [~,I]=sort(c(:,10));
            I(:,2)=ceil(5*[1:length(I)]/length(I));
            b1=[];
            for i=1:5
                ind=I(:,2)==i;
                b1(i,:,:)=squeeze(nanmean(b(ind,:,3:10,2),1));
            end
            b2=[];
            for i=1:8
                b2=[b2,b1(:,:,i)'];
            end
            b3=[b3;b2];
        end
    end
    subplot(4,4,in);imagesc([1 40],[0 720],b3);set(gca,'YDir','normal');
    xline([5.5:5:40],'k');yline(360,'k');
    xticks([3:5:40]);xticklabels(unittype10(3:10));ylabel('Theta phase: out | in');
end
cd(homedir);figure_title='thetaseq_spike_local_run_outin';save_current_figure(figure_title,8);
figure;k=1;
for uk=1:8
    b3=[];
    for dir=1:2
        a=thetaseq_align1{dir,k};
        b=thetaspike_align1{dir,k};
        thetacycle=thetacycle_align1{dir,k};
        b1=[];
        for in=1:16
            ind=thetacycle(:,42)==in;
            c=thetacycle(ind,:);
            a_=a(ind,:,:);
            b_=b(ind,:,:,:);
            [~,I]=sort(c(:,10));
            I(:,2)=ceil(5*[1:length(I)]/length(I));
            for i=1:5
                ind=I(:,2)==i;
                b1(i,:,in)=squeeze(nanmean(b_(ind,:,2+uk,2),1));
            end
        end
        b2=[];
        for in=1:16
            b2=[b2,b1(:,:,in)'];
        end
        b3=[b3;b2];
    end
    subplot(2,4,uk);imagesc([1 80],[0 720],b3);set(gca,'YDir','normal');
    xline([5.5:5:80],'k');yline(360,'k');title(unittype10(2+uk));
    ylabel('Theta phase: out | in');xlabel('16 sessions: sorted by vel');
end
cd(homedir);figure_title='thetaseq_spike_local_run_outin2';save_current_figure(figure_title);

figure;
for in=1:16
    for k=1
        as=[];
        for dir=1:2
            a=thetaseq_align1{dir,k};
            b=thetaspike_align1{dir,k};
            thetacycle=thetacycle_align1{dir,k};
            ind=thetacycle(:,42)==in;
            c=thetacycle(ind,:);
            a=a(ind,:,:);
            b=b(ind,:,:,:);
            [~,I]=sort(c(:,10));
            I(:,2)=ceil(5*[1:length(I)]/length(I));
            a1=[];
            for i=1:5
                ind=I(:,2)==i;
                a1=[a1,squeeze(nanmean(a(ind,:,:),1))'];
            end
            as=[as;a1(31:71,:)];
        end
    end
    subplot(4,4,in);imagesc([ ],[-1 1]*40,as);set(gca,'YDir','normal');
    xline([36.5:36:36*5],'k');yline(0,'k');xlabel('Theta phase: sorted by vel');ylabel('Location: out | in');
end
cd(homedir);figure_title='thetaseq_local_run_outin_vel';save_current_figure(figure_title,8);

indv=thetacycle0(:,10)>=10;
indin=thetacycle0(:,11)==-1;
indout=thetacycle0(:,11)==1;
thetaseq_align0=[];
for dir=1:2
    if dir==1
        indd=indout;
    else
        indd=indin;
    end
    ind=indv & indd;
    thetacycle=thetacycle0(ind,:);
    thetaseq=thetaseq0(ind,:,:);
    thetaseq_align=align_thetaseq(thetacycle,thetaseq);
    thetaseq_align0{dir}=thetaseq_align;
    thetacycle_align0{dir}=thetacycle;
end
figure;
for in=1:16
    as=[];
    for dir=1:2
        a=thetaseq_align0{dir};
        thetacycle=thetacycle_align0{dir};
        ind=thetacycle(:,42)==in;
        a1=squeeze(nanmean(a(ind,:,:),1))';
        as=[as,a1,a1];
    end
    subplot(4,4,in);imagesc([0 4]*360,[-100 100],as);set(gca,'YDir','normal');
    xline(360*2,'k');ylim([-1 1]*40);
end
figure;
for i=1:size(thetacycle2,1)
    if (thetacycle2(i,69)==thetacycle2(i,26) | thetacycle2(i,71)==thetacycle2(i,26)) & abs(thetacycle2(i,11))==1 & thetacycle2(i,10)>=10
        a=squeeze(thetaseq2(i,:,:))';
        imagesc([a,a]);set(gca,'YDir','normal');yline([15,50],'c');
        if thetacycle2(i,69)==thetacycle2(i,26)
            yline(thetacycle2(i,27)/2,'r');
        elseif thetacycle2(i,71)==thetacycle2(i,26)
            yline(thetacycle2(i,27)/2+35,'r');
        end
        title(num2str(thetacycle2(i,11)));
        pause;clf;
    end
end
%----------------------------------------------------------------------------
% 2. reward/inter pause
indv=thetacycle1(:,10)<10;
indinter=abs(thetacycle1(:,11))==0.5;
indrew=thetacycle1(:,11)==0;
indlocal=thetacycle1(:,69)==thetacycle1(:,26);
indremote=thetacycle1(:,69)~=thetacycle1(:,26);
thetaseq_align1=[];thetaseq_noalign1=[];
for dir=1:2
    if dir==1
        indd=indrew;
    else
        indd=indinter;
    end
    for k=1:2
        if k==1
            ind=indv & indd & indlocal;% 
        else
            ind=indv & indd & indremote;% 
        end
        thetacycle=thetacycle1(ind,:);
        thetaseq=thetaseq1(ind,:,:);
        thetaseq_align=align_thetaseq(thetacycle,thetaseq);
        thetaseq_align1{dir,k}=thetaseq_align;
        thetacycle_align1{dir,k}=thetacycle;
        thetaseq_noalign1{dir,k}=thetaseq;
    end
end
figure;
for in=1:16
    ass=[];
    for k=1:2
        as=[];
        for dir=1:2
            a=thetaseq_noalign1{dir,k};
            thetacycle=thetacycle_align1{dir,k};
            ind=thetacycle(:,42)==in;
            a1=squeeze(nanmean(a(ind,:,:),1))';
            as=[as,a1,a1];
        end
        ass=[ass,as];
    end
    subplot(4,4,in);imagesc([0 8]*360,[-100 100],ass);set(gca,'YDir','normal');
    xline(360*2*[1:4],'k');ylim([-1 1]*40);
end

group_local_ratio=[];
for i=1:9
    ind=thetacycleall(:,39)==i-1;
    group_local_ratio(i,1)=nanmean(thetacycleall(ind,26)==thetacycleall(ind,61));
    group_local_ratio(i,2)=nanmean(thetacycleall(ind,26)==thetacycleall(ind,61) | (thetacycleall(ind,26)==thetacycleall(ind,63) & thetacycleall(ind,64)>=0.1) );
    group_local_ratio(i,3)=nanmean(thetacycleall(ind,26)==thetacycleall(ind,61) | (thetacycleall(ind,26)==thetacycleall(ind,63) & thetacycleall(ind,64)>=0.1) | (thetacycleall(ind,26)==thetacycleall(ind,65) & thetacycleall(ind,66)>=0.1) );
end
figure;
yyaxis left;histogram(thetacycleall(:,39),[-0.5:8.5]);
xticks([0:8]);xticklabels({'None','Reward Pause','Reward Run','Inter Pause','Inter Run','In Pause','In Run','Out Pause','Out run'});
xlabel('Thetacycle Group');ylabel('Theta Cycle N');
yyaxis right;
for i=1:3
    hold on;plot([0:8],group_local_ratio(:,i),'o-');
end
ylabel('Local Encoding Percent');yline(0.5,'k--');

figure
subplot(2,2,1);histogram(thetacycle0(:,26),[-0.5:8.5]);xlabel('Current arm or center');title('thetacycle0');
subplot(2,2,2);
a=histcounts2(thetacycle1(:,69),thetacycle1(:,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm');title('thetacycle1');
subplot(2,2,3);
a=histcounts2(thetacycle2(:,69),thetacycle2(:,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm 1');title('thetacycle2');
subplot(2,2,4);
a=histcounts2(thetacycle2(:,71),thetacycle2(:,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm 2');title('thetacycle2');
%---------------------------------------------------------------------------------
% below average thetaseq specifically for different groups

% 1. thetaseq in/out, run, local
indv=thetacycleall(:,10)>=10;
indin=thetacycleall(:,11)==-1;
indout=thetacycleall(:,11)==1;
thetaseq1=[];phasebin=10;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_2','linear_maze','deseqall');
    deseqall(:,27)=ceil(deseqall(:,27)/phasebin);
    ind=thetacycleall(:,42)==in ; 
    for k=1%:2
        if k==1
            indk=thetacycleall(:,60)==k & indout & indv & ind & thetacycleall(:,26)==thetacycleall(:,69);
        elseif k==2
            indk=thetacycleall(:,60)==k & indout & indv & ind & (thetacycleall(:,26)==thetacycleall(:,69) | thetacycleall(:,26)==thetacycleall(:,71));
        end
        seq=[];
        for pos=16:30
            indp=indk & ceil(thetacycleall(:,27)/2)==pos;
            cycle=thetacycleall(indp,:);
            lia=ismember(deseqall(:,40),cycle(:,21));
            for phase=1:360/phasebin
                indph=deseqall(:,27)==phase;
                if k==1
                    seq(pos-15,phase,1:15)=squeeze(nanmean(nansum(linear_maze( lia & indph, 1:15, :),1),3));
                    seq(pos-15,phase,16:50)=squeeze(nansum(linear_maze( lia & indph, 16:50, cycle(i,69)),1));
                elseif k==2
                    seq(i,phase,1:15)=squeeze(nanmean(nanmean(linear_maze( lia & indph, 1:15, :),1),3));
                    seq(i,phase,16:50)=squeeze(nanmean(linear_maze( lia & indph, 16:50, cycle(i,69)),1));
                    seq(i,phase,51:85)=squeeze(nanmean(linear_maze( lia & indph, 16:50, cycle(i,71)),1));
                end
            end
        end
        thetaseq1(:,:,:,in)=seq;
    end
    in
end

figure;centerbin=51;thetaseq1_align=[];
for in=1:16
    as=zeros(101,36,15);as2=[];
    for pos=1:15
        a=squeeze(thetaseq1(pos,:,:,in))';
        p1=centerbin-(pos+15)+1;
        p2=centerbin-(pos+15)+50;
        as(p1:p2,:,pos)=a;
        as2=[as2,as(:,:,pos)];
    end
    thetaseq1_align(:,:,in)=nansum(as,3);
    figure(1);subplot(4,4,in);a=nansum(as,3);imagesc([a,a]);set(gca,'YDir','normal');ylim([20 80]);yline(centerbin,'r');
    figure(2);subplot(8,2,in);imagesc(as2);set(gca,'YDir','normal');ylim([20 80]);yline(centerbin,'r');
end
%----------------------------------------------------------------------------
% below focus on reward pause when only 1 arm is decoded
ind=thetacycleall(:,39)==1 & thetacycleall(:,40)==1;n=sum(ind);
thetacycle1=thetacycleall(ind,:);
depos_dist1=depos_distall(ind,:,:);
thetacycle1(:,53:54)=nan;
[row,col]=find(squeeze(depos_dist1(:,1,:)));
[~,I]=sort(row);
lia=ismember([1:n],row);
if sum(~lia)==0 & length(row)==n
    thetacycle1(:,53)=col(I)-1; % encoded arm
else
    disp(['Error with encoded arm!']);
end
[row,col]=find(squeeze(depos_dist1(:,2,:)));
[~,I]=sort(row);row=row(I);col=col(I);
if sum(diff(row)==0)>0
    disp(['Encoded Visit ID has duplicates!']);
else
    lia=ismember([1:n],row);
    thetacycle1(lia,54)=col-1; % encoded visit id
end

indcorrect=thetacycle1(:,56)==1;
figure;
subplot(2,2,1);a=histcounts2(thetacycle1(indcorrect,53),thetacycle1(indcorrect,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm');
title('Reward Pause : Correct Trials');colorbar;
subplot(2,2,2);a=histcounts2(thetacycle1(indcorrect,54),round(thetacycle1(indcorrect,28)),[-0.5:8.5],[-.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current visit ID');ylabel('Encoded Visit ID');
title('Reward Pause : Correct Trials');colorbar;
subplot(2,2,3);a=histcounts2(thetacycle1(~indcorrect,53),thetacycle1(~indcorrect,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm');title('Reward Pause');
title('Reward Pause : Error Trials');colorbar;
subplot(2,2,4);a=histcounts2(thetacycle1(~indcorrect,54),round(thetacycle1(~indcorrect,28)),[-0.5:8.5],[-7.5:8.5]);
imagesc([-7 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current visit ID');ylabel('Encoded Visit ID');
title('Reward Pause : Error Trials');colorbar;

%  below for theta sequences
thetaseq1=[];phasebin=20;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info','linear_maze','deseqall');
    ind=thetacycle1(:,42)==in;
    cycle=thetacycle1(ind,:);
    deseqall(:,27)=ceil(deseqall(:,27)/phasebin); 
    seq=[];
    for phase=1:360/phasebin
        indph=deseqall(:,27)==phase;
        for i=1:size(cycle,1)
            lia=ismember(deseqall(:,40),cycle(i,21));
            if cycle(i,53)>0
                seq(i,phase,:)=squeeze(nanmean(linear_maze( lia & indph, :, cycle(i,53)),1));
            else
                seq(i,phase,:)=squeeze(nanmean(nanmean(linear_maze( lia & indph, :, :),1),3));
            end
        end
    end
    thetaseq1=cat(1,thetaseq1,seq);
    in
end
thetamax=[];
for i=1:18
    [M,I]=max(squeeze(thetaseq1(:,i,:)),[],2);
    thetamax(:,i,1)=I;
    thetamax(:,i,2)=M;
end
thetaseqcorr=[];
for i=1:size(thetaseq1,1)
    a=thetamax(i,:,1);
    [r,p]=corrcoef([1:18],a);
    thetaseqcorr(i,:)=[r(2),p(2)];
end
indlocal=thetacycle1(:,53)==thetacycle1(:,26);
figure;rcutoff=0.3;indpos= thetacycle1(:,27)>=60 & thetacycle1(:,27)<=100 & thetacycle1(:,4)>2;
for in=1:16
for k=1:2
    if k==1
        indk=indlocal & thetacycle1(:,42)==in & indpos;
    else
        indk=~indlocal & thetacycle1(:,42)==in & indpos;
    end
    subplot(2,3,1+(k-1)*3);
    ind=thetaseqcorr(:,1)>rcutoff & indk;
    a=squeeze(nanmean(thetaseq1(ind,:,:),1))';
    imagesc([0 720],[0 100],[a,a]);set(gca,'YDir','normal');title(['+ : N = ',num2str(sum(ind))]);
    
    subplot(2,3,2+(k-1)*3);
    ind=thetaseqcorr(:,1)<-rcutoff & indk;
    a=squeeze(nanmean(thetaseq1(ind,:,:),1))';
    imagesc([0 720],[0 100],[a,a]);set(gca,'YDir','normal');title(['- : N = ',num2str(sum(ind))]);
    
    subplot(2,3,3+(k-1)*3);
    ind=abs(thetaseqcorr(:,1))<=rcutoff & indk;
    a=squeeze(nanmean(thetaseq1(ind,:,:),1))';
    imagesc([0 720],[0 100],[a,a]);set(gca,'YDir','normal');title(['0 : N = ',num2str(sum(ind))]);
end
pause;clf;
end

figure;rcutoff=0.3;indpos= thetacycle1(:,27)>=60 & thetacycle1(:,27)<=100 & thetacycle1(:,4)>1;
for in=1:16
    indk=thetacycle1(:,42)==in & indpos;
    subplot(1,3,1);
    ind=thetaseqcorr(:,1)>rcutoff & indk;
    a=squeeze(nanmean(thetaseq1(ind,:,:),1))';
    imagesc([0 720],[0 100],[a,a]);set(gca,'YDir','normal');title(['+ : N = ',num2str(sum(ind))]);
    
    subplot(1,3,2);
    ind=thetaseqcorr(:,1)<-rcutoff & indk;
    a=squeeze(nanmean(thetaseq1(ind,:,:),1))';
    imagesc([0 720],[0 100],[a,a]);set(gca,'YDir','normal');title(['- : N = ',num2str(sum(ind))]);
    
    subplot(1,3,3);
    ind=abs(thetaseqcorr(:,1))<=rcutoff & indk;
    a=squeeze(nanmean(thetaseq1(ind,:,:),1))';
    imagesc([0 720],[0 100],[a,a]);set(gca,'YDir','normal');title(['0 : N = ',num2str(sum(ind))]);
    pause;clf;
end

figure;[~,I]=sort(thetacycle1(:,4),'descend');
for ii=1:size(thetaseq1,1)
    i=I(ii);
    a=squeeze(thetaseq1(i,:,:))';
    a=[a,a];
    imagesc([0 720],[0 100],a);set(gca,'YDir','normal');
    a=thetamax(i,:,1);ind=~isnan(squeeze(thetamax(i,:,2)));
    [r,p]=corrcoef([1:18],a);%find(ind),
    xlabel('Theta Phase');ylabel('Position (cm)');title(num2str([r(2),p(2)],2));
    pause;clf;
end
%-----------------------------------------------------------------------------
% below use a new way to organize data
thetacycleall(:,57)=sum(squeeze(depos_distall(:,1,2:9))>0,2);
thetacycleall(:,58)=max(squeeze(depos_distall(:,1,:)),[],2)./thetacycleall(:,24); % dominance index
ind=thetacycleall(:,39)==1 ;n=sum(ind);
thetacycle_rewp=thetacycleall(ind,:);
depos_dist_rewp=depos_distall(ind,:,:);
thetacycle_rewp(:,59)=nan;

ind1=thetacycle_rewp(:,58)>=0.7;
thetacycle_rewp(ind1,59)=1;
thetacycle1=thetacycle_rewp(ind1,:);
depos_dist1=depos_dist_rewp(ind1,:,:);
thetacycle_rewp(:,53:54)=nan;
[M1,I]=max(squeeze(depos_dist1(:,1,:)),[],2);
thetacycle_rewp(ind1,53)=I-1;
[M2,I]=max(squeeze(depos_dist1(:,2,:)),[],2);
I(M1>M2)=nan;
thetacycle_rewp(ind1,54)=I-1;

ind2=thetacycle_rewp(:,58)<0.7;

[~,I]=sort(row);
lia=ismember([1:n],row);
if sum(~lia)==0 & length(row)==n
    thetacycle1(:,53)=col(I)-1; % encoded arm
else
    disp(['Error with encoded arm!']);
end
[row,col]=find(squeeze(depos_dist1(:,2,:)));
[~,I]=sort(row);row=row(I);col=col(I);
if sum(diff(row)==0)>0
    disp(['Encoded Visit ID has duplicates!']);
else
    lia=ismember([1:n],row);
    thetacycle1(lia,54)=col-1; % encoded visit id
end


cyclegroup_N=[];
for den=0:3
    ind=thetacycleall(:,39)==1 & thetacycleall(:,57)==den;
    n=sum(ind);cyclegroup_N(den+1)=n;
    thetacycle=thetacycleall(ind,:);
    depos_dist=depos_distall(ind,:,:);

ind=thetacycleall(:,39)==1 & thetacycleall(:,40)==2;n=sum(ind);
thetacycle2=thetacycleall(ind,:);
depos_dist2=depos_distall(ind,:,:);
decycle=nan*ones(n,2,2);
[row,col]=find(squeeze(depos_dist2(:,1,:)));
[~,I]=sort(row);
lia=ismember([1:n],row);
if sum(~lia)==0 & length(row)==2*n
    decycle(:,:,1)=reshape(col(I)-1,[2,n])'; % encoded arm
else
    disp(['Error with encoded arm!']);
end
[row,col]=find(squeeze(depos_dist1(:,2,:)));
[~,I]=sort(row);row=row(I);col=col(I);
if sum(diff(row)==0)>0
    disp(['Encoded Visit ID has duplicates!']);
else
    lia=ismember([1:n],row);
    thetacycle1(lia,54)=col-1; % encoded visit id
end

indcorrect=thetacycle1(:,56)==1;
figure;
subplot(2,2,1);a=histcounts2(thetacycle1(indcorrect,53),thetacycle1(indcorrect,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm');
title('Reward Pause : Correct Trials');colorbar;
subplot(2,2,2);a=histcounts2(thetacycle1(indcorrect,54),round(thetacycle1(indcorrect,28)),[-0.5:8.5],[-.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current visit ID');ylabel('Encoded Visit ID');
title('Reward Pause : Correct Trials');colorbar;
subplot(2,2,3);a=histcounts2(thetacycle1(~indcorrect,53),thetacycle1(~indcorrect,26),[-0.5:8.5],[-0.5:8.5]);
imagesc([0 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current Arm');ylabel('Encoded Arm');title('Reward Pause');
title('Reward Pause : Error Trials');colorbar;
subplot(2,2,4);a=histcounts2(thetacycle1(~indcorrect,54),round(thetacycle1(~indcorrect,28)),[-0.5:8.5],[-7.5:8.5]);
imagesc([-7 8],[0 8],a);set(gca,'YDir','normal');xlabel('Current visit ID');ylabel('Encoded Visit ID');
title('Reward Pause : Error Trials');colorbar;
  
%--------------------------------------------------------------------------------
% below compare properties of thetacyles btw pause vs run
thetacycleall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Analyze_Theta_Cycle_Broad','thetacycle','propname','propset');
    ind=thetacycle(:,4)>0 & thetacycle(:,3)<0.2;thetacycle=thetacycle(ind,:);
    thetacycle(:,21)=[1:size(thetacycle,1)];thetacycle(:,22)=in;
    thetacycleall=[thetacycleall;thetacycle];
end
propset=[3,4,7,5,8,15,14,16:20,10];
propname={'Duration','Theta Power','Peak Percentile','Trough Percentile','Theta Prom','Ascending Half Dur',...
    'Descending Half Dur','Descend/Ascend','MML Sink','MML Sink Percentile','MML Source','MML Source Percentile','Velocity'};
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG IML/GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};

indv=thetacycleall(:,10)<10;binN=80;
figure;
for p=1:13
    pid=propset(p);
    subplot(3,5,p);histogram(thetacycleall(indv,pid),binN);
    hold on;histogram(thetacycleall(~indv,pid),binN);
    m1=max(thetacycleall(:,pid));m2=min(thetacycleall(:,pid));
    xline([m1,m2],'r');[h,pv]=ttest2(thetacycleall(indv,pid),thetacycleall(~indv,pid));
    title([propname{p},': P = ',num2str(pv,2)]);
end
legend('Pause','Run');
end
end