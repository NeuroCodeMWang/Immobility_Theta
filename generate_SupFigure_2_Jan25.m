function generate_SupFigure_2_Jan25(ifsavedata)

if nargin<1
    ifsavedata=0;
end

if 0
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
% below confirmed that pksets saved from two places are identical except
% for visit id 
pksetsall0=[];pksetsall=[];interrunall0=[];interrunall1=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('ThetaCycle_Decode_Info_dHPC2','pksets','inter_run');
    pksets(:,end+1)=in;pksetsall=[pksetsall;pksets];
    inter_run(:,end+1)=in;interrunall1=[interrunall1;inter_run];
    load('InOutBound_Behavior_Analysis.mat','pksets','inter_run');
    pksets(:,end+1)=in;pksetsall0=[pksetsall0;pksets];
    inter_run(:,end+1)=in;interrunall0=[interrunall0;inter_run];
end
% below recalculate the next visit id for inter_run
interrunall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    %load('ThetaCycle_Decode_Info_dHPC2','pksets','inter_run');
    load('InOutBound_Behavior_Analysis.mat','pksets','inter_run');
    inter_run(:,9)=pksets(inter_run(:,2)+1,7); % next visit id
    a=pksets(inter_run(:,2)+1,1);
    ind=a~=inter_run(:,1);
    inter_run(ind,9)=0;
    inter_run(:,end+1)=in;interrunall=[interrunall;inter_run];
end
% only consider pk with good behaviors
indexcludepk= ( round(pksetsall(:,7))>-4 & round(pksetsall(:,7))<=0 ) | ( round(pksetsall(:,7))~=pksetsall(:,7) );
indexcludeinter= ( round(interrunall(:,9))>-4 & round(interrunall(:,9))<=0 ) | abs(round(interrunall(:,9))-interrunall(:,9))>0.1 ; 
pksetsall(:,43)=1;pksetsall(indexcludepk,43)=0;
interrunall(:,22)=1;interrunall(indexcludeinter,22)=0;
cd(homedir);save('pksets_interrun_new','pksetsall','interrunall');
figure;
subplot(1,2,1);
a=pksetsall(~indexcludepk,12)-pksetsall(~indexcludepk,10);histogram(a,[0:2:80]);xline(30,'r');
xlabel('Time Spent at Reward Site (s)');ylabel('Count');title(num2str(mean(a<=30),2));
subplot(1,2,2);
a=interrunall(~indexcludeinter,7)-interrunall(~indexcludeinter,5);histogram(a,[0:1:30]);xline(20,'r');
xlabel('Time Spent at Choice Point (s)');ylabel('Count');title(num2str(mean(a<=20),2));
end
%----------------------------------------------------------------------------------
if 0
% below to quantify velocity & theta power during behaviors
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
load('pksets_interrun_new');
velocity_cutoff=0;binum=10;posbin=round(40/binum);region_num=11;
theta_outset=[];theta_inset=[];theta_interset=[];theta_rewardset=[];
vel_outall=[];vel_inall=[];vel_rewardall=[];vel_interall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_delta_LFP_HPC.mat', 'lfptime','thetapower','deltapower');
    load('new_reward_center_in_out','pos_behave_marker');pos_behave_marker=pos_behave_marker(:,3:4);
    %load('InOutBound_Behavior_Analysis.mat', 'pos_behave_marker');
    load('Position_Data_Maze.mat'); 
    ind=interrunall(:,21)==in;inter_run=interrunall(ind,:);
    ind=pksetsall(:,42)==in;pksets=pksetsall(ind,:);
    [velocity,vel_angle,locdis,loc_angle,locdis_change,~,vel_angle_accelerate]=calculate_location_movement(Position_Data);
    thetadelta=thetapower./deltapower;
    thetapower=zscore(thetapower,0,1);
    ind=lfptime(:,2)>0;lfptime_trial=lfptime(ind,:);thetapower=thetapower(ind,:);thetadelta=thetadelta(ind,:);
    %for i=1:11
        %thetadelta(:,i) = smooth(lfptime_trial(:,1),thetadelta(:,i),0.01,'rloess');
    %end
    theta_out=zeros(size(pksets,1),binum,region_num,2,2);
    theta_in=zeros(size(pksets,1),binum,region_num,2,2);
    theta_reward=zeros(size(pksets,1),binum,region_num,2);
    theta_inter=zeros(size(inter_run,1),binum,region_num,2);
    vel_out=zeros(size(pksets,1),binum,2);
    vel_in=zeros(size(pksets,1),binum,2);
    vel_reward=zeros(size(pksets,1),binum);
    vel_inter=zeros(size(inter_run,1),binum);
    for i=1:size(pksets,1)
    % outbound
    ind= pos_behave_marker(:,1)==4 & pos_behave_marker(:,2)==i & velocity>=velocity_cutoff &...
        locdis_change>0;
    if sum(ind)>0
        p1=pksets(i,17);
        p2=pksets(i,18);
        a=ceil((locdis-p1)/(p2-p1)*binum);
        a1=ceil((locdis-30)/posbin);
        for j=1:binum
            ind0=ind & a==j;
            lia=ismember(lfptime_trial(:,4),find(ind0==1));
            theta_out(i,j,:,1,1)=nanmean(thetapower(lia,:),1);
            theta_out(i,j,:,1,2)=nanmean(thetadelta(lia,:),1);
            vel_out(i,j,1)=nanmean(velocity(ind0));
            
            ind0=ind & a1==j;
            lia=ismember(lfptime_trial(:,4),find(ind0==1));
            theta_out(i,j,:,2,1)=nanmean(thetapower(lia,:),1);
            theta_out(i,j,:,2,2)=nanmean(thetadelta(lia,:),1);
            vel_out(i,j,2)=nanmean(velocity(ind0));
        end
    end

    % inbound
    ind= pos_behave_marker(:,1)==3 & pos_behave_marker(:,2)==i & ...
        velocity>=velocity_cutoff & locdis_change<0;
    if sum(ind)>0
        p1=pksets(i,19);
        p2=pksets(i,20);
        a=ceil((locdis-p1)/(p2-p1)*binum);
        a1=ceil((locdis-30)/posbin);
        for j=1:binum
            ind0=ind & a==j;
            lia=ismember(lfptime_trial(:,4),find(ind0==1));
            theta_in(i,j,:,1,1)=nanmean(thetapower(lia,:),1);
            theta_in(i,j,:,1,2)=nanmean(thetadelta(lia,:),1);
            vel_in(i,j,1)=nanmean(velocity(ind0));

            ind0=ind & a1==j;
            lia=ismember(lfptime_trial(:,4),find(ind0==1));
            theta_in(i,j,:,2,1)=nanmean(thetapower(lia,:),1);
            theta_in(i,j,:,2,2)=nanmean(thetadelta(lia,:),1);
            vel_in(i,j,2)=nanmean(velocity(ind0));
        end
    end

    % reward
    ind= pos_behave_marker(:,1)==1 & pos_behave_marker(:,2)==i ;
    if sum(ind)>0
        posid=find(ind==1);
        posid=posid(:);
        posid(:,2)=ceil([1:size(posid,1)]*binum/size(posid,1));
        for j=1:binum
            ind0=posid(:,2)==j;
            lia=ismember(lfptime_trial(:,4),posid(ind0,1));
            theta_reward(i,j,:,1)=nanmean(thetapower(lia,:),1);
            theta_reward(i,j,:,2)=nanmean(thetadelta(lia,:),1);
            vel_reward(i,j)=nanmean(velocity(posid(ind0,1)));
        end
    end
    end

    for i=1:size(inter_run,1)
        % inter run
    ind= abs(pos_behave_marker(:,1))==2 & pos_behave_marker(:,2)==i ;
    if sum(ind)>0
        posid=find(ind==1);
        posid=posid(:);
        posid(:,2)=ceil([1:size(posid,1)]*binum/size(posid,1));
        for j=1:binum
            ind0=posid(:,2)==j;
            lia=ismember(lfptime_trial(:,4),posid(ind0,1));
            theta_inter(i,j,:,1)=nanmean(thetapower(lia,:),1);
            theta_inter(i,j,:,2)=nanmean(thetadelta(lia,:),1);
            vel_inter(i,j)=nanmean(velocity(posid(ind0,1)));
        end    
    end
    end
    theta_outset{in}=theta_out;
    theta_inset{in}=theta_in;
    theta_rewardset{in}=theta_reward;
    theta_interset{in}=theta_inter;
    vel_outall=cat(1,vel_outall,vel_out);
    vel_inall=cat(1,vel_inall,vel_in);
    vel_rewardall=[vel_rewardall;vel_reward];
    vel_interall=[vel_interall;vel_inter];
    disp(['Session Completed: ',num2str(in)]);
end

theta_out_mean=[];theta_in_mean=[];theta_reward_mean=[];theta_inter_mean=[];
vel_out_mean=[];vel_in_mean=[];vel_reward_mean=[];vel_inter_mean=[];
runduration_cutoff=4;
for in=1:16
    ind=interrunall(:,21)==in;
    inter_run=interrunall(ind,:);
    vel_inter=vel_interall(ind,:);
    ind=pksetsall(:,42)==in;
    pksets=pksetsall(ind,:);
    vel_out=vel_outall(ind,:,:);
    vel_in=vel_inall(ind,:,:);
    vel_reward=vel_rewardall(ind,:);
    theta_out=theta_outset{in};
    theta_in=theta_inset{in};
    theta_reward=theta_rewardset{in};
    theta_inter=theta_interset{in};
    pksets(:,7)=round(pksets(:,7));
    inter_run(:,4)=round(inter_run(:,4));
    for i=1:3
    if i==1 % correct forced visit
        indout= pksets(:,7)>0 & pksets(:,7)<=4 & pksets(:,43)==1;%& pksets(:,14)<=runduration_cutoff & pksets(:,17)<=30 & pksets(:,18)>=60 
        indin= pksets(:,7)>0 & pksets(:,7)<=4 & pksets(:,43)==1;%& pksets(:,15)<=runduration_cutoff & pksets(:,19)<=30 & pksets(:,20)>=60 
        indreward= pksets(:,7)>0 & pksets(:,7)<=4 & pksets(:,43)==1;%& pksets(:,11)>=60 & pksets(:,13)>=60 
        indinter= inter_run(:,4)>0 & inter_run(:,4)<4 & inter_run(:,11)<=60 & inter_run(:,22)==1;
    elseif i==2 % correct free visit
        indout= pksets(:,7)>4 & pksets(:,7)<=8 & pksets(:,43)==1;%& pksets(:,14)<=runduration_cutoff & pksets(:,17)<=30 & pksets(:,18)>=60 
        indin= pksets(:,7)>4 & pksets(:,7)<=8 & pksets(:,43)==1;%& pksets(:,15)<=runduration_cutoff & pksets(:,19)<=30 & pksets(:,20)>=60 
        indreward= pksets(:,7)>4 & pksets(:,7)<=8 & pksets(:,43)==1;%& pksets(:,11)>=60 & pksets(:,13)>=60 
        indinter= inter_run(:,4)>=4 & inter_run(:,4)<8 & inter_run(:,22)==1;%& inter_run(:,11)<=40 
    elseif i==3 % error visit
        indout= pksets(:,7)<-3 & pksets(:,43)==1;%& pksets(:,14)<=runduration_cutoff & pksets(:,17)<=30 & pksets(:,18)>=60 
        indin= pksets(:,7)<-3 & pksets(:,43)==1;%& pksets(:,15)<=runduration_cutoff & pksets(:,19)<=30 & pksets(:,20)>=60 
        indreward= pksets(:,7)<-3 & pksets(:,43)==1;%& pksets(:,11)>=60 & pksets(:,13)>=60 
        indinter= inter_run(:,4)<-3 & inter_run(:,22)==1;%& inter_run(:,11)<=40 
    end
    theta_out_mean(i,:,:,:,:,in)=squeeze(nanmean(theta_out(indout,:,:,:,:),1));
    theta_in_mean(i,:,:,:,:,in)=squeeze(nanmean(theta_in(indin,:,:,:,:),1));
    theta_reward_mean(i,:,:,:,in)=squeeze(nanmean(theta_reward(indreward,:,:,:),1));
    theta_inter_mean(i,:,:,:,in)=squeeze(nanmean(theta_inter(indinter,:,:,:),1));

    vel_out_mean(i,:,:,in)=squeeze(nanmean(vel_out(indout,:,:),1));
    vel_in_mean(i,:,:,in)=squeeze(nanmean(vel_in(indin,:,:),1));
    vel_reward_mean(i,:,in)=squeeze(nanmean(vel_reward(indreward,:),1));
    vel_inter_mean(i,:,in)=squeeze(nanmean(vel_inter(indinter,:),1));
    end
end
cd(homedir);save('theta_velocity_inoutbound_inter_reward_new',"theta_interset",'theta_rewardset','theta_inset',...
    'theta_outset','theta_inter_mean',"theta_reward_mean","theta_in_mean","theta_out_mean","vel_interall",...
    'vel_rewardall','vel_inall','vel_outall','vel_inter_mean',"vel_reward_mean","vel_in_mean","vel_out_mean",'-v7.3');
end

homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);
load('SessionSet16');
load('theta_velocity_inoutbound_inter_reward_new');
HPC_layer_name={'CA1 sp','CA1 sr','CA1 slm','DG oml','DG mml','DG gcl','CA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
powertypename={'Theta Power (z)','Theta / Delta Power Ratio'};
plotcolors=[0 0.4470 0.7410];
plotcolors(2,:)=[0.8500 0.3250 0.0980];
plotcolors(3,:)=0;
powerid=2;ns=[];
for postype=1:2
figure;
for i=1:2
    subplot(3,4,10+i);
    plot(0,0,'-','Color',plotcolors(1,:));hold on;
    plot(0,0,'-','Color',plotcolors(2,:));hold on;plot(0,0,'k-');
end
for h=1:11
    subplot(3,4,h);
    for i=1:3
        a1=squeeze(theta_out_mean(i,:,h,postype,powerid,:));
        a2=squeeze(theta_in_mean(i,:,h,postype,powerid,:));
        a3=squeeze(theta_reward_mean(i,:,h,powerid,:));
        a4=squeeze(theta_inter_mean(i,:,h,powerid,:));
        %% 
        hold on;shaded_errbar([0.1:0.1:1]+1,a1,plotcolors(i,:));
        hold on;shaded_errbar(-[0.1:0.1:1]+4,a2,plotcolors(i,:));
        hold on;shaded_errbar([0.1:0.1:1],a4,plotcolors(i,:));
        hold on;shaded_errbar(2+[0.1:0.1:1],a3,plotcolors(i,:));
    end
    b1=squeeze(nanmean(theta_out_mean(3,:,h,postype,powerid,:),2));
    b2=squeeze(nanmean(nanmean(theta_out_mean(1:2,:,h,postype,powerid,:),2),1));
    p1=signrank(b1-b2);

    b1=squeeze(nanmean(theta_in_mean(3,:,h,postype,powerid,:),2));
    b2=squeeze(nanmean(nanmean(theta_in_mean(1:2,:,h,postype,powerid,:),2),1));
    p2=signrank(b1-b2);

    b1=squeeze(nanmean(theta_reward_mean(3,:,h,powerid,:),2));
    b2=squeeze(nanmean(nanmean(theta_reward_mean(1:2,:,h,powerid,:),2),1));
    p3=signrank(b1-b2);

    b1=squeeze(nanmean(theta_inter_mean(3,:,h,powerid,:),2));
    b2=squeeze(nanmean(nanmean(theta_inter_mean(1:2,:,h,powerid,:),2),1));
    p4=signrank(b1-b2);
    xticks([0.5:1:3.5]);xticklabels({'Center','Outbound','Reward','Inbound'});
    xlabel('Position or Time Index : Choice Point | Out | Reward | In');ylim([0 21]);
    ylabel(powertypename{powerid});title([HPC_layer_name{h} '; P = ',num2str([p4,p1,p3,p2],2)]);
end
legend('Forced Visit','Free Visit','Error Visit');
% below plot velocity profile
subplot(3,4,12);
for i=1:3
    a1=squeeze(vel_out_mean(i,:,postype,:));
    a2=squeeze(vel_in_mean(i,:,postype,:));
    a3=squeeze(vel_reward_mean(i,:,:));
    a4=squeeze(vel_inter_mean(i,:,:));
    hold on;shaded_errbar([0.1:0.1:1]+1,a1,plotcolors(i,:));
    hold on;shaded_errbar(-[0.1:0.1:1]+4,a2,plotcolors(i,:));
    hold on;shaded_errbar([0.1:0.1:1],a4,plotcolors(i,:));
    hold on;shaded_errbar(2+[0.1:0.1:1],a3,plotcolors(i,:));
end
b1=squeeze(nanmean(vel_out_mean(3,:,postype,:),2));
b2=squeeze(nanmean(nanmean(vel_out_mean(1:2,:,postype,:),2),1));
p1=signrank(b1-b2);

b1=squeeze(nanmean(vel_in_mean(3,:,postype,:),2));
b2=squeeze(nanmean(nanmean(vel_in_mean(1:2,:,postype,:),2),1));
p2=signrank(b1-b2);

b1=squeeze(nanmean(vel_reward_mean(3,:,:),2));
b2=squeeze(nanmean(nanmean(vel_reward_mean(1:2,:,:),2),1));
p3=signrank(b1-b2);

b1=squeeze(nanmean(vel_inter_mean(3,:,:),2));
b2=squeeze(nanmean(nanmean(vel_inter_mean(1:2,:,:),2),1));
p4=signrank(b1-b2);
xticks([0.5:1:3.5]);xticklabels({'Center','Outbound','Reward','Inbound'});
xlabel('Position or Time Index : Choice Point | Out | Reward | In');
ylabel('Velocity (cm/s)');title(['Velocity : P = ',num2str([p4,p1,p3,p2],2)]);
legend('Forced Visit','Free Visit','Error Visit');
figure_title=['SupFigure_2_Theta_Velocity_InOutInterReward_',num2str(postype)];save_current_figure(figure_title);
end

if ifsavedata
% below to save data
postype=2;powerid=2;hname={'dCA1 pyr_LFP','dCA1 st rad_LFP','dCA1 slm_LFP','DG OML_LFP','DG MML_LFP','DG GCL_LFP','dCA3_LFP'};
iname={'Forced','Choice','Error'};letternum={'a','b','c','d','e','f','g','h'};
colname={'Center 1','Center 2','Center 3','Center 4','Center 5','Center 6','Center 7','Center 8','Center 9','Center 10',...
    'Outbound 1','Outbound 2','Outbound 3','Outbound 4','Outbound 5','Outbound 6','Outbound 7','Outbound 8','Outbound 9','Outbound 10',...
    'Reward 1','Reward 2','Reward 3','Reward 4','Reward 5','Reward 6','Reward 7','Reward 8','Reward 9','Reward 10',...
    'Inbound 1','Inbound 2','Inbound 3','Inbound 4','Inbound 5','Inbound 6','Inbound 7','Inbound 8','Inbound 9','Inbound 10'};
rowname=[];
for i=1:3
    for in=1:16
        rowname{in+(i-1)*16}=[iname{i},' ',num2str(in)];
    end
end
for h=1:8
    figdata=[];
    for i=1:3
        if h==1
            a1=squeeze(vel_out_mean(i,:,postype,:));
            a2=squeeze(vel_in_mean(i,10:-1:1,postype,:));
            a3=squeeze(vel_reward_mean(i,:,:));
            a4=squeeze(vel_inter_mean(i,:,:));
            figdata=[figdata;[a4',a1',a3',a2']];
        else
            a1=squeeze(theta_out_mean(i,:,h-1,postype,powerid,:));
            a2=squeeze(theta_in_mean(i,10:-1:1,h-1,postype,powerid,:));
            a3=squeeze(theta_reward_mean(i,:,h-1,powerid,:));
            a4=squeeze(theta_inter_mean(i,:,h-1,powerid,:));
            figdata=[figdata;[a4',a1',a3',a2']];
        end
    end
    T = array2table(figdata,'RowNames',rowname,'VariableNames',colname);
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Sup Figure 2',letternum{h}],'WriteRowNames',true);
end
end