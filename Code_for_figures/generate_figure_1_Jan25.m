function generate_figure_1_Jan25(ifsavedata)
if nargin<1
    ifsavedata=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
plotcolors5=slanCM('binary',7);plotcolors5=plotcolors5([7:-2:1],:);
defaultcolor7=slanCM('deep',7);defaultcolor7(1,:)=[0.8828 0.7305 0.2656];
plotcolors3=[0 0.4470 0.7410];
plotcolors3(2,:)=[0.8500 0.3250 0.0980];
plotcolors3(3,:)=0;
plotcolors3c=plotcolors3;plotcolors3c(3,:)=[0.8828 0.7305 0.2656];
%------------------------------------------------------------------------------
% 2. figure 1
vname={'Run','Pause','Immobile'};velocity_cutoff=10;
if 0
    velmat=[];
    for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('new_reward_center_in_out','pos_behave_marker');pos_behave_marker=pos_behave_marker(:,3:4);
    %load('InOutBound_Behavior_Analysis.mat', 'pos_behave_marker');
    load('Position_Data_Maze.mat');
    for b=1:4
        if b==1
            indv=pos_behave_marker(:,1)==0;
        elseif b==2
            indv=abs(pos_behave_marker(:,1))==0.5;
        elseif b==3
            indv=pos_behave_marker(:,1)==1;
        elseif b==4
            indv=pos_behave_marker(:,1)==-1;
        end
        velmat(in,b)=[nanmean(Position_Data(indv,5))];
    end
end

    b=3;specv=[];thedelv=[];
    for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('spectrum_HPC_layer_trial_sleep_2s','Etimeset','specset','theta_delta_power_set','f');
    Etime=Etimeset{b};
    spec=specset{b};
    theta_delta_power=theta_delta_power_set{b};
    indstill0=Etime(:,4)<1;
    indstill=Etime(:,4)<velocity_cutoff;
    indrun=Etime(:,4)>=velocity_cutoff;
    specv(:,:,1,in)=squeeze(nanmean(spec(:,indstill,:),2));
    specv(:,:,2,in)=squeeze(nanmean(spec(:,indrun,:),2));
    specv(:,:,3,in)=squeeze(nanmean(spec(:,indstill0,:),2));
    thedelv(:,:,1,in)=squeeze(nanmean(theta_delta_power(indstill,:,:),1))';
    thedelv(:,:,2,in)=squeeze(nanmean(theta_delta_power(indrun,:,:),1))';
    thedelv(:,:,3,in)=squeeze(nanmean(theta_delta_power(indstill0,:,:),1))';
    %subplot(4,4,in);histogram(Etime(:,4));title([num2str([mean(Etime(:,4)<5),mean(Etime(:,4)>10)])]);
end
    cd(homedir);save('Figure1_Data','thedelv','specv','velmat','f');
else
    cd(homedir);load('Figure1_Data');
    load('Rem_SW_sleep_analysis_2s','spec_sleep','thedel_sleep','f');
end
%----------------------------------------------------------------------------------
% below generate figure 1 c & e
pr=2;pc=2;regionnum=7;fign=1;
figure(fign);subplot(pr,pc,4); % theta/delta ratio
comp=[];stats=[];
for i=1:regionnum
    if 0
    [p,~,stat]=signrank(squeeze(thedelv(3,i,2,:)),1,"tail","right",'method','exact'); % movement
    comp(i,1)=p;stats(i,1)=stat.signedrank;
    [p,~,stat]=signrank(squeeze(thedelv(3,i,3,:)),1,"tail","right",'method','exact'); % pause
    comp(i,2)=p;stats(i,2)=stat.signedrank;
    [p,~,stat]=signrank(squeeze(thedel_sleep(3,i,3,:)),1,"tail","right",'method','exact'); % sws
    comp(i,3)=p;stats(i,3)=stat.signedrank;
    else
    [p,~,stat]=signrank(squeeze(thedelv(3,i,2,:)),1,"tail","right",'method','approximate'); % movement
    comp(i,1)=p;stats(i,1)=stat.zval;
    [p,~,stat]=signrank(squeeze(thedelv(3,i,3,:)),1,"tail","right",'method','approximate'); % pause
    comp(i,2)=p;stats(i,2)=stat.zval;
    [p,~,stat]=signrank(squeeze(thedel_sleep(3,i,3,:)),1,"tail","right",'method','approximate'); % sws
    comp(i,3)=p;stats(i,3)=stat.zval;
    end
end
hold on;shaded_errbar([1:regionnum],squeeze(thedelv(3,1:regionnum,2,:)),plotcolors3(1,:));
hold on;shaded_errbar([1:regionnum],squeeze(thedelv(3,1:regionnum,3,:)),plotcolors3(2,:));
hold on;shaded_errbar([1:regionnum],squeeze(thedel_sleep(3,1:regionnum,3,:)),plotcolors3c(3,:));
xticks([1:regionnum]);xticklabels(HPC_layer_name(1:regionnum));yline(1,'k--');title(num2str(comp',2));xlabel(num2str(stats',2));
ylabel('Theta / Delta Power Ratio');legend('Maze: run','','Maze: pause','','Slow wave sleep','');xlim([0.5 regionnum+.5]);
if ifsavedata
figdata=[squeeze(thedelv(3,1:regionnum,2,:))';squeeze(thedelv(3,1:regionnum,3,:))';squeeze(thedel_sleep(3,1:regionnum,3,:))'];
rowname=[];
for i=1:3
    for in=1:16
        if i==1
            rowname{in+(i-1)*16}=['Run_',num2str(in)];
        elseif i==2
            rowname{in+(i-1)*16}=['Pause_',num2str(in)];
        elseif i==3
            rowname{in+(i-1)*16}=['SWS_',num2str(in)];
        end
    end
end
T = array2table(figdata,'VariableNames',{'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3'},'RowNames',rowname);
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet','Figure 1h','WriteRowNames',true);
end

indf=f<=50;letternum={'e','f','h'};fname=[];
for i=1:sum(indf)
    fname{i}=num2str(f(i));
end
% run & pause spectrum 
% below quantify theta peaks
thetapeaks=[];layer=4;thetapks=[];
for in=1:16
    for v=1:2
        [pks,locs] = findpeaks(squeeze(specv(:,layer,v+1,in)));
        ind=f(locs)>4 & f(locs)<10;
        thetapks(in,v,:)=[sum(ind),f(locs(ind))];
        thetapeaks=[thetapeaks;[pks(ind),locs(ind),v*ones(sum(ind),1)]];
    end
end
thetamean=mean(thetapks(:,:,2),1);
thetastd=std(thetapks(:,:,2),1);
for v=1:2
    figure(fign);subplot(pr,pc,v);figdata=[];rowname=[];
    for layer=1:7
        hold on;shaded_errbar(f,squeeze(specv(:,8-layer,v+1,:)),defaultcolor7(layer,:));
        figdata=[figdata;squeeze(specv(indf,layer,v+1,:))'];
        for in=1:16
            rowname{in+(layer-1)*16}=[HPC_layer_name{layer},'_',num2str(in)];
        end
    end
    xlim([1 50]);xlabel('Frequency (Hz)');ylabel('Power (dB)');title([vname{v},': ',num2str(thetamean(v),2),' +- ',num2str(thetastd(v),2)]);
    legend('CA3','','DG GCL','','DG MML','','DG OML','','CA1 slm','','CA1 st rad','','CA1 pyr','');
    if ifsavedata
        T = array2table(figdata,'VariableNames',fname,'RowNames',rowname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 1',letternum{v}],'WriteRowNames',true);
    end
end

% slow wave sleep
figure(fign);subplot(pr,pc,3);i=1;figdata=[];
for layer=1:7
    hold on;shaded_errbar(f,squeeze(spec_sleep(:,8-layer,4-i,:)),defaultcolor7(layer,:));
    figdata=[figdata;squeeze(spec_sleep(indf,layer,4-i,:))'];
end
n=sum(~isnan(squeeze(nanmean(spec_sleep(:,1,4-i,:),1))));
xlim([1 50]);xlabel('Frequency (Hz)');ylabel('Power (dB)');title(['slow wave sleep, N = ',num2str(n)]);
legend('CA3','','DG GCL','','DG MML','','DG OML','','CA1 slm','','CA1 st rad','','CA1 pyr','');
if ifsavedata
    T = array2table(figdata,'VariableNames',fname,'RowNames',rowname);
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Sup Figure 3',letternum{i}],'WriteRowNames',true);
end
figure(fign);figure_title='Figure1_Jan2025';save_current_figure(figure_title);
%--------------------------------------------------------------------------------------------
if 0
% 3. sup figure 4
cd(homedir);load('SlowWaveSleep_Example',"timeset");
ins=find(timeset(:,1)>0);
for i=1:length(ins)
    in2=ins(i);
    generate_SupFigure_4_sws_lfp_example(1,in2);
end
%----------------------------------------------------------------------------
% below for sup vel theta/delta ratio
cd(homedir);
load('theta_velocity_inoutbound_inter_reward','theta_inter_mean',"theta_reward_mean","theta_in_mean","theta_out_mean", ...
    'vel_inter_mean',"vel_reward_mean","vel_in_mean","vel_out_mean");
powertypename={'Theta Power (z)','Theta / Delta Power Ratio'};
powerid=2;postype=1;
% below plot for figure1
figure(2);pr=3;pc=3;
subplot(pr,pc,1);% below plot velocity profile
for i=1:3
    a1=squeeze(vel_out_mean(i,:,postype,:));
    a2=squeeze(vel_in_mean(i,:,postype,:));
    a3=squeeze(vel_reward_mean(i,:,:));
    a4=squeeze(vel_inter_mean(i,:,:));
    hold on;shaded_errbar([0.1:0.1:1]-0.05+1,a1,plotcolors3(i,:));
    hold on;shaded_errbar(-[0.1:0.1:1]+4.05,a2,plotcolors3(i,:));
    hold on;shaded_errbar([0.1:0.1:1]-0.05,a4,plotcolors3(i,:));
    hold on;shaded_errbar(2+[0.1:0.1:1]-0.05,a3,plotcolors3(i,:));
end
xlabel('Position or Time Index : Choice Point | Out | Reward | In');
ylabel('Velocity (cm/s)');legend('Forced visit','','Correct choice visit','','Error choice visit','');
FIG_INDEX=['fig1_1'];save_fig(FIG_INDEX,ifsavefig);
figure(2);subplot(pr,pc,2);
velout=squeeze(nanmean(vel_out_mean(:,:,postype,:),2));
velin=squeeze(nanmean(vel_in_mean(:,:,postype,:),2));
velrew=squeeze(nanmean(vel_reward_mean,2));
velcen=squeeze(nanmean(vel_inter_mean,2));
velmean=[nanmean(velcen,2),nanmean(velout,2),nanmean(velrew,2),nanmean(velin,2)]';
b=bar(velmean);
for i=1:2
    b(i).FaceColor = plotcolors3(i,:);
end
b(3).FaceColor =0.6*[1 1 1];
for in=1:16
    hold on;plot([0.78 1 1.22],velcen(:,in),'k.');
    hold on;plot([0.78 1 1.22]+1,velout(:,in),'k.');
    hold on;plot([0.78 1 1.22]+2,velrew(:,in),'k.');
    hold on;plot([0.78 1 1.22]+3,velin(:,in),'k.');
end
xticks([1:4]);xticklabels({'Reward','Center','Outbound','Inbound'});title('Velocity');%yline(4,'k--');
ylabel('Mean Velocity (cm/s)');legend('Forced visit','Correct choice visit','Error choice visit');
FIG_INDEX=['fig1_2'];save_fig(FIG_INDEX,ifsavefig);

figure(2);subplot(pr,pc,3);
for h=4;
    for i=1:3
        a1=squeeze(theta_out_mean(i,:,h,postype,powerid,:));
        a2=squeeze(theta_in_mean(i,:,h,postype,powerid,:));
        a3=squeeze(theta_reward_mean(i,:,h,powerid,:));
        a4=squeeze(theta_inter_mean(i,:,h,powerid,:));
        hold on;shaded_errbar([0.1:0.1:1]-0.05+1,a1,plotcolors3(i,:));
        hold on;shaded_errbar(-[0.1:0.1:1]+4.05,a2,plotcolors3(i,:));
        hold on;shaded_errbar([0.1:0.1:1]-0.05,a4,plotcolors3(i,:));
        hold on;shaded_errbar(2+[0.1:0.1:1]-0.05,a3,plotcolors3(i,:));
    end
end
xlabel('Position or Time Index : Choice Point | Out | Reward | In');
ylabel(powertypename{powerid});title([HPC_layer_name{h}]);
legend('Forced visit','','Correct choice visit','','Error choice visit','');
ylim([0 inf]);%yline(1,'k--');
FIG_INDEX=['fig1_3'];save_fig(FIG_INDEX,ifsavefig);

figure(2);subplot(pr,pc,4);
thetaout=squeeze(nanmean(theta_out_mean(:,:,h,postype,powerid,:),2));
thetain=squeeze(nanmean(theta_in_mean(:,:,h,postype,powerid,:),2));
thetarew=squeeze(nanmean(theta_reward_mean(:,:,h,powerid,:),2));
thetacen=squeeze(nanmean(theta_inter_mean(:,:,h,powerid,:),2));
thetamean=[nanmean(thetacen,2),nanmean(thetaout,2),nanmean(thetarew,2),nanmean(thetain,2)]';
b=bar(thetamean);
for i=1:2
    b(i).FaceColor = plotcolors3(i,:);
end
b(3).FaceColor =0.6*[1 1 1];
for in=1:16
    hold on;plot([0.78 1 1.22],thetacen(:,in),'k.');
    hold on;plot([0.78 1 1.22]+1,thetaout(:,in),'k.');
    hold on;plot([0.78 1 1.22]+2,thetarew(:,in),'k.');
    hold on;plot([0.78 1 1.22]+3,thetain(:,in),'k.');
end
xticks([1:4]);xticklabels({'Reward','Center','Outbound','Inbound'});
ylabel(powertypename{powerid});title([HPC_layer_name{h}]);legend('Forced visit','Correct choice visit','Error choice visit');
FIG_INDEX=['fig1_4'];save_fig(FIG_INDEX,ifsavefig);

vname={'Run','Pause','Immobile'};
figure(2);subplot(pr,pc,5);
mu=nanmean(velmat,1);bg=[ones(16,1),2*ones(16,1),3*ones(16,1),4*ones(16,1)];
swarmchart(bg(:),velmat(:),'k.');xlim([0 5]);ylim([0 60]);
for b=1:4
    hold on;plot([0.7 1.3]+(b-1),mu(b)*[1 1],'r');
end
xticks([1:4]);xticklabels({'Reward','Center','Outbound','Inbound'});
ylabel('Mean Velocity (cm/s)');legend('One session','Mean');
FIG_INDEX=['fig1_5'];save_fig(FIG_INDEX,ifsavefig);

figure(2);subplot(pr,pc,6);regionnum=7;
hold on;shaded_errbar([1:regionnum],squeeze(thedelv(3,1:regionnum,2,:)),plotcolors3(1,:));
hold on;shaded_errbar([1:regionnum],squeeze(thedelv(3,1:regionnum,3,:)),plotcolors3(2,:));
hold on;shaded_errbar([1:regionnum],squeeze(thedel_sleep(3,1:regionnum,3,:)),plotcolors3c(3,:));
xticks([1:regionnum]);xticklabels(HPC_layer_name(1:regionnum));yline(1,'k--');
ylabel('Theta / Delta Power Ratio');legend('Maze: run','','Maze: pause','','Slow wave sleep','');xlim([0.5 regionnum+.5]);
FIG_INDEX=['fig1_6'];save_fig(FIG_INDEX,ifsavefig);

for v=1:2
    figure(2);subplot(pr,pc,v+6);
    for layer=1:7
        hold on;shaded_errbar(f,squeeze(specv(:,8-layer,v+1,:)),defaultcolor7(layer,:));
    end
    xlim([1 50]);xlabel('Frequency (Hz)');ylabel('Power (dB)');title(vname{v});
    %legend('dCA1 pyr: Run','dCA1 pyr: Pause','dCA1 st rad: Run','dCA1 st rad: Pause','dCA1 slm: Run','dCA1 slm: Pause','DG OML: Run','DG OML: Pause','DG MML: Run','DG MML: Pause','DG GCL: Run','DG GCL: Pause','dCA3: Run','dCA3: Pause');
    FIG_INDEX=['fig1_',num2str(v+6)];save_fig(FIG_INDEX,ifsavefig);
end

figure(2);subplot(pr,pc,9);
for layer=1:7
    hold on;shaded_errbar(f,squeeze(spec_sleep(:,8-layer,3,:)),defaultcolor7(layer,:));
end
xlim([1 50]);xlabel('Frequency (Hz)');ylabel('Power (dB)');title('Slow Wave Sleep');
legend('dCA3','','DG GCL','','DG MML','','DG OML','','dCA1 slm','','dCA1 st rad','','dCA1 pyr','');
FIG_INDEX=['fig1_9'];save_fig(FIG_INDEX,ifsavefig);
cd(figdir);figure_title='fig1_thetapower_run_pause';save_current_figure(figure_title);
end
if 0
% below for sup vel theta/delta ratio
cd(homedir);
load('theta_velocity_inoutbound_inter_reward','theta_inter_mean',"theta_reward_mean","theta_in_mean","theta_out_mean", ...
    'vel_inter_mean',"vel_reward_mean","vel_in_mean","vel_out_mean");
figdir='X:\Mengni\Data_Analysis\Paper_Figures';
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
powertypename={'Theta Power (z)','Theta / Delta Power Ratio'};
plotcolors=[0 0.4470 0.7410];
plotcolors(2,:)=[0.8500 0.3250 0.0980];
plotcolors(3,:)=0;
powerid=2;
for postype=1
figure(3);
for i=1:2
subplot(4,2,1+(i-1)*7);
plot(0,0,'-','Color',plotcolors(1,:));hold on;
plot(0,0,'-','Color',plotcolors(2,:));hold on;plot(0,0,'k-');
end
for h=1:7
    figure(3);subplot(4,2,h+1);
    for i=1:3
        a1=squeeze(theta_out_mean(i,:,h,postype,powerid,:));
        a2=squeeze(theta_in_mean(i,:,h,postype,powerid,:));
        a3=squeeze(theta_reward_mean(i,:,h,powerid,:));
        a4=squeeze(theta_inter_mean(i,:,h,powerid,:));
        hold on;shaded_errbar([0.1:0.1:1]-0.05+1,a1,plotcolors(i,:));
        hold on;shaded_errbar(-[0.1:0.1:1]+4.05,a2,plotcolors(i,:));
        hold on;shaded_errbar([0.1:0.1:1]-0.05,a4,plotcolors(i,:));
        hold on;shaded_errbar(2+[0.1:0.1:1]-0.05,a3,plotcolors(i,:));
    end
    yline(1,'k--');
    if powerid==2
        ylim([0 20]);
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
    xlabel('Position or Time Index : Choice Point | Out | Reward | In');
    ylabel(powertypename{powerid});title([HPC_layer_name{h} '; P = ',num2str([p4,p1,p3,p2],2)]);
    FIG_INDEX=['fig1sup2_',num2str(h+1)];save_fig(FIG_INDEX,ifsavefig);
end
legend('Forced Visit','Free Visit','Error Visit');
% below plot velocity profile
figure(3);subplot(4,2,1);
for i=1:3
    a1=squeeze(vel_out_mean(i,:,postype,:));
    a2=squeeze(vel_in_mean(i,:,postype,:));
    a3=squeeze(vel_reward_mean(i,:,:));
    a4=squeeze(vel_inter_mean(i,:,:));
    hold on;shaded_errbar([0.1:0.1:1]-0.05+1,a1,plotcolors(i,:));
    hold on;shaded_errbar(-[0.1:0.1:1]+4.05,a2,plotcolors(i,:));
    hold on;shaded_errbar([0.1:0.1:1]-0.05,a4,plotcolors(i,:));
    hold on;shaded_errbar(2+[0.1:0.1:1]-0.05,a3,plotcolors(i,:));
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
xlabel('Position or Time Index : Choice Point | Out | Reward | In');
ylabel('Velocity (cm/s)');title(['Velocity : P = ',num2str([p4,p1,p3,p2],2)]);
legend('Forced Visit','Free Visit','Error Visit');
FIG_INDEX=['fig1sup2_1'];save_fig(FIG_INDEX,ifsavefig);
figure_title=['ThetaDelta_Velocity_InOutInterReward_',num2str(postype)];save_current_figure(figure_title);
end
end

if 0
typenum=16;typeindex=[typenum:-1:1]*floor(256/typenum);
plotcolors=colormap(jet);
plotcolors=plotcolors(typeindex,:);
figure;
for layer=1:7
    subplot(2,4,layer);
    for in=1:16
        hold on;plot(log10(f),nanmean(specv(:,layer,1,in),4),'color',plotcolors(in,:));
        hold on;plot(log10(f),nanmean(specv(:,layer,2,in),4),'--','color',plotcolors(in,:));
    end
end

figure;
subplot(2,2,1);regionnum=7;
err1=nanstd(thedelv(3,1:regionnum,1,:),0,4)/4;
err2=nanstd(thedelv(3,1:regionnum,2,:),0,4)/4;
hold on;errorbar([1:regionnum],nanmean(thedelv(3,1:regionnum,2,:),4),err2,'bo-','LineWidth',2);
hold on;errorbar([1:regionnum],nanmean(thedelv(3,1:regionnum,1,:),4),err1,'ro-','LineWidth',2);
xticks([1:regionnum]);xticklabels(HPC_layer_name(1:regionnum));
ylabel('Theta / Delta Power Ratio');legend('Run','Pause');xlim([0.5 regionnum+.5]);
subplot(2,2,2);
for layer=1:7
    hold on;plot(log10(f),nanmean(specv(:,layer,2,:),4),'--','color',plotcolors(layer,:));
    hold on;plot(log10(f),nanmean(specv(:,layer,1,:),4),'color',plotcolors(layer,:));
end
xticks(log10([1,4,8,16,50,100,250]));xticklabels({'1','4','8','16','50','100','250'});
xlabel('Frequency');ylabel('Power (dB)');
legend('dCA1 pyr: Run','dCA1 pyr: Pause','dCA1 st rad: Run','dCA1 st rad: Pause','dCA1 slm: Run','dCA1 slm: Pause','DG OML: Run','DG OML: Pause','DG MML: Run','DG MML: Pause','DG GCL: Run','DG GCL: Pause','dCA3: Run','dCA3: Pause');
subplot(2,2,3);
shaded_errbar(log10(f)',squeeze(specv(:,1,2,:)),plotcolors(1,:));
hold on;shaded_errbar(log10(f)',squeeze(specv(:,1,1,:)),plotcolors(3,:));
hold on;shaded_errbar(log10(f)',squeeze(specv(:,4,2,:)),plotcolors(5,:));
hold on;shaded_errbar(log10(f)',squeeze(specv(:,4,1,:)),plotcolors(7,:));
xticks(log10([1,4,8,16,50,100,250]));xticklabels({'1','4','8','16','50','100','250'});
xlabel('Frequency');ylabel('Power (dB)');
legend('dCA1 pyr: Run','','dCA1 pyr: Pause','','DG OML: Run','','DG OML: Pause','');
cd(homedir);figure_title='fig1_thetapower_run_pause';save_current_figure(figure_title);
end