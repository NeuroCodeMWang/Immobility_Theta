function figure_1_spectrum_HPC_layer_trial_sleep(in)

homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
% 1. below calculate theta & delta filtered LFP in dca1 channels
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('HPC_Layer_Channel_LFP');
    [thetapower,thetaphase,thetafilted]=thetafilt_npx(LFP_HPC,LFP_Frequency);
    [deltapower,deltaphase,deltafilted]=deltafilt_npx(LFP_HPC,LFP_Frequency);
    cd(savedir);save('theta_delta_LFP_HPC','lfptime','thetapower','thetaphase','thetafilted',...
        "deltafilted",'deltaphase','deltapower','LFP_Frequency','-v7.3');
    disp(['Session Completed: ',num2str(in)]);
end

homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
p2range=[250 420];p3range=[100 250];segave=0;LFP_Frequency=625;
params.Fs=LFP_Frequency;params.fpass=[1,300];win=2;
blankpos=[];winE=[];
figure;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    load('Position_Data_Raw.mat');
    load('HPC_Layer_Channel_LFP');
    load('theta_delta_LFP_HPC','thetapower','deltapower');
    
    ind=Position_Data_Raw(:,4)<0;
    position=Position_Data_Raw(ind,:);
    if (in~=13) & (in~=14)
        p2range1=p2range;p3range1=p3range;
    else
        p2range1=p2range+150;p3range1=p3range+50;
    end
    %subplot(4,5,in);plot(position(:,2),position(:,3),'.');
    %xline(p2range1,'r');yline(p3range1,'r');

    ind= (position(:,2)<=p2range1(2) & position(:,2)>=p2range1(1)) &...
        (position(:,3)<=p3range1(2) & position(:,3)>=p3range1(1)) ;
    position(~ind,2:3)=0;
    blankpos(in,:)=mean(position(:,2:3)==0,1);
    
    position1=position(position(:,4)==-1,:);
    position2=position(position(:,4)==-2,:);
    [Velocity1]=CALCULATE_VELOCITY(position1);position1(:,5)=Velocity1;
    [Velocity2]=CALCULATE_VELOCITY(position2);position2(:,5)=Velocity2;

    LFP_Frequency=625;
    inds1=lfptime(:,2)==-1;
    inds2=lfptime(:,2)==-2;
    indv=lfptime(:,2)>0;
    specset=[];Etimeset=[];theta_delta_power_set=[];
    for b=1:3
        if b==1
            data=LFP_HPC(inds1,:);datatime=lfptime(inds1,1);
            thetapower_=thetapower(inds1,:);deltapower_=deltapower(inds1,:);
            position=position1;
        elseif b==2
            data=LFP_HPC(inds2,:);datatime=lfptime(inds2,1);
            thetapower_=thetapower(inds2,:);deltapower_=deltapower(inds2,:);
            position=position2;
        elseif b==3
            data=LFP_HPC(indv,:);datatime=lfptime(indv,1);
            thetapower_=thetapower(indv,:);deltapower_=deltapower(indv,:);
            position=Position_Data;
        end
        N=size(datatime,1); % length of segmented data
        Eind1=1:LFP_Frequency*win:N-LFP_Frequency*win;
        Eind2=Eind1+LFP_Frequency*win;
        Etime=[datatime(Eind1),datatime(Eind2)];
        Etime(:,3)=Etime(:,2)-Etime(:,1);
        Etime(:,5:6)=[Eind1(:),Eind2(:)];
        if b==3
            subplot(4,4,in);histogram(Etime(:,3));title(num2str(mean(Etime(:,3)==win),2));
        end
        indE=Etime(:,3)==win;
        Etime=Etime(indE,:);
        winE(in,:,b)=[size(Etime,1),mean(indE)];

        % below to calculate spectrum
        spec=[];
        for c=1:11
            [S,f]=mtspectrumsegc(data(:,c),win,params,segave);
            if size(S,2)~=length(indE)
                disp(['Error S size : ',num2str(in)]);
            else
                S=10*log10(S(:,indE));
                spec(:,:,c)=S;
            end
        end
        specset{b}=spec;

        % below to measure rat's movement & theta, delta power
        theta_delta_power=[];
        for i=1:size(Etime,1)
            ind= position(:,1)>=Etime(i,1) & position(:,1)<=Etime(i,2) ;
            if b==3
                Etime(i,4)=mean(position(ind,5));
            else
                Etime(i,4)=mean(position(ind,5)==0 & position(ind,2)>0); % still percent
            end

            indlfp=[Etime(i,5):Etime(i,6)];
            theta_delta_power(i,:,1)=nanmean(thetapower_(indlfp,:),1);
            theta_delta_power(i,:,2)=nanmean(deltapower_(indlfp,:),1);
            theta_delta_power(i,:,3)=theta_delta_power(i,:,1)./theta_delta_power(i,:,2);
        end
        Etimeset{b}=Etime;
        theta_delta_power_set{b}=theta_delta_power;
    end
    cd(savedir);save('spectrum_HPC_layer_trial_sleep_2s','specset','Etimeset','theta_delta_power_set','winE','position1',"position2",'blankpos','f','win','-v7.3');
    disp(['Session Completed: ',num2str(in)]);
    clear thetapower deltapower LFP_HPC
end
%----------------------------------------------------------------------------------
% below plot for figure1
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
b=3;specv=[];thedelv=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('spectrum_HPC_layer_trial_sleep','Etimeset','specset','theta_delta_power_set','f');
    Etime=Etimeset{b};
    spec=specset{b};
    theta_delta_power=theta_delta_power_set{b};
    indstill=Etime(:,4)<10;
    indrun=Etime(:,4)>=10;
    specv(:,:,1,in)=squeeze(nanmean(spec(:,indstill,:),2));
    specv(:,:,2,in)=squeeze(nanmean(spec(:,indrun,:),2));
    thedelv(:,:,1,in)=squeeze(nanmean(theta_delta_power(indstill,:,:),1))';
    thedelv(:,:,2,in)=squeeze(nanmean(theta_delta_power(indrun,:,:),1))';
    %subplot(4,4,in);histogram(Etime(:,4));title([num2str([mean(Etime(:,4)<5),mean(Etime(:,4)>10)])]);
end
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
typenum=7;typeindex=[typenum:-1:1]*floor(256/typenum);
plotcolors=colormap(jet);
plotcolors=plotcolors(typeindex,:);
figure;
subplot(2,2,1);regionnum=11;
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
%------------------------------------------------------------------------------------------------------
% below plot for not figure related 
spec_mean=[];thetadelta_mean=[];binum=10;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('spectrum_HPC_layer_trial_sleep');
    for b=1:3
        spec=specset{b};
        Etime=Etimeset{b};
        theta_delta_power=theta_delta_power_set{b};
        if b==3
            [~,I]=sort(Etime(:,4));
            Etime(I,4)=[1:length(I)]/length(I);
        else
            Etime(:,4)=1-Etime(:,4);
        end
        for i=1:binum
            ind=ceil(Etime(:,4)*binum)==i;
            spec_mean(:,i,:,b,in)=squeeze(nanmean(spec(:,ind,:),2));
            thetadelta_mean(:,i,b,in)=squeeze(nanmean(theta_delta_power(ind,:,3),1));
        end
    end
end
behavename={'Sleepbox1','Sleepbox2','Trials'};
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG IML/GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
figure;
for b=1:3
    as=[];
    for in=1:16
        a=squeeze(thetadelta_mean(:,:,b,in));
        as=[as,a];
    end
    subplot(3,1,b);
    imagesc([0 16],[0 11],as);set(gca,'YDir','normal');
    xline(1:15,'k');title(behavename{b});xlabel('Activity: Low - High : Session ID');
    yticks([.5:10.5]);yticklabels(HPC_layer_name);colorbar;
end
cd(homedir);figure_title='Theta_Delta_Ratio_HPC_Layer_Trial_Sleepbox';save_current_figure(figure_title);

figure;
for b=1:3
    as=[];
    for r=1:11
        a=squeeze(nanmean(spec_mean(:,:,r,b,:),5));
        as=[as,a];
    end
    subplot(3,1,b);
    imagesc([0 11],[f],as);set(gca,'YDir','normal');
    xline(1:10,'k');title(behavename{b});xlabel('Activity: Low - High');
    ylabel('LFP Frequency (Hz)');%ylim([1 40]);
    xticks([.5:10.5]);xticklabels(HPC_layer_name);
end
cd(homedir);figure_title='Spectrum_HPC_Layer_Trial_Sleepbox_2';save_current_figure(figure_title);
