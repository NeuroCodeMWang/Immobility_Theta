% 1. detect sleep periods and sleep stages
% 2. find good slow wave sleep examples
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
%figdir='C:\Users\Mengni\Documents\Paper\Figure_Nov2024';
if 0
% 1. below detect sleep periods and sleep stages
stillratiocut=0.95;min_sleep_dur=50;min_rem_dur=10;min_nonrem_dur=0;
remall=[];segall=[];nonremall=[];s1all=[];s2all=[];spec_sleep=[];thedel_sleep=[];figure;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('spectrum_HPC_layer_trial_sleep_2s','Etimeset','specset','f','theta_delta_power_set');
    spec1=specset{1};
    spec2=specset{2};
    s1=Etimeset{1};
    s2=Etimeset{2};
    theta_delta1=theta_delta_power_set{1};
    theta_delta2=theta_delta_power_set{2};
    s1(:,7:13)=squeeze(theta_delta1(:,1:7,3));
    s2(:,7:13)=squeeze(theta_delta2(:,1:7,3));
    s1(:,14)=[1:size(s1,1)];
    s2(:,14)=[1:size(s2,1)]; 
    s1(:,15)=0;
    s2(:,15)=0;
    blanksegment1=get_segment(s1(:,4)>=stillratiocut);
    blanksegment2=get_segment(s2(:,4)>=stillratiocut);
    blanksegment1(:,4)=s1(blanksegment1(:,1),1);
    blanksegment1(:,5)=s1(blanksegment1(:,2),2);
    blanksegment1(:,6)=blanksegment1(:,5)-blanksegment1(:,4);
    blanksegment1(:,7)=1;
    blanksegment1(:,8)=in;
    sleep1=blanksegment1(blanksegment1(:,6)>=min_sleep_dur,:);
    rem1s=[];nonrem1s=[];
    for i=1:size(sleep1,1)
        ss1=s1(sleep1(i,1):sleep1(i,2),:);
        rem1=get_segment(ss1(:,7)>1);
        if ~isempty(rem1)
            rem1(:,4)=ss1(rem1(:,1),1);
            rem1(:,5)=ss1(rem1(:,2),2);
            rem1(:,6)=rem1(:,5)-rem1(:,4);
            ind=rem1(:,6)>=min_rem_dur;
            if sum(ind)>0
                rem1=rem1(ind,:);
                for j=1:size(rem1,1)
                    rem1(j,7)=nanmean(ss1(rem1(j,1):rem1(j,2),7));
                    indj=ss1(rem1(j,1),14):ss1(rem1(j,2),14);
                    s1(indj,15)=1; % rem sleep
                end
                rem1(:,8)=in;
                rem1s=[rem1s;rem1];
            end
        end
        nonrem1=get_segment(ss1(:,7)<1);
        if ~isempty(nonrem1)
            nonrem1(:,4)=ss1(nonrem1(:,1),1);
            nonrem1(:,5)=ss1(nonrem1(:,2),2);
            nonrem1(:,6)=nonrem1(:,5)-nonrem1(:,4);
            ind=nonrem1(:,6)>=min_nonrem_dur;
            if sum(ind)>0
                nonrem1=nonrem1(ind,:);
                for j=1:size(nonrem1,1)
                    nonrem1(j,7)=nanmean(ss1(nonrem1(j,1):nonrem1(j,2),7));
                    indj=ss1(nonrem1(j,1),14):ss1(nonrem1(j,2),14);
                    s1(indj,15)=2; % slow wave sleep
                end
                nonrem1(:,8)=in;
                nonrem1s=[nonrem1s;nonrem1];
            end
        end
    end

    blanksegment2(:,4)=s2(blanksegment2(:,1),1);% start time
    blanksegment2(:,5)=s2(blanksegment2(:,2),2);% end time
    blanksegment2(:,6)=blanksegment2(:,5)-blanksegment2(:,4);% duration
    blanksegment2(:,7)=2;
    blanksegment2(:,8)=in;
    sleep2=blanksegment2(blanksegment2(:,6)>=min_sleep_dur,:);
    rem2s=[];nonrem2s=[];
    for i=1:size(sleep2,1)
        ss2=s2(sleep2(i,1):sleep2(i,2),:);
        rem2=get_segment(ss2(:,7)>1);
        if ~isempty(rem2)
            rem2(:,4)=ss2(rem2(:,1),1);
            rem2(:,5)=ss2(rem2(:,2),2);
            rem2(:,6)=rem2(:,5)-rem2(:,4);
            ind=rem2(:,6)>=min_rem_dur;
            if sum(ind)>0
                rem2=rem2(ind,:);
                for j=1:size(rem2,1)
                    rem2(j,7)=nanmean(ss2(rem2(j,1):rem2(j,2),7));
                    indj=ss2(rem2(j,1),14):ss2(rem2(j,2),14);
                    s2(indj,15)=1; % rem sleep
                end
                rem2(:,8)=in;
                rem2s=[rem2s;rem2];
            end
        end
        nonrem2=get_segment(ss2(:,7)<1);
        if ~isempty(nonrem2)
            nonrem2(:,4)=ss2(nonrem2(:,1),1);
            nonrem2(:,5)=ss2(nonrem2(:,2),2);
            nonrem2(:,6)=nonrem2(:,5)-nonrem2(:,4);
            ind=nonrem2(:,6)>=min_nonrem_dur;
            if sum(ind)>0
                nonrem2=nonrem2(ind,:);
                for j=1:size(nonrem2,1)
                    nonrem2(j,7)=nanmean(ss2(nonrem2(j,1):nonrem2(j,2),7));
                    indj=ss2(nonrem2(j,1),14):ss2(nonrem2(j,2),14);
                    s2(indj,15)=2; % slow wave sleep
                end
                nonrem2(:,8)=in;
                nonrem2s=[nonrem2s;nonrem2];
            end
        end
    end
    s1(:,16)=in;
    s1all=[s1all;s1];
    s2(:,16)=in;
    s2all=[s2all;s2];
    remall=[remall;rem1s;rem2s];
    nonremall=[nonremall;nonrem1s;nonrem2s];
    segall=[segall;blanksegment1;blanksegment2];
    subplot(4,4,in);histogram(blanksegment1(:,6));
    hold on;histogram(blanksegment2(:,6));

    spec=cat(2,spec1,spec2);
    s12=[s1;s2];
    theta_delta_power=cat(1,theta_delta1,theta_delta2);
    for i=1:3
        ind=s12(:,15)==i-1; % awake, rem, slow wave sleep
        spec_sleep(:,:,i,in)=squeeze(nanmean(spec(:,ind,:),2));
        thedel_sleep(:,:,i,in)=squeeze(nanmean(theta_delta_power(ind,:,:),1))';
    end
end
sleepsession=[];
for in=1:16
    rem=remall(remall(:,8)==in,:);
    nonrem=nonremall(nonremall(:,8)==in,:);
    sleepsession(in,:)=[sum(rem(:,6)),size(rem,1),sum(nonrem(:,6)),size(nonrem,1)];
end
cd(homedir);save('Rem_SW_sleep_analysis_2s','segall','remall','stillratiocut','s1all','s2all','f',...
    'nonremall','min_sleep_dur','min_rem_dur','min_nonrem_dur','sleepsession','spec_sleep','thedel_sleep','-v7.3');
end
cd(homedir);load('Rem_SW_sleep_analysis_2s');
savedir=SessionSet16{1};cd(savedir);load('spectrum_HPC_layer_trial_sleep_2s','f');
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
defaultcolor4=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560]];
defaultcolor7=[[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];[0.4940 0.1840 0.5560];[0 0 0];[0 0.4470 0.7410];[0.9290 0.6940 0.1250];[1 1 1]*0.6];
sleepname={'Awake Sleepbox','REM','Slow Wave'};
figure;
subplot(2,2,1);regionnum=7;
for i=1:3
    hold on;shaded_errbar([1:regionnum],squeeze(thedel_sleep(3,1:regionnum,i,:)),defaultcolor4(i,:));
end
xticks([1:regionnum]);xticklabels(HPC_layer_name(1:regionnum));
ylabel('Theta / Delta Power Ratio');legend('Awake Sleepbox','','REM','','Slow Wave','');xlim([0.5 regionnum+.5]);
for i=1:3
    subplot(2,2,i+1);
    for layer=1:7
        hold on;shaded_errbar(f,squeeze(spec_sleep(:,8-layer,i,:)),defaultcolor7(8-layer,:));
    end
    xlim([1 50]);xlabel('Frequency');ylabel('Power (dB)');title(sleepname{i});
    legend('dCA3','','DG GCL','','DG MML','','DG OML','','dCA1 slm','','dCA1 st rad','','dCA1 pyr','');
end
figure_title='fig1_thetapower_sleep';save_current_figure(figure_title);
%------------------------------------------------------------------------------------
% 2. below identify good slow wave sleep examples for individual sessions
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
load('Rem_SW_sleep_analysis','remall','nonremall');close all;
factor=90;
timeset=zeros(16,2);timeset(in,:)=[950.3 955.3];
for in=1:16
    rem=remall(remall(:,8)==in,:);
    nonrem=nonremall(nonremall(:,8)==in,:);
    %if (~isempty(nonrem)) %| (~isempty(rem))
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_delta_LFP_HPC.mat', 'thetaphase', 'lfptime', 'thetapower');
    load('HPC_Layer_Channel_LFP.mat', 'LFP_HPC','LFP_Frequency');
    load('CSD_LFP_Visualization.mat', 'CSDall', 'dHPC_layer7channel3', 'channels');

    dHPC_layer7channel3(:,4)=0;
    for i=1:7
        dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2));
    end

    % 0. define reftheta as from OML DG & identify individual theta cycles
    filted=zscore(LFPfilt(LFP_HPC(:,4),LFP_Frequency));
    phase=thetaphase(:,4)-180;phase=phase+(phase<0)*360; % make thetapeak as 180; trough as 0
    reftheta=[lfptime(:,1),filted,zscore(thetapower(:,4)),phase];
    indsleep=lfptime(:,2)<0;
    reftheta=reftheta(indsleep,:);
    lfpz_sleep=zscore(LFP_HPC(indsleep,:),0,1);
    lfptime_sleep=lfptime(indsleep,:);
    CSDall_sleep=CSDall(:,indsleep);
    lfp7=LFP_HPC(indsleep,1:7);
    reftheta(:,5)=CSDall_sleep(dHPC_layer7channel3(5,4),:)';
    reftheta(:,6)=CSDall_sleep(dHPC_layer7channel3(5,4)+1,:)';
    if 0
    thetacycle=identify_thetacycles(reftheta);
    thetacycle(:,17:20)=thetacycle(:,10:13);
    thetacycle(:,13)=[1:size(thetacycle,1)];
    %clear LFP_HPC thetapower thetaphase 
    lfptime_sleep(:,6)=0;
    for i=1:size(nonrem,1)
        indlfp= lfptime_sleep(:,1)>=nonrem(i,4) & lfptime_sleep(:,1)<=nonrem(i,5);
        lfptime_sleep(indlfp,6)=2;
    end
    for i=1:size(rem,1)
        indlfp= lfptime_sleep(:,1)>=rem(i,4) & lfptime_sleep(:,1)<=rem(i,5);
        lfptime_sleep(indlfp,6)=1;
    end
    lfptime_sleep(:,5)=0; % associated thetacycle id
    thetacycle(:,21)=0;
    for i=1:size(thetacycle,1)
        ind=thetacycle(i,1):thetacycle(i,2);
        lfptime_sleep(ind,5)=i;
        thetacycle(i,21)=nanmean(lfptime_sleep(ind,6));
    end
    end
    figure;
    for i=1:size(nonrem,1)
        indlfp= lfptime_sleep(:,1)>=nonrem(i,4) & lfptime_sleep(:,1)<=nonrem(i,5);
        subplot(2,1,1);        
        for j=1:7
            if j==4
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'r');
            else
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
            end
        end
        xline(nonrem(i,4:5),'k');xlim([min(lfptime_sleep(indlfp,1)),max(lfptime_sleep(indlfp,1))]);
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        
        subplot(2,1,2);
        imagesc(lfptime_sleep(indlfp,1),[channels(1,3) channels(end,3)],CSDall_sleep(:,indlfp));set(gca,'YDir','normal');

        xline(nonrem(i,4:5),'k');xlabel('Time in Slow Wave Sleep');xlim([min(lfptime_sleep(indlfp,1)),max(lfptime_sleep(indlfp,1))]);
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        pause;clf;
    end

    timeset(in,:)=[513.5 518.5];%[15360.3 15365.3];%[830 836];%[1910 1916];%[12945 12951];%[1897.8 1904.2];%[1666.5 1672.5];%[1673.5 1678.5];%[1092 1098];%[950.3 955.3];
    figure;
    for i=1%:size(timeset,1)
        indlfp= lfptime_sleep(:,1)>=timeset(in,1) & lfptime_sleep(:,1)<=timeset(in,2);
        subplot(2,1,1);        
        for j=1:7
            if j==1
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor*2+dHPC_layer7channel3(j,2),'k');
            elseif j==4
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'r');
            else
                hold on;plot(lfptime_sleep(indlfp,1),lfp7(indlfp,j)*factor+dHPC_layer7channel3(j,2),'k');
            end
        end
        %xline(nonrem(i,4:5),'k');
        xlim([min(lfptime_sleep(indlfp,1)),max(lfptime_sleep(indlfp,1))]);title(num2str(in));
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        
        subplot(2,1,2);
        imagesc(lfptime_sleep(indlfp,1),[channels(1,3) channels(end,3)],CSDall_sleep(:,indlfp));set(gca,'YDir','normal');
        %xline(nonrem(i,4:5),'k');
        xlabel('Time in Slow Wave Sleep (s)');xlim([min(lfptime_sleep(indlfp,1)),max(lfptime_sleep(indlfp,1))]);
        yticks(dHPC_layer7channel3(end:-1:1,2));yticklabels({'dCA3','DG GCL','DG MML','DG OML','dCA1 slm','dCA1 st rad','dCA1 pyr'});
        %pause;clf;
    end
    set(findall(gcf, 'Type', 'Line'),'LineWidth',1.5,'MarkerSize', 6);set(gcf,'color','w');set(findobj(gcf,'type','axes'),'FontSize',15);
   
end
cd(homedir);save('SlowWaveSleep_Example',"timeset");