%function generate_Figure_2_neuron_Jan25(ifsavedata)

%if nargin<1
    ifsavedata=0;
%end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
neuron_colname={'1.Neuron ID','2.Depth on Probe','3. Quality Class','4.Unit Type'};
spike_colname={'1.Spike Time','2.Neuron ID','3.Theta phase','4.Thetacycle id 13','5.Velocity level','6.Behavior marker'};
if 0
thetacycle_N=[];neuronall=[];neuron_firerateall=[];timedur=[];LFP_Frequency=625;
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('theta_trough_trough','thetacycle');
    load('Analyze_Theta_Cycle_Broad','reftheta','lfptime_trial');
    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    ind=unit_decode(:,10)<=3;
    neuron=unit_decode(ind,[1,5,10,16,17]);
    neuron(:,1)=[1:size(neuron,1)];
    ind=spike_trial(:,4)<=3;
    spike=spike_trial(ind,1:2);spike(:,3)=0;
    for i=1:size(neuron,1)
        ind=spike(:,2)==neuron(i,5);
        spike(ind,3)=neuron(i,1);
    end
    spike=spike(:,[1,3]);
    neuron=neuron(:,1:4);
    ruler=lfptime_trial(:,1);template=spike(:,1);[outputindex,error]=match(template,ruler,0);
    spike(:,3:4)=[reftheta(outputindex,4),lfptime_trial(outputindex,5)];
    
    indstill0=thetacycle(:,10)<1 & thetacycle(:,10)>=0;
    indstill1=thetacycle(:,10)<5 & thetacycle(:,10)>=1;
    indstill2= thetacycle(:,10)<10 & thetacycle(:,10)>=5;
    indlowv= thetacycle(:,10)>=10 ;
    thetacycle_N(in,:)=[sum(indstill0),sum(indstill1),sum(indstill2),sum(indlowv)];
    for v=1:4
        if v==1
            indv=thetacycle(:,10)<1 & thetacycle(:,10)>=0;
        elseif v==2
            indv=thetacycle(:,10)<5 & thetacycle(:,10)>=1;
        elseif v==3
            indv=thetacycle(:,10)<10 & thetacycle(:,10)>=5;
        elseif v==4
            indv=thetacycle(:,10)>=10 ;
        end
        for b=1:5
            if b<5
                indb=indv & thetacycle(:,79)==0 & thetacycle(:,62)==b;
            else
                indb=indv & thetacycle(:,79)==0;
            end
            lia=ismember(lfptime_trial(:,5),thetacycle(indb,13));
            timedur(in,b,v,:)=histcounts(reftheta(lia,4),[0:10:360])/LFP_Frequency;
        end
    end
    neuron(:,5)=0;spike(:,5:7)=0;neuron_firerate=nan*ones(size(neuron,1),5,4,36,4);neuron_fireindex=nan*ones(size(neuron,1),5,4,8);
    for n=1:size(neuron,1)
        indn=spike(:,2)==neuron(n,1);
        neuron(n,5)=sum(indn);
        if sum(indn)>=10
            for v=1:4
                if v==1
                    indv=thetacycle(:,10)<1 & thetacycle(:,10)>=0;
                elseif v==2
                    indv=thetacycle(:,10)<5 & thetacycle(:,10)>=1;
                elseif v==3
                    indv=thetacycle(:,10)<10 & thetacycle(:,10)>=5;
                elseif v==4
                    indv=thetacycle(:,10)>=10 ;
                end
                for b=1:5
                    if b<5
                        indb=indv & thetacycle(:,79)==0 & thetacycle(:,62)==b ;
                    else
                        indb=indv & thetacycle(:,79)==0 ;
                    end
                    lia=ismember(spike(:,4),thetacycle(indb,13));
                    ind= lia & indn;
                    if sum(ind)>=10
                        spike(ind,5)=v;
                        if b<5
                            spike(ind,6)=b;
                        end
                        neuron_firerate(n,b,v,:,1)=histcounts(spike(ind,3),[0:10:360]); % raw spike n
                        neuron_firerate(n,b,v,:,2)=squeeze(neuron_firerate(n,b,v,:,1))./squeeze(timedur(in,b,v,:)); % fire rate
                        a=gaussian_smooth(squeeze(neuron_firerate(n,b,v,:,2)));
                        neuron_firerate(n,b,v,:,3)=a; % smoothed fire rate
                        neuron_firerate(n,b,v,:,4)=(a-min(a))/(max(a)-min(a)); % fire index
                        [M,I]=max(a);
                        [mu,r]=circ_meand(spike(ind,3));
                        spi=spike(ind,3)*pi/180;
                        [pval] = circ_rtest(spi);
                        neuron_fireindex(n,b,v,:)=[sum(ind),(max(a)-min(a))/max(a),I*10-5,max(a),mean(a),mu,r,pval];
                        % 1.spike n | 2.mod index | 3.peak phase | 4.peak rate | 5.mean rate | 6.mean phase | 7.resultant length | 8. pval phase lock
                    end
                end
            end
        end
    end
    neuron(:,6)=in;
    neuronall=[neuronall;neuron];
    neuron_firerateall=cat(1,neuron_firerateall,neuron_firerate);
    cd(savedir);save('Decode_Spike',"neuron","spike",'spike_colname','neuron_colname','neuron_firerate','neuron_fireindex','timedur','thetacycle','-v7.3');
    in
end
types=[1,2,5,6];neuronall=[];neuron_firerateall=[];neuron_fireindexall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Decode_Spike');
    lia=ismember(neuron(:,4),types);
    neuronall=[neuronall;neuron(lia,:)];
    neuron_firerateall=cat(1,neuron_firerateall,neuron_firerate(lia,:,:,:,:));
    neuron_fireindexall=cat(1,neuron_fireindexall,neuron_fireindex(lia,:,:,:));
end
% below quantify major, minor peak for ca1 cells
ca1ind=neuronall(:,4)==1;
ca1_fire=neuron_firerateall(ca1ind,:,:,:,:);
ca1_index=neuron_fireindexall(ca1ind,:,:,:);
ca1=neuronall(ca1ind,:);
pkall=[];spike_N_cutoff=30;ca1_bimodality=[];
for i=1:size(ca1_index,1)
    for b=1:5
        for v=1:4
            if  ca1_index(i,b,v,1)>=spike_N_cutoff
                a=squeeze(ca1_fire(i,b,v,:,3))';
                a=[a,a];
                [M,I]=min(a);
                a1=a(I:I+35);
                [pks,locs]=findpeaks(a1(:),'MinPeakDistance',1);%,'MinPeakWidth',2
                pk1=[locs+I-1,pks,locs];
                lc=pk1(:,1)-36*(pk1(:,1)>36);
                pk1(:,5)=lc;
                pk1(:,6)=b;
                pk1(:,7)=v;
                pk1(:,8)=ca1(i,3);
                pk1(:,6)=0;
                ind= lc>=7 & lc<=28 ;pk1(ind,6)=1; % major peak
                pk1(~ind,6)=2; % minor peak
                ind=pk1(:,6)>0;pk1=pk1(ind,:);
                for k=1:2
                    ind=pk1(:,6)==k;
                    if sum(ind)==1
                        ca1_bimodality(i,b,v,k,1:5)=pk1(ind,1:5);
                    elseif sum(ind)==0
                        ca1_bimodality(i,b,v,k,1:5)=nan;
                    elseif sum(ind)>1
                        pk2=pk1(ind,:);
                        [M,I]=max(pk2(:,2));
                        ca1_bimodality(i,b,v,k,1:5)=pk2(I,1:5);
                    end
                end
                ca1_bimodality(i,b,v,:,6)=[nanmean(a(7:28)),nanmean(a(29:42))];
                ca1_bimodality(i,b,v,:,7:8)=ca1_bimodality(i,b,v,:,[2,6])./ca1_bimodality(i,b,v,1,[2,6]);
                pkall=[pkall;pk1];
            end
        end
    end
end
ca1_bimodality(:,:,:,:,9)=ca1_bimodality(:,:,:,:,5)*10-5;

    cd(datadir);save('Figure2_unit_fire_data','ca1_bimodality','ca1',"ca1_fire",'ca1_index','neuronall','neuron_firerateall', ...
    'neuron_fireindexall','pkall','spike_N_cutoff','-v7.3');
else
    cd(datadir);load('Figure2_unit_fire_data');
end

unittype4 = {'dCA1 Principal','dCA3DG Principal','dCA1 Interneuron','dCA3DG Interneuron'};
vgroup={'0-1 cm/s','1-5 cm/s','5-10 cm/s','> 10 cm/s'};
plotcolorscat=slanCM('YlGnBu',5);plotcolorscat=plotcolorscat([5:-1:1],:);plotcolorscat(5,:)=[0.8828 0.7305 0.2656];plotcolorscat(1,:)=0;
plotcolors4=slanCM('matter',8);plotcolors4=plotcolors4([8:-2:1],:);plotcolors4(4,:)=[0.8828 0.7305 0.2656];
types=[1,2,5,6];ylimset=[6,6,27,27];
colname=[];
for v=1:4
    for i=1:36
        colname{i+(v-1)*36}=[vgroup{v},'_',num2str(i)];
    end
end
fign=2;pn1=2;pn2=4;b=5;spike_N_cutoff=30;
% mean firing rate
for i=1:4
    figure(fign);subplot(pn1,pn2,i);
    ind=min(neuron_fireindexall(:,b,:,1),[],3)>=spike_N_cutoff & neuronall(:,4)==types(i);
    for v=1:4
        y=squeeze(neuron_firerateall(ind,b,v,:,2))';y=[y;y];
        hold on;shaded_errbar([5:10:715],y,plotcolors4(v,:));
    end
    xlim([0 720]);xlabel('Theta Phase');title([unittype4{i},': N = ',num2str(sum(ind))]);
    ylabel('Mean Firing Rate');ylim([0 ylimset(i)]);xticks([0:360:720]);
    legend('0-1 cm/s','','1-5 cm/s','','5-10 cm/s','','> 10 cm/s','');
    if ifsavedata
        figdata=[squeeze(neuron_firerateall(ind,b,1,:,2)),squeeze(neuron_firerateall(ind,b,2,:,2)),squeeze(neuron_firerateall(ind,b,3,:,2)),squeeze(neuron_firerateall(ind,b,4,:,2))];
        T = array2table(figdata,'VariableNames',colname);
        writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2j_',num2str(i)]);
    end
end

% sig neuron percent
figure(fign);subplot(pn1,pn2,5);figdata=[];comp=[];ns=[];
for i=1:4
    ind=min(neuron_fireindexall(:,b,:,1),[],3)>=spike_N_cutoff & neuronall(:,4)==types(i);
    y=mean(squeeze(neuron_fireindexall(ind,b,:,8))<=0.05,1);
    hold on;plot([1:4],y,'Color',plotcolorscat(i+1,:));
    figdata(i,:)=y;
    n1=sum(ind);n2=sum(squeeze(neuron_fireindexall(ind,b,:,8))<=0.05,1);
    ns(i)=n1;
    for v=1:4
        comp(i,v)=binopdf(n2(v),n1,0.05);
    end
end    
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Theta of Velocity Ranges');
ylabel('Neuron Sig Percent');title('dHPC Unit');ylim([0 1]);legend(unittype4);title({num2str(comp,2);num2str(ns)});
if ifsavedata
T = array2table(figdata,'VariableNames',{'0-1 cm/s','1-5 cm/s','5-10 cm/s','> 10 cm/s'},'RowNames',{'dCA1 Principal','dCA3DG Principal','dCA1 Interneuron','dCA3DG Interneuron'});
writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2k'],'WriteRowNames',true);
end

% theta mod index
lett={'A','F','K','P'};comp=[];ns=[];stats=[];
for i=1:4
    figure(fign);subplot(pn1,pn2,6);
    ind=min(neuron_fireindexall(:,b,:,1),[],3)>=spike_N_cutoff & neuronall(:,4)==types(i);
    y=squeeze(neuron_fireindexall(ind,b,:,2));
    hold on;shaded_errbar([1:4],y',plotcolorscat(i+1,:));
    [pv,tbl]= anova1(y);
    comp(i)=pv;ns(i)=size(y,1);stats(i)=tbl{2,5};
    %comp(i)=signrank(y(:,4)-y(:,1));
    if ifsavedata
    figdata=y;
    T = array2table(figdata,'VariableNames',{[unittype4{i},'_0-1 cm/s'],[unittype4{i},'_1-5 cm/s'],[unittype4{i},'_5-10 cm/s'],[unittype4{i},'_> 10 cm/s']});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2l'],'WriteRowNames',true,'Range',[lett{i},'1']);
    end
end    
figure(fign);subplot(pn1,pn2,6);
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Theta of Velocity Ranges');
ylabel('Neuron Theta Modulation Index');title({num2str(stats);num2str(comp,2);num2str(ns)});ylim([0 1]);
legend({unittype4{1},'',unittype4{2},'',unittype4{3},'',unittype4{4},''});

% peak rate
figure(fign);subplot(pn1,pn2,7);ca1name={'CA1 Major','CA1 Minor'};lett={'A','F','K','P','U'};comp=[];ns=[];stats=[];
for i=1%:2
    ind=min(ca1_index(:,b,:,1),[],3)>=spike_N_cutoff ;
    y=squeeze(ca1_bimodality(ind,b,:,i,6));y=y./y(:,1);indy=sum(isinf(y),2)==0;y=y(indy,:);
    hold on;shaded_errbar([1:4],y',plotcolorscat(i,:));
    [pv,~,stat]=signrank(y(:,4),1,"tail","right",'method','approximate');
    comp(i)=pv;stats(i)=stat.zval;
    ns(i)=size(y,1);
    if ifsavedata
    figdata=y;
    T = array2table(figdata,'VariableNames',{[ca1name{i},'_0-1 cm/s'],[ca1name{i},'_1-5 cm/s'],[ca1name{i},'_5-10 cm/s'],[ca1name{i},'_> 10 cm/s']});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2m'],'WriteRowNames',true,'Range',[lett{i},'1']);
    end
end    
for i=2:4
    ind=min(neuron_fireindexall(:,b,:,1),[],3)>=spike_N_cutoff & neuronall(:,4)==types(i);
    y=squeeze(neuron_fireindexall(ind,b,:,4));y=y./y(:,1);
    hold on;shaded_errbar([1:4],y',plotcolorscat(i+1,:));
    [pv,~,stat]=signrank(y(:,4),1,"tail","right",'method','approximate');
    comp(i)=pv;stats(i)=stat.zval;
    ns(i)=size(y,1);
    if ifsavedata
    figdata=y;
    T = array2table(figdata,'VariableNames',{[unittype4{i},'_0-1 cm/s'],[unittype4{i},'_1-5 cm/s'],[unittype4{i},'_5-10 cm/s'],[unittype4{i},'_> 10 cm/s']});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2m'],'WriteRowNames',true,'Range',[lett{i+1},'1']);
    end
end    
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Theta of Velocity Ranges');
ylabel('Neuron peak rate/ rate at 0-1');title({num2str(stats);num2str(comp,2);num2str(ns)});%ylim([0 1]);
legend({'CA1 P major','',unittype4{2},'',unittype4{3},'',unittype4{4},''});%,'CA1 P minor',''

% peak phase
figure(fign);subplot(pn1,pn2,8);comp=[];ns=[];stats=[];
for i=1%:2
    figure(fign);subplot(pn1,pn2,8);
    ind=min(ca1_index(:,b,:,1),[],3)>=spike_N_cutoff ;
    y=squeeze(ca1_bimodality(ind,b,:,i,9));
    indy=sum(isinf(y) | isnan(y),2)==0;y=y(indy,:);
    hold on;er=shaded_errbar([1:4],y',plotcolorscat(i,:),1);
    n1=size(y,1);idp=ones(n1,1)*[1:4];
    [pval table] = circ_wwtest(y(:)*pi/180,idp(:));
    comp(i)=pval;stats(i)=table{2,5};
    ns(i)=size(y,1);
    if ifsavedata
    figdata=y;
    T = array2table(figdata,'VariableNames',{[ca1name{i},'_0-1 cm/s'],[ca1name{i},'_1-5 cm/s'],[ca1name{i},'_5-10 cm/s'],[ca1name{i},'_> 10 cm/s']});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2n'],'WriteRowNames',true,'Range',[lett{i},'1']);
    end
end    
for i=2:4
    figure(fign);subplot(pn1,pn2,8);
    ind=min(neuron_fireindexall(:,b,:,1),[],3)>=spike_N_cutoff & neuronall(:,4)==types(i);
    y=squeeze(neuron_fireindexall(ind,b,:,3));
    indy=sum(isinf(y) | isnan(y),2)==0;y=y(indy,:);
    hold on;shaded_errbar([1:4],y',plotcolorscat(i+1,:),1);
    n1=size(y,1);idp=ones(n1,1)*[1:4];
    [pval table] = circ_wwtest(y(:)*pi/180,idp(:));
    comp(i)=pval;stats(i)=table{2,5};
    ns(i)=size(y,1);
    if ifsavedata
    figdata=y;
    T = array2table(figdata,'VariableNames',{[unittype4{i},'_0-1 cm/s'],[unittype4{i},'_1-5 cm/s'],[unittype4{i},'_5-10 cm/s'],[unittype4{i},'_> 10 cm/s']});
    writetable(T,'C:\Users\Mengni\Documents\Paper\Figure_Nov2024\Figure_Data.xlsx','Sheet',['Figure 2n'],'WriteRowNames',true,'Range',[lett{i+1},'1']);
    end
end    
xticks([1:4]);xticklabels(vgroup);xlim([.5 4.5]);xlabel('Theta of Velocity Ranges');
ylabel('Neuron peak phase');title('dHPC Unit');ylim([0 380]);title({num2str(stats);num2str(comp,2);num2str(ns)});
legend({'CA1 P major','',unittype4{2},'',unittype4{3},'',unittype4{4},''});%,'CA1 P minor',''
figure(fign);figure_title='Figure_2_unit_fire_Jan2025';save_current_figure(figure_title);