%function analyze_SessionSet16

% 1. save sessionset16 dir
clear all;close all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';
%cd(homedir);load('SessionSet');
%SessionSet16=SessionSet(1:14);%SessionSet16{15}=SessionSet{17};SessionSet16{16}=SessionSet{18};sessionnum=16;
%cd(homedir);save('SessionSet16','SessionSet16',"basedir",'homedir',"sessionnum");
cd(homedir);load('SessionSet16');
unittype8 = {'dCA1 Pyramidal','dCA3DG Pyramidal','vDGCA3 Pyramidal','v subiculum Pyramidal','dCA1 Interneuron','dCA3DG Interneuron','vDGCA3 Interneuron','v subiculum Interneuron'};
% 2. save spike data
spikeratio=[];inter_pyr_cutoff=0.5; % ms
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    load('Spike_Data.mat','spike','unit_metric_decode','cluster_columnname');
    load('LFP_CSD.mat','channelcsd');
    unit=unit_metric_decode;
    unit(:,13)=unit(:,6)>=inter_pyr_cutoff;% if excitatory (1)
    ind=unit(:,6)>1;unit(ind,13)=-1;% spike width too wide; mark as abnormal
    unit(:,17)=unit(:,12);ind=unit(:,12)>=5;unit(ind,17)=-10;
    vHPC_layer=struct;indvca3=[];indvdg=[];indvca1=[];indvdgca3=[];indnotHPC=[];
    % below curate for individual sessions
    if ceil(in/2)==1 % vCA3-vDG-vCA1
        vHPC_layer.vCA1=[10,16];
        vHPC_layer.MoDG1=[24,31];
        vHPC_layer.DG_Unit=29;
        vHPC_layer.vCA3=49;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    elseif in==3 % vDG-vCA3-vDG-vCA1
        vHPC_layer.vCA1=[1,5];
        vHPC_layer.MoDG2=[14,17,21,43,48];
        vHPC_layer.vCA3=34;
        indvdgca3= unit(:,5)>=channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
    elseif in==4 % vDG-vCA3-vDG-vCA1
        vHPC_layer.vCA1=[1,6];
        vHPC_layer.MoDG2=[15,18,22,44,49];
        vHPC_layer.vCA3=35;
        indvdgca3= unit(:,5)>=channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
    elseif in==5 % vDG-vCA3-vDG-vCA1
        vHPC_layer.vCA1=[21,26];
        vHPC_layer.MoDG2=[35,41,67,72];
        vHPC_layer.vCA3=49;
        vHPC_layer.HPCCUT=10;
        indvdgca3= unit(:,5)>=channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)>channelcsd(vHPC_layer.HPCCUT,3) & unit(:,5)<channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
        indnotHPC=unit(:,5)<=channelcsd(vHPC_layer.HPCCUT,3) & unit(:,12)>=5;
    elseif in==6 % vDG-vCA3-vDG-vCA1
        vHPC_layer.vCA1=[21,26];
        vHPC_layer.MoDG2=[34,41,61,66,70];
        vHPC_layer.vCA3=50;
        vHPC_layer.HPCCUT=10;
        indvdgca3= unit(:,5)>=channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)>channelcsd(vHPC_layer.HPCCUT,3) & unit(:,5)<channelcsd(vHPC_layer.MoDG2(1),3) & unit(:,12)>=5;
        indnotHPC=unit(:,5)<=channelcsd(vHPC_layer.HPCCUT,3) & unit(:,12)>=5;
    elseif ceil(in/2)==4 % vCA3-vDG-vCA1
        vHPC_layer.vCA1=[33,39];
        vHPC_layer.MoDG1=[48,57,62];
        vHPC_layer.DG_Unit=53;
        vHPC_layer.vCA3=67;
        vHPC_layer.HPCCUT=25;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)>channelcsd(vHPC_layer.HPCCUT,3) & unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
        indnotHPC=unit(:,5)<=channelcsd(vHPC_layer.HPCCUT,3) & unit(:,12)>=5;
    elseif ceil(in/2)==5 % vCA3-vDG-vCA1
        vHPC_layer.vCA1=[9,15];
        vHPC_layer.MoDG1=[24,31,37,43];
        vHPC_layer.DG_Unit=33;
        vHPC_layer.vCA3=49;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    elseif ceil(in/2)==6 % vCA3-vDG-vCA1
        vHPC_layer.vCA1=[1,4];
        vHPC_layer.MoDG1=[15,21,26,32];
        vHPC_layer.DG_Unit=26;
        vHPC_layer.vCA3=39;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    elseif in==13 % vCA3-vDG-vCA1
        vHPC_layer.vCA1=[3,9];
        vHPC_layer.MoDG1=[17,23];
        vHPC_layer.vCA3=35;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    elseif in==14 % vCA3-vDG-vCA1
        vHPC_layer.vCA1=[3,10];
        vHPC_layer.MoDG1=[18,23];
        vHPC_layer.vCA3=36;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    elseif in==15 % vCA3-'vDG'-vCA1
        vHPC_layer.vCA1=[3,10];
        vHPC_layer.MoDG1=[34,30,38];
        vHPC_layer.vCA3=43;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    elseif in==16 % vCA3-'vDG'-vCA1
        vHPC_layer.vCA1=[3,10];
        vHPC_layer.MoDG1=[33,29,37];
        vHPC_layer.vCA3=42;
        indvca3= unit(:,5)>channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvdg= unit(:,5)>=channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,5)<=channelcsd(max(vHPC_layer.MoDG1),3) & unit(:,12)>=5;
        indvca1= unit(:,5)<channelcsd(vHPC_layer.MoDG1(1),3) & unit(:,12)>=5;
    end
    if ~isempty('indvca1')
        unit(indvca1,[12,17])=7;
    end
    if ~isempty('indvdg')
        unit(indvdg,[12,17])=6;
    end
    if ~isempty('indvca3')
        unit(indvca3,[12,17])=5;
    end
    if ~isempty('indvdgca3')
        unit(indvdgca3,[12,17])=5.5;
    end
    if ~isempty('indnotHPC')
        unit(indnotHPC,[12,17])=-1;
    end
    indinter=unit(:,6)<inter_pyr_cutoff;
    layer=zeros(size(unit,1),1);
    ind=unit(:,12)==1;layer(ind)=1; %dca1
    ind=unit(:,12)>=2 & unit(:,12)<=4;layer(ind)=2; % dca3dg
    ind=unit(:,12)>=5 & unit(:,12)<=6;layer(ind)=3; % vca3dg
    ind=unit(:,12)==7;layer(ind)=4; % ventral subiculum
    unit(:,16)=double(indinter)*4+layer;
    ind=unit(:,6)>1 | unit(:,12)==-1;unit(ind,16)=9;
    unit_layer(in,:)=histcounts(unit(:,16),[.5:1:9.5]);
    
    ind=unit(:,16)<9 & unit(:,16)>0;unit_decode=unit(ind,1:16);
    [~,I]=sort(unit_decode(:,5),'descend'); % sorted by recording depth from dorsal to ventral
    unit_decode=unit_decode(I,:);
    unit_decode(:,17)=[1:size(unit_decode,1)]; % new unit id
    cluster_columnname{16}='Unit_type_region_code(1-9)';
    cluster_columnname{17}='new unit id';

    % below for spike data
    spike(:,9)=0; % unit type
    for k=1:9
        ind=unit(:,16)==k;
        lia=ismember(spike(:,6),unit(ind,11));
        spike(lia,9)=k;
    end
    spike(:,10)=0; % sleepbox or trial id
    for s=1:2
        ind= spike(:,1)>=sleep(s,1) & spike(:,1)<=sleep(s,2) ;
        spike(ind,10)=-s;
    end
    for i=1:size(trialepoch,1)
        ind= spike(:,1)>=trialepoch(i,1) & spike(:,1)<=trialepoch(i,2) ;
        spike(ind,10)=i;
    end
    spike(:,11)=0; % new id that corresponds to unit_decode
    for i=1:size(unit_decode,1)
        ind=spike(:,6)==unit_decode(i,11);
        spike(ind,11)=unit_decode(i,17);
    end
    spikedata=spike(:,[1,11,6,5,9,10]); % 1.time | 2.unit id | 3.old goodokay unit id | 4.quality level | 5.unittype(1-9) | 6. sleepbox or trialid | 7. position match index
    [~,I]=sort(spikedata(:,1));
    spikedata=spikedata(I,:);

    ind= spikedata(:,6)>0 & spikedata(:,2)>0;
    spike_trial=spikedata(ind,:);
    ind= (spikedata(:,6)==-1 | spikedata(:,6)==-2) & spikedata(:,2)>0;
    spike_sleepbox=spikedata(ind,:);

    ruler=Position_Data(:,1);template=spike_trial(:,1);[outputindex,error]=match(template,ruler,0);
    spike_trial(:,7)=outputindex;
    spikedata_columnname='1.time | 2.unit id | 3.old goodokay unit id | 4.quality level | 5.unittype(1-9) | 6. sleepbox or trialid | 7. position match index';
    cd(savedir);save('Spike_Data_decode_trial_sleepbox','unit','unit_decode','spikedata','spike_trial','spike_sleepbox',...
        'cluster_columnname','spikedata_columnname','vHPC_layer','-v7.3');
  
    spikeratio(in,:)=[size(spike,1),size(spikedata,1),size(spike_trial,1),size(spike_sleepbox,1)];
    figure(1);subplot(4,4,in);histogram(unit(:,17));title(sum(unit(:,17)==-10));
    figure(2);subplot(4,4,in);histogram(error); 
end

% below save hpc layer channel data as a separate dateset
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG IML/GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
hpcall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('CA1_CA3DG_Ripples.mat', 'dHPC_layer7channel3');
    load('Spike_Data_decode_layer_defined.mat', 'vHPC_layer');
    load('LFP_CSD.mat', 'channelcsd','lfpcsd_ref','LFP_Frequency','lfptime');
    if isfield(vHPC_layer,'MoDG1') & isfield(vHPC_layer,'MoDG2')
        'error!'
        in
    end
    if isfield(vHPC_layer,'MoDG1')
        dgch=vHPC_layer.MoDG1(1);
    elseif isfield(vHPC_layer,'MoDG2')
        dgch=vHPC_layer.MoDG2(1);
    else
        dgch=0;
    end
    vHPC_layer_Channel=[vHPC_layer.vCA3,dgch,vHPC_layer.vCA1([2,1])]';
    vHPC_layer_Channel(:,2:3)=channelcsd(vHPC_layer_Channel(:,1),3:4);
    HPC_Layer_Channel=[dHPC_layer7channel3;vHPC_layer_Channel];
    LFP_HPC=lfpcsd_ref(:,HPC_Layer_Channel(:,1));
    cd(savedir);save('HPC_Layer_Channel_LFP','HPC_Layer_Channel','LFP_HPC','lfptime','LFP_Frequency','dHPC_layer7channel3','vHPC_layer');
    hpcall(:,:,in)=HPC_Layer_Channel;
    clear dHPC_layer7channel3 vHPC_layer channelcsd lfpcsd_ref lfptime dgch vHPC_layer_Channel HPC_Layer_Channel LFP_HPC
    in
end

decode_unit_N=[];unitall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Spike_Data_decode_trial_sleepbox','unit_decode');
    decode_unit_N(in,:)=histcounts(unit_decode(:,16),[0.5:9.5]);
    unit_decode(:,18)=in;
    unitall=[unitall;unit_decode];
end
%---------------------------------------------------------------------------------
% 3. calculate place field
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    CALCULATE_PLACE_FIELDS_8ARMS_NPX(savedir);
    disp(['Completed for ',num2str(in)]);
end
%---------------------------------------------------------------------------------
% 4. decode with 250ms for decoding accuracy
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.25;WinAdvance=DecodeWin;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');

    indv=Position_Data(:,5)>=Velocity_Cutoff;
    PositionTrain=Position_Data(indv,2:3);
    PosProbDist=zeros(size(PositionBasis,1),1);
    for i=1:size(PositionBasis,1)
        PosProbDist(i)=sum(mvnpdf(PositionBasis(i,:),PositionTrain,Sigma_Pos));
    end
    PosProbDist=PosProbDist/sum(PosProbDist);
    PosProbDist0=PosProbDist;
    indnan1d=log10(PosProbDist0)<-maze_cutoff;PosProbDist0(indnan1d)=nan; 
    maze1d=[~indnan1d,PositionBasis];
    maze1d_position=maze1d(maze1d(:,1)>0,2:3);
    Position_Transition=zeros(size(maze1d_position,1));
    for i=1:size(maze1d_position,1)
        indpos= find( Position_Data(1:end-1,5)>=10 & abs(Position_Data(1:end-1,2)-maze1d_position(i,1))<=binsize/2 & abs(Position_Data(1:end-1,3)-maze1d_position(i,2))<=binsize/2 );
        pospost=Position_Data(indpos+1,2:3);
        for j=1:size(maze1d_position,1)
            indpospost= (abs(pospost(:,1)-maze1d_position(j,1))<=binsize/2 & abs(pospost(:,2)-maze1d_position(j,2))<=binsize/2);
            Position_Transition(i,j)=mean(indpospost);
        end
    end
    Position_Transition_Smooth=[];
    for i=1:size(Position_Transition,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=Position_Transition_SmoothFactor*binsize;
        Position_Transition_Smooth(i,:)=nanmean(Position_Transition(ind,:),1);
    end

    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    spikedata=spike_trial; % select for hpc cells
    spikedata(:,7:10)=Position_Data(spikedata(:,7),[5,2,3,4]);
    TrialSet=unique(Position_Data(:,4));
    TrialNum=length(TrialSet); 
    decode_result=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7,5]);
        [Decoded_Sequence,~,dedata_smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        ind=Decoded_Sequence(:,9)>=10;
        decode_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8,13]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8,13]),1)];
        disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    decode_session(in,:)=mean(decode_result,1);
    PosProbDist_set{in}=PosProbDist;
    maze1d_set{in}=maze1d;
    Position_Transition_Smooth_set{in}=Position_Transition_Smooth;
    Position_Transition_set{in}=Position_Transition;
    
    % below to shuffle unit id
    for shuffle=1:0%100
    p=randperm(size(Field_Data,3));
    shuffle_result=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7]);
        [Decoded_Sequence,~]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data(:,:,p),maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        ind=Decoded_Sequence(:,9)>=10;
        shuffle_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8]),1)];
        %disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    shuffle_result_session(shuffle,:,in)=mean(shuffle_result,1);
    end
    disp(['Finish Session  ',num2str(in),' ; Performance = ',num2str(decode_session(in,:),2)]);
    %disp(['Shuffle : ',num2str(mean(shuffle_result_session(:,:,in),1),2)]);
end
cd(homedir);save('sorted_spike_decoding_8arm_setup_trial_sleepbox','Position_Transition_set','Position_Transition_Smooth_set','maze1d_set',...
    "PosProbDist_set",'decode_session','BinSigma','maze_cutoff','DecodeWin','WinAdvance','Position_Transition_SmoothFactor',...
    "shuffle_result_session",'SessionSet16');
%------------------------------------------------------------------------------------
% 5. decoding with 20 ms window for info content
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.02;WinAdvance=0.01;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    indv=Position_Data(:,5)>=Velocity_Cutoff;
    PositionTrain=Position_Data(indv,2:3);
    PosProbDist=zeros(size(PositionBasis,1),1);
    for i=1:size(PositionBasis,1)
        PosProbDist(i)=sum(mvnpdf(PositionBasis(i,:),PositionTrain,Sigma_Pos));
    end
    PosProbDist=PosProbDist/sum(PosProbDist);
    PosProbDist0=PosProbDist;
    indnan1d=log10(PosProbDist0)<-maze_cutoff;PosProbDist0(indnan1d)=nan; 
    maze1d=[~indnan1d,PositionBasis];
    maze1d_position=maze1d(maze1d(:,1)>0,2:3);
    Position_Transition=zeros(size(maze1d_position,1));
    for i=1:size(maze1d_position,1)
        indpos= find( Position_Data(1:end-1,5)>=10 & abs(Position_Data(1:end-1,2)-maze1d_position(i,1))<=binsize/2 & abs(Position_Data(1:end-1,3)-maze1d_position(i,2))<=binsize/2 );
        pospost=Position_Data(indpos+1,2:3);
        for j=1:size(maze1d_position,1)
            indpospost= (abs(pospost(:,1)-maze1d_position(j,1))<=binsize/2 & abs(pospost(:,2)-maze1d_position(j,2))<=binsize/2);
            Position_Transition(i,j)=mean(indpospost);
        end
    end
    Position_Transition_Smooth=[];
    for i=1:size(Position_Transition,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=Position_Transition_SmoothFactor*binsize;
        Position_Transition_Smooth(i,:)=nanmean(Position_Transition(ind,:),1);
    end
    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    spikedata=spike_trial; % select for hpc cells
    spikedata(:,7:10)=Position_Data(spikedata(:,7),[5,2,3,4]);
    TrialSet=unique(Position_Data(:,4));
    TrialNum=length(TrialSet); 
    decode_result=[];deseqall=[];dedataall=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7,5]);
        [Decoded_Sequence,~,Decoded_Data_Smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        Decoded_Sequence(:,14)=trialid;
        deseqall=[deseqall;Decoded_Sequence];
        dedataall=[dedataall,Decoded_Data_Smooth];
        ind=Decoded_Sequence(:,9)>=10;
        decode_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8,13]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8,13]),1)];
        disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    decode_session(in,:)=mean(decode_result,1);
    PosProbDist_set{in}=PosProbDist;
    maze1d_set{in}=maze1d;
    Position_Transition_Smooth_set{in}=Position_Transition_Smooth;
    Position_Transition_set{in}=Position_Transition;
    cd(savedir);save(['sorted_spike_decoding_8arm5'],'decode_result','deseqall','dedataall','PosProbDist',...
        'maze1d','Position_Transition_Smooth','Position_Transition','BinSigma','maze_cutoff','DecodeWin','WinAdvance',...
        'savedir','Position_Transition_SmoothFactor','-v7.3');
    disp(['Finish Session  ',num2str(in),' ; Performance = ',num2str(decode_session(in,:),2)]);
end
%----------------------------------------------------------------------------------
% 6. below add spec power to deseqall
clear all;close all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('local_remote_decoding.mat');
    load('Spectrum_InOutInterReward','spec_trial_wavelet','lfptime_trial','fre');
    ca3spec=spec_trial_wavelet{7};clear spec_trial_wavelet
    indhighfre=fre>=150;
    indfastgamma=fre>=100 & fre<150;
    indmidgamma=fre>=50 & fre<100;
    indlowgamma=fre>=30 & fre<50;
    ca3power=[nanmean(ca3spec(indhighfre,:),1)',nanmean(ca3spec(indlowgamma,:),1)',nanmean(ca3spec(indmidgamma,:),1)',nanmean(ca3spec(indfastgamma,:),1)'];
    ruler=lfptime_trial(:,1);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
    deseqall(:,15:18)=ca3power(outputindex,:);
    cd(savedir);save('local_remote_decoding','deseqall','linear_pos','pksets','inter_run','trialset','-v7.3');
    disp(['Finish Session : ',num2str(in)]);
end
%----------------------------------------------------------------------------------
% 7. below add CSD value to deseqall
clear all;close all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('local_remote_decoding.mat');
    load('CSD_LFP_Visualization.mat','CSDall','lfptime','channels','dHPC_layer7channel3');
    dHPC_layer7channel3(:,4)=0;
    for i=1:7
        dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2));
    end
    CSDall7=CSDall(dHPC_layer7channel3(:,4),:)';
    ruler=lfptime(:,1);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
    deseqall(:,40:46)=CSDall7(outputindex,:);

    load('Ripple_DS_Spike.mat','dspeaks_d4');
    ind=dspeaks_d4(:,18)>=2;
    dsmatch=dspeaks_d4(ind,:);
    for ds=1:3
        indds=dsmatch(:,6)==ds;
        ruler=dsmatch(indds,9);template=deseqall(:,1);[outputindex,error]=match(template,ruler,0);
        deseqall(:,46+ds)=error;
    end
    cd(savedir);save('local_remote_decoding','deseqall','linear_pos','pksets','inter_run','trialset','-v7.3');
    clear CSDall lfptime dHPC_layer7channel3 channels deseqall linear_pos dspeaks_d4
    disp(['Finish Session : ',num2str(in)]);
end
%----------------------------------------------------------------------------------
% 8. decoding with 20 ms window for info content during sleepbox
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.02;WinAdvance=0.01;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('sorted_spike_decoding_8arm5','PosProbDist','maze1d','Position_Transition_Smooth','Position_Transition');
    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_sleepbox');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    deseqall=[];dedataall=[];
    for i=1:2
        ind=spike_sleepbox(:,6)==-i; 
        Spike_Data=spike_sleepbox(ind,[1:4,6,5]);
        [Decoded_Sequence,~,Decoded_Data_Smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,[],DecodeWin,WinAdvance);
        Decoded_Sequence(:,14)=-i;
        deseqall=[deseqall;Decoded_Sequence];
        dedataall=[dedataall,Decoded_Data_Smooth];
    end
    cd(savedir);save(['sorted_spike_decoding_8arm5_sleepbox'],'deseqall','dedataall','PosProbDist',...
        'maze1d','Position_Transition_Smooth','Position_Transition','BinSigma','maze_cutoff','DecodeWin','WinAdvance',...
        'savedir','Position_Transition_SmoothFactor','-v7.3');
    disp(['Finish Session  ',num2str(in)]);
end
%----------------------------------------------------------------------------------
% 9. decoding with 20 ms window 5ms overlap for info content during trials
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.02;WinAdvance=0.005;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    indv=Position_Data(:,5)>=Velocity_Cutoff;
    PositionTrain=Position_Data(indv,2:3);
    PosProbDist=zeros(size(PositionBasis,1),1);
    for i=1:size(PositionBasis,1)
        PosProbDist(i)=sum(mvnpdf(PositionBasis(i,:),PositionTrain,Sigma_Pos));
    end
    PosProbDist=PosProbDist/sum(PosProbDist);
    PosProbDist0=PosProbDist;
    indnan1d=log10(PosProbDist0)<-maze_cutoff;PosProbDist0(indnan1d)=nan; 
    maze1d=[~indnan1d,PositionBasis];
    maze1d_position=maze1d(maze1d(:,1)>0,2:3);
    Position_Transition=zeros(size(maze1d_position,1));
    for i=1:size(maze1d_position,1)
        indpos= find( Position_Data(1:end-1,5)>=10 & abs(Position_Data(1:end-1,2)-maze1d_position(i,1))<=binsize/2 & abs(Position_Data(1:end-1,3)-maze1d_position(i,2))<=binsize/2 );
        pospost=Position_Data(indpos+1,2:3);
        for j=1:size(maze1d_position,1)
            indpospost= (abs(pospost(:,1)-maze1d_position(j,1))<=binsize/2 & abs(pospost(:,2)-maze1d_position(j,2))<=binsize/2);
            Position_Transition(i,j)=mean(indpospost);
        end
    end
    Position_Transition_Smooth=[];
    for i=1:size(Position_Transition,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=Position_Transition_SmoothFactor*binsize;
        Position_Transition_Smooth(i,:)=nanmean(Position_Transition(ind,:),1);
    end
    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    spikedata=spike_trial; % select for hpc cells
    spikedata(:,7:10)=Position_Data(spikedata(:,7),[5,2,3,4]);
    TrialSet=unique(Position_Data(:,4));
    TrialNum=length(TrialSet); 
    decode_result=[];deseqall=[];dedataall=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7,5]);
        [Decoded_Sequence,~,Decoded_Data_Smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        Decoded_Sequence(:,14)=TrialSet(trialid);
        deseqall=[deseqall;Decoded_Sequence];
        dedataall{trialid}=Decoded_Data_Smooth;
        ind=Decoded_Sequence(:,9)>=10;
        decode_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8,13]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8,13]),1)];
        disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    decode_session(in,:)=mean(decode_result,1);
    PosProbDist_set{in}=PosProbDist;
    maze1d_set{in}=maze1d;
    Position_Transition_Smooth_set{in}=Position_Transition_Smooth;
    Position_Transition_set{in}=Position_Transition;
    cd(savedir);save(['sorted_spike_decoding_8arm6'],'decode_result','deseqall','dedataall','PosProbDist',...
        'maze1d','Position_Transition_Smooth','Position_Transition','BinSigma','maze_cutoff','DecodeWin','WinAdvance',...
        'savedir','Position_Transition_SmoothFactor','-v7.3');
    disp(['Finish Session  ',num2str(in),' ; Performance = ',num2str(decode_session(in,:),2)]);
end
%----------------------------------------------------------------------------------
% 10. generate layer estimation plot
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
layerall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('HPC_Layer_Channel_LFP.mat', 'HPC_Layer_Channel');
    layerall(:,in)=HPC_Layer_Channel(:,2);
end
HPC_layer_name={'dCA1 pyr','dCA1 st rad','dCA1 slm','DG OML','DG MML','DG GCL','dCA3','vCA3','vDG','vCA1 rad','v CA1/Sub'};
figure;
for in=1:8
    a=layerall(:,in*2-1)-layerall(1,in*2-1);
    a(:,2)=layerall(:,in*2)-layerall(1,in*2);
    a=nanmean(a,2);
    if in<=4
        hold on;plot([1:11],a,'ro-');
    else
        hold on;plot([1:11],a,'bo-');
    end
end
xticks([1:11]);xticklabels(HPC_layer_name);xlim([0.5 11.5]);title('Estimation of Dorsal & Ventral HPC Layers on Probe');
ylabel('Depth on Probe Relative to dCA1 Pyr (um)');xlabel('HPC Layers');
%----------------------------------------------------------------------------------
% 11. decode with only dhpc units at 250ms for decoding accuracy
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.25;WinAdvance=DecodeWin;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');

    indv=Position_Data(:,5)>=Velocity_Cutoff;
    PositionTrain=Position_Data(indv,2:3);
    PosProbDist=zeros(size(PositionBasis,1),1);
    for i=1:size(PositionBasis,1)
        PosProbDist(i)=sum(mvnpdf(PositionBasis(i,:),PositionTrain,Sigma_Pos));
    end
    PosProbDist=PosProbDist/sum(PosProbDist);
    PosProbDist0=PosProbDist;
    indnan1d=log10(PosProbDist0)<-maze_cutoff;PosProbDist0(indnan1d)=nan; 
    maze1d=[~indnan1d,PositionBasis];
    maze1d_position=maze1d(maze1d(:,1)>0,2:3);
    Position_Transition=zeros(size(maze1d_position,1));
    for i=1:size(maze1d_position,1)
        indpos= find( Position_Data(1:end-1,5)>=10 & abs(Position_Data(1:end-1,2)-maze1d_position(i,1))<=binsize/2 & abs(Position_Data(1:end-1,3)-maze1d_position(i,2))<=binsize/2 );
        pospost=Position_Data(indpos+1,2:3);
        for j=1:size(maze1d_position,1)
            indpospost= (abs(pospost(:,1)-maze1d_position(j,1))<=binsize/2 & abs(pospost(:,2)-maze1d_position(j,2))<=binsize/2);
            Position_Transition(i,j)=mean(indpospost);
        end
    end
    Position_Transition_Smooth=[];
    for i=1:size(Position_Transition,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=Position_Transition_SmoothFactor*binsize;
        Position_Transition_Smooth(i,:)=nanmean(Position_Transition(ind,:),1);
    end

    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    spikedata=spike_trial; % select for dhpc cells
    lia=ismember(unit_decode(:,16),[1,2,5,6]);
    unit_dhpc=unit_decode(lia,:);
    unit_dhpc(:,18)=[1:size(unit_dhpc,1)];
    Field_Data=Field_Data(:,:,lia);
    lia=ismember(spikedata(:,2),unit_dhpc(:,17));
    spikedata=spikedata(lia,:);spikedata(:,8)=0;
    for i=1:size(unit_dhpc,1)
        ind=spikedata(:,2)==unit_dhpc(i,17);
        spikedata(ind,8)=unit_dhpc(i,18);
    end
    spikedata(:,2)=spikedata(:,8);
    spikedata(:,7:10)=Position_Data(spikedata(:,7),[5,2,3,4]);
    TrialSet=unique(Position_Data(:,4));
    TrialNum=length(TrialSet); 
    decode_result=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7,5]);
        [Decoded_Sequence,~,dedata_smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        ind=Decoded_Sequence(:,9)>=10;
        decode_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8,13]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8,13]),1)];
        disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    decode_session(in,:)=mean(decode_result,1);
    PosProbDist_set{in}=PosProbDist;
    maze1d_set{in}=maze1d;
    Position_Transition_Smooth_set{in}=Position_Transition_Smooth;
    Position_Transition_set{in}=Position_Transition;
    
    % below to shuffle unit id
    for shuffle=1:0%100
    p=randperm(size(Field_Data,3));
    shuffle_result=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7]);
        [Decoded_Sequence,~]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data(:,:,p),maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        ind=Decoded_Sequence(:,9)>=10;
        shuffle_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8]),1)];
        %disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    shuffle_result_session(shuffle,:,in)=mean(shuffle_result,1);
    end
    disp(['Finish Session  ',num2str(in),' ; Performance = ',num2str(decode_session(in,:),2)]);
    %disp(['Shuffle : ',num2str(mean(shuffle_result_session(:,:,in),1),2)]);
end
cd(homedir);save('sorted_spike_decoding_8arm_dhpc_trial_sleepbox','Position_Transition_set','Position_Transition_Smooth_set','maze1d_set',...
    "PosProbDist_set",'decode_session','BinSigma','maze_cutoff','DecodeWin','WinAdvance','Position_Transition_SmoothFactor',...
    "shuffle_result_session",'SessionSet16');
%-----------------------------------------------------------------------------------
% 12. decode with only vhpc units at 250ms for decoding accuracy
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.25;WinAdvance=DecodeWin;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');

    indv=Position_Data(:,5)>=Velocity_Cutoff;
    PositionTrain=Position_Data(indv,2:3);
    PosProbDist=zeros(size(PositionBasis,1),1);
    for i=1:size(PositionBasis,1)
        PosProbDist(i)=sum(mvnpdf(PositionBasis(i,:),PositionTrain,Sigma_Pos));
    end
    PosProbDist=PosProbDist/sum(PosProbDist);
    PosProbDist0=PosProbDist;
    indnan1d=log10(PosProbDist0)<-maze_cutoff;PosProbDist0(indnan1d)=nan; 
    maze1d=[~indnan1d,PositionBasis];
    maze1d_position=maze1d(maze1d(:,1)>0,2:3);
    Position_Transition=zeros(size(maze1d_position,1));
    for i=1:size(maze1d_position,1)
        indpos= find( Position_Data(1:end-1,5)>=10 & abs(Position_Data(1:end-1,2)-maze1d_position(i,1))<=binsize/2 & abs(Position_Data(1:end-1,3)-maze1d_position(i,2))<=binsize/2 );
        pospost=Position_Data(indpos+1,2:3);
        for j=1:size(maze1d_position,1)
            indpospost= (abs(pospost(:,1)-maze1d_position(j,1))<=binsize/2 & abs(pospost(:,2)-maze1d_position(j,2))<=binsize/2);
            Position_Transition(i,j)=mean(indpospost);
        end
    end
    Position_Transition_Smooth=[];
    for i=1:size(Position_Transition,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=Position_Transition_SmoothFactor*binsize;
        Position_Transition_Smooth(i,:)=nanmean(Position_Transition(ind,:),1);
    end

    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    spikedata=spike_trial; % select for vhpc cells
    lia=ismember(unit_decode(:,16),[3,4,7,8]);
    unit_dhpc=unit_decode(lia,:);
    unit_dhpc(:,18)=[1:size(unit_dhpc,1)];
    Field_Data=Field_Data(:,:,lia);
    lia=ismember(spikedata(:,2),unit_dhpc(:,17));
    spikedata=spikedata(lia,:);spikedata(:,8)=0;
    for i=1:size(unit_dhpc,1)
        ind=spikedata(:,2)==unit_dhpc(i,17);
        spikedata(ind,8)=unit_dhpc(i,18);
    end
    spikedata(:,2)=spikedata(:,8);
    spikedata(:,7:10)=Position_Data(spikedata(:,7),[5,2,3,4]);
    TrialSet=unique(Position_Data(:,4));
    TrialNum=length(TrialSet); 
    decode_result=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7,5]);
        [Decoded_Sequence,~,dedata_smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        ind=Decoded_Sequence(:,9)>=10;
        decode_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8,13]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8,13]),1)];
        disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    decode_session(in,:)=mean(decode_result,1);
    PosProbDist_set{in}=PosProbDist;
    maze1d_set{in}=maze1d;
    Position_Transition_Smooth_set{in}=Position_Transition_Smooth;
    Position_Transition_set{in}=Position_Transition;
    
    % below to shuffle unit id
    for shuffle=1:0%100
    p=randperm(size(Field_Data,3));
    shuffle_result=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7]);
        [Decoded_Sequence,~]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data(:,:,p),maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        ind=Decoded_Sequence(:,9)>=10;
        shuffle_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8]),1)];
        %disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    shuffle_result_session(shuffle,:,in)=mean(shuffle_result,1);
    end
    disp(['Finish Session  ',num2str(in),' ; Performance = ',num2str(decode_session(in,:),2)]);
    %disp(['Shuffle : ',num2str(mean(shuffle_result_session(:,:,in),1),2)]);
end
cd(homedir);save('sorted_spike_decoding_8arm_vhpc_trial_sleepbox','Position_Transition_set','Position_Transition_Smooth_set','maze1d_set',...
    "PosProbDist_set",'decode_session','BinSigma','maze_cutoff','DecodeWin','WinAdvance','Position_Transition_SmoothFactor',...
    "shuffle_result_session",'SessionSet16');
%-----------------------------------------------------------------------------------
% 13. dHPC decoding with 20 ms window 5ms overlap for info content during trials
close all;clear all;set(0,'DefaultFigureColormap',feval('turbo'));
% 0. common parameters
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
binsize=2;binmax=109;bin_edge=binmax+binsize/2;
maze_cutoff=4;Velocity_Cutoff=10;
x = -binmax:binsize:binmax;
y = -binmax:binsize:binmax;
[X,Y] = meshgrid(x,y);
PositionBasis=[X(:),Y(:)]; % x,y
x_edge=[-bin_edge:binsize:bin_edge];
y_edge=[-bin_edge:binsize:bin_edge];
BinSigma=40;Sigma_Pos=[BinSigma,0;0,BinSigma]; 
DecodeWin=0.02;WinAdvance=0.005;Position_Transition_SmoothFactor=8;shuffle_result_session=[];
decode_session=[];maze1d_set=[];Position_Transition_Smooth_set=[];Position_Transition_set=[];PosProbDist_set=[];
for in=14:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    indv=Position_Data(:,5)>=Velocity_Cutoff;
    PositionTrain=Position_Data(indv,2:3);
    PosProbDist=zeros(size(PositionBasis,1),1);
    for i=1:size(PositionBasis,1)
        PosProbDist(i)=sum(mvnpdf(PositionBasis(i,:),PositionTrain,Sigma_Pos));
    end
    PosProbDist=PosProbDist/sum(PosProbDist);
    PosProbDist0=PosProbDist;
    indnan1d=log10(PosProbDist0)<-maze_cutoff;PosProbDist0(indnan1d)=nan; 
    maze1d=[~indnan1d,PositionBasis];
    maze1d_position=maze1d(maze1d(:,1)>0,2:3);
    Position_Transition=zeros(size(maze1d_position,1));
    for i=1:size(maze1d_position,1)
        indpos= find( Position_Data(1:end-1,5)>=10 & abs(Position_Data(1:end-1,2)-maze1d_position(i,1))<=binsize/2 & abs(Position_Data(1:end-1,3)-maze1d_position(i,2))<=binsize/2 );
        pospost=Position_Data(indpos+1,2:3);
        for j=1:size(maze1d_position,1)
            indpospost= (abs(pospost(:,1)-maze1d_position(j,1))<=binsize/2 & abs(pospost(:,2)-maze1d_position(j,2))<=binsize/2);
            Position_Transition(i,j)=mean(indpospost);
        end
    end
    Position_Transition_Smooth=[];
    for i=1:size(Position_Transition,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=Position_Transition_SmoothFactor*binsize;
        Position_Transition_Smooth(i,:)=nanmean(Position_Transition(ind,:),1);
    end
    load('Spike_Data_decode_trial_sleepbox.mat','unit_decode','spike_trial');
    load('Field_Data_trial_sleepbox.mat','Field_Data');
    spikedata=spike_trial; % select for dhpc cells
    lia=ismember(unit_decode(:,16),[1,2,5,6]);
    unit_dhpc=unit_decode(lia,:);
    unit_dhpc(:,18)=[1:size(unit_dhpc,1)];
    Field_Data=Field_Data(:,:,lia);
    lia=ismember(spikedata(:,2),unit_dhpc(:,17));
    spikedata=spikedata(lia,:);spikedata(:,8)=0;
    for i=1:size(unit_dhpc,1)
        ind=spikedata(:,2)==unit_dhpc(i,17);
        spikedata(ind,8)=unit_dhpc(i,18);
    end
    spikedata(:,2)=spikedata(:,8);
    spikedata(:,7:10)=Position_Data(spikedata(:,7),[5,2,3,4]);
    TrialSet=unique(Position_Data(:,4));
    TrialNum=length(TrialSet); 
    decode_result=[];deseqall=[];dedataall=[];
    for trialid=1:length(TrialSet)
        lia=ismember(spikedata(:,10),TrialSet(trialid));% trial id
        Spike_Data=spikedata(lia,[1,2,8,9,7,5]);
        [Decoded_Sequence,~,Decoded_Data_Smooth]=BAYESIAN_DECODING_8ARM1D_Continuous(Spike_Data,Field_Data,maze1d,Position_Transition_Smooth,Position_Data,DecodeWin,WinAdvance);
        Decoded_Sequence(:,14)=TrialSet(trialid);
        deseqall=[deseqall;Decoded_Sequence];
        dedataall{trialid}=Decoded_Data_Smooth;
        ind=Decoded_Sequence(:,9)>=10;
        decode_result(trialid,:)=[size(Decoded_Sequence,1),nanmean(Decoded_Sequence(:,[8,13]),1),mean(ind),nanmean(Decoded_Sequence(ind,[8,13]),1)];
        disp(['Finish Trial  ',num2str(trialid),' ; Error = ',num2str(decode_result(trialid,:),2)]);
    end
    decode_session(in,:)=mean(decode_result,1);
    PosProbDist_set{in}=PosProbDist;
    maze1d_set{in}=maze1d;
    Position_Transition_Smooth_set{in}=Position_Transition_Smooth;
    Position_Transition_set{in}=Position_Transition;
    cd(savedir);save(['sorted_spike_decoding_8arm_dHPC'],'decode_result','deseqall','dedataall','PosProbDist',...
        'maze1d','Position_Transition_Smooth','Position_Transition','BinSigma','maze_cutoff','DecodeWin','WinAdvance',...
        'savedir','Position_Transition_SmoothFactor','-v7.3');
    disp(['Finish Session  ',num2str(in),' ; Performance = ',num2str(decode_session(in,:),2)]);
end
%-------------------------------------------------------------------------------
% 14. adjust dca1 channel according to visual inspection of
% 'dHPC_layer7channel3', 08/17/2024
clear all;close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
dHPC_layer7channel4=[];sessionchange=[-1,0,0,0,-1,-1,-1,0,0,0,-1,-1,-1,-1,-1,0];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('HPC_Layer_Channel_LFP.mat','dHPC_layer7channel3');
    load('LFP_CSD.mat', 'channelcsd', 'lfpcsd_ref','lfptime');
    load('CSD_LFP_Visualization.mat', 'channels');
    dHPC_layer7channel3(1,1)=dHPC_layer7channel3(1,1)+sessionchange(in);
    dHPC_layer7channel3(:,2:3)=channelcsd(dHPC_layer7channel3(:,1),3:4);
    dHPC_layer7channel3(:,4)=0;
    for i=1:7
        dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2) & channels(:,4)==dHPC_layer7channel3(i,3));
    end
    LFP_HPC=lfpcsd_ref(:,dHPC_layer7channel3(:,1));
    dHPC_layer7channel4(:,:,in)=dHPC_layer7channel3;
    save('dHPC_layer7channel3_ca1_adjusted_LFP','LFP_HPC','lfptime','dHPC_layer7channel3','-v7.3');
end
cd(homedir);save('dHPC_layer7channel3_ca1_adjusted',"dHPC_layer7channel4",'sessionchange');

clear all;close all;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
dHPC_layer7channel5=[];sessionchange=[0,0,-1,0,0,0,1,0,0,-1,1,0,1,1,0,0];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('dHPC_layer7channel3_ca1_adjusted_LFP','dHPC_layer7channel3');
    load('LFP_CSD.mat', 'channelcsd', 'lfpcsd_ref','lfptime');
    load('CSD_LFP_Visualization.mat', 'channels');
    dHPC_layer7channel3(1,1)=dHPC_layer7channel3(1,1)+sessionchange(in);
    dHPC_layer7channel3(:,2:3)=channelcsd(dHPC_layer7channel3(:,1),3:4);
    dHPC_layer7channel3(:,4)=0;
    for i=1:7
        dHPC_layer7channel3(i,4)=find(channels(:,3)==dHPC_layer7channel3(i,2) & channels(:,4)==dHPC_layer7channel3(i,3));
    end
    LFP_HPC=lfpcsd_ref(:,dHPC_layer7channel3(:,1));
    dHPC_layer7channel5(:,:,in)=dHPC_layer7channel3;
    save('dHPC_layer7channel3_ca1_adjusted_LFP_v2','LFP_HPC','lfptime','dHPC_layer7channel3','-v7.3');
end
cd(homedir);save('dHPC_layer7channel3_ca1_adjusted_v2',"dHPC_layer7channel5",'sessionchange');
%------------------------------------------------------------------------------------
% 15. new identification of behavior epochs, simply location to center
reward_dis_cut=70;pos_behave_marker_colname={'1.old marker','2.old pkid','3.new marker','4.new pkid'};as=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);  
    load('InOutBound_Behavior_Analysis.mat','pos_behave_marker','pksets','inter_run');
    load('Position_Data_Maze.mat');
    Position_Data(:,6)=vecnorm(Position_Data(:,2:3),2,2);
    Position_Data(1:end-1,7)=diff(Position_Data(:,6));
    Position_Data(end,7)=Position_Data(end-1,7);
    ind=abs(pos_behave_marker(:,1))==0.5;
    pos_behave_marker(ind,2)=inter_run(pos_behave_marker(ind,2),2);
    pos_behave_marker(:,3)=0;
    ind=Position_Data(:,6)>reward_dis_cut & pos_behave_marker(:,1)~=-100;pos_behave_marker(ind,3)=1;
    ind=Position_Data(:,6)<=30 & pos_behave_marker(:,1)~=-100;pos_behave_marker(ind,3)=2;
    ind=Position_Data(:,6)>30 & Position_Data(:,6)<=reward_dis_cut & pos_behave_marker(:,1)==-1;pos_behave_marker(ind,3)=3;
    ind=Position_Data(:,6)>30 & Position_Data(:,6)<=reward_dis_cut & pos_behave_marker(:,1)==1;pos_behave_marker(ind,3)=4;
    pos_behave_marker(:,4)=pos_behave_marker(:,2);
    for i=1:size(pksets,1)-1
        if pksets(i,1)==pksets(i+1,1)
            ind=Position_Data(:,1)>pksets(i,4)&Position_Data(:,1)<pksets(i+1,4)&pos_behave_marker(:,3)==2;
            pos_behave_marker(ind,4)=i;
        else
            ind=Position_Data(:,1)>pksets(i,4) & Position_Data(:,1)<pksets(i+1,4) & pos_behave_marker(:,3)==2 & Position_Data(:,4)==pksets(i,1);
            pos_behave_marker(ind,4)=-i;
            ind=Position_Data(:,1)>pksets(i,4) & Position_Data(:,1)<pksets(i+1,4) & pos_behave_marker(:,3)==2 & Position_Data(:,4)==pksets(i+1,1);
            pos_behave_marker(ind,4)=-i-1;
        end
    end
    ind=Position_Data(:,1)<pksets(1,4)&pos_behave_marker(:,3)==2;pos_behave_marker(ind,4)=-1;
    ind=Position_Data(:,1)>pksets(end,4)&pos_behave_marker(:,3)==2;pos_behave_marker(ind,4)=-size(pksets,1);
    pos_behave_marker(:,6:7)=Position_Data(:,6:7);
    cd(savedir);save('new_reward_center_in_out','pos_behave_marker','pksets','inter_run','reward_dis_cut','pos_behave_marker_colname');
    pos_behave_marker(:,8)=in;
    as=[as;pos_behave_marker];
    if 0
    pos_behave_marker(:,5)=pos_behave_marker(:,2)-pos_behave_marker(:,4);
    
    ind1=pos_behave_marker(:,3)==2&pos_behave_marker(:,1)==-1;
    ind2=pos_behave_marker(:,3)==2&abs(pos_behave_marker(:,1))==0.5;
    ind3=pos_behave_marker(:,3)==2&pos_behave_marker(:,1)==1;
    subplot(4,4,in);histogram(pos_behave_marker(ind1,5));
    hold on;histogram(pos_behave_marker(ind2,5));
    hold on;histogram(pos_behave_marker(ind3,5));
    end
end
%-----------------------------------------------------------------------------------