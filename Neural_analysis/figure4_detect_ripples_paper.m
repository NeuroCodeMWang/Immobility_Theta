function figure4_detect_ripples_paper(in)

homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
load('dHPC_layer7channel3_ca1_adjusted',"dHPC_layer7channel4");
savedir=SessionSet16{in};
cd(savedir);
load('HPC_Layer_Channel_LFP.mat', 'LFP_HPC','lfptime');
load('LFP_CSD.mat', 'lfpcsd_ref');
dHPC_layer7channel3=squeeze(dHPC_layer7channel4(:,:,in));
% 0. Detection of SWR
LFP_Frequency=625;
LFP_Ripple1=[lfptime(:,1),lfpcsd_ref(:,dHPC_layer7channel3(1,1)-1:dHPC_layer7channel3(1,1)+1)];%
LFP_SPW=[lfptime(:,1),LFP_HPC(:,2)];
LFP_DGCA3=[lfptime(:,1),LFP_HPC(:,6:7)];

[CA1_Ripple_Events,Ripple_Amplitude_CA1,Ripple_Filtered_CA1,Ripple_Amplitude_LFP_CA1,Ripple_Filtered_LFP_CA1]=ripple_filter_v3(LFP_Ripple1,LFP_Frequency,[],3);% >15ms
[SPW_Amplitude,SPW_Events,SPW_Filtered_LFP_mean]=SharpWave_filter(LFP_SPW,LFP_Frequency); % 2 sd cutoff
[CA3_Ripple_Events,Ripple_Amplitude_DGCA3,Ripple_Filtered_DGCA3,Ripple_Amplitude_LFP_DGCA3,Ripple_Filtered_LFP_DGCA3]=ripple_filter_v3(LFP_DGCA3,LFP_Frequency,[],3); % 3 std

CA1_Ripple_Events(:,6)=0;timewin=0;
for i=1:size(CA1_Ripple_Events,1)
    ind= SPW_Events(:,1)<=CA1_Ripple_Events(i,2)+timewin & SPW_Events(:,1)>=CA1_Ripple_Events(i,1)-timewin;
    CA1_Ripple_Events(i,6)=sum(ind);
end
indr= CA1_Ripple_Events(:,6)>=1 ;
disp(['Total CA1 Ripples : ',num2str(size(CA1_Ripple_Events,1)),' ; with SPW : ',num2str(mean(indr)*100,2),' % ']);
CA1ripple=CA1_Ripple_Events(indr,:);
CA1ripple(:,6)=1;CA3_Ripple_Events(:,6)=2;
ripples_sleep_maze=[CA1ripple;CA3_Ripple_Events];
[~,I]=sort(ripples_sleep_maze(:,3));
ripples_sleep_maze=ripples_sleep_maze(I,:);
ruler=lfptime(:,1);template=ripples_sleep_maze(:,3);[outputindex,error]=match(template,ruler,0);
ripples_sleep_maze(:,9)=lfptime(outputindex,3);% sleepbox or velocity
ripples_sleep_maze=ripples_sleep_maze(ripples_sleep_maze(:,9)<10,:);
ripples=ripples_sleep_maze( ripples_sleep_maze(:,9)>=0,:);

load('InOutBound_Behavior_Analysis.mat','pos_behave_marker');
load('Position_Data_Maze.mat');
load('Analyze_Theta_Cycle_Broad','reftheta');
ruler=Position_Data(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
ripples(:,[19:20,9:12])=[Position_Data(outputindex,2:5),pos_behave_marker(outputindex,:)];
ripples(:,13)=vecnorm(ripples(:,19:20),2,2);
ripples(:,14)=in;
ind=ripples(:,10)<10 & ripples(:,11)~=-100 & ripples(:,12)~=-100;ripples=ripples(ind,:);
ruler=reftheta(:,1);template=ripples(:,3);[outputindex,error]=match(template,ruler,0);
ripples(:,16:17)=reftheta(outputindex,3:4); % theta power | theta phase
cd(savedir);save('Ripple_Events_PAPER5.mat','ripples','ripples_sleep_maze','CA3_Ripple_Events', ...
    'Ripple_Amplitude_DGCA3',"CA1_Ripple_Events","Ripple_Amplitude_CA1",'SPW_Events',...
    'Ripple_Filtered_CA1','Ripple_Amplitude_LFP_CA1','Ripple_Filtered_LFP_CA1',...
    'Ripple_Filtered_DGCA3','Ripple_Amplitude_LFP_DGCA3','Ripple_Filtered_LFP_DGCA3','-v7.3');
[in size(ripples,1) ]