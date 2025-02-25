function CALCULATE_PLACE_FIELDS_8ARMS_NPX_move_immobile

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 
% This function calculates the locations within the behavioral arena or
% track where the neurons fire.  It will first calculate 2-dimensional
% place fields.  If the track is linear (Track_Type==1), the program will
% also calculate linear (uni-directional and bi-directional) place fields.
% 
% Bin_Size is the size (in cm) of the bins that are used to make the place
% fields.  Velocity_Cutoff is the velocity (in cm/sec) that the rat must be
% moving in order for the spikes to count toward place fields.
% Firing_Rate_Cutoff is the minimum peak firing rate (in Hz) that a cell
% must have to be considered to have a place field.  Otherwise, its place
% field is eliminated from future analysis.  Analyze_Linear (1 or 0) tells
% the program whether or not to calculate linear place fields (1=yes,0=no).
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% If the user did not fully define the input variables, the most likely
% values are either searched for or automatically provided.
position0=110;Bin_Size=2;Velocity_Cutoff=10;Firing_Rate_Cutoff=0.001;dt=0.0335;
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
% This section calculates the 2-dimensional place fields
for in=1:16
savedir=SessionSet16{in};
cd(savedir);
load('Ripple_Events_PAPER5.mat','ripples');
load('Spike_Data_decode_trial_sleepbox.mat', 'spike_trial','unit_decode');
spike_trial(:,8)=1;
for i=1:size(ripples,1)
    ind= spike_trial(:,1)>=ripples(i,1) & spike_trial(:,1)<=ripples(i,2) ;
    spike_trial(ind,8)=0;
end
spike_trial=spike_trial(spike_trial(:,8)==1,:);
spike_trial(:,5)=spike_trial(:,7);spikedata=spike_trial;
unit_metric_decode=unit_decode;
load Position_Data_Maze.mat
ind1=ripples(:,10)<=10 & ripples(:,2)-ripples(:,1)>=0.035 & ripples(:,6)==1;
ripples=ripples(ind1,:);

x=Position_Data(:,2);
y=Position_Data(:,3);
[Velocity,vel_angle,locdis,loc_angle,locdis_change]=calculate_location_movement(Position_Data);
locdis=smoothdata(locdis,'gaussian',50);

% This eliminates all position and spike data where the rat was moving
% slower than the Velocity_Cutoff. I then bin the position information.
positionbin=ceil(([x,y]+position0)/Bin_Size);
binnum=ceil(position0*2/Bin_Size);
positionbin_ripple=ceil((ripples(:,19:20)+position0)/Bin_Size);

% Here, the total time spent moving in each position bin is calculated
Time_In_Position=zeros(binnum,binnum,4);
for i=1:binnum
    for j=1:binnum
        ind= positionbin(:,1)==i & positionbin(:,2)==j & Velocity>=Velocity_Cutoff ;
        Time_In_Position(j,i,1)=sum(ind)*dt;
        ind= positionbin(:,1)==i & positionbin(:,2)==j & Velocity<Velocity_Cutoff ;
        indr= positionbin_ripple(:,1)==i & positionbin_ripple(:,2)==j & ripples(:,10)<Velocity_Cutoff ;
        Time_In_Position(j,i,2)=sum(ind)*dt-nansum(ripples(indr,2)-ripples(indr,1));
        ind= positionbin(:,1)==i & positionbin(:,2)==j & Velocity<5 ;
        indr= positionbin_ripple(:,1)==i & positionbin_ripple(:,2)==j & ripples(:,10)<5 ;
        Time_In_Position(j,i,3)=sum(ind)*dt-nansum(ripples(indr,2)-ripples(indr,1));
        ind= positionbin(:,1)==i & positionbin(:,2)==j & Velocity<=1 ;
        indr= positionbin_ripple(:,1)==i & positionbin_ripple(:,2)==j & ripples(:,10)<=1 ;
        Time_In_Position(j,i,4)=sum(ind)*dt-nansum(ripples(indr,2)-ripples(indr,1));
    end
end

% Here, the number of spikes in each position bin is calculated for each cell
spikedata(:,7)=positionbin(spikedata(:,5),1);
spikedata(:,8)=positionbin(spikedata(:,5),2);
ind=min(spikedata(:,7:8),[],2)>0;
spikedata=spikedata(ind,:);
spikedata(:,9)=Velocity(spikedata(:,5));
Spikes_In_Position=zeros(binnum,binnum,size(unit_metric_decode,1),4);
for N=1:size(spikedata,1)
    if spikedata(N,9)>=Velocity_Cutoff
        Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),1)=Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),1)+1;
    end
    if spikedata(N,9)<Velocity_Cutoff
        Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),2)=Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),2)+1;
    end
    if spikedata(N,9)<5
        Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),3)=Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),3)+1;
    end
     if spikedata(N,9)<=1
        Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),4)=Spikes_In_Position(spikedata(N,8),spikedata(N,7),spikedata(N,2),4)+1;
    end
end

% Here, I calculate the actual firing rate (spikes/time) and save it as Field_Data
Firing_Rate_In_Position=zeros(size(Spikes_In_Position));
for N=1:size(Spikes_In_Position,3)
    for v=1:4
        Firing_Rate_In_Position(:,:,N,v)=Spikes_In_Position(:,:,N,v)./Time_In_Position(:,:,v);
    end
end
Firing_Rate_In_Position( isnan(Firing_Rate_In_Position) )=0;
Firing_Rate_In_Position( isinf(Firing_Rate_In_Position) )=0;
Smoothed_Firing_Rate=Firing_Rate_In_Position;
Field_Data_raw=Firing_Rate_In_Position;
Two_D_Filter=fspecial('gaussian',[20 20]*3,2);  %This is a gaussian filter with a kernel St. Dev. of 2 bins (4 cm) that filters out to 10 bins in all directions (20 cm)
for N=1:size(Firing_Rate_In_Position,3)
    for v=1:4
        Smoothed_Firing_Rate(:,:,N,v)=filter2(Two_D_Filter,Firing_Rate_In_Position(:,:,N,v));
    end
end
%Smoothed_Firing_Rate(find(isnan(Smoothed_Firing_Rate)))=0;%Smoothed_Firing_Rate(find(Smoothed_Firing_Rate<0))=0;
Field_Data=Smoothed_Firing_Rate;

% below calculate basic properties of place fields
field_property=[];
for N=1:size(Firing_Rate_In_Position,3)
    for v=1:4
        a=Field_Data(:,:,N,v);
        field_property(N,1,v)=max(a(:)); % peak field rate
        field_property(N,2,v)=sum(a(:)>max(a(:))*0.1)*4; % field size: cm^2
        field_property(N,3,v)=1-mean(a(:)>=0); % proportion of abnormal value
    end
end
% peak field rate| field size: cm^2 | proportion of abnormal value | info
% rate
% below calculate information rate for each cell
% Here, the total time spent moving in each position bin is calculated
spatial_info=[];
for v=1:4
probx =Time_In_Position(:,:,v)/sum(sum(Time_In_Position(:,:,v),1),2);
for i=1:size(Field_Data,3)
    A=squeeze(Field_Data(:,:,i,v));
    clear F;F=sum(sum(probx.*A));
    if F>0
        a=log2(A/F);a((A==0))=0;
        spatial_info(i,v)=sum(sum(probx.*A.*a))/F;
    else
        spatial_info(i,v)=nan;
    end
end
end
field_property(:,4,:)=spatial_info;
cd(savedir);save('Field_Data_move_immobile','Field_Data','Field_Data_raw','field_property','unit_decode','spikedata','-v7.3');
disp(['Place field analysis completed for ',num2str([in,size(Field_Data,3)]),' Units!']);
clear spikedata unit_decode Field_Data Field_Data_raw field_property spike_trial ripples Position_Data
end
if 0
    load("Spike_Data.mat",'spike');
    figure;
for i=1:size(Field_Data,3)
    subplot(1,3,1);
    imagesc(Field_Data(:,:,i));
    title(num2str([i,field_property(i,:)],2));
    subplot(1,3,2);plot(unit_waveform_decode(i,:));
    ind=spike(:,6)==i;
    subplot(1,3,3);histogram(spike(ind,4));title(unit_metric_decode(i,10));
    pause(0.2);clf;
end
end