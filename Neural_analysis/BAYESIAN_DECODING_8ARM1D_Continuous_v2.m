function [Decoded_Sequence,Decoded_Data,Decoded_Data_Smooth]=BAYESIAN_DECODING_8ARM1D_Continuous_v2(Spike_Data,Field_Data,maze1d,Position_Transition,Position_Data,Time_Window,Time_Advance,Start_Time,End_Time)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 
% This function uses the Bayesian decoding algorithm from [Davidson et al,
% 2009] and decodes the information contained within a spike train.
% 
% Spike_Data is the list of spike times and cell IDs that you want to
% decode.  The Spike_Data data needs to be a two-dimensional array where
% the first column is the time of each spike and the second column is the
% cell ID.  The function LOAD_CLUSTERED_SPIKE_DATA generates spike data
% that is appropriate for the current function.
% 
% Field_Data is the firing rate fields of the cells included in Spike_Train.
% Field_Data is assumed to have either one-dimensional or two-dimensional
% fields (if you are decoding three-dimensional or higher, you'll need to
% re-write this code).  If the fields are one-dimensional (such as linear
% track place fields or head direction fields), Field_Data should be a
% two-dimensional array, with the first dimension (rows) being the firing
% rate at every bin and the second dimension (columns) being the
% corresponding cell ID.  If the fields are two-dimensional (such as open
% arena place fields), Field_Data should be a three-dimensional array, with the
% first two dimensions being the X- and Y-positions and the third
% dimensions being the corresponding cell ID.  The function
% CALCULATE_PLACE_FIELDS generates fields that are appropriate for the
% current function.
% 
% Time_Window is a single value that lists the time (in seconds) that will
% be decoded on every pass.  Time_Advance is a single value that lists the
% time (in seconds) that the Time_Window will be advanced.  So if you
% wanted your window of decoding to be 20 ms, and you wanted that to
% advance at 5 ms through the Spike_Train, you'd make Time_Window==0.02 and
% Time_Advance==0.005.  
% 
% Decoded_Sequence lists the starting timepoints (in seconds) and the peak
% decoded locations (in either one or two dimensions) for each window of
% decoding.  Decoded_Data is the actual complete posterior probabilities
% for all of the decoding windows. 
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% If the user did not fully define the input variables, the most likely
% values are either searched for or automatically provided.
if nargin<2
    error('ERROR! User must provide both a spike train and fields for the decoding analysis.')
elseif nargin<3
    Time_Window=0.02;
    Time_Advance=0.005;
elseif nargin<4
    Time_Advance=0.005;
end
if (max(Spike_Data(:,1))-min(Spike_Data(:,1)))<Time_Window,
    warning('The window of decoding is larger than the entirety of the provided spike data.')
end
if Time_Window<=0 || Time_Advance<=0 || length(Time_Window)>1 || length(Time_Advance)>1,
    error('The third and fourth inputs must not be a matrix and must have positive values.')
end
if max(Spike_Data(:,1))-min(Spike_Data(:,1))<0,
    error('The provided spike data does not span any time.')
end

% First, because the decoding algorithm multiplies fields together for a
% part of the calculation, any field with a value of 0 prevents the
% algorithm from ever representing that location.  Therefore, I need to go
% through and replace all values of 0 with a very low, non-zero number.
Field_Data2=Field_Data;
if size(Field_Data,3)>1 %for two-dimensional fields
    for N=1:size(Field_Data2,3)
        Field=Field_Data2(:,:,N);
        Minimum=min(min(Field((Field>0))));
        if ~isempty(Minimum)
            if Minimum/10>0
                Field((Field<=0))=Minimum/10;
            else
                Field((Field<=0))=Minimum;
            end
        else
            Field(:)=1;  % If a neuron has a firing rate of 0 throughout the entire environment, I change it to 1 to effectively eliminate it from the product analysis.
        end
        Field_Data2(:,:,N)=Field;
        clear Minimum;
        clear Field;
    end
elseif size(Field_Data,3)==1 %for one-dimensional fields
    for N=1:size(Field_Data2,2)
        Field=Field_Data2(:,N);
        Minimum=min(Field((Field>0)));
        if ~isempty(Minimum)
            if Minimum/10>0
                Field((Field<=0))=Minimum/10;
            else
                Field((Field<=0))=Minimum;
            end
        else
            Field(:)=1;  % If a neuron has a firing rate of 0 throughout the entire environment, I change it to 1 to effectively eliminate it from the product analysis.
        end
        Field_Data2(:,N)=Field;
        clear Minimum;
        clear Field;
    end
end

if size(Field_Data,3)>1
    Field_Data20=Field_Data2;Field_Data2=[];
    for N=1:size(Field_Data,3)
        A=Field_Data20(:,:,N);
        Field_Data2(:,N)=A(:);
    end
    Field_Data2=Field_Data2(maze1d(:,1)>0,:);
end
maze1d_position=maze1d(maze1d(:,1)>0,2:3);
pospre0=ones(size(Position_Transition,1),1)/size(Position_Transition,1);
% Here, I pre-allocate a variable to store the decoded data
t=[Start_Time:Time_Advance:End_Time];
% Here is where the decoding actual happens
Decoded_Data=zeros(size(Position_Transition,1),length(t));
Decoded_Sequence=zeros(length(t),13);
decodespikenum=zeros(length(t),9);
for Line=1:length(t)
    Start_Time=t(Line);
    Subset_Spike_Data=Spike_Data((Spike_Data(:,1)>=Start_Time & Spike_Data(:,1)<(Start_Time+Time_Window)),:);
    if size(Field_Data,3)>1
        if ~isempty(Subset_Spike_Data)
            decodespikenum(Line,1:10)=[size(Subset_Spike_Data,1),length(unique(Subset_Spike_Data(:,2))),histcounts(Subset_Spike_Data(:,3),[0.5:8.5])]; 
            Decoded_Matrix=prod(Field_Data2(:,Subset_Spike_Data(:,2)),2).*exp(-Time_Window*sum(Field_Data2,2));      
            Divider=1;
            while isinf(max(Decoded_Matrix(:)))
                Divider=Divider*2;
                Decoded_Matrix=prod((Field_Data2(:,Subset_Spike_Data(:,2))/Divider),2).*exp(-Time_Window*sum((Field_Data2/Divider),2));
            end
            clear Divider;
        else
            Decoded_Matrix=exp(-Time_Window*sum(Field_Data2,2));
        end
        Decoded_Matrix((Decoded_Matrix<0))=0;
        if max(Decoded_Matrix)>0
            Decoded_Matrix=Decoded_Matrix/sum(Decoded_Matrix);
        end
        if Line==1
            pospre=pospre0;
        else
            pospre=Decoded_Data(:,Line-1);
        end
        Decoded_Matrix=Decoded_Matrix.*(Position_Transition'*pospre);
        if max(Decoded_Matrix)>0
            Decoded_Matrix=Decoded_Matrix/sum(Decoded_Matrix);
            [~,I]=max(Decoded_Matrix);Max_X_Position=maze1d_position(I,1);Max_Y_Position=maze1d_position(I,2);
        else
            Max_Y_Position=0;
            Max_X_Position=0;
        end
        Decoded_Sequence(Line,[1:3])=[Start_Time+Time_Window/2,Max_X_Position,Max_Y_Position];%,mean(Subset_Spike_Data(:,3:5),1)];
        Decoded_Data(:,Line)=Decoded_Matrix;
    elseif size(Field_Data,3)==1
            Decoded_Matrix=prod(Field_Data2(:,Subset_Spike_Data(:,2)),2).*exp(-Time_Window*sum(Field_Data,2));
            Divider=1;
            while isinf(max(Decoded_Matrix))
                Divider=Divider*2;
                Decoded_Matrix=prod((Field_Data2(:,Subset_Spike_Data(:,2))/Divider),2).*exp(-Time_Window*sum((Field_Data/Divider),2));
            end
            clear Divider;
            Decoded_Matrix((Decoded_Matrix<0))=0;
            if max(max(Decoded_Matrix))>0
                Decoded_Matrix=Decoded_Matrix/sum(Decoded_Matrix);
            end
            if max(max(Decoded_Matrix))>0
                Max_Position=find(Decoded_Matrix==max(Decoded_Matrix),1,'first');
            else
                Max_Position=0;
            end
            Decoded_Sequence(Line,:)=[Start_Time+Time_Window/2,Max_Position];
            if nargout==2
                Decoded_Data(:,Line)=Decoded_Matrix;
            end
    end
end
ind=decodespikenum(:,1)==0;
Decoded_Sequence(ind,2:end)=nan;
Decoded_Data(:,ind)=nan;
Decoded_Sequence(:,[6:7,19:26])=decodespikenum;
if ~isempty(Position_Data)
    ruler=Position_Data(:,1);template=Decoded_Sequence(:,1);[outputindex,error]=match(template,ruler,0);
    Decoded_Sequence(:,[4:5,9:10])=[Position_Data(outputindex,[2:3,5]),error];
    %below calculate decode error
    Decoded_Sequence(:,8)=vecnorm(Decoded_Sequence(:,2:3)-Decoded_Sequence(:,4:5),2,2); % decode error
end
if nargout==3
    DecodeSmoothFactor=2;binsize=2;
    Decoded_Data_Smooth=zeros(size(Decoded_Data));
    for i=1:size(Decoded_Data,1)
        ind=vecnorm(maze1d_position(:,1:2)-maze1d_position(i,1:2),2,2)<=DecodeSmoothFactor*binsize;
        Decoded_Data_Smooth(i,:)=nanmean(Decoded_Data(ind,:),1);
    end
    [~,I]=max(Decoded_Data_Smooth,[],1);
    Decoded_Sequence(:,11:12)=maze1d_position(I,:); % decoded position after smooth
    if ~isempty(Position_Data)
        Decoded_Sequence(:,13)=vecnorm(Decoded_Sequence(:,11:12)-Decoded_Sequence(:,4:5),2,2); % smooth decode error
    end
end