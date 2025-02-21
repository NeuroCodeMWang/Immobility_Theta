function generate_SupFigure_6_Jan25(ifsavedata)
if nargin<1
    ifsavedata=0;
end
homedir='X:\Mengni\Data_Analysis\Session_combined_0324';cd(homedir);load('SessionSet16');
datadir='X:\Mengni\Data_Analysis\Paper_Figures\Figure_Data_Nov2024';
if 0
CALCULATE_PLACE_FIELDS_8ARMS_NPX_move_immobile_theta;
unitall=[];property=[];fieldall=[];rewcenfire=[];spikeall=[];
for in=1:16
    savedir=SessionSet16{in};
    cd(savedir);
    load('Field_Data_move_immobile_theta','field_property','unit_decode','Field_Data','spikedata');
    a=[];
    for v=1:7
        [linearized_maze]=linearize_placefield(Field_Data(:,:,:,v));
        a(:,:,v)=[squeeze(nanmean(linearized_maze(:,36:50,:),2)),squeeze(nanmean(nanmean(linearized_maze(:,1:15,:),2),3))];
    end
    lia=ismember(unit_decode(:,16),[1,2,5,6]);
    unit=unit_decode(lia,:);
    prop=field_property(lia,:,:);
    fieldall=cat(3,fieldall,Field_Data(:,:,lia,:));
    lia2=ismember(spikedata(:,2),unit(:,17));
    spikedata(:,11)=in;
    spikeall=[spikeall;spikedata(lia2,:)];
    unit(:,18)=in;
    unitall=[unitall;unit];
    property=cat(1,property,prop);
    rewcenfire=cat(1,rewcenfire,a(lia,:,:));
end
rewcenfire(:,10,:)=nanstd(rewcenfire(:,1:8,:),0,2)./nanmean(rewcenfire(:,1:8,:),2);
rewcenfire(:,11,:)=rewcenfire(:,9,:)./nanmean(rewcenfire(:,1:8,:),2);
cd(homedir);save('SupFigure_6','rewcenfire','property','unitall','spikeall',"fieldall",'-v7.3');
else
    cd(homedir);load('SupFigure_6');
end
unitid=[71,13,55,15];v1=1;vid=4;position0=110;Bin_Size=2;ns=[];
unitcode=[1,2,5,6];unitname={'ca1 pyr','ca3 pyr','ca1 int','ca3 int'};
propname={'peak field rate','field size: cm^2','proportion of abnormal value','info rate'};
vname={'>10 cm/s','<10 cm/s','<5 cm/s','< 1cm/s'};
figure;
for u=1:4
    indu=unitall(:,16)==unitcode(u);
    unit=unitall(indu,:);
    prop=property(indu,:,:);ns(u,:)=[size(unit,1),sum(sum(isnan(squeeze(prop(:,1,[1,4]))),2)==0)];
    rew=squeeze(rewcenfire(indu,10,:));
    id=unitid(u);
    spi=spikeall(spikeall(:,11)==unit(id,18) & spikeall(:,2)==unit(id,17),:);
    in=unit(id,18);
    savedir=SessionSet16{in};
    cd(savedir);
    load('Position_Data_Maze.mat');
    positionbin=ceil((Position_Data(:,2:3)+position0)/Bin_Size);
    if 0
    subplot(5,4,u);
    ind=Position_Data(:,5)>=10;
    plot(positionbin(ind,1),positionbin(ind,2),'.','Color',[1 1 1]*0.5);hold on;
    ind=spi(:,9)>=10;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);
    title([unitname{u},vname{1}]);
    subplot(5,4,u+4);
    ind=Position_Data(:,5)<5;
    plot(positionbin(ind,1),positionbin(ind,2),'.','Color',[1 1 1]*0.5);hold on;
    ind=spi(:,9)<5 & spi(:,10)==1;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);
    title([unitname{u},vname{3}]);
    end
    subplot(5,4,u+4*2);
    plot(prop(:,1,v1),prop(:,1,vid),'k.');
    a=ceil(max([prop(:,1,v1);prop(:,1,vid)])/10)*10;hold on;plot([0,a],[0,a],'r--');xlim([0 a]);ylim([0 a]);
    %b=prop(:,1,v1)-prop(:,1,vid);[~,p]=ttest(b(isnan(b)),'tail','right');%p=signrank(b,0,'tail','right');
    [~,p,~,stat]=ttest(prop(:,1,v1),prop(:,1,vid),'tail','right');
    title([unitname{u},' : field peak, p=',num2str(p,2),'; t=',num2str(stat.tstat,4)]);xlabel(vname{v1});ylabel(vname{vid});
    subplot(5,4,u+4*3);
    plot(prop(:,4,v1),prop(:,4,vid),'k.');
    a=ceil(max([prop(:,4,v1);prop(:,4,vid)]));hold on;plot([0,a],[0,a],'r--');xlim([0 a]);ylim([0 a]);
    %b=prop(:,4,v1)-prop(:,4,vid);[~,p]=ttest(b,'tail','right');%p=signrank(b,0,'tail','right');
    [~,p,~,stat]=ttest(prop(:,4,v1),prop(:,4,vid),'tail','right');
    title([unitname{u},' : info rate, p=',num2str(p,2),'; t=',num2str(stat.tstat,4)]);xlabel(vname{v1});ylabel(vname{vid});
    subplot(5,4,u+4*4);
    plot(rew(:,v1),rew(:,vid),'k.');
    a=ceil(max([rew(:,v1);rew(:,vid)]));hold on;plot([0,a],[0,a],'r--');xlim([0 a]);ylim([0 a]);
    %b=rew(:,v1)-rew(:,vid);[~,p]=ttest(b,'tail','right');%p=signrank(b,0,'tail','right');
    [~,p,~,stat]=ttest(rew(:,v1),rew(:,vid),'tail','right');
    title([unitname{u},' : std/mean, p=',num2str(p,2),'; t=',num2str(stat.tstat,4)]);xlabel(vname{v1});ylabel(vname{vid});
end
figure_title='SupFigure6_Jan2025';save_current_figure(figure_title);

if 0
% below to look for good examples
figure;
for u=2%:4
    indu=unitall(:,16)==unitcode(u);
    unit=unitall(indu,:);
    for i=1:size(unit,1)
        if in~=unit(i,18)
            in=unit(i,18);
            savedir=SessionSet16{in};
            cd(savedir);
            load('Position_Data_Maze.mat');
            positionbin=ceil((Position_Data(:,2:3)+position0)/Bin_Size);
        end
        spi=spikeall(spikeall(:,11)==unit(i,18) & spikeall(:,2)==unit(i,17),:);
        subplot(1,2,1);
        ind=Position_Data(:,5)>=10;
        plot(positionbin(ind,1),positionbin(ind,2),'.','Color',[1 1 1]*0.5);hold on;
        ind=spi(:,9)>=10;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);
        %title([unitname{u},vname{1}]);
        subplot(1,2,2);
        ind=Position_Data(:,5)<5;
        plot(positionbin(ind,1),positionbin(ind,2),'.','Color',[1 1 1]*0.5);hold on;
        ind=spi(:,9)<5 & spi(:,10)==1;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);
        title(num2str([in i]));
        pause;clf;
    end
end
vname={'>= 10cm/s','theta <10','nontheta <10','theta <5','nontheta <5','theta <1','nontheta <1'};
figure;pname={'reward','center'};
for u=1:4
    indu=unitall(:,16)==unitcode(u);
    rew=rewcenfire(indu,10:11,:);
    for p=1:2
        subplot(2,4,u+(p-1)*4);
        plot(squeeze(rew(:,p,1)),squeeze(rew(:,p,2)),'.');
        hold on;plot([0 3],[0 3],'r');
        title([pname{p},' + ',unitname{u}]);
    end
end
figure;
for u=3%:4
    indu=unitall(:,16)==unitcode(u);
    unit=unitall(indu,:);
    for i=1:size(unit,1)
        spi=spikeall(spikeall(:,11)==unit(i,18) & spikeall(:,2)==unit(i,17),:);
        subplot(2,4,1);
        ind=spi(:,9)>=10;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{1});
        subplot(2,4,2);
        ind=spi(:,9)<10 & spi(:,10)==1;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{2});
        subplot(2,4,3);
        ind=spi(:,9)<10 & spi(:,10)==0;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{3});
        subplot(2,4,4);
        ind=spi(:,9)<5 & spi(:,10)==1;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{4});
        subplot(2,4,5);
        ind=spi(:,9)<5 & spi(:,10)==0;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{5});
        subplot(2,4,6);
        ind=spi(:,9)<=1 & spi(:,10)==1;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{6});
        subplot(2,4,7);
        ind=spi(:,9)<=1 & spi(:,10)==0;plot(spi(ind,7),spi(ind,8),'r.');xlim([0 110]);ylim([0 110]);title(vname{7});
        title(num2str([i,unit(i,10)]));
        pause;clf;
    end
end

figure;vid=4;v1=1;
for u=1:4
    indu=unitall(:,16)==unitcode(u);
    field=fieldall(:,:,indu,:);
    prop=property(indu,:,:);
    rew=squeeze(rewcenfire(indu,10,:));
    subplot(5,4,u);
    imagesc(nanmean(field(:,:,:,v1),3));set(gca,'YDir','normal');
    title([unitname{u},vname{v1}]);
    subplot(5,4,u+4);
    imagesc(nanmean(field(:,:,:,vid),3));set(gca,'YDir','normal');
    title([unitname{u},vname{vid}]);
    subplot(5,4,u+4*2);
    plot(prop(:,1,v1),prop(:,1,vid),'k.');
    a=ceil(max([prop(:,1,v1);prop(:,1,vid)])/10)*10;hold on;plot([0,a],[0,a],'r');
    %b=prop(:,1,v1)-prop(:,1,vid);[~,p]=ttest(b(isnan(b)),'tail','right');%p=signrank(b,0,'tail','right');
    [~,p]=ttest(prop(:,1,v1),prop(:,1,vid),'tail','right');
    title([unitname{u},' : field peak, p=',num2str(p,2)]);xlabel(vname{v1});ylabel(vname{vid});
    subplot(5,4,u+4*3);
    plot(prop(:,4,v1),prop(:,4,vid),'k.');
    a=ceil(max([prop(:,4,v1);prop(:,4,vid)]));hold on;plot([0,a],[0,a],'r');
    %b=prop(:,4,v1)-prop(:,4,vid);[~,p]=ttest(b,'tail','right');%p=signrank(b,0,'tail','right');
    [~,p]=ttest(prop(:,4,v1),prop(:,4,vid),'tail','right');
    title([unitname{u},' : info rate, p=',num2str(p,2)]);xlabel(vname{v1});ylabel(vname{vid});
    subplot(5,4,u+4*4);
    plot(rew(:,v1),rew(:,vid),'k.');
    a=ceil(max([rew(:,v1);rew(:,vid)]));hold on;plot([0,a],[0,a],'r');
    %b=rew(:,v1)-rew(:,vid);[~,p]=ttest(b,'tail','right');%p=signrank(b,0,'tail','right');
    [~,p]=ttest(rew(:,v1),rew(:,vid),'tail','right');
    title([unitname{u},' : std/mean, p=',num2str(p,2)]);xlabel(vname{v1});ylabel(vname{vid});
end
%axis equal 
figure_title='SupFigure14_unit_spatial_move_immobile';save_current_figure(figure_title);

unitcode=[1,2,5,6];unitname={'ca1 pyr','ca3 pyr','ca1 int','ca3 int'};
propname={'peak field rate','field size: cm^2','proportion of abnormal value','info rate'};
vname={'>10 cm/s','<10 cm/s','<5 cm/s','< 1cm/s'};
figure;
for u=1:4
    indu=unitall(:,16)==unitcode(u);
    prop=property(indu,:,:);
    for p=1:4
        subplot(4,4,u+(p-1)*4);
        for i=1:size(prop,1)
            hold on;plot([1:3],squeeze(prop(i,p,:)),'o-');
        end
        title([propname{p},' + ',unitname{u}]);
    end
end
figure;
for u=1:4
    indu=unitall(:,16)==unitcode(u);
    prop=property(indu,:,:);
    for p=1:4
        subplot(4,4,u+(p-1)*4);
        plot(prop(:,p,1),prop(:,p,3),'.'); 
        hold on;plot([0,4],[0,4],'r');
        title([propname{p},' + ',unitname{u}]);
    end
end
figure;
for u=2%:4
    indu=unitall(:,16)==unitcode(u);
    field=fieldall(:,:,indu,:);
    for i=1:size(field,3)
        for v=1:3
            subplot(1,3,v);imagesc(field(:,:,i,v));
            set(gca,'YDir','normal');
        end
        pause;clf;
    end
end
figure;
for u=1:4
    indu=unitall(:,16)==unitcode(u);
    field=squeeze(nanmean(fieldall(:,:,indu,:),3));
    for v=1:4
        subplot(4,4,u+(v-1)*4);
        imagesc(field(:,:,v));
        set(gca,'YDir','normal');
        title([vname{v},' + ',unitname{u}]);
    end
end
end