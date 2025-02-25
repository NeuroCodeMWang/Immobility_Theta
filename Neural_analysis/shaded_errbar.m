function err=shaded_errbar(x,y,colorspec,ifphase,marker,legendstring)

if nargin<3 | isempty(colorspec)
    colorspec=[0 0.4470 0.7410];
end
if nargin<4 | isempty(ifphase)
    ifphase=0;
end
if nargin<5 | isempty(marker)
    marker='-';
end
n=sum(~isnan(y),2);
x=x(:);
if ifphase==0
    err=nanstd(y,0,2).*((n).^(-0.5));
    y=nanmean(y,2);
else
    ym=[];err=[];
    for i=1:length(x)
        [mu,~,err1]=circ_meand(y(i,:));
        ym(i)=mu;
        err(i)=err1;%*((n(i))^(-0.5));
    end
    y=ym(:);
    err=err(:);
    if max(y)-min(y)>180
        y=y+(y<180)*360;
    end
end
if isstring(colorspec)
    hold on;plot(x,y,marker,colorspec);
else
    hold on;plot(x,y,marker,'color',colorspec);
end
hold on;patch([x;flip(x)], [y-err;flip(y+err)], colorspec, 'FaceAlpha',0.3, 'EdgeColor','none');
if 0
legendstring1=[];
for i=1:length(legendstring)
    legendstring1{2*i-1}=legendstring{i};
    legendstring1{2*i}=[];
end
legend(legendstring1);
end