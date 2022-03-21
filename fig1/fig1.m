clc;clear;load coast;myp=othercolor('BrBG11');n1=1:1999;addpath('../utilities');path='../data/';
%%
load '../data/vol.txt';
Erup_Type=1;Erup_Forcing= -0;
key_dates=find(vol(:,1)>1470  & vol(:,2)==Erup_Type & vol(:,3)<Erup_Forcing);
events_real=vol(key_dates,1);
events_real=flipud(events_real);
events_noel=events_real-1470+1;
n_events=length(events_noel);before = 5;after = 5;wj=before+after+1;
%% May¨CSeptember precipitation for the Asian monsoon region (Feng et al. 2013)--------------------------1
lats=ncread([path 'Feng.P.recon.1470-1999.nc'],'lat');lons=ncread([path 'Feng.P.recon.1470-1999.nc'],'lon');
pr0=ncread([path 'Feng.P.recon.1470-1999.nc'],'pre')/30; [m,n,ll]=size(pr0);
Xevents3d = nan(m,n,wj,n_events); Xevents3ds = nan(m,n,wj,n_events); 
%%
for i=1:n_events
    Xevents3d(:,:,:,i) = pr0(:,:,events_noel(i)-before:events_noel(i)+after);
end
%% anomaly
for i=1:wj
     Xevents3ds(:,:,i,:) = Xevents3d(:,:,i,:)-nanmean(Xevents3d(:,:,1:before,:),3);
end
isam=500; sig= 10 ;cor=zeros(m,n,2);
for i=1:m
    for j=1:n
            sam=zeros( isam,1);
            for ia=1:isam
                for k=1:n_events
                    ii=max(1,floor( rand()*(wj) )) ;
                    jj=max(1,floor( rand()*n_events ));
                    sam(ia,1)=sam(ia,1)+  Xevents3ds(i,j,ii,jj)/n_events;
                end 
            end 
            [nn,xout]=hist(sam(:,1),100);
            for ia=1:100
                per=sum(nn(1,1:ia)/isam*100);
                if(per<=sig)
                    cor(i,j,1)=xout(1,ia);
                end 
                per=sum(nn(1,100-ia+1:100 )/isam*100);
                if(per<=sig)
                    cor(i,j,2)=xout(1,100-ia+1);
                end 
            end 
        end;end 
year=2;
Xcomps3d = squeeze(nanmean(Xevents3ds(:,:,before+year,:),4)); 

figure1 =figure('WindowState','maximized');fs=28;
subplot(1,2,1)
sx0=nan(m*n,1); sy0=nan(m*n,1);k0=0;dom=zeros(m,n);
%% sig for abs(ans)> 0.1
for i=1:m for j=1:n
        k0=k0+1;
        if isnan(Xcomps3d(i,j))
        else
        if(Xcomps3d(i,j)>0.1 && Xcomps3d(i,j)> cor(i,j,2) || Xcomps3d(i,j)<-0.1 && Xcomps3d(i,j)< cor(i,j,1))
            sx0(k0,1)=NaN; sy0(k0,1)=NaN;
        else
            sx0(k0,1)=lons(i); sy0(k0,1)=lats(j);
        end;end;end;end
contourf(lons,lats,Xcomps3d',[-10:0.01:0.0001 0.0001:0.01:10],'linestyle','none');caxis([-0.5 0.5]); colormap(myp);hold on;
hc=colorbar('hori','Position',[0.133203125 0.0811320754716981 0.24 0.02],'FontSize',fs); %   [x0 y0 Width Height] 
hc.Label.String = 'mm day^-^1';
hc.Label.Position= [0.713355533074867 2.68461551116063 0];

scatter(sx0,sy0,2, 'MarkerEdgeColor','none','MarkerFaceColor','k');hold on;
text(100,55.74,'a','FontWeight','bold','FontSize',fs+4,'FontName','Helvetica');

set(gca,'xlim',[100 134.85], 'xtick',[100:10:140], 'XTickLabels',{'100¡ãE','110¡ãE','120¡ãE','130¡ãE','140¡ãE'},.....
    'ylim',[5 54],'ytick',[10:10:50], 'yticklabels',{'10¡ãN','20¡ãN' ,'30¡ãN','40¡ãN','50¡ãN' },'TickDir','out','FontSize',fs,......
    'FontName','Helvetica','PlotBoxAspectRatio',[1 1 1],'FontWeight','bold' );
hold on;h1=plot(long,lat,'k',long+360,lat,'k'  );set(h1,'color',[0 0 0]); hold on; 
h=line([0 356.5],[0 0]);set(h,'linestyle','--','color',[0.6 0.6 0.6]); box on;hold on;
%%
S = shaperead('worldrivers.shp','UseGeoCoords',true);
x3=S(82).Lon;y3=S(82).Lat;x4=S(84).Lon;y4=S(84).Lat;   % ³¤½­
hold on;plot(x3,y3,'r',x4,y4,'r','linewidth',1)

x1=106;x2=122;y1=24.5;y2=34;
 
rectangle('position',...
    [lons(lons==x1+0.25) lats(lats==y1+0.25) lons(lons==x2+0.25)-lons(lons==x1+0.25) lats(lats==y2+0.25)-lats(lats==y1+0.25)],...
    'LineWidth',2,'LineStyle','--','EdgeColor',[1 0 0]);hold on;
%%-------------------------------------------------------------------------------------------------------1
ps=zeros(ll,1);
kk=0;
for i=1:m   
    for j=1:n 
        if isnan(pr0(i,j,:)) else
            if (lats(j) > y1 && lats(j) < y2 && lons(i) > x1 && lons(i) < x2)
            cc=pr0(i,j,:)*cos(lats(j)/180*pi);
            ps=ps+cc(:);
            kk=kk+cos(lats(j)/180*pi);
        end;end;end
end;ps=ps/kk;
save fengEasm.mat ps;
%%
Xeventfeng = nan(wj,n_events);  
%%
for i=1:n_events
    Xeventfeng(:,i) = ps(events_noel(i)-before:events_noel(i)+after);
    Xeventfeng(:,i) = Xeventfeng(:,i)-nanmean(Xeventfeng(1:before,i)); 
end
Xcomp = nanmean(Xeventfeng,2); 
%% -----------------------------------------------------------------------------------Plot
subplot(1,2,2);
corf=repeat_cc3(Xeventfeng,n_events,wj);

te=(-before:1:after-2);
h=bar(-before:-1,Xcomp(1:before),'FaceColor','b','EdgeColor','none');hold on; h.FaceAlpha = 0.5;
h=bar(0:after-2,Xcomp(after+1:wj-2),'FaceColor','r','EdgeColor','none');hold on; h.FaceAlpha = 0.5;

line([-before-0.4 after-2+0.4],[corf,corf],'Color',0.5*[1 1 1],'linestyle','-');hold on;

line([0 0],[-1.1 1.1],'linewidth',2.0,'color',0.7*[1 1 1]),hold on,

xlabel('Year','FontWeight','bold','FontSize',fs,'FontName','Helvetica');
ylabel('Pr (mm day^-^1)','FontWeight','bold','FontSize',fs,'FontName','Helvetica');
axis([-before-1 after+1-2 -0.21 0.21])
set(gca,'xtick',[-before:after-2],'XTickLabel',{'-5','-4','-3','-2','-1','0','+1','+2','+3'},...
    'FontName','Helvetica','fontsize',fs,'TickDir','out' ,'PlotBoxAspectRatio',[1 1 1] ,'FontWeight','bold')
text(-6,0.2237,'b','FontWeight','bold','FontSize',fs+4,'FontName','Helvetica');