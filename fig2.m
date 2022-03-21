clc;clear;load coast;myp=othercolor('BrBG11');addpath ../utilities;
path='../data/';

%% bandpass define
[b,a]=butter(4,[1/9,1/2]/0.5001,'bandpass');
ninm=nan(2007,11);
%% 
ila= -1 ;
for ik=[1,2,3,5,8]
    n1=load([path 'n' num2str(ik) '.txt']);n1(:,1)=n1(:,1)+ila; 
    n1(:,2)=n1(:,2)-mean(n1(:,2));n1(:,2)=n1(:,2)/std(n1(:,2));
    n1s=filtfilt(b,a,n1(:,2));ninm(n1(1,1):n1(end,1),ik)=n1s; clear n1 n1s;
end
for ik=[6,9,10]
    n1=load([path 'n' num2str(ik) '.txt']);n1(:,2)=-n1(:,2);n1(:,1)=n1(:,1)+ila;
    n1(:,2)=n1(:,2)-mean(n1(:,2));n1(:,2)=n1(:,2)/std(n1(:,2));
    n1s=filtfilt(b,a,n1(:,2));ninm(n1(1,1):n1(end,1),ik)=n1s; clear n1 n1s;
end

ik=4; load  ([path 'n4.txt']);n4=flipud(n4); n4(:,1)=n4(:,1)+ila ;
n4(:,2)=n4(:,2)-mean(n4(:,2));n4(:,2)=n4(:,2)/std(n4(:,2));
n4s=filtfilt(b,a,n4(:,2));ninm(n4(1,1):n4(end,1),ik)=n4s; clear n4 n4s;

ik=7; load ([path 'n7.txt']); n7(:,1)=n7(:,1)+ila+1; 
n7(:,2)=n7(:,2)-mean(n7(:,2));n7(:,2)=n7(:,2)/std(n7(:,2));
n7s=filtfilt(b,a,n7(:,2));ninm(n7(1,1):n7(end,1),ik)=n7s; clear n7 n7s;
%% add Dee -- Coral
load ([path 'July_Corals.mat']);

n1c=CDT;
n1c(:,2)=-n1c(:,2);

n1c(:,3)=(n1c(:,2)-nanmean(n1c(:,2)))/std(n1c(~isnan(n1c(:,2)),2));
n1cs=nan(861,1);
  n1cs(1:319)=filtfilt(b,a,n1c(1:319,3));
n1cs(491:559)=filtfilt(b,a,n1c(491:559,3));
n1cs(741:861)=filtfilt(b,a,n1c(741:861,3));
ninm(n1c(1,1):n1c(end,1),11)=n1cs; 

ninm_mean=nanmean(ninm,2);
ninostd=std(ninm_mean(~isnan(ninm_mean)));
%%
load '../data/vol.txt';
Erup_Type=1;Erup_Forcing= -0;
key_dates=find(vol(:,1)>1470  & vol(:,2)==Erup_Type & vol(:,3)<Erup_Forcing);
events_real=vol(key_dates,1);
events_real=flipud(events_real);before = 5;after = 5;wj=before+after+1;

kk=0;ks=0;%events_nino=nan(22,2);
for iv=1:length(events_real)
    if ninm_mean(events_real(iv)) > 0.5 * ninostd
        kk=kk+1;
        events_nino(kk,1)=events_real(iv);
    else
        ks=ks+1;
        events_nino(ks,2)=events_real(iv);
    end
end

evenes_all=(1470+5:1999-5)';ka=0;
for iv=1:length(evenes_all)
    if ninm_mean(evenes_all(iv)) > 0.5 * ninostd & evenes_all(iv)~=events_nino(:,1)
        ka=ka+1;
       events_nino(ka,3)=evenes_all(iv);
    end
end
events_nino(events_nino==0)=nan;

figure1 =figure('WindowState','maximized');fs=20;

for is=1:3
    
events_noel=events_nino(~isnan(events_nino(:,is)),is)-1470+1;
n_events=length(events_noel);

subplot(1,3,is)
%% May–September precipitation for the Asian monsoon region (Feng et al. 2013)---------------------------------
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
isam= 500;sig= 10;
cor=repeat_3d(Xevents3ds,n_events,wj,isam,sig);
year=2;
Xcomps3d = squeeze(nanmean(Xevents3ds(:,:,before+year,:),4)); 

sx0=zeros(m*n,1); sy0=zeros(m*n,1);k0=0;

for i=1:m for j=1:n
        k0=k0+1;
        if isnan(Xcomps3d(i,j))
        else
            if(Xcomps3d(i,j)> 0.1 && Xcomps3d(i,j)> cor(i,j,2) || Xcomps3d(i,j)< -0.1 && Xcomps3d(i,j)< cor(i,j,1))
                sx0(k0,1)=NaN; sy0(k0,1)=NaN;
            else
                sx0(k0,1)=lons(i); sy0(k0,1)=lats(j);
            end;end;end;end

contourf(lons,lats,Xcomps3d',[-10:0.02:0.0001 0.0001:0.02:10],'linestyle','none');caxis([-0.5 0.5]); colormap(myp);hold on;
scatter(sx0,sy0,2, 'MarkerEdgeColor','none','MarkerFaceColor','k');hold on;

set(gca,'xlim',[100 134.85], 'xtick',[100:10:140], 'XTickLabels',{'100°E','110°E','120°E','130°E','140°E'},.....
    'ylim',[5 54],'ytick',[10:10:50], 'yticklabels',{'10°N','20°N' ,'30°N','40°N','50°N' },'TickDir','out','FontSize',fs,......
    'FontName','Helvetica','PlotBoxAspectRatio',[1 1 1],'FontWeight','bold' );
hold on;h1=plot(long,lat,'k',long+360,lat,'k'  );set(h1,'color',[0 0 0]); hold on; 
h=line([0 356.5],[0 0]);set(h,'linestyle','--','color',[0.6 0.6 0.6]); box on;hold on;
%%
S = shaperead('worldrivers.shp','UseGeoCoords',true);
x3=S(82).Lon;y3=S(82).Lat;x4=S(84).Lon;y4=S(84).Lat;   % 长江
hold on;plot(x3,y3,'r',x4,y4,'r','linewidth',1)

x1=106;x2=122;y1=24.5;y2=34;
 
rectangle('position',...
    [lons(lons==x1+0.25) lats(lats==y1+0.25) lons(lons==x2+0.25)-lons(lons==x1+0.25) lats(lats==y2+0.25)-lats(lats==y1+0.25)],...
    'LineWidth',2,'LineStyle','--','EdgeColor',[1 0 0]);hold on;
end
subplot(1,3,1)
% title({'with El Ni$\ \tilde{n}$o'},'FontSize',18,'Interpreter','latex')
title(['vol with El Nino (' num2str(sum(~isnan(events_nino(:,1)))) ' cases)'],'FontSize',fs,'FontName','Helvetica','FontWeight','bold')
text(100,56.5,'a','FontWeight','bold','FontSize',fs+4,'FontName','Helvetica');
subplot(1,3,2)
title(['vol without El Nino (' num2str(sum(~isnan(events_nino(:,2)))) ' cases)'],'FontSize',fs,'FontName','Helvetica','FontWeight','bold')
text(100,56.5,'b','FontWeight','bold','FontSize',fs+4,'FontName','Helvetica');
subplot(1,3,3)
title(['El Nino without vol (' num2str(sum(~isnan(events_nino(:,3)))) ' cases)'],'FontSize',fs,'FontName','Helvetica','FontWeight','bold')
text(100,56.5,'c','FontWeight','bold','FontSize',fs+4,'FontName','Helvetica');
hc=colorbar('Position',[0.91 0.34 0.01 0.36]); %   [x0 y0 Width Height] 
hc.Label.String = 'mm day^-^1';
hc.Label.Rotation= 0;
hc.Label.Position= [2.89919996261599 -0.747597490335656 0];
annotation(figure1,'textbox',[0.235765625 0.769654088050314 0.00954687499999995 0.0117924528301887],...
    'String',{'~'},'LineStyle','none','FontWeight','bold','FontSize',fs,'FitBoxToText','off');
annotation(figure1,'textbox',[0.522874999999998 0.769654088050313 0.00954687499999995 0.0117924528301887],...
    'String',{'~'},'LineStyle','none','FontWeight','bold','FontSize',fs,'FitBoxToText','off');
annotation(figure1,'textbox',[0.730296875000001 0.769654088050313 0.00954687500000007 0.0117924528301887],...
    'String',{'~'},'LineStyle','none','FontWeight','bold','FontSize',fs,'FitBoxToText','off');