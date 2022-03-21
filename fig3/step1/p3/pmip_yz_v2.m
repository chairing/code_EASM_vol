clear;clc;
%%
path='../../../../data/pmip3/';
%% GAO
q00=ncread([path '2pr_Amon_CCSM4_past1000_r1i1p1_085001-185012.nc'],'pr');
st_clim(:,:,:,1)=q00(:,:,13:(1849-850+1)*12);clear q00;
q00=ncread([path '2pr_Amon_IPSL-CM5A-LR_past1000_r1i1p1_085001-185012.nc'],'pr');
st_clim(:,:,:,2)=q00(:,:,13:(1849-850+1)*12);clear q00;
q00=ncread([path '2pr_Amon_FGOALS-s2_past1000_r1i1p1_085001-185012.nc'],'pr');
st_clim(:,:,:,3)=q00(:,:,13:(1849-850+1)*12);clear q00;
q00=ncread([path '2pr_Amon_MRI-CGCM3_past1000_r1i1p1_085001-185012.nc'],'pr');
st_clim(:,:,:,4)=q00(:,:,13:(1849-850+1)*12);clear q00;
q00=ncread([path '2pr_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc'],'pr');
st_clim(:,:,:,5)=q00(:,:,13:(1849-850+1)*12);clear q00;
%% CROWLEY
% q00=ncread([path '2pr_Amon_CSIRO-Mk3L-1-2_past1000_r1i1p1_085101-185012.nc'],'pr');
% st_clim(:,:,:,1)=q00(:,:,1:(1849-851+1)*12);clear q00;%% start 851
% q00=ncread([path '2pr_Amon_MPI-ESM-P_past1000_r1i1p1_085001-184912.nc'],'pr');
% st_clim(:,:,:,2)=q00(:,:,13:(1849-850+1)*12);clear q00;
% q00=ncread([path '2pr_Amon_GISS-E2-R_past1000_r1i1p121_085001-185012.nc'],'pr');
% st_clim(:,:,:,3)=q00(:,:,13:(1849-850+1)*12);clear q00;
% q00=ncread([path '2pr_Amon_HadCM3_past1000_r1i1p1_085001-185012.nc'],'pr');
% st_clim(:,:,:,4)=q00(:,:,13:(1849-850+1)*12);clear q00;
% q00=ncread([path '2pr_Amon_MIROC-ESM_past1000_r1i1p1_085001-184912.nc'],'pr');
% st_clim(:,:,:,5)=q00(:,:,13:(1849-850+1)*12);clear q00;

[m,n,ll,modelnum]=size(st_clim);lv=ll/12;st_annual=zeros(m,n,lv,modelnum);n1=851:1849;
%% JJA mean
for k=1:lv
    st_annual(:,:,k,:)=mean(st_clim(:,:,(k-1)*12+5:(k-1)*12+9,:),3);
end

key_dates_Gao=[1600;1641;1673;1809;1815;1835];

[erup_years,events]=intersect(n1,key_dates_Gao); 
% n_events=length(events);
n_events=7 ; % due to Crowley has Seven events;
%% Superposed Epoch Analysis
before = 5;after = 5;wj=before+after+1;
Xevents = nan(m,n,wj,n_events,modelnum); Xeventss = nan(m,n,wj,n_events,modelnum); 
%% Gao
for i=1:n_events-1
    Xevents(:,:,:,i,[1,3,4,5]) = st_annual(:,:,events(i)-before:events(i)+after,[1,3,4,5]);
end
%% anomaly
for i=1:wj
     Xeventss(:,:,i,:,[1,3,4,5]) = Xevents(:,:,i,:,[1,3,4,5])-nanmean(Xevents(:,:,1:before,:,[1,3,4,5]),3);
end
% %% Crowley
% key_dates_Cro=[1576;1594;1674;1696;1809;1816;1835]; 
% [erup_years_Cro,events_Cro]=intersect(n1,key_dates_Cro); 
% n_events_Cro=length(events_Cro);
% %% 
% for i=1:n_events_Cro
%     Xevents(:,:,:,i,1:5) = st_annual(:,:,events_Cro(i)-before:events_Cro(i)+after,1:5);
% end
% %% anomaly
% for i=1:wj
%      Xeventss(:,:,i,:,1:5) = Xevents(:,:,i,:,1:5)-nanmean(Xevents(:,:,1:before,:,1:5),3);
% end
%% IPSL
key_dates_ipsl=[1600;1641;1673;1809;1815;1835];
[erup_years_ipsl,events_ipsl]=intersect(n1,key_dates_ipsl); 
n_events_ipsl=length(events_ipsl);
%% 
for i=1:n_events_ipsl
    Xevents(:,:,:,i,2) = st_annual(:,:,events_ipsl(i)-before:events_ipsl(i)+after,2);
end
%% anomaly
for i=1:wj
     Xeventss(:,:,i,:,2) = Xevents(:,:,i,:,2)-nanmean(Xevents(:,:,1:before,:,2),3);
end
%% m,n,wj,n_events,modelnum

pmip_pr_year1=Xeventss(:,:,before+2,:,:);

save pmip_pr_year1_6mjjasGao.mat pmip_pr_year1;
%save pmip_pr_wj_10djf.mat Xeventss;
