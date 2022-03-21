clear all; clc;
%%
path='E:\1_work_moreVols\data\pmip3_q\';
%% GAO
% q00=ncread([path '850_2hus_ccsm4_085101-184912.nc'],'hus');
% st_clim(:,:,:,1)=q00;clear q00;
% q00=ncread([path '850_2hus_ipsl_085101-184912.nc'],'hus'); 
% st_clim(:,:,:,2)=q00;clear q00;
% q00=ncread([path '850_2hus_fgoasl_085101-184912.nc'],'hus');
% st_clim(:,:,:,3)=q00;clear q00;
% q00=ncread([path '850_2hus_mri_085101-184912.nc'],'hus');
% st_clim(:,:,:,4)=q00;clear q00;
% q00=ncread([path '850_2hus_bcc_085101-184912.nc'],'hus');
% st_clim(:,:,:,5)=q00;clear q00;
%% CROWLEY
q00=ncread([path '850_2hus_csiro_085101-184912.nc'],'hus');
st_clim(:,:,:,1)=q00;clear q00;
q00=ncread([path '850_2hus_mpi_085101-184912.nc'],'hus');
st_clim(:,:,:,2)=q00;clear q00;
q00=ncread([path '850_2hus_giss_085101-184912.nc'],'hus');
st_clim(:,:,:,3)=q00;clear q00;
%q00=ncread([path '850_2hus_hadcm_085101-184912.nc'],'hus');
%st_clim(:,:,:,9)=q00;clear q00;
q00=ncread([path '850_2hus_miroc_085101-184912.nc'],'hus');
st_clim(:,:,:,4)=q00;clear q00;
%%
[m,n,ll,modelnum]=size(st_clim);lv=ll/12;st_annual=zeros(m,n,lv,modelnum);n1=851:1849;
%% MJJAS mean
for k=1:lv
    st_annual(:,:,k,:)=mean(st_clim(:,:,(k-1)*12+5:(k-1)*12+9,:),3);
end
key_dates_Gao=[1600;1641;1673;1809;1815;1835];

[erup_years,events]=intersect(n1,key_dates_Gao); 
% n_events=length(events);
n_events=7; % due to Crowley has seven events;
%% Superposed Epoch Analysis
before = 5;after = 5;wj=before+after+1;
Xevents = nan(m,n,wj,n_events,modelnum); Xeventss = nan(m,n,wj,n_events,modelnum); 
% %% Gao
% for i=1:n_events-1
%     Xevents(:,:,:,i,[1,3,4,5]) = st_annual(:,:,events(i)-before:events(i)+after,[1,3,4,5]);
% end
% %% anomaly
% for i=1:wj
%      Xeventss(:,:,i,:,[1,3,4,5]) = Xevents(:,:,i,:,[1,3,4,5])-nanmean(Xevents(:,:,1:before,:,[1,3,4,5]),3);
% end
%% Crowley
key_dates_Cro=[1576;1594;1674;1696;1809;1816;1835]; 
[erup_years_Cro,events_Cro]=intersect(n1,key_dates_Cro); 
n_events_Cro=length(events_Cro);
%% 
for i=1:n_events_Cro
    Xevents(:,:,:,i,1:4) = st_annual(:,:,events_Cro(i)-before:events_Cro(i)+after,1:4);
end
%% anomaly
for i=1:wj
     Xeventss(:,:,i,:,1:4) = Xevents(:,:,i,:,1:4)-nanmean(Xevents(:,:,1:before,:,1:4),3);
end
% %% IPSL
% key_dates_ipsl=[1600;1641;1673;1809;1815;1835];
% [erup_years_ipsl,events_ipsl]=intersect(n1,key_dates_ipsl); 
% n_events_ipsl=length(events_ipsl);
% %% 
% for i=1:n_events_ipsl
%     Xevents(:,:,:,i,2) = st_annual(:,:,events_ipsl(i)-before:events_ipsl(i)+after,2);
% end
% %% anomaly
% for i=1:wj
%      Xeventss(:,:,i,:,2) = Xevents(:,:,i,:,2)-nanmean(Xevents(:,:,1:before,:,2),3);
% end
%% m,n,wj,n_events,modelnum

% pmip_hus_year1_clim=Xevents(:,:,before+2,:,:);
% pmip_hus_year1=Xevents(:,:,before+2,:,:);

% save pmip_hus_year1_10mjjas_climCrow.mat pmip_hus_year1_clim;
% save pmip_hus_year1_10mjjasCrow.mat pmip_hus_year1;

save pmip_hus_wj_7mjjas_climCrow.mat Xevents;
save pmip_hus_wj_7mjjasCrow.mat Xeventss;