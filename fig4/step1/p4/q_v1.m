clc;clear;
%% common info
path='E:\1_work_moreVols\data\pmip4\';
before = 5;after = 5;wj=before+after+1;
n1=850:1849;
%% Volcanic Data
toohey(:,1)=ncread([path 'eVolv2k_v2.1_ds_1.nc'],'yearCE');
toohey(:,2)=ncread([path 'eVolv2k_v2.1_ds_1.nc'],'lat');
toohey(:,3)=ncread([path 'eVolv2k_v2.1_ds_1.nc'],'hemi');
toohey(:,4)=ncread([path 'eVolv2k_v2.1_ds_1.nc'],'ssi');
%% find vol year in Toohey
volhemi= -1;
key_dates=find(toohey(:,3)~= volhemi & toohey(:,1) > 1470+after & toohey(:,1) < 1850-before & toohey(:,4) > 7 ); % & abs(toohey(:,2)) < 10 
toohey_trop=toohey(key_dates,:);
%% 
ncfile=dir([path '*hus*.nc']);
lons=ncread([path ncfile(1).name],'lon');
lats=ncread([path ncfile(1).name],'lat');
%%
ncfilenumber=numel(ncfile)
%%
[m]=size(lons,1);[n]=size(lats,1);
%%
[erup_years,events]=intersect(n1,toohey(key_dates,1));
n_events=length(events)
pr_year1=nan(m,n,n_events,ncfilenumber);

for f=1:ncfilenumber
    st_clim=ncread([path ncfile(f).name],'hus');

    [m,n,ll]=size(st_clim);lv=ll/12;st_annual=zeros(m,n,lv);
    %% MJJAS mean
    for k=1:lv
        st_annual(:,:,k)=mean(st_clim(:,:,(k-1)*12+5:(k-1)*12+9),3);
    end
    %% Superposed Epoch Analysis
    Xevents = nan(m,n,wj,n_events); Xeventss = nan(m,n,wj,n_events);
    %% Gao
    for i=1:n_events
        Xevents(:,:,:,i) = st_annual(:,:,events(i)-before:events(i)+after);
    end
    %% anomaly
    for i=1:wj
        Xeventss(:,:,i,:) = Xevents(:,:,i,:)-nanmean(Xevents(:,:,1:before,:),3);
    end
    %% m,n,n_events,modelnum
%     q_year1_clim(:,:,:,f)=squeeze(Xevents(:,:,before+2,:));
%          q_year1(:,:,:,f)=squeeze(Xeventss(:,:,before+2,:));
    Xeventss_all(:,:,:,:,f)=Xeventss;
    Xevents_all(:,:,:,:,f)=Xevents;
end

% save q_year1_10mjjas_clim.mat q_year1_clim;
% save q_year1_10mjjas.mat q_year1;

save q_wj_20mjjas_clim.mat Xevents_all;
save q_wj_20mjjas.mat Xeventss_all;