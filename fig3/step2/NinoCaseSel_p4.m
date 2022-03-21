clc;clear;
addpath E:\1_work_moreVols\data\pmip4;
%%
toohey(:,1)=ncread('eVolv2k_v2.1_ds_1.nc','yearCE');
toohey(:,2)=ncread('eVolv2k_v2.1_ds_1.nc','lat');
toohey(:,3)=ncread('eVolv2k_v2.1_ds_1.nc','hemi');
toohey(:,4)=ncread('eVolv2k_v2.1_ds_1.nc','ssi');

%% find vol year in Toohey
volhemi= -1;
key_dates=find(toohey(:,3)~= volhemi & toohey(:,1) > 1470+5 & toohey(:,1) < 1850-5 & toohey(:,4) > 7); % & abs(toohey(:,2)) < 10 

toohey_trop=toohey(key_dates,:);

%%
load nino34_pmip4.mat;glmrst=glmrst-273.15;
glmrst(:,2:3)=[];load nino34_pmip4_tos.mat;glmrst(:,2:3)=glmtos;

[ll,modelnum]=size(glmrst);lv=ll/12;
%% DJF mean
glm_naa=nan(lv,modelnum);glm_naa_a=nan(lv,modelnum);

for k=1:lv-1
    glm_naa(k,:)=nanmean(glmrst((k-1)*12+12:(k-1)*12+12+2,:));
end
%% events number
n1=850:1849;
[erup_years,events]=intersect(n1,toohey(key_dates,1));

before = 5;after = 5;wj=before+after+1;
n_events = length(erup_years);
%% Superposed Epoch Analysis
Xevents = nan(wj,n_events,3);    
%% Perform SEA:

for i=1:n_events
    Xevents(:,i,:) = glm_naa(events(i)-before:events(i)+after,:);
    Xevents(:,i,:) = (Xevents(:,i,:)-nanmean(Xevents(1:before,i,:)));%/max(abs(Xevents(:,i,iv))); % remove mean over "before" of window
end

nino=squeeze(Xevents(before+1,:,:));

save nino_tos.mat Xevents;
%% nino std check

[ll,vv]=size(glmrst);lv=ll/12;
%%
glm_ano_ce=nan(ll,vv);
for k=1:ll
    glm_ano_ce(k,:)=glmrst(k,:)-nanmean(glmrst);
end
%% DJF mean
glm_naa_cesm=nan(lv-1,vv);
for k=1:lv-1
    glm_naa_cesm(k,:)=mean(glm_ano_ce((k-1)*12+12:(k-1)*12+12+2,:));
end

std_vol=mean(squeeze(std(Xevents)));
std_all=std(glm_naa_cesm);

0.5*(std_all-std_vol)