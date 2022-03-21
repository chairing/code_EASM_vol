clc;clear;load coast;my=othercolor('BuDRd_18'); myp=othercolor('BrBG11');before=5; 
addpath(genpath('F:/1_work_moreVols/utilities/')); % import auxiliary functions
lons=ncread('F:/1_work_moreVols/data/pmip3/2ts_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc','lon');
lats=ncread('F:/1_work_moreVols/data/pmip3/2ts_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc','lat');
land=ncread('F:/1_work_moreVols/data/2sftlf_fx_past1000_r0i0p0.nc','sftlf');
%% PMIP3
load ../step1/p3/pmip_pr_year1_7mjjasGao.mat;          
pmip_pr_year1(:,:,:,:,6:10)=cell2mat(struct2cell(load('../step1/p3/pmip_pr_year1_7mjjasCrow.mat')));
pr0=squeeze(pmip_pr_year1)*86400;[m,n,n_events,modelnum]=size(pr0); 

load ../step1/p3/pmip_sst_year1_7mjjasGao.mat;
pmip_sst_year1(:,:,:,:,6:10)=cell2mat(struct2cell(load('../step1/p3/pmip_sst_year1_7mjjasCrow.mat')));
pmip_sst_year1=squeeze(pmip_sst_year1);

load ../step1/p3/pmip_ua_year1_7mjjasGao.mat; 
pmip_ua_year1(:,:,:,:,6:10)=cell2mat(struct2cell(load('../step1/p3/pmip_ua_year1_7mjjasCrow.mat')));
pmip_ua_year1=squeeze(pmip_ua_year1); 

load ../step1/p3/pmip_va_year1_7mjjasGao.mat; 
pmip_va_year1(:,:,:,:,6:10)=cell2mat(struct2cell(load('../step1/p3/pmip_va_year1_7mjjasCrow.mat')));
pmip_va_year1=squeeze(pmip_va_year1);

%% PMIP4
p4_ts_year1=cell2mat(struct2cell(load('../step1/p4/ts_year1_9mjjas.mat')));[m4,n4,n_events4,modelnum4]=size(p4_ts_year1);
p4_pr_year1=cell2mat(struct2cell(load('../step1/p4/pr_year1_9mjjas.mat')));
p4_ua_year1=cell2mat(struct2cell(load('../step1/p4/ua_year1_9mjjas.mat')));
p4_va_year1=cell2mat(struct2cell(load('../step1/p4/va_year1_9mjjas.mat')));
%% remove land or ocean info. due to some model using ST instead SST
for i=1:m
    for j=1:n
        if land(i,j) >0
 pmip_sst_year1(i,j,:,:)=nan;
    p4_ts_year1(i,j,:,:)=nan;
        else
            pr0(i,j,:,:)=nan;
    p4_pr_year1(i,j,:,:)=nan;
        end
    end
end
%--------------------------------------------------------------------------vol
%%
pre_m(:,:,1:10)=squeeze(nanmean(pr0,3));%--PMIP3
pre_m(:,:,11:13)=squeeze(nanmean(p4_pr_year1,3));%--PMIP4
prem=nanmean(pre_m,3);

sst_m(:,:,1:10)=squeeze(nanmean(pmip_sst_year1,3));%--PMIP3
sst_m(:,:,11:13)=squeeze(nanmean(p4_ts_year1,3));%--PMIP4
sstm=nanmean(sst_m,3);

ua_m(:,:,1:10)=squeeze(nanmean(pmip_ua_year1,3));%--PMIP3
ua_m(:,:,11:13)=squeeze(nanmean(p4_ua_year1,3));%--PMIP4
uam=nanmean(ua_m,3);

va_m(:,:,1:10)=squeeze(nanmean(pmip_va_year1,3));%--PMIP3
va_m(:,:,11:13)=squeeze(nanmean(p4_va_year1,3));%--PMIP4
vam=nanmean(va_m,3);
%% check
% figure;for iv=1:13; subplot(2,7,iv);contourf(pre_m(:,:,iv)');end
% figure;for iv=1:13; subplot(2,7,iv);contourf(sst_m(:,:,iv)');end
% figure;for iv=1:13; subplot(2,7,iv);contourf(ua_m(:,:,iv)'); end
% figure;for iv=1:13;subplot(2,7,iv);contourf(va_m(:,:,iv)'); end
%  stop;
%--------------------------------------------------------------------------volelnino
load ../step2/pmip3_cesmNE_nino7.mat; % E:\fig81600-971\FigN.m
degree=0.5;

nino_std=squeeze(std(Xcomp_10(:,:,1:10)));%% for PMIP
s2=squeeze(Xcomp_10(before+1,:,1:10));%% for PMIP
dsst3=nan(m,n,n_events,modelnum);dpre3=nan(m,n,n_events,modelnum);dua3=nan(m,n,n_events,modelnum);dva3=nan(m,n,n_events,modelnum);
dsst3p=nan(m,n,n_events,modelnum);dpre3p=nan(m,n,n_events,modelnum);dua3p=nan(m,n,n_events,modelnum);dva3p=nan(m,n,n_events,modelnum);

kinP=0;kinPp=0;
for iv=1:modelnum
    for ik=1:n_events
        if s2(ik,iv,1) > degree*nino_std(ik,iv)
            kinP=kinP+1;
            dsst3(:,:,ik,iv)=pmip_sst_year1(:,:,ik,iv);
            dpre3(:,:,ik,iv)=pr0(:,:,ik,iv);
            dua3(:,:,ik,iv)=pmip_ua_year1(:,:,ik,iv);
            dva3(:,:,ik,iv)=pmip_va_year1(:,:,ik,iv);
        end
        if s2(ik,iv,1) < degree*nino_std(ik,iv)
            kinPp=kinPp+1;
            dsst3p(:,:,ik,iv)=pmip_sst_year1(:,:,ik,iv);
            dpre3p(:,:,ik,iv)=pr0(:,:,ik,iv);
            dua3p(:,:,ik,iv)=pmip_ua_year1(:,:,ik,iv);
            dva3p(:,:,ik,iv)=pmip_va_year1(:,:,ik,iv);
        end
        
    end
end

dsst_m(:,:,1:10)=squeeze(nanmean(dsst3,3));
dpre_m(:,:,1:10)=squeeze(nanmean(dpre3,3));
 dua_m(:,:,1:10)=squeeze(nanmean(dua3,3));
 dva_m(:,:,1:10)=squeeze(nanmean(dva3,3));
 
dsstp_m(:,:,1:10)=squeeze(nanmean(dsst3p,3));
dprep_m(:,:,1:10)=squeeze(nanmean(dpre3p,3));
 duap_m(:,:,1:10)=squeeze(nanmean(dua3p,3));
 dvap_m(:,:,1:10)=squeeze(nanmean(dva3p,3)); 
 
%% PMIP4
p4nin340=cell2mat(struct2cell(load('../step2/nino_tos.mat')));

p4nin34=squeeze(p4nin340(before+1,:,:));
p4nin34_std=squeeze(std(p4nin340));

dsstp4=nan(m4,n4,n_events4,modelnum4);dprep4=nan(m4,n4,n_events4,modelnum4);duap4=nan(m4,n4,n_events4,modelnum4);dvap4=nan(m4,n4,n_events4,modelnum4);
dsstp4p=nan(m4,n4,n_events4,modelnum4);dprep4p=nan(m4,n4,n_events4,modelnum4);duap4p=nan(m4,n4,n_events4,modelnum4);dvap4p=nan(m4,n4,n_events4,modelnum4);

kin4=0;kin4p=0;
for iv=1:modelnum4
    for ik=1:n_events4
        if p4nin34(ik,iv) > degree*p4nin34_std(ik,iv)
            kin4=kin4+1;
            dprep4(:,:,ik,iv)=p4_pr_year1(:,:,ik,iv);
            dsstp4(:,:,ik,iv)=p4_ts_year1(:,:,ik,iv);
            duap4(:,:,ik,iv)=p4_ua_year1(:,:,ik,iv);
            dvap4(:,:,ik,iv)=p4_va_year1(:,:,ik,iv);
        end
        if p4nin34(ik,iv) < degree*p4nin34_std(ik,iv)
            kin4p=kin4p+1;
            dprep4p(:,:,ik,iv)=p4_pr_year1(:,:,ik,iv);
            dsstp4p(:,:,ik,iv)=p4_ts_year1(:,:,ik,iv);
            duap4p(:,:,ik,iv)=p4_ua_year1(:,:,ik,iv);
            dvap4p(:,:,ik,iv)=p4_va_year1(:,:,ik,iv);
        end
    end
end

dsst_m(:,:,11:13)=squeeze(nanmean(dsstp4,3));
dpre_m(:,:,11:13)=squeeze(nanmean(dprep4,3));
dua_m(:,:,11:13)=squeeze(nanmean(duap4,3));
dva_m(:,:,11:13)=squeeze(nanmean(dvap4,3));

dsstp_m(:,:,11:13)=squeeze(nanmean(dsstp4p,3));
dprep_m(:,:,11:13)=squeeze(nanmean(dprep4p,3));
duap_m(:,:,11:13)=squeeze(nanmean(duap4p,3));
dvap_m(:,:,11:13)=squeeze(nanmean(dvap4p,3));
%%
dsstm=squeeze(nanmean(dsst_m,3));
dprem=squeeze(nanmean(dpre_m,3));
duam=squeeze(nanmean(dua_m,3));
dvam=squeeze(nanmean(dva_m,3));

dsstmp=squeeze(nanmean(dsstp_m,3));
dpremp=squeeze(nanmean(dprep_m,3));
duamp=squeeze(nanmean(duap_m,3));
dvamp=squeeze(nanmean(dvap_m,3));
%% check
% figure;for iv=1:13; subplot(2,7,iv);contourf(dpre_m(:,:,iv)');end
% figure;for iv=1:13; subplot(2,7,iv);contourf(dsst_m(:,:,iv)');end
% figure;for iv=1:13; subplot(2,7,iv);contourf(dua_m(:,:,iv)'); end
% figure;for iv=1:13;subplot(2,7,iv);contourf(dva_m(:,:,iv)'); end
% 
% figure;for iv=1:13; subplot(2,7,iv);contourf(dprep_m(:,:,iv)');end
% figure;for iv=1:13; subplot(2,7,iv);contourf(dsstp_m(:,:,iv)');end
% figure;for iv=1:13; subplot(2,7,iv);contourf(duap_m(:,:,iv)'); end
% figure;for iv=1:13;subplot(2,7,iv);contourf(dvap_m(:,:,iv)'); end

casenum=kin4+kinP

casenump=kin4p+kinPp

%% T-test
  pre_m0=pre_m;  pre_m0(abs(pre_m0)<0.1)=nan;
 dpre_m0=dpre_m;dpre_m0(abs(dpre_m0)<0.1)=nan;
dprep_m0=dprep_m;dprep_m0(abs(dprep_m0)<0.1)=nan;

  sst_m0=sst_m;  sst_m0(abs(sst_m0)<0.1)=nan;
 dsst_m0=dsst_m;dsst_m0(abs(dsst_m0)<0.1)=nan;
dsstp_m0=dsstp_m;dsstp_m0(abs(dsstp_m0)<0.1)=nan;

figure;
subplot(1,2,1)
contourf(lons,lats,pre_m(:,:,1)');
subplot(1,2,2)
contourf(lons,lats,pre_m0(:,:,1)');

alpha=0.1;
presig=ttest(pre_m0,0,'Alpha',alpha,'Dim',3);
dpresig=ttest(dpre_m0,0,'Alpha',alpha,'Dim',3);
dprepsig=ttest(dprep_m0,0,'Alpha',alpha,'Dim',3);

sstsig=ttest(sst_m0,0,'Alpha',alpha,'Dim',3);
dsstsig=ttest(dsst_m0,0,'Alpha',alpha,'Dim',3);
dsstpsig=ttest(dsstp_m0,0,'Alpha',alpha,'Dim',3);

uasig=ttest(ua_m,0,'Alpha',alpha,'Dim',3);
duasig=ttest(dua_m,0,'Alpha',alpha,'Dim',3);
duapsig=ttest(duap_m,0,'Alpha',alpha,'Dim',3);
%
preall(:,:,1)=prem;
preall(:,:,2)=dprem;
preall(:,:,3)=presig;
preall(:,:,4)=dpresig;
preall(:,:,5)=dpremp;
preall(:,:,6)=dprepsig;
nc_out(preall,lons,lats,'pre','1.all vol and 2.vol with elnino case in year 1 pre; 3 and 4 are ttest for 1 and 2; 5 without elnino and 6test')
%%
sstall(:,:,1)=sstm;
sstall(:,:,2)=dsstm;
sstall(:,:,3)=sstsig;
sstall(:,:,4)=dsstsig;
sstall(:,:,5)=dsstmp;
sstall(:,:,6)=dsstpsig;
nc_out(sstall,lons,lats,'sst1','1.all vol and 2.vol with elnino case in year 1 sst; 3 and 4 are ttest for 1 and 2; 5 without elnino and 6test')
%%
uall(:,:,1)=uam;
uall(:,:,2)=duam;
uall(:,:,3)=uasig;
uall(:,:,4)=duasig;
uall(:,:,5)=duamp;
uall(:,:,6)=duapsig;
nc_out(uall,lons,lats,'u','1.all vol and 2.vol with elnino case in year 1 u; 3 and 4 are ttest for 1 and 2; 5 without elnino and 6test')
%%
vall(:,:,1)=vam;
vall(:,:,2)=dvam;
vall(:,:,3)=dvamp;

nc_out(vall,lons,lats,'v','1.all vol and 2.vol with elnino case in year 1 v; 3 without elnino')