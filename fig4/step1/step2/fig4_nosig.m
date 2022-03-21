clc;clear;
addpath(genpath('F:/1_work_moreVols/utilities/')); % import auxiliary functions
%% common info
lons=ncread('F:/1_work_moreVols/data/pmip3/2ts_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc','lon');
lats=ncread('F:/1_work_moreVols/data/pmip3/2ts_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc','lat');
land0=ncread('F:/1_work_moreVols/data/2sftlf_fx_past1000_r0i0p0.nc','sftlf');
x1=106;x2=122;y1=25;y2=34; % same as Recon. region
degree=0.5;before=5;
%% m,n,wj,n_events,modelnum
%% pmip
                  pmip_q=cell2mat(struct2cell(load('../step1/p3/pmip_hus_wj_7mjjasGao.mat')));
     pmip_q(:,:,:,:,6:9)=cell2mat(struct2cell(load('../step1/p3/pmip_hus_wj_7mjjasCrow.mat')));
       pmip_q(:,:,:,:,6)=[]; % remving CSIRO due to missing wap data
             clim_pmip_q=cell2mat(struct2cell(load('../step1/p3/pmip_hus_wj_7mjjas_climGao.mat')));
clim_pmip_q(:,:,:,:,6:9)=cell2mat(struct2cell(load('../step1/p3/pmip_hus_wj_7mjjas_climCrow.mat')));
  clim_pmip_q(:,:,:,:,6)=[]; % remving CSIRO due to missing wap data

     q00(:,:,:,:,1:8)=pmip_q;     clear pmip_q;
clim_q00(:,:,:,:,1:8)=clim_pmip_q;clear clim_pmip_q;    

     w00(:,:,:,:,1:5)=cell2mat(struct2cell(load('../step1/p3/pmip_wap_wj_7mjjasGao.mat')));
     w00(:,:,:,:,6:8)=cell2mat(struct2cell(load('../step1/p3/pmip_wap_wj_7mjjasCrow.mat')));
clim_w00(:,:,:,:,1:5)=cell2mat(struct2cell(load('../step1/p3/pmip_wap_wj_7mjjas_climGao.mat')));     
clim_w00(:,:,:,:,6:8)=cell2mat(struct2cell(load('../step1/p3/pmip_wap_wj_7mjjas_climCrow.mat')));     
%% pmip4
q0=nan(181,91,11,9,11);clim_q0=nan(181,91,11,9,11);
w0=nan(181,91,11,9,11);clim_w0=nan(181,91,11,9,11);
     q0(:,:,:,1:7,1:8)=q00;
clim_q0(:,:,:,1:7,1:8)=clim_q00;
     w0(:,:,:,1:7,1:8)=w00;
clim_w0(:,:,:,1:7,1:8)=clim_w00;

     q0(:,:,:,:,9:11)=cell2mat(struct2cell(load('../step1/p4/q_wj_9mjjas.mat'))); 
clim_q0(:,:,:,:,9:11)=cell2mat(struct2cell(load('../step1/p4/q_wj_9mjjas_clim.mat'))); 
     w0(:,:,:,:,9:11)=cell2mat(struct2cell(load('../step1/p4/w_wj_9mjjas.mat')));
clim_w0(:,:,:,:,9:11)=cell2mat(struct2cell(load('../step1/p4/w_wj_9mjjas_clim.mat')));
%%
     w0=-w0;
clim_w0=-clim_w0;
%%
[m,n,wj,n_events,modelnum]=size(q0);nino_std0=nan(9,11);s20=nan(9,11);
%% criterion
load 'F:\1_work_moreVols\fig3\step2\pmip3_cesmNE_nino7.mat'; 

nino_std0(1:7,1:8)=squeeze(std(Xcomp_10(:,:,[1:5,7,8,10])));%% for PMIP3
      s20(1:7,1:8)=squeeze(Xcomp_10(before+1,:,[1:5,7,8,10]));%% for PMIP3

p4nin340=cell2mat(struct2cell(load('F:\1_work_moreVols\fig3\step2\nino_tos.mat')));
nino_std0(:,9:11)=squeeze(std(p4nin340));
      s20(:,9:11)=squeeze(p4nin340(before+1,:,:));
%% calculation start
x_a0=nan(modelnum,9);y_a0=nan(modelnum,9);all_x_a0=nan(modelnum,9);all_y_a0=nan(modelnum,9);
cor_x=nan(6,modelnum);cor_y=nan(6,modelnum);all_cor_x=nan(6,modelnum);all_cor_y=nan(6,modelnum);
for im=1:modelnum
     q=q0(:,:,:,:,im);
     w=w0(:,:,:,:,im);
clim_q=clim_q0(:,:,:,:,im);
clim_w=clim_w0(:,:,:,:,im);
nino_std=nino_std0(:,im);
      s2=s20(:,im);
      
      run figs5cal3.m
      
      x_a0(im,:)=glmq_climw;cor_x(:,im)=cor_glmq_climw;
      y_a0(im,:)=glmclimq_w;cor_y(:,im)=cor_glmclimq_w0;
      
      all_x_a0(im,:)=lani_glmq_climw;all_cor_x(:,im)=cor_lani_glmq_climw0;
      all_y_a0(im,:)=lani_glmclimq_w;all_cor_y(:,im)=cor_lani_glmclimq_w0;       
end
all_x_a0(5,[2,3,6])=nan;
all_y_a0(5,[3,6])=nan;
%% -----------------------------------------------------------------------------------------------FigureOut
moden={'CCSM4','IPSL','FGOALS-s2','MRI','BCC','MPI','GISS','MIROC','EC-Earth3','MIROC-ES2L','MRI-ESM2'};

figure
x_int=4/90*10e-6;y_int=x_int*8;x_intm=6/90*10e-6;y_intm=x_intm*8;ass=6;
%% ALL
color=[0 0 0];
load all.mat
x0=all_x(~isnan(all_x));
y0=all_y(~isnan(all_y));
x_amme=mean(x0);
y_amme=mean(y0);

%% T-test
alpha=0.1;
x_sig=ttest(x0,0,'Alpha',alpha,'Dim',1); 
y_sig=ttest(y0,0,'Alpha',alpha,'Dim',1);

if x_sig == 1   
    line([x_amme+x_intm x_amme],[y_amme y_amme],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
end

if y_sig == 1 
    line([x_amme x_amme],[y_amme-y_intm y_amme],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
end

hal=scatter(x_amme,y_amme,'MarkerFaceColor',color,'MarkerEdgeColor',color,'SizeData',150,'Marker','o');hold on;
text(x_amme+0.000001,y_amme,'ALL','FontName','Helvetica','fontsize',10,'FontWeight','bold','color',color);hold on;

clear x0 y0 x_amme y_amme;

color=[1 99/255 71/255];
x_a=nanmean(x_a0,2);y_a=nanmean(y_a0,2);
%% MME
x0=x_a0(~isnan(x_a0));y0=y_a0(~isnan(y_a0));
x_amme=mean(x0);y_amme=mean(y0);

%% T-test
alpha=0.1;
x_sig=ttest(x0(~isnan(x0)),0,'Alpha',alpha,'Dim',1);clear x0;
y_sig=ttest(y0(~isnan(y0)),0,'Alpha',alpha,'Dim',1);clear y0;

if x_sig == 1
    line([x_amme+x_intm x_amme],[y_amme y_amme],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
end
 
if y_sig == 1
    line([x_amme x_amme],[y_amme-y_intm y_amme],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
end
%%
hv=scatter(x_a(:),y_a(:),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',color,'SizeData',50);hold on;
text(x_a(:)-0.0000002,y_a(:)+0.000003,moden,'FontName','Helvetica','fontsize',10,'FontWeight','bold','color',color);hold on;

scatter(x_amme,y_amme,'MarkerFaceColor',color,'MarkerEdgeColor',color,'SizeData',150,'Marker','o');hold on;
text(x_amme+0.000001,y_amme,'EL','FontName','Helvetica','fontsize',10,'FontWeight','bold','color',color);hold on;
%%
color=[0 0.75 0.75];
clear x_a y_a cor_x cor_y;
x_a=nanmean(all_x_a0,2);
y_a=nanmean(all_y_a0,2);
cor_x=all_cor_x;
cor_y=all_cor_y;

%% MME
x0=all_x_a0(~isnan(all_x_a0));
y0=all_y_a0(~isnan(all_y_a0));
x_amme=mean(x0);
y_amme=mean(y0);

%% T-test
alpha=0.1;
x_sig=ttest(x0(~isnan(x0)),0,'Alpha',alpha,'Dim',1);clear x0;
y_sig=ttest(y0(~isnan(y0)),0,'Alpha',alpha,'Dim',1);clear y0;

if x_sig == 1  
    line([x_amme+x_intm x_amme],[y_amme y_amme],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
end
 
if y_sig == 1  
    line([x_amme x_amme],[y_amme-y_intm y_amme],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
end

for iv=1:11
    if(x_a(iv) > 0 && x_a(iv) > cor_x(6,iv)) || (x_a(iv) < 0 && x_a(iv) < cor_x(3,iv))
        line([x_a(iv)+x_int x_a(iv)],[y_a(iv) y_a(iv)],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
    end
end
for iv=1:11
    if(y_a(iv) > 0 && y_a(iv) > cor_y(6,iv)) || (y_a(iv) < 0 && y_a(iv) < cor_y(3,iv))
    line([x_a(iv) x_a(iv)],[y_a(iv)-y_int y_a(iv)],'Marker','*','Color',color,'MarkerSize',ass,'LineStyle','none','LineWidth',1);hold on;
    end
end

h=scatter(x_a(:),y_a(:),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',color,'SizeData',50);hold on;
text(x_a(:)-0.0000002,y_a(:)+0.000002,moden,'FontName','Helvetica','fontsize',10,'FontWeight','bold','color',color);hold on;

scatter(x_amme,y_amme,'MarkerFaceColor',color,'MarkerEdgeColor',color,'SizeData',150,'Marker','o');hold on;
text(x_amme+0.000001,y_amme,'nEL','FontName','Helvetica','fontsize',10,'FontWeight','bold','color',color);hold on;

line([0 0],[-1 1],'color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','--'); hold on;
line([-1 1],[0 0],'color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','--'); hold on;
xlabel('thermodynamical (q'' w, Pa s^-^1)')
ylabel('dynamical (q w'', Pa s^-^1)')
le=legend([hal,hv(1),h(1)],'all vol cases','vol with El Nino','vol without El Nino');set(le,'EdgeColor',[1 1 1],'Location','NorthWest'); 

axis([-1.9 0.86 -7.2 16]*10e-6)
axis square;
set(gca,'FontName','Helvetica','fontsize',14,'TickDir','out','FontWeight','bold');grid on;box on;