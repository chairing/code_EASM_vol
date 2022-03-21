clc;clear;
%% read   index
path='F:\fig81600-971\glmnino3.4\';
glmslpmatfile=dir([path '*.mat']);
for i=1:numel(glmslpmatfile)
    glmslpmatfile(i).name
    bb=load([path glmslpmatfile(i).name]);
    glmrst0(:,i) = bb.nin34';
end
glmslp=zeros(11988,10);
k=0;
% for i=[1,10,3,8,9,7,5,4,2,6] % !!!
for i=[1,10,3,8,9,2,7,4,5,6] % !!!
    k=k+1;
    glmrst(:,k)=glmrst0(:,i);
end
clearvars -except glmrst
[ll,vv]=size(glmrst);lv=ll/12;before = 5;after = 5;wj=before+after+1;
%%
glm_ano=nan(ll,vv);
for k=1:ll
    glm_ano(k,:)=glmrst(k,:)-nanmean(glmrst);
end
%% DJF mean
glm_naa=nan(lv,vv);
for k=1:lv-1
    glm_naa(k,:)=nanmean(glm_ano((k-1)*12+12:(k-1)*12+12+2,:));
end
load F:/SEA/tropvolgaoandcrowley.mat  % G:\3\PMIP\ensemblevol\code\Step1selectTropicalVol.m
%% 
n1(:,1)=851:1849;
n1(:,2:11)=glm_naa;
%% find vol year in Gao
volsize= 0;
key_dates_Gao=find(tropvolfirstyear(:,3)>volsize & tropvolfirstyear(:,1)>1470+before); 
[erup_years,events]=intersect(n1(:,1),tropvolfirstyear(key_dates_Gao,1)); 
%% events number
n_events = length(erup_years);
%% Superposed Epoch Analysis
Xevents = nan(wj,n_events,11);  Xcomp_10=nan(wj,n_events+1,10); % due to Crowley vols = 7
%% Perform SEA:
for iv=[2,4,5,6]
    for i=1:n_events
        Xevents(:,i,iv) = n1(events(i)-before:events(i)+after,iv);
        Xevents(:,i,iv) = (Xevents(:,i,iv)-nanmean(Xevents(1:before,i,iv)));%/max(abs(Xevents(:,i,iv))); % remove mean over "before" of window
    end
end
Xcomp(:,[1,3,4,5]) = squeeze(mean(Xevents(:,:,[2,4,5,6]),2)); 

Xcomp_10(:,1:6,[1,3,4,5]) = Xevents(:,:,[2,4,5,6]); 

%% find vol year in Crowley
% volsize_Cro= 0.10;
% key_dates_Cro=find(tropvolfirstyear(:,2)>volsize_Cro & tropvolfirstyear(:,1)>851+before & tropvolfirstyear(:,1)~=1600); 
% [erup_years_Cro,events_Cro]=intersect(n1(:,1),tropvolfirstyear(key_dates_Cro,1)); 
key_dates_Cro=[1576;1594;1674;1696;1809;1816;1835]; 

[erup_years_Cro,events_Cro]=intersect(n1(:,1),key_dates_Cro); 

%% events number
n_events_Cro = length(erup_years_Cro);
%% Superposed Epoch Analysis
Xevents_Cor = nan(wj,n_events_Cro,11);          
%% Perform SEA:
for iv=7:11
    for i=1:n_events_Cro
        Xevents_Cor(:,i,iv) = n1(events_Cro(i)-before:events_Cro(i)+after,iv);
        Xevents_Cor(:,i,iv) = (Xevents_Cor(:,i,iv)-nanmean(Xevents_Cor(1:before,i,iv)));%/max(abs(Xevents_Cor(:,i,iv))); % remove mean over "before" of window
    end
end
Xcomp(:,6:10) = squeeze(mean(Xevents_Cor(:,:,7:11),2)); 

Xcomp_10(:,:,6:10) = Xevents_Cor(:,:,7:11); 
%% find vol year in ipsl 
%erup_years_ipsl=[1213,1258,1278,1286,1452,1600,1641,1809,1815,1835];
erup_years_ipsl=[1600; 1641; 1673; 1809; 1815; 1835 ];
[erup_years_ipsl,events_ipsl]=intersect(n1(:,1),erup_years_ipsl);
%% events number
n_events_ipsl = length(erup_years_ipsl);
%% Superposed Epoch Analysis
Xevents_Cor = nan(wj,n_events_ipsl,11);          
%% Perform SEA:
for iv=3
    for i=1:n_events_ipsl
        Xevents_Cor(:,i,iv) = n1(events_ipsl(i)-before:events_ipsl(i)+after,iv);
        Xevents_Cor(:,i,iv) = (Xevents_Cor(:,i,iv)-nanmean(Xevents_Cor(1:before,i,iv)));%/max(abs(Xevents_Cor(:,i,iv))); % remove mean over "before" of window
    end
end
Xcomp(:,2) = squeeze(mean(Xevents_Cor(:,:,3),2)); 

Xcomp_10(:,1:6,2) = Xevents_Cor(:,:,3); 


%% cesm
% load nino34_all.mat;

glmcesm=cell2mat(struct2cell(load('F:/fig81600-971/newnino34_cesm.mat')));

[ll,vv]=size(glmcesm);lv=ll/12; 
%%
glm_ano_ce=nan(ll,vv);
for k=1:ll
    glm_ano_ce(k,:)=glmcesm(k,:)-nanmean(glmcesm);
end
%% DJF mean
glm_naa_cesm=nan(lv,vv);
for k=1:lv-1
    glm_naa_cesm(k,:)=mean(glm_ano_ce((k-1)*12+12:(k-1)*12+12+2,:));
end
%% 
n1(:,1)=851:1849;
n1(:,12:26)=glm_naa_cesm;
%% Superposed Epoch Analysis
Xevents = nan(wj,n_events,26);          
%% Perform SEA:
for iv=12:26
    for i=1:n_events
        Xevents(:,i,iv) = n1(events(i)-before:events(i)+after,iv);
        Xevents(:,i,iv) = (Xevents(:,i,iv)-nanmean(Xevents(1:before,i,iv)));%/max(abs(Xevents(:,i,iv))); % remove mean over "before" of window
    end
end
Xcomp(:,11:25) = squeeze(mean(Xevents(:,:,12:26),2)); 

Xcomp_10(:,1:6,11:25) = Xevents(:,:,12:26); 

save pmip3_cesmNE_nino7.mat Xcomp_10;

% Xcomp_out=cell2mat(struct2cell(load('E:\1_work_dir\1_work_dir\work_dir\work dir\pmip\recon\model\data\pmip3_cesmNE_nino10.mat')));
% 
% Xcomp_10-Xcomp_out
% eq(Xcomp_10,Xcomp_out);for iv=1:25;plot(ans(:,:,iv));end
stop;
figure;
subplot(1,2,1);
%% color design
colors=[87/255 1 226/255;
    99/255 132/255 123/255;
    144/255 1 111/255;
    159/255 87/255 1;
    1 238/255 11/255; 
    255/255 0 0; 
    1 6/255 1; 
    1 96/255 99/255; 
    1 194/255 130/255; 
    0 0 1];

rr=1:11;
for iv=1:10
hm(iv)=plot(rr,Xcomp(rr,iv),'LineWidth',2,'Color',colors(iv,:));hold on;
end;
plot(rr,mean(Xcomp(rr,1:5),2),'LineWidth',3,'Color','r');hold on;
plot(rr,mean(Xcomp(rr,6:10),2),'LineWidth',3,'Color','b');hold on;

set(gca,'xlim',[1 10],'xtick',[1:1:10],'XTickLabel',[-5:1:5-1],'ylim',[-0.3 1],'ytick',[-0.2:0.2:1] ,'FontSize',12,'FontName',......
     'Times New Roman','FontWeight','bold','XMinorTick','on','linewidth',2,'TickDir','out','PlotBoxAspectRatio',[1 1 1]);
line([1 11],[0 0],'linestyle',':','color',[0 0 0]); box on;hold on;
line([6 6],[-3 2],'linestyle',':','color',[0 0 0]); box on;hold on;
legend1=legend(hm,'CCSM4','IPSL','FGOALS','MRI','BCC','MPI','HadCM3','GISS','CSIRO','MIROC');
set(legend1,'EdgeColor',[1 1 1],'Location','NorthWest');
xlabel('Year from eruption','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
title('PMIP','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
subplot(1,2,2);

he=plot(rr,Xcomp(rr,11:25),'LineWidth',1,'Color',[0.5 0.5 0.5]);hold on;
plot(rr,mean(Xcomp(rr,11:25),2),'LineWidth',2,'Color','r');hold on;

set(gca,'xlim',[1 10],'xtick',[1:1:10],'XTickLabel',[-5:1:5-1],'ylim',[-0.3 1],'ytick',[-0.2:0.2:1] ,'FontSize',12,'FontName',......
     'Times New Roman','FontWeight','bold','XMinorTick','on','linewidth',2,'TickDir','out','PlotBoxAspectRatio',[1 1 1]);
line([1 11],[0 0],'linestyle',':','color',[0 0 0]); box on;hold on;
line([6 6],[-3 2],'linestyle',':','color',[0 0 0]); box on;hold on;
xlabel('Year from eruption','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
title('CESM','FontWeight','bold','FontSize',12,'FontName','Times New Roman');

