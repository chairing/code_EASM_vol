%%  
[m,n,wj,n_events]=size(q);
glmq=zeros(wj,n_events);clim_glmq=zeros(wj,n_events);
glmw=zeros(wj,n_events);clim_glmw=zeros(wj,n_events);

for iv=1:n_events
    for ij=1:wj
        ks=0;
        for i=1:m
            for j=1:n
                if (~isnan(q(i,j,ij,iv)) && land0(i,j)>0 && lons(i)>=x1 && lons(i)<=x2 && lats(j)>=y1 && lats(j) <=y2)
                    glmq(ij,iv)=glmq(ij,iv)+q(i,j,ij,iv)*cos(lats(j)/180*pi);
               clim_glmq(ij,iv)=clim_glmq(ij,iv)+clim_q(i,j,ij,iv)*cos(lats(j)/180*pi);
                    ks=ks+cos(lats(j)/180*pi);
                end
            end
        end
        glmq(ij,iv)=glmq(ij,iv)/ks;clim_glmq(ij,iv)=clim_glmq(ij,iv)/ks;
    end
end

for iv=1:n_events
    for ij=1:wj
        ks=0;
        for i=1:m
            for j=1:n
                if (~isnan(w(i,j,ij,iv)) && land0(i,j)>0 && lons(i)>=x1 && lons(i)<=x2 && lats(j)>=y1 && lats(j) <=y2)
                    glmw(ij,iv)=glmw(ij,iv)+w(i,j,ij,iv)*cos(lats(j)/180*pi);
                    clim_glmw(ij,iv)=clim_glmw(ij,iv)+clim_w(i,j,ij,iv)*cos(lats(j)/180*pi);                  
                    ks=ks+cos(lats(j)/180*pi);
                end
            end
        end
        glmw(ij,iv)=glmw(ij,iv)/ks;clim_glmw(ij,iv)=clim_glmw(ij,iv)/ks;  
    end
end

kin=0;kde=0;
qnino=nan(wj,n_events);qlani=nan(wj,n_events);
wnino=nan(wj,n_events);wlani=nan(wj,n_events);
clim_qnino=nan(wj,n_events);clim_qlani=nan(wj,n_events);
clim_wnino=nan(wj,n_events);clim_wlani=nan(wj,n_events);

for ik=1:n_events
    if s2(ik,1) > degree*nino_std(ik)
        kin=kin+1;
        qnino(:,ik)=glmq(:,ik);
        wnino(:,ik)=glmw(:,ik);
        
        clim_qnino(:,ik)=clim_glmq(:,ik);
        clim_wnino(:,ik)=clim_glmw(:,ik);
    end
     if s2(ik,1) < degree*nino_std(ik)
        kde=kde+1;
        qlani(:,ik)=glmq(:,ik);
        wlani(:,ik)=glmw(:,ik);
        
        clim_qlani(:,ik)=clim_glmq(:,ik);
        clim_wlani(:,ik)=clim_glmw(:,ik);
     end
end; clear nino_std s2;
kin
kde
%%
nino_q=squeeze(qnino(before+2,:));
nino_w=squeeze(wnino(before+2,:));
nino_clim_q=squeeze(clim_qnino(before+2,:));
nino_clim_w=squeeze(clim_wnino(before+2,:));

glmq_climw=nino_q.*nino_clim_w;
glmclimq_w=nino_clim_q.*nino_w;

lani_q=squeeze(qlani(before+2,:));
lani_w=squeeze(wlani(before+2,:));
lani_clim_q=squeeze(clim_qlani(before+2,:));
lani_clim_w=squeeze(clim_wlani(before+2,:));

lani_glmq_climw=lani_q.*lani_clim_w;
lani_glmclimq_w=lani_clim_q.*lani_w;
%%
     glmq_climw0=qnino.*clim_wnino;
glmclimq_w0=clim_qnino.*wnino;

     lani_glmq_climw0=qlani.*clim_wlani;
lani_glmclimq_w0=clim_qlani.*wlani;
%% sig
aa=glmq_climw0(:,:);en1=aa(~isnan(aa));en1=reshape(en1,[wj,length(en1)/wj]);[wj,nsn_events]=size(en1);
cor_glmq_climw=repeat_cc3(en1,nsn_events,wj);clear en1 aa;
%
aa=glmclimq_w0(:,:);en1=aa(~isnan(aa));en1=reshape(en1,[wj,length(en1)/wj]);[wj,nsn_events]=size(en1);
cor_glmclimq_w0=repeat_cc3(en1,nsn_events,wj);clear en1 aa;
%
aa=lani_glmq_climw0(:,:);en1=aa(~isnan(aa));en1=reshape(en1,[wj,length(en1)/wj]);[wj,nsn_events]=size(en1);
cor_lani_glmq_climw0=repeat_cc3(en1,nsn_events,wj);clear en1 aa;
%
aa=lani_glmclimq_w0(:,:);en1=aa(~isnan(aa));en1=reshape(en1,[wj,length(en1)/wj]);[wj,nsn_events]=size(en1);
cor_lani_glmclimq_w0=repeat_cc3(en1,nsn_events,wj);clear en1 aa;
% %%
% glmq_climw=squeeze(nanmean(glmq_climw0(before+2,:),2));
% glmclimq_w=squeeze(nanmean(glmclimq_w0(before+2,:),2));
% lani_glmq_climw=squeeze(nanmean(lani_glmq_climw0(before+2,:),2));
% lani_glmclimq_w=squeeze(nanmean(lani_glmclimq_w0(before+2,:),2));