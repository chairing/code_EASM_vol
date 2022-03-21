load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
begin
;---Read in netCDF file
  aa = addfile("./step3/pre.nc","r")
  pre = aa->pre(:,:,:)  
  lon= aa->lon(:)
  lat= aa->lat(:)
; 
  lon@units = "degrees_east"
  lat@units = "degrees_north"
  pre!1 = "lat"
  pre&lat= lat
  pre!2 ="lon"
  pre&lon= lon
    if (any(isnan_ieee(pre))) then
      value = 1.e20
      replace_ieeenan (pre, value, 0)
      pre@_FillValue = value
  end if
;  
  bb = addfile("./step3/sst1.nc","r")
  ts = bb->sst1(:,:,:)
  ts!1 = "lat"
  ts&lat= lat
  ts!2 ="lon"
  ts&lon= lon

; define missing value
  if (any(isnan_ieee(ts))) then
      value = 1.e20
      replace_ieeenan (ts, value, 0)
      ts@_FillValue = value
  end if

  cc = addfile("./step3/u.nc","r")
  dd = addfile("./step3/v.nc","r")

  uwind = cc->u(:,:,:)
  vwind = dd->v(:,:,:)

  ;---Read in netCDF file
  aa2 = addfile("../figNR2_pmip/step4/pre.nc","r")
  pre2 = aa2->pre(:,:,:)  
  pre2!1 = "lat"
  pre2&lat= lat
  pre2!2 ="lon"
  pre2&lon= lon
    if (any(isnan_ieee(pre2))) then
      value = 1.e20
      replace_ieeenan (pre2, value, 0)
      pre2@_FillValue = value
  end if
;  
  bb2 = addfile("../figNR2_pmip/step4/sst.nc","r")
  ts2 = bb2->sst(:,:,:)
  ts2!1 = "lat"
  ts2&lat= lat
  ts2!2 ="lon"
  ts2&lon= lon
; define missing value
  if (any(isnan_ieee(ts2))) then
      value = 1.e20
      replace_ieeenan (ts2, value, 0)
      ts2@_FillValue = value
  end if

  cc2 = addfile("../figNR2_pmip/step4/u.nc","r")
  dd2 = addfile("../figNR2_pmip/step4/v.nc","r")

  uwind2 = cc2->u(:,:,:)
  vwind2 = dd2->v(:,:,:)
;;************************************************
;; plot parameters that do not change
;;************************************************  

wks  = gsn_open_wks("ps","figs")         ; Opens a ps file

  resources=True
  resources@gsnDraw =False
  resources@gsnFrame = False
  resources@gsnLeftString=""
  resources@gsnRightString="" 
  resources@gsnAddCyclic = False
  resources@tiXAxisString=""

  resources@vpWidthF  = 0.65
  resources@vpHeightF = 0.2
  vcres=resources
  res2Q=resources

  resources@cnFillOn = True
  resources@cnRasterSmoothingOn =True
  resources@cnLinesOn = False
  resources@cnLineLabelsOn  =False
  resources@lbBoxLinesOn=False

  resources@gsnLeftStringFontHeightF = 0.02
  resources@gsnRightStringFontHeightF =0.02
  resources@cnLevelSelectionMode="ManualLevels"
  resources@pmLabelBarWidthF   =0.26     ; The width of the label
  resources@pmLabelBarHeightF =0.05
  resources@pmLabelBarOrthogonalPosF = 0.16
  resources@pmLabelBarParallelPosF = -0.28
  resources@lbOrientation="Horizontal"
  resources@lbLabelFontHeightF =0.0175 
  resources@lbLabelStride=5
  respre=resources
  respre2=resources
  respre3=resources
  respre4=resources

  resources@mpShapeMode = "FreeAspect"
  resources@mpMinLonF            =  50                ;-- min longitude
  resources@mpCenterLonF         =  180
  resources@mpMaxLonF            =  270               ;-- max longitude
  resources@mpMinLatF            =  -10                 ;-- min latitude
  resources@mpMaxLatF            =  45

  resources@tmYLMode ="Explicit"
  resources@tmYLValues=(/-10,0,15,30,45/)
  resources@tmYLLabels=(/"10~S~o~N~S","EQ","15~S~o~N~N","30~S~o~N~N","45~S~o~N~N"/)
  resources@tmYLLabelFontHeightF = 0.015
  resources@tmYLLabelDeltaF      =-0.6
  resources@tmXBMode             = "Explicit"       ; Set labels for X axis
  resources@tmXBValues           = (/60,120,180,240/)
  resources@tmXBLabels           = (/"60~S~o~N~E","120~S~o~N~E","180","120~S~o~N~W"/)
  resources@tmXBLabelFontHeightF = 0.015
  resources@tmXBLabelDeltaF      =-0.6
  
  resources@mpOutlineOn                 = True
  resources@mpGeophysicalLineColor    ="gray50"         ;地图边界线颜色
  resources@mpGeophysicalLineThicknessF      = 1

;resources@mpDataSetName = "Earth..4.1"
;resources@mpDataBaseVersion = "MediumRes"
;resources@mpOutlineSpecifiers = (/"changjiang"/)
;resources@mpCountyLineColor = (/"red"/)
;resources@mpCountyLineThicknessF = 1

; for t test 
  res2Q@cnFillOn = True ; 打开
  res2Q@cnLinesOn = False ; 关闭线条
  res2Q@cnLevelSelectionMode = "ManualLevels"
  res2Q@cnMinLevelValF = -0.5 ; 最小值
  res2Q@cnMaxLevelValF = 0.5 ; 最大值
  res2Q@cnLevelSpacingF = 0.5 ; 间距
  res2Q@lbLabelBarOn   = False
  res2Q@cnMonoFillColor = True               
  res2Q@cnMonoFillPattern = False              
  res2Q@cnMonoFillScale   = False       
  res2Q@cnFillPatterns    = (/ -1, 17,17, -1/)
  res2Q@cnFillDotSizeF=0.001
  
  res2Q@cnFillScales      = (/0.45,0.9,0.9,0.45/)
  res2Q@cnInfoLabelOn        = False          ; turn off info label
  res2Q@cnLineLabelsOn   = False
;************************************************
; resource list for first data array
;************************************************
  
  resources@gsnLeftString = "a all vol cases" 
  resources@gsnRightString= "PMIP"
  resources@cnFillPalette   =   "BlueDarkRed18"
  resources@cnMinLevelValF=-1
  resources@cnMaxLevelValF=1
  resources@cnLevelSpacingF=0.04
  resources@lbLabelBarOn   = False
  plot_sst = gsn_csm_contour_map(wks,ts(0,:,:),resources) 
  
  plot_sstt =gsn_csm_contour(wks,ts(2,:,:),res2Q)
  overlay(plot_sst,plot_sstt)
;************************************************
; resource list for second data array
;************************************************

  respre@cnFillPalette   =   "precip_diff_12lev"
  respre@cnMinLevelValF=-0.6
  respre@cnMaxLevelValF=0.6
  respre@cnLevelSpacingF=0.02
  respre@lbOrientation="Vertical"
  respre@lbLabelBarOn   = False
  respre@cnInfoLabelOn     = False

  plot_pre = gsn_csm_contour(wks,pre(0,:,:),respre)
  plot_pret =gsn_csm_contour(wks,pre(2,:,:),res2Q)
  overlay(plot_pre,plot_pret)
  overlay(plot_sst,plot_pre )

;************************************************
; resource list for wind data array
;************************************************
  vcres@vcGlyphStyle            = "Fillarrow"    ; turn on curly vectors
  vcres@vcMinFracLengthF        = 0.2 
  vcres@vcMinDistanceF          = 0.018
  vcres@vcRefMagnitudeF         = 2.             ; define vector ref mag
  vcres@vcRefLengthF            = 0.04             ; define length of vec ref

  vcres@vcRefAnnoOn             = False
  ;vcres@vcRefAnnoPerimOn          = False          ;参考矢量的边框
  ;vcres@vcRefAnnoString1          = "2 m/s"            ;not draw the word"reference vector"
  ;vcres@vcRefAnnoOrthogonalPosF = -.15            ; move ref vector up
 ; vcres@vcRefAnnoParallelPosF = 0.1  ; 水平位置
  ;vcres@vcRefAnnoString2On =True
  vcres@lbLabelBarOn            = False  
  vcres@vcRefAnnoFontHeightF =0.012
  vcres@vcMonoFillArrowFillColor = False ; in different colors
  ;u_t=where(ttua.ne.0,2,0) 
  
  colors=(/"grey60","black"/)
  vcres@vcLevelPalette   = colors; "vegetation_ClarkUdrought_severity" ;"MPL_cool" ;"BlueDarkRed18"             ; set color map
  
 
  vector= gsn_csm_vector_scalar(wks,uwind(0,:,:),vwind(0,:,:),uwind(2,:,:),vcres)         ;20
  ;vector= gsn_csm_vector(wks,uwind(:,:),vwind(:,:),vcres)         ;20
  overlay(plot_sst,vector)
;************************************************
; resource list for second data array
;************************************************
  resources@gsnLeftString = "c vol with El Nin~H-13V2F35~D~FV-2H3~o (24 cases)" 
  resources@gsnRightString= "PMIP"
  resources@cnMinLevelValF=-1
  resources@cnMaxLevelValF=1
  resources@cnLevelSpacingF=0.04
  resources@lbLabelBarOn   = False
  plot_sst2 = gsn_csm_contour_map_ce(wks,ts(1,:,:),resources) 
  
  plot_sstt2 =gsn_csm_contour(wks,ts(3,:,:),res2Q)
  overlay(plot_sst2,plot_sstt2)
  respre2@cnFillPalette   =   "precip_diff_12lev"
  respre2@cnMinLevelValF=-0.6
  respre2@cnMaxLevelValF=0.6
  respre2@cnLevelSpacingF=0.02
  respre2@lbLabelBarOn   = False
  respre2@cnInfoLabelOn     = False

  plot_pre2 = gsn_csm_contour(wks,pre(1,:,:),respre2)
  plot_pret2 =gsn_csm_contour(wks,pre(3,:,:),res2Q)
  overlay(plot_pre2,plot_pret2)
  overlay(plot_sst2,plot_pre2 )

  vcres@vcRefMagnitudeF         = 2.             ; define vector ref mag
  vector2= gsn_csm_vector_scalar(wks,uwind(1,:,:),vwind(1,:,:),uwind(3,:,:),vcres)         ;20
  ;vector= gsn_csm_vector(wks,uwind(:,:),vwind(:,:),vcres)         ;20
  overlay(plot_sst2,vector2)

  resources@gsnLeftString = "b vol without El Nin~H-13V2F35~D~FV-2H3~o (68 cases)" 
  resources@gsnRightString= "PMIP"
  resources@cnMinLevelValF=-1
  resources@cnMaxLevelValF=1
  resources@cnLevelSpacingF=0.04
  resources@lbLabelBarOn   = False
  plot_sst3 = gsn_csm_contour_map_ce(wks,ts(4,:,:),resources) 
  
  plot_sstt3 =gsn_csm_contour(wks,ts(5,:,:),res2Q)
  overlay(plot_sst3,plot_sstt3)
  respre3@cnFillPalette   =   "precip_diff_12lev"
  respre3@cnMinLevelValF=-0.6
  respre3@cnMaxLevelValF=0.6
  respre3@cnLevelSpacingF=0.02
  respre3@lbLabelBarOn   = False
  respre3@cnInfoLabelOn     = False

  plot_pre3 = gsn_csm_contour(wks,pre(4,:,:),respre3)

  plot_pret3 =gsn_csm_contour(wks,pre(5,:,:),res2Q)
  overlay(plot_pre3,plot_pret3)
  overlay(plot_sst3,plot_pre3 )

  vcres@vcRefMagnitudeF         = 2.             ; define vector ref mag
  vector3= gsn_csm_vector_scalar(wks,uwind(4,:,:),vwind(2,:,:),uwind(5,:,:),vcres)         ;20
  ;vector= gsn_csm_vector(wks,uwind(:,:),vwind(:,:),vcres)         ;20
  overlay(plot_sst3,vector3)

  resources@gsnLeftString = "d El Nin~H-13V2F35~D~FV-2H3~o without vol (1218 cases)" 
  resources@gsnRightString= "PMIP"
  resources@cnMinLevelValF=-1
  resources@cnMaxLevelValF=1
  resources@cnLevelSpacingF=0.04
  resources@lbLabelBarOn   = True
  plot_sst4 = gsn_csm_contour_map_ce(wks,ts2(0,:,:),resources) 
  
  plot_sstt4 =gsn_csm_contour(wks,ts2(1,:,:),res2Q)
  overlay(plot_sst4,plot_sstt4)

  respre4@cnFillPalette   =   "precip_diff_12lev"
  respre4@cnMinLevelValF=-0.6
  respre4@cnMaxLevelValF=0.6
  respre4@cnLevelSpacingF=0.02
  respre4@lbLabelBarOn   = True
  respre4@pmLabelBarOrthogonalPosF = -0.25
  respre4@pmLabelBarParallelPosF = 0.17

  plot_pre4 = gsn_csm_contour(wks,pre2(0,:,:),respre4)

  plot_pret4 =gsn_csm_contour(wks,pre2(1,:,:),res2Q)
  overlay(plot_pre4,plot_pret4)
  overlay(plot_sst4,plot_pre4 )

  vcres@vcRefMagnitudeF         = 2.             ; define vector ref mag
  vector4= gsn_csm_vector_scalar(wks,uwind2(0,:,:),vwind2(0,:,:),uwind2(1,:,:),vcres)         ;20
  ;vector= gsn_csm_vector(wks,uwind(:,:),vwind(:,:),vcres)         ;20
  overlay(plot_sst4,vector4)

;>============================================================<
;                     add box 
;>------------------------------------------------------------<
  ypts = (/ 19.0,  19.0,  9.0,   9.0, 19.0/)
  xpts = (/120.0, 150.0,150.0, 120.0,120.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "red"                     ; color of lines
  resp@gsLineThicknessF = 2.0                       ; thickness of lines
  resp@gsLineLabelString= "box 1"                    ; adds a line label string
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.
  duma = new(4,graphic)
   dum = new(4,graphic)  
  dumpp = new(4,graphic) 
  dumpp4 = new(4,graphic) 
; draw each line separately. Each line must contain two points.
  do i = 0 , 3
   duma(i)=gsn_add_polyline(wks,plot_sst,xpts(i:i+1),ypts(i:i+1),resp)
    dum(i)=gsn_add_polyline(wks,plot_sst2,xpts(i:i+1),ypts(i:i+1),resp)   
   dumpp(i)=gsn_add_polyline(wks,plot_sst3,xpts(i:i+1),ypts(i:i+1),resp)  
   dumpp4(i)=gsn_add_polyline(wks,plot_sst4,xpts(i:i+1),ypts(i:i+1),resp)  
       
  end do 

  ypts = (/ 32.0,  32.0,  22.0,   22.0, 32.0/)
  xpts = (/130.0, 160.0,160.0, 130.0,130.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "red"                     ; color of lines
  resp@gsLineThicknessF = 2.0                       ; thickness of lines
  resp@gsLineLabelString= "box 2"                    ; adds a line label string
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.
  duma2 = new(4,graphic) 
   dum2 = new(4,graphic)  
  dump2 = new(4,graphic) 
  dump4 = new(4,graphic) 
; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    duma2(i)=gsn_add_polyline(wks,plot_sst,xpts(i:i+1),ypts(i:i+1),resp)
     dum2(i)=gsn_add_polyline(wks,plot_sst2,xpts(i:i+1),ypts(i:i+1),resp)  
    dump2(i)=gsn_add_polyline(wks,plot_sst3,xpts(i:i+1),ypts(i:i+1),resp)  
    dump4(i)=gsn_add_polyline(wks,plot_sst4,xpts(i:i+1),ypts(i:i+1),resp)  
  end do 
;************************************************
; resource list for common data array
;************************************************
  xres = True
  xres@txFontHeightF = 0.012
  gsn_text_ndc(wks,"K",0.515,0.307,xres)
  gsn_text_ndc(wks,"mm day~S~-1~N~",0.745,0.307,xres)

  res  =True

  res@gsnFrame              = False
  res@amOrthogonalPosF      = -0.65
  res@vpWidthF             = 1
  res@vpHeightF            = 1
  res@gsnPanelYF           =(/0.9,0.7/)
  res@gsnPanelXF           =(/0.9,0.7/)

  gsn_panel(wks,(/plot_sst,plot_sst2,plot_sst3,plot_sst4/),(/2,2/),res)
  end
