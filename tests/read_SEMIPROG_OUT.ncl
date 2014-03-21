load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

begin
 ni = 32
 nj = 1
 nk = 27
 fin = "../SEMIPROG_OUT"

 u_crm = fbindirread(fin,0,(/nk,nj,ni/),"double")
 v_crm = fbindirread(fin,1,(/ni,nj,nk/),"double")
 w_crm = fbindirread(fin,2,(/ni,nj,nk/),"double")
 t_crm = fbindirread(fin,3,(/ni,nj,nk/),"double")
 q_crm = fbindirread(fin,4,(/ni,nj,nk/),"double")
 qn_crm = fbindirread(fin,5,(/ni,nj,nk/),"double")
 qp_crm = fbindirread(fin,6,(/ni,nj,nk/),"double")
 qc_crm = fbindirread(fin,7,(/ni,nj,nk/),"double")
 qi_crm = fbindirread(fin,8,(/ni,nj,nk/),"double")
 qpc_crm = fbindirread(fin,9,(/ni,nj,nk/),"double")
 qpi_crm = fbindirread(fin,10,(/ni,nj,nk/),"double")
 qrs_crm = fbindirread(fin,11,(/ni,nj,nk/),"double")
 qrl_crm = fbindirread(fin,12,(/ni,nj,nk/),"double")

 wks = gsn_open_wks("png","t_crm")
 gsn_define_colormap(wks,"gui_default")

 res           = True
 res@cnFillOn  = True
 res@vpWidthF  = 0.5
 res@vpHeightF = 0.3

 plot = gsn_csm_contour(wks,t_crm(:,0,:),res)

 wks = gsn_open_wks("png","u_crm")
 plot = gsn_csm_contour(wks,u_crm(:,0,:),res)

 wks = gsn_open_wks("png","v_crm")
 plot = gsn_csm_contour(wks,v_crm(:,0,:),res)

 wks = gsn_open_wks("png","w_crm")
 plot = gsn_csm_contour(wks,w_crm(:,0,:),res)

 wks = gsn_open_wks("png","q_crm")
 plot = gsn_csm_contour(wks,q_crm(:,0,:),res)

;...

 wks = gsn_open_wks("png","qrs_crm")
 plot = gsn_csm_contour(wks,qrs_crm(:,0,:),res)

 wks = gsn_open_wks("png","qrl_crm")
 plot = gsn_csm_contour(wks,qrl_crm(:,0,:),res)
end
