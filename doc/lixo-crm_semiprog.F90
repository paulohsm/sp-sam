PROGRAM crm_semiprog
IMPLICIT NONE
!-----------------------------------------------------------------------
!
! Purpose:
! A kind of dummy AGCM to call subroutine CRM as a semi-prognostic test.
! 
! Method: 
! Read prescribed fields needed by CRM from file, use it as input
! 
! Author: CPTEC/INPE, Ocean Group contact: Paulo Santiago
! 
!-----------------------------------------------------------------------

 integer :: i, lchnk, ncol, ppcnst, pver
 real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

! The 'state' type, as found in CAM3.1 source code, specifically in file
! cam1/models/atm/cam/src/physics/cam1/physics_types.F90
! The strategy is to provide values for derived types bellow from values
! stored in files
 TYPE state
     integer                                :: &
         lchnk,     &! chunk index
         ncol        ! number of active columns
     real(r8), dimension(pcols)             :: &
!        lat,       &! latitude (radians)
!        lon,       &! longitude (radians)
         ps,        &! surface pressure
!        psdry,     &! dry surface pressure
         phis        ! surface geopotential
     real(r8), dimension(pcols,pver)        :: &
         t,         &! temperature (k)
         u,         &! zonal wind (m/s)
         v,         &! meridional wind (m/s)
!        s,         &! dry static energy
!        omega,     &! vertical pressure velocity (Pa/s)
         pmid,      &! midpoint pressure (Pa)
!        pmiddry,   &! midpoint pressure dry (Pa)
         pdel,      &! layer thickness (Pa)
!        pdeldry,   &! layer thickness dry (Pa)
!        lnpmid,    &! ln(pmid)
!        lnpmiddry, &! log midpoint pressure dry (Pa)
!        exner,     &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
         zm          ! geopotential height above surface at midpoints (m)
     real(r8), dimension(pcols,pver,ppcnst) :: &
         q           ! constituent mixing ratio (kg/kg moist or dry air depending on type)
     real(r8), dimension(pcols,pver+1)      :: &
!        pint,      &! interface pressure (Pa)
!        pintdry,   &! interface pressure dry (Pa)
!        lnpint,    &! ln(pint)
!        lnpintdry, &! log interface pressure dry (Pa)
         zi          ! geopotential height above surface at interfaces (m)
!    real(r8), dimension(pcols)             :: &
!        te_ini,    &! vertically integrated total (kinetic + static) energy of initial state
!        te_cur,    &! vertically integrated total (kinetic + static) energy of current state
!        tw_ini,    &! vertically integrated total water of initial state
!        tw_cur     &! vertically integrated total water of new state
!    integer :: count! count of values with significant energy or water imbalances
 END TYPE state
! 
 TYPE ptend
     real(r8), dimension(pcols,pver)        :: &
         s, &! heating rate (J/kg/s)
         u, &! u momentum tendency (m/s/s)
         v, &! v momentum tendency (m/s/s)
     real(r8), dimension(pcols,pver,ppcnst) :: &
         q  &! constituent tendencies (kg/kg/s)
 END TYPE ptend
! Initializing crm subroutine's arguments
 state%lchnk = 1
 state%ncol = 1
 lchnk = state%lchnk
 ncol  = state%ncol 
! <read 'state' variables from file>
 ptend%q(:,:,1) = 0.
 ptend%q(:,:,ixcldliq) = 0.
 ptend%q(:,:,ixcldice) = 0.
 ptend%s(:,:) = 0.
 precc(:) = 0.
 precl(:) = 0.
 precsc(:) = 0.
 precsl(:) = 0.
 cltot(:) = 0.
 clhgh(:) = 0.
 clmed(:) = 0.
 cllow(:) = 0.
 cld(:,:) = 0.
 cldr(:,:) = 0.
 cldtop(:,:) = 0.
 gicewp(:,:) = 0.
 gliqwp(:,:) = 0.
 stat_buffer(:,:) = 0.
 mctot(:,:) = 0.
 mcup(:,:) = 0.
 mcdn(:,:) = 0.
 mcuup(:,:) = 0.
 mcudn(:,:) = 0.
 crm_qc(:,:) = 0.
 crm_qi(:,:) = 0.
 crm_qs(:,:) = 0.
 crm_qg(:,:) = 0.
 crm_qr(:,:) = 0.
 flux_qt(:,:) = 0.
 flux_u(:,:) = 0.
 flux_v(:,:) = 0.
 fluxsgs_qt(:,:) = 0.
 tkez(:,:) = 0.
 tkesgsz(:,:) = 0.
 flux_qp(:,:) = 0.
 precflux(:,:) = 0.
 qt_ls(:,:) = 0.
 qt_trans(:,:) = 0.
 qp_trans(:,:) = 0.
 qp_fall(:,:) = 0.
 qp_evp(:,:) = 0.
 qp_src(:,:) = 0.
 z0m(:) = 0.
 taux_crm(:) = 0.
 tauy_crm(:) = 0.
 t_ls(:,:) = 0.

 DO i = 1,ncol ! run crm for every agcm column
!-----------------------------------------------------------------------
! Where do each argument in crm subroutine come from?
! lchnk -> physics_types module, through physics_state type
! ncol -> physics_types module, through physics_state type
! state%* -> physics_types module
! ztodt -> given as argument
! pver -> given as 
! ptend%* -> physics_types module, through physics_ptend type
! *_crm ->
!-----------------------------------------------------------------------
! The following subroutine call was extracted from 'tphysbc.F90'
! Variable's descriptions from 'crm/crm.F'
 CALL crm (lchnk,                 &! chunk identifier
           i,                     &! column identifier
           state%t(i,:),          &!t global grid temperature (K)
           state%q(i,:,1),        &!t global grid water vapor (g/g)
           state%q(i,:,ixcldliq), &!t global grid cloud liquid water (g/g)
!#ifdef SP_DIR_NS
!          state%q(i,:,ixcldice), &!t global grid cloud ice (g/g)
!          state%v(i,:),          &!t global grid u (m/s)
!          state%u(i,:),          &!t global grid v (m/s)
!#ELSE
           state%q(i,:,ixcldice), &! global grid cloud ice (g/g)
           state%u(i,:),          &! global grid u (m/s)
           state%v(i,:),          &! global grid v (m/s)
!#ENDIF
           state%ps(i),           &!t global grid surface pressure (Pa)
           state%pmid(i,:),       &!t global grid pressure (Pa)
           state%pdel(i,:),       &!t layer's pressure thickness (Pa)
           state%phis(i),         &!t global grid surface geopotential (m2/s2)
           state%zm(i,:),         &!t global grid height (m)
           state%zi(i,:),         &!t global grid interface height (m)
           ztodt,                 &! global model's time step
           pver,                  &! number of levels
           ptend%u(i,:),          &! tendency of zonal wind
           ptend%v(i,:),          &! tendency of merid wind
           ptend%q(i,:,1),        &!* tendency of water vapor
           ptend%q(i,:,ixcldliq), &!* tendency of cloud liquid water
           ptend%q(i,:,ixcldice), &!* tendency of cloud ice
           ptend%s(i,:),          &!* tendency of static energy
           u_crm(i,:,:,:),        &! crm u-wind component
           v_crm(i,:,:,:),        &! crm v-wind component
           w_crm(i,:,:,:),        &! crm w-wind component
           t_crm(i,:,:,:),        &! crm temperature
           q_crm(i,:,:,:),        &! crm total water
           qn_crm(i,:,:,:),       &! crm cloud water/ice
           qp_crm(i,:,:,:),       &! crm precipitation
           qc_crm(i,:,:,:),       &! crm cloud water
           qi_crm(i,:,:,:),       &! crm cloud ice
           qpc_crm(i,:,:,:),      &! crm precip water
           qpi_crm(i,:,:,:),      &! crm precip ice
           prec_crm(i,:,:),       &! crm precipitation rate
           qrs_crm(i,:,:,:),      &! crm sw rad. heating
           qrl_crm(i,:,:,:),      &! crm lw rad. heating
           fsds_crm(i,:,:),       &! downward solar flux at surface
           fsns_crm(i,:,:),       &! surface solar absorbed flux
           fsntoa_crm(i,:,:),     &! net downward solar flux at top of atmosphere
           fsutoa_crm(i,:,:),     &! upwelling shortwave flux at toa
           flwds_crm(i,:,:),      &! surface longwave down flux
           flns_crm(i,:,:),       &! srf longwave cooling (up-down) flux
           flut_crm(i,:,:),       &! outgoing lw flux at model top
           fsntoac_crm(i,:,:),    &! clearsky net solar flux at toa
           fsdsc_crm(i,:,:),      &! clearsky downward solar flux at surface
           flutc_crm(i,:,:),      &! clearsky outgoing lw flux at model top
           flnsc_crm(i,:,:),      &! clearsky srf longwave cooling (up-down) flux
           t_rad(i,:,:,:),        &! rad temperature
           qv_rad(i,:,:,:),       &! rad vapor 
           qc_rad(i,:,:,:),       &! rad cloud water
           qi_rad(i,:,:,:),       &! rad cloud ice
           precc(i),              &!* convective precip rate (m/s)
           precl(i),              &!* stratiform precip rate (m/s)
           precsc(i),             &!* convective snow rate (m/s)
           precsl(i),             &!* stratiform snow rate (m/s)
           cltot(i),              &!* shaded cloud fraction
           clhgh(i),              &!* shaded cloud fraction
           clmed(i),              &!* shaded cloud fraction
           cllow(i),              &!* shaded cloud fraction
           stat_buffer(i,:),      &!* buffer for on-column statistics
           cld(i,:),              &!* cloud fraction
           cldr(i,:),             &!* cloud fraction based on -30dBz radar reflectivity
           cldtop(i,:),           &!* cloud top pdf
           gicewp(i,:),           &!* ice water path
           gliqwp(i,:),           &!* ice water path
           mctot(i,:),            &!* cloud mass flux
           mcup(i,:),             &!* updraft cloud mass flux
           mcdn(i,:),             &!* downdraft cloud mass flux
           mcuup(i,:),            &!* unsat updraft cloud mass flux
           mcudn(i,:),            &!* unsat downdraft cloud mass flux
           crm_qc(i,:),           &!* mean cloud water
           crm_qi(i,:),           &!* mean cloud ice
           crm_qs(i,:),           &!* mean snow
           crm_qg(i,:),           &!* mean graupel
           crm_qr(i,:),           &!* mean rain
           tkez(i,:),             &!* tke profile
           tkesgsz(i,:),          &!* sgs tke profile
           flux_u(i,:),           &!*
           flux_v(i,:),           &!*
           flux_qt(i,:),          &!* nonprecipitation water flux
           fluxsgs_qt(i,:),       &!* sgs nonprecipitation water flux
           flux_qp(i,:),          &!* precipitation water flux
           precflux(i,:),         &!* precipitation flux
           qt_ls(i,:),            &!* tendency of nonprec water due to large-scale
           qt_trans(i,:),         &!* tendency of nonprec water due to transport
           qp_trans(i,:),         &!* tendency of prec water due to transport
           qp_fall(i,:),          &!* tendency of prec water due to fall-out
           qp_evp(i,:),           &!* tendency of prec water due to evp
           qp_src(i,:),           &!* tendency of prec water due to conversion
           t_ls(i,:),             &!* tendency of lwse due to large-scale
           prectend(i),           &! column integrated tendency in precipitating water+ice (kg/m2/s)
           precstend(i),          &! column integrated tendency in precipitating ice (kg/m2/s)
           ocnfrac(i),            &! area fraction of the ocean
           wnd,                   &! large-scale surface wind (m/s)
           tau00,                 &! large-scale surface stress (N/m2)
           bflx,                  &! large-scale surface buoyance flux (K m/s)
           taux_crm(i),           &!* zonal crm surface stress perturbation (N/m2)
           tauy_crm(i),           &!* merid crm surface stress perturbation (N/m2)
           z0m(i),                &!* surface stress (N/m2)
           timing_factor(i)       &! crm cpu efficiency
 )

 END DO ! DO i = 1,ncol

END PROGRAM crm_semiprog
