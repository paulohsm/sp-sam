PROGRAM crm_semiprog

 USE crmdims, ONLY: crm_nx, crm_ny, crm_nz

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
#include 'semiprog_params.inc'
 integer :: i, k, m, lchnk, ppcnst, ixcldliq, ixcldice
! exclusive variables for 'radctl'
 real(r8), intent(in) :: lwup(pcols)        ! Surface longwave up flux
 real(r8) :: emis(pcols,pver)                  ! Cloud longwave emissivity
 real(r8), pointer, dimension(:,:) :: qcwat, tcwat, lcwat, cld ! physics buffer fields to compute tendencies for cloud condensation package
 real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
 real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
 real(r8) :: coszrs(pcols)                     ! Cosine solar zenith angle
 real(r8), intent(in) :: asdir(pcols)       ! Albedo: shortwave, direct
 real(r8), intent(in) :: asdif(pcols)       ! Albedo: shortwave, diffuse
 real(r8), intent(in) :: aldir(pcols)       ! Albedo: longwave, direct
 real(r8), intent(in) :: aldif(pcols)       ! Albedo: longwave, diffuse
 real(r8) :: pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
 integer :: nmxrgn(pcols)                      ! Number of maximally overlapped regions
 real(r8), intent(inout) :: fsns(pcols)     ! Surface solar absorbed flux
 real(r8), intent(inout) :: fsnt(pcols)     ! Net column abs solar flux at model top
 real(r8), intent(inout) :: flns(pcols)     ! Srf longwave cooling (up-down) flux
 real(r8), intent(inout) :: flnt(pcols)     ! Net outgoing lw flux at model top
 real(r8), intent(inout) :: qrs(pcols,pver) ! Shortwave heating rate
 real(r8), intent(inout) :: qrl(pcols,pver) ! Longwave  heating rate
 real(r8), intent(inout) :: flwds(pcols)    ! Surface longwave down flux
 real(r8) :: rel(pcols,pver)                   ! Liquid cloud particle effective radius
 real(r8) :: rei(pcols,pver)                   ! Ice effective drop size (microns)
 real(r8), intent(inout) :: sols(pcols)     ! Direct beam solar rad. onto srf (sw)
 real(r8), intent(inout) :: soll(pcols)     ! Direct beam solar rad. onto srf (lw)
 real(r8), intent(inout) :: solsd(pcols)    ! Diffuse solar radiation onto srf (sw)
 real(r8), intent(inout) :: solld(pcols)    ! Diffuse solar radiation onto srf (lw
 real(r8), intent(in) :: landfrac(pcols)    ! land fraction
 real(r8), intent(out) :: fsds(pcols)       ! Surface solar down flux
! radctl's ifdef CRM
 real(r8) :: fsntoa(pcols)  ! Net solar flux at TOA
 real(r8) :: fsntoac(pcols) ! Clear sky net solar flux at TOA
 real(r8) :: fsdsc(pcols)   ! Clear sky flux Shortwave Downwelling Surface
 real(r8) :: flwdsc(pcols)  ! Clear-sky Surface longwave down flux
 real(r8) :: fsntc(pcols)   ! Clear sky total column abs solar flux
 real(r8) :: fsnsc(pcols)   ! Clear sky surface abs solar flux
 real(r8) :: fsutoa(pcols)  ! Flux Shortwave Upwelling TOA
 real(r8) :: fsutoac(pcols) ! Clear sky Flux Shortwave Upwelling TOA
 real(r8) :: flut(pcols)    ! Upward flux at top of model
 real(r8) :: flutc(pcols)   ! Upward Clear Sky flux at top of model
 real(r8) :: flntc(pcols)   ! Clear sky lw flux at model top
 real(r8) :: flnsc(pcols)   ! Clear sky lw flux at srf (up-down)
 real(r8) :: solin(pcols)   ! Solar incident flux
! other variables
 real(r8), intent(inout) :: u_crm  (pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: v_crm  (pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: w_crm  (pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: t_crm  (pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: q_crm  (pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: qn_crm  (pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: qp_crm (pcols, crm_nx, crm_ny, crm_nz)
 real(r8) :: qc_crm (pcols, crm_nx, crm_ny, crm_nz)
 real(r8) :: qi_crm (pcols, crm_nx, crm_ny, crm_nz)
 real(r8) :: qpc_crm(pcols, crm_nx, crm_ny, crm_nz)
 real(r8) :: qpi_crm(pcols, crm_nx, crm_ny, crm_nz)
 real(r8) :: prec_crm(pcols, crm_nx, crm_ny)
 real(r8), intent(inout) :: qrs_crm(pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: qrl_crm(pcols, crm_nx, crm_ny, crm_nz)
 real(r8), intent(inout) :: fsds_crm(pcols,crm_nx,crm_ny)   ! Flux Shortwave Downwelling Surface
 real(r8), intent(inout) :: fsns_crm(pcols,crm_nx,crm_ny)   ! Surface solar absorbed flux
 real(r8), intent(inout) :: fsntoa_crm(pcols,crm_nx,crm_ny) ! Net column abs solar flux at model top
 real(r8), intent(inout) :: fsutoa_crm(pcols,crm_nx,crm_ny) ! Flux Shortwave Upwelling TOA
 real(r8), intent(inout) :: flwds_crm(pcols,crm_nx,crm_ny)  ! Surface longwave down flux
 real(r8), intent(inout) :: flns_crm(pcols,crm_nx,crm_ny)   ! Srf longwave cooling (up-down) flux
 real(r8), intent(inout) :: flut_crm(pcols,crm_nx,crm_ny)   ! Outgoing lw flux at model top
 real(r8), intent(inout) :: fsntoac_crm(pcols,crm_nx,crm_ny)! Clear sky total column abs solar flux
 real(r8), intent(inout) :: fsdsc_crm(pcols,crm_nx,crm_ny)  ! Clear sky downard solar flux surface
 real(r8), intent(inout) :: flutc_crm(pcols,crm_nx,crm_ny)  ! Clear sky outgoing lw flux at model top
 real(r8), intent(inout) :: flnsc_crm(pcols,crm_nx,crm_ny)  ! Clear sky lw flux at srf (up-down)
 real(r8) :: t_rad(pcols, crm_nx, crm_ny, crm_nz) ! rad temperuture
 real(r8) :: qv_rad(pcols, crm_nx, crm_ny, crm_nz) ! rad vapor
 real(r8) :: qc_rad(pcols, crm_nx, crm_ny, crm_nz) ! rad cloud water
 real(r8) :: qi_rad(pcols, crm_nx, crm_ny, crm_nz) ! rad cloud ice
 real(r8), intent(out) :: precc(pcols)                  ! Convective-scale preciptn rate
 real(r8), intent(out) :: precl(pcols)                  ! Large-scale precipitation rate
 real(r8), intent(out) :: precsc(pcols)                 ! C.S. snowfall rate
 real(r8), intent(out) :: precsl(pcols)                 ! L.S. snowfall rate
 real(r8) :: cltot(pcols)                      ! Diagnostic total cloud cover
 real(r8) :: clhgh(pcols)                      !       "     hgh  cloud cover
 real(r8) :: clmed(pcols)                      !       "     mid  cloud cover
 real(r8) :: cllow(pcols)                      !       "     low  cloud cover
 real(r8) :: stat_buffer(pcols,19) ! One-column CRM statistics for the ARM diagnostics
 real(r8) :: cldn(pcols,pver) ! cloud top pdf
 real(r8) :: cldr(pcols,pver) ! cloud fraction based on -30dBZ radar reflectivity
 real(r8) :: cldtop(pcols,pver)
 real(r8) :: gicewp(pcols,pver)      ! grid-box cloud ice water path
 real(r8) :: gliqwp(pcols,pver)      ! grid-box cloud liquid water path
 real(r8) :: mctot(pcols,pver)          ! total cloud mass flux
 real(r8) :: mcup(pcols,pver)           ! cloud updraft mass flux
 real(r8) :: mcdn(pcols,pver)           ! cloud downdraft mass flux
 real(r8) :: mcuup(pcols,pver)          ! unsaturated updraft mass flux
 real(r8) :: mcudn(pcols,pver)          ! unsaturated downdraft mass flux
 real(r8) :: crm_qc(pcols,pver)         ! cloud water
 real(r8) :: crm_qi(pcols,pver)         ! cloud ice
 real(r8) :: crm_qs(pcols,pver)         ! snow
 real(r8) :: crm_qg(pcols,pver)         ! graupel
 real(r8) :: crm_qr(pcols,pver)         ! rain
 real(r8) :: tkez(pcols,pver)     ! tke profile
 real(r8) :: tkesgsz(pcols,pver)     ! sgs tke profile
 real(r8) :: flux_u(pcols,pver)        ! x-momentum flux
 real(r8) :: flux_v(pcols,pver)        ! y-momentum flux
 real(r8) :: flux_qt(pcols,pver)        ! nonprecipitating water flux
 real(r8) :: fluxsgs_qt(pcols,pver)     ! sgs nonprecipitating water flux
 real(r8) :: flux_qp(pcols,pver)        ! precipitating water flux
 real(r8) :: precflux(pcols,pver)       ! precipitation flux
 real(r8) :: qt_ls(pcols,pver)        ! water tendency due to large-scale
 real(r8) :: qt_trans(pcols,pver)     ! nonprecip water tendency due to transport
 real(r8) :: qp_trans(pcols,pver)     ! precip water tendency due to transport
 real(r8) :: qp_fall(pcols,pver)      ! precip water tendency due to fall-out
 real(r8) :: qp_evp(pcols,pver)       ! precip water tendency due to evap
 real(r8) :: qp_src(pcols,pver)       ! precip water tendency due to conversion
 real(r8) :: t_ls(pcols,pver)        ! tendency of crm's liwse due to large-scale
 real(r8) :: prectend(pcols) ! tendency in precipitating water and ice
 real(r8) :: precstend(pcols) ! tendency in precipitating ice
 real(r8), intent(IN) :: ocnfrac(pcols)                ! land fraction
 real(r8) :: wnd  ! surface wnd
 real(r8) :: tau00  ! surface stress
 real(r8) :: bflx   ! surface buoyancy flux (Km/s)
 real(r8) :: taux_crm(pcols)  ! zonal CRM surface stress perturbation
 real(r8) :: tauy_crm(pcols)  ! merid CRM surface stress perturbation
 real(r8) :: z0m(pcols)
 real :: timing_factor(pcols) ! factor for crm cpu-usage: 1 means no subcycling

! The 'state' type, as found in CAM3.1 source code, specifically in file
! cam1/models/atm/cam/src/physics/cam1/physics_types.F90
! The strategy is to provide values for derived types bellow from values
! stored in files
 TYPE state
     integer                                :: &
         lchnk,     &! chunk index
         ncol        ! number of active columns
     real(r8), dimension(pcols)             :: &
         ps,        &! surface pressure
         phis        ! surface geopotential
     real(r8), dimension(pcols,pver)        :: &
         t,         &! temperature (k)
         u,         &! zonal wind (m/s)
         v,         &! meridional wind (m/s)
         pmid,      &! midpoint pressure (Pa)
         pdel,      &! layer thickness (Pa)
         lnpmid,    &! ln(pmid)
         zm          ! geopotential height above surface at midpoints (m)
     real(r8), dimension(pcols,pver,ppcnst) :: &
         q           ! constituent mixing ratio (kg/kg moist or dry air depending on type)
     real(r8), dimension(pcols,pver+1)      :: &
         lnpint,  &! ln(pint)
         zi          ! geopotential height above surface at interfaces (m)
 END TYPE state
! The 'ptend' type, as found in CAM3.1 source code specifically in file
! cam1/models/atm/cam/src/physics/cam1/physics_types.F90
 TYPE ptend
     real(r8), dimension(pcols,pver)        :: &
         s, &! heating rate (J/kg/s)
         u, &! u momentum tendency (m/s/s)
         v, &! v momentum tendency (m/s/s)
     real(r8), dimension(pcols,pver,ppcnst) :: &
         q  &! constituent tendencies (kg/kg/s)
 END TYPE ptend
!-----------------------------------------------------------------------
! Some CRM arguments are calculated by radiation driver. Thus, it's
!+necessary to call radctl before calling crm itself. 2013-04-21
 CALL radctl (lchnk,        & ! chunk identifier, comes from...
              ncol,         & ! number of atmospheric columns
              lwup,         & ! Longwave up flux at surface
              emis,         & ! Cloud emissivity
              state%pmid,   & ! Model level pressures, comes from CAM3's physics_types
              state%pint,   & ! Model interface pressures, ''
              state%lnpmid, & ! Natural log of pmid, ''
              state%lnpint, & ! Natural log of pint, ''
              state%t,      &
              state%q,      &
              cld,          & ! Fractional cloud cover, comes from...
              cicewp,       & ! in-cloud cloud ice water path, comes from...
              cliqwp,       & ! in-cloud cloud liquid water path, comes from...
              coszrs,       & ! Cosine solar zenith angle, comes from...
              asdir,        & ! albedo shortwave direct, comes from...
              asdif,        & ! albedo shortwave diffuse, comes from...
              aldir,        & ! albedo longwave direct, comes from...
              aldif,        & ! albedo longwave diffuse, comes from...
              pmxrgn,       & ! Maximum values of pmid for each, comes from...
              nmxrgn,       & ! Number of maximally overlapped regions, comes from...
              fsns,         & ! Surface absorbed solar flux, comes from...
              fsnt,         & ! Net column abs solar flux at model top, comes from...
              flns,         & ! Srf longwave cooling (up-down) flux, comes from...
              flnt,         & ! Net outgoing lw flux at model top, comes from...
              qrs,          & ! Solar heating rate, comes from...
              qrl,          & ! Longwave cooling rate, comes from...
              flwds,        & ! Surface down longwave flux, comes from...
              rel,          & ! liquid effective drop size (microns), comes from...
              rei,          & ! ice effective drop size (microns), comes from...
              sols,         & ! Downward solar rad onto surface (sw direct), comes from...
              soll,         & ! Downward solar rad onto surface (lw direct), comes from...
              solsd,        & ! Downward solar rad onto surface (sw diffuse), comes from...
              solld,        & ! Downward solar rad onto surface (lw diffuse), comes from...
              landfrac,     & ! land fraction, comes from...
              state%zm,     &
              state,        & ! A derived type of kind physics_state, comes from CAM3's physics_types
              fsds,         & ! Flux Shortwave Downwelling Surface, comes from...
!#ifdef CRM
              fsntoa,       &
              fsntoac,      &
              fsdsc,        &
              flwdsc,       &
              fsntc,        &
              fsnsc,        &
              fsutoa,       &
              fsutoac,      &
              flut,         &
              flutc,        &
              flntc,        &
              flnsc,        &
              solin,        &
              .TRUE.,       &
              doabsems,     &
              dosw,         &
              dolw          &
!#endif
             )

! Initializing crm subroutine's arguments
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

 DO i=1,ncol
     prec_crm(i,:,:) = 0.
 END DO

 DO k=1,crm_nz
     DO i=1,ncol
         m = pver-k+1
         u_crm(i,:,:,k) = state%u(i,m)
         v_crm(i,:,:,k) = state%v(i,m)
         w_crm(i,:,:,k) = 0.
         t_crm(i,:,:,k) = state%t(i,m)
         q_crm(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
         qn_crm(i,:,:,k) = state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
         qp_crm(i,:,:,k) = 0.
         qc_crm(i,:,:,k) = 0.
         qi_crm(i,:,:,k) = 0.
         qpc_crm(i,:,:,k) = 0.
         qpi_crm(i,:,:,k) = 0.
         qrs_crm(i,:,:,k) = (qrs(i,m))/cpair
         qrl_crm(i,:,:,k) = (qrl(i,m))/cpair
     END DO
 END DO

 CRM_NCOL: DO i = 1,ncol ! run crm for every agcm column
     tau00 = SQRT(taux(i)**2 + tauy(i)**2)
     wnd = SQRT(state%u(i,pver)**2 + state%v(i,pver)**2)
     bflx = shf(i)/cpair + 0.61*state%t(i,pver)*lhf(i)/latvap
!-----------------------------------------------------------------------
! The following subroutine call was extracted from 'tphysbc.F90'
! Variable's descriptions from 'crm/crm.F'
    CALL crm (lchnk,           &! chunk identifier
         i,                     &! column identifier
         state%t(i,:),          &!t global grid temperature (K)
         state%q(i,:,1),        &!t global grid water vapor (g/g)
         state%q(i,:,ixcldliq), &!t global grid cloud liquid water (g/g)
         state%q(i,:,ixcldice), &! global grid cloud ice (g/g)
         state%u(i,:),          &! global grid u (m/s)
         state%v(i,:),          &! global grid v (m/s)
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
         u_crm(i,:,:,:),        &!* crm u-wind component
         v_crm(i,:,:,:),        &!* crm v-wind component
         w_crm(i,:,:,:),        &!* crm w-wind component
         t_crm(i,:,:,:),        &!* crm temperature
         q_crm(i,:,:,:),        &!* crm total water
         qn_crm(i,:,:,:),       &!* crm cloud water/ice
         qp_crm(i,:,:,:),       &!* crm precipitation
         qc_crm(i,:,:,:),       &!* crm cloud water
         qi_crm(i,:,:,:),       &!* crm cloud ice
         qpc_crm(i,:,:,:),      &!* crm precip water
         qpi_crm(i,:,:,:),      &!* crm precip ice
         prec_crm(i,:,:),       &!* crm precipitation rate
         qrs_crm(i,:,:,:),      &!* crm sw rad. heating
         qrl_crm(i,:,:,:),      &!* crm lw rad. heating
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
         prectend(i),           &!c column integrated tendency in precipitating water+ice (kg/m2/s)
         precstend(i),          &!c column integrated tendency in precipitating ice (kg/m2/s)
         ocnfrac(i),            &!* area fraction of the ocean
         wnd,                   &!* large-scale surface wind (m/s)
         tau00,                 &!* large-scale surface stress (N/m2)
         bflx,                  &!* large-scale surface buoyance flux (K m/s)
         taux_crm(i),           &!* zonal crm surface stress perturbation (N/m2)
         tauy_crm(i),           &!* merid crm surface stress perturbation (N/m2)
         z0m(i),                &!* surface stress (N/m2)
         timing_factor(i)       &!¹ crm cpu efficiency
     )
 END DO CRM_NCOL ! DO i = 1,ncol

END PROGRAM crm_semiprog
