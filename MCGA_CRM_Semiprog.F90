PROGRAM MCGA_CRM_SEMIPROG

!*************************************************************************
!
! Purpose:
! A kind of dummy AGCM to call subroutine CRM as a semi-prognostic test.
! 
! Method: 
! Read prescribed fields needed by CRM from file, use it as input
! 
! Author: CPTEC/INPE, Ocean Group contact: Paulo Santiago
! 
!*************************************************************************

USE shr_kind_mod, ONLY: r8 => shr_kind_r8, i4 => shr_kind_i4 !, r4 => shr_kind_r4
USE physconst,    ONLY: gravit, latvap, latice, cpair, rair, epsilo
USE crmdims,   ONLY: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy, crm_dt, YES3D

IMPLICIT NONE

! Variables used to store data read from file SEMIPROG_IN
INTEGER(i4) :: k, l, m, nz, nt !, plev, pplev
REAL(r8) :: iotmp, tstep, long, lati, topo, lsmk
CHARACTER(12) :: tlev
CHARACTER(12), ALLOCATABLE, DIMENSION(:) :: timestamp
REAL, ALLOCATABLE, DIMENSION(:) :: pslc, usst, vsst, cssf, clsf, ocis, &
        oces, iswf, roce, olis, oles, role, swtc, ocic, lwtc, lwbc, pres
REAL, ALLOCATABLE, DIMENSION(:,:) :: temp, umes, liqm, icem, &
        uvel, vvel, swrh, lwrh

! AGCM dummy parameters
INTEGER, PARAMETER :: lchnk = 1 ! Chunk identifier
INTEGER, PARAMETER :: icol  = 1 ! AGCM column identifier

! Variables passed from dummy AGCM (read from file)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: &
   tl,   & ! AGCM's grid temperature (k)
   ql,   & ! AGCM's grid water vapor (g/g)
   qcl,  & ! AGCM's grid cloud liquid water (g/g)
   qci,  & ! AGCM's grid cloud ice (g/g)
   ul,   & ! AGCM's grid u (m/s)
   vl,   & ! AGCM's grid v (m/s)
   pmid, & ! Global grid pressure (Pa)
   pdel, & ! Layer's pressure thicknes (Pa)
   zmid, & ! Global grid height (m)
   zint, & ! [plev+1] Global grid interface height (m)
   qrs,  & ! Shortwave heating profile
   qrl     ! Longwave heating profile
REAL(r8) :: ps   ! Global grid surface pressure (Pa)
REAL(r8) :: phis ! Global grid surface geopotential (m2/s2)
!REAL(r8) :: shf  ! Surface sensible heat flux (W/m2)
!REAL(r8) :: lhf  ! Surface latent heat flux (W/m2)
!REAL(r8) :: taux ! Surface stress zonal (N/m2)
!REAL(r8) :: tauy ! Surface stress merid (N/m2)

! AGCM parameters, read from file
REAL(r8) :: dt_gl ! Global model's time step
INTEGER :: plev, pplev ! Number of AGCM vertical levels
                       ! plev  : CAM3's middle layer
                       ! pplev : CAM3's interface layer

! Global grid tendencies output by CRM parameterization
REAL(r8), ALLOCATABLE, DIMENSION(:) :: &
   ultend,  & ! Tendency of zonal wind
   vltend,  & ! Tendency of merid wind
   qltend,  & ! Tendency of water vapor
   qcltend, & ! Tendency of cloud liquid water
   qcitend, & ! Tendency of cloud ice
   sltend     ! Tendency of static energy

! CRM variables (inout); initialized locally
REAL(KIND=r8) :: u_crm    (crm_nx,crm_ny,crm_nz) ! CRM u-wind component
REAL(r8) :: v_crm    (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
REAL(r8) :: w_crm    (crm_nx,crm_ny,crm_nz) ! CRM w-wind component
REAL(r8) :: t_crm    (crm_nx,crm_ny,crm_nz) ! CRM temperature
REAL(r8) :: q_crm    (crm_nx,crm_ny,crm_nz) ! CRM total water
REAL(r8) :: qn_crm   (crm_nx,crm_ny,crm_nz) ! CRM cloud water/ice
REAL(r8) :: qp_crm   (crm_nx,crm_ny,crm_nz) ! CRM precipitation
REAL(r8) :: qc_crm   (crm_nx,crm_ny,crm_nz) ! CRM cloud water
REAL(r8) :: qi_crm   (crm_nx,crm_ny,crm_nz) ! CRM cloud ice
REAL(r8) :: qpc_crm  (crm_nx,crm_ny,crm_nz) ! CRM precip water
REAL(r8) :: qpi_crm  (crm_nx,crm_ny,crm_nz) ! CRM precip ice
REAL(r8) :: prec_crm (crm_nx,crm_ny)        ! CRM precipitation rate
REAL(r8) :: qrs_crm  (crm_nx,crm_ny,crm_nz) ! CRM SW rad. heating
REAL(r8) :: qrl_crm  (crm_nx,crm_ny,crm_nz) ! CRM LW rad. heating

! Variables provided by AGCM's radiation parameterization (read from file)
!REAL(r8) :: fsds  ! Flux SW downwelling surface
!REAL(r8) :: fsns  ! Surface solar absorbed flux
!REAL(r8) :: fsnt  ! Net column abs solar flx at model top
!REAL(r8) :: fsut  ! Flux SW upwelling TOA
!REAL(r8) :: flwds ! Surface LW down flux
!REAL(r8) :: flns  ! Surface LW cooling (up-down) flux
!REAL(r8) :: flut  ! Outgoing LW flux at model top
!REAL(r8) :: fsntc ! Clear sky total column abs solar flux
!REAL(r8) :: fsdsc ! Clr. sky downward solar flux surface
!REAL(r8) :: flutc ! Clr. sky outgoing LW flx at model top
!REAL(r8) :: flnsc ! Clr. sky LW flux at srf (up-down)

! Radiation variables handled by CRM; initialized locally
REAL(r8) :: fsds_crm    (crm_nx,crm_ny) ! Flux SW downwelling surface
REAL(r8) :: fsns_crm    (crm_nx,crm_ny) ! Surface solar absorbed flux
REAL(r8) :: fsntoa_crm  (crm_nx,crm_ny) ! Net column abs solar flx at model top
REAL(r8) :: fsutoa_crm  (crm_nx,crm_ny) ! Flux SW upwelling TOA
REAL(r8) :: flwds_crm   (crm_nx,crm_ny) ! Surface LW down flux
REAL(r8) :: flns_crm    (crm_nx,crm_ny) ! Surface LW cooling (up-down) flux
REAL(r8) :: flut_crm    (crm_nx,crm_ny) ! Outgoing LW flux at model top
REAL(r8) :: fsntoac_crm (crm_nx,crm_ny) ! Clear sky total column abs solar flux
REAL(r8) :: fsdsc_crm   (crm_nx,crm_ny) ! Clr. sky downward solar flux surface
REAL(r8) :: flutc_crm   (crm_nx,crm_ny) ! Clr. sky outgoing LW flx at model top
REAL(r8) :: flnsc_crm   (crm_nx,crm_ny) ! Clr. sky LW flux at srf (up-down)

! CRM radiation variables (out)
REAL(r8) :: t_rad  (crm_nx,crm_ny,crm_nz) ! Rad. temperature
REAL(r8) :: qv_rad (crm_nx,crm_ny,crm_nz) ! Rad. vapor
REAL(r8) :: qc_rad (crm_nx,crm_ny,crm_nz) ! Rad. cloud water
REAL(r8) :: qi_rad (crm_nx,crm_ny,crm_nz) ! Rad. cloud ice

! Output variables from CRM
REAL(r8) :: precc  ! Convective precip rate (m/s)
REAL(r8) :: precl  ! Stratiform precip rate (m/s)
REAL(r8) :: precsc ! Convective snow rate (m/s)
REAL(r8) :: precsl ! Stratiform snow rate (m/s)

! CRM input/output - originally provided by AGCM's cloud diagnostics 
! package, but vanished before superparameterization CRM subroutine call. 
! That is, locally initialized ( all = 0. )
REAL(r8) :: cltot ! Shaded cloud fraction / diagnostic total cloud cover
REAL(r8) :: clhgh ! Shaded cloud fraction
REAL(r8) :: clmed ! Shaded cloud fraction
REAL(r8) :: cllow ! Shaded cloud fraction

! Buffer for one-column statistics (inout)
REAL(r8) :: stat_buffer(19)

! Output from CRM
REAL(r8), ALLOCATABLE, DIMENSION(:) :: &
   cld,        & ! Cloud fraction
   cldr,       & ! Cloud frac. based on -30dBZ radar reflectivity
   cldtop,     & ! Cloud top PDF
   gicewp,     & ! Ice water path
   gliqwp,     & ! Ice water path
   mc,         & ! Cloud mass flux
   mcup,       & ! Updraft cloud mass flux
   mcdn,       & ! Downdraft cloud mass flux
   mcuup,      & ! Unsat. updraft cloud mass flux
   mcudn,      & ! Unsat. downdraft cloud mass flux
   crm_qc,     & ! Mean cloud water
   crm_qi,     & ! Mean cloud ice
   crm_qs,     & ! Mean snow
   crm_qg,     & ! Mean graupel
   crm_qr,     & ! Mean rain
   tkez,       & ! TKE profile
   tkesgsz,    & ! SGS TKE profile
   flux_u,     & 
   flux_v,     & 
   flux_qt,    & ! Non-precipitating water flux
   fluxsgs_qt, & ! SGS non-precipitating water flux
   flux_qp,    & ! Precipitating water flux
   pflx,       & ! Precipitation flux
   qt_ls,      & ! Tendency of non-prec. water due to large-scale
   qt_trans,   & ! Tendency of non-prec. water due to transport
   qp_trans,   & ! Tendency of prec. water due to transport
   qp_fall,    & ! Tendency of prec. water due to fall out
   qp_evp,     & ! Tendency of prec. water due to evaporation
   qp_src,     & ! Tendency of prec. water due to conversion
   t_ls          ! Tendency of lwse due to large-scale
REAL(r8) :: &
   prectend,   & ! Col. integrated tend. in prec. water+ice (kg/m2/s)
   precstend     ! Col. integrated tend. in prec. ice (kg/m2/s) 

! Variables for surface prescription
REAL(r8) :: ocnfrac ! Area fraction of the ocean
REAL(r8) :: wnd     ! Large-scale surface wind (m/s)
REAL(r8) :: tau00   ! Large-scale surface stress (N/m2)
REAL(r8) :: bflxls  ! Large-scale surface buoyancy flux (K m/s)

! Output from CRM parameterization
REAL(r8) :: taux_crm      ! Zonal CRM surface stress perturbation (N/m2)
REAL(r8) :: tauy_crm      ! Merid. CRM surface stress perturbation (N/m2)
REAL(r8) :: z0m           ! Surface stress (N/m2)
REAL(r8) :: timing_factor ! CRM CPU efficiency: 1 means no subcycling

! Local variables
INTEGER io_out
REAL wnd_r4 ! due to incompability in z0_est function, crm subroutine
REAL(r8) :: tave ! used to calculate aver. temperature of pressure lauers
REAL(r8), ALLOCATABLE, DIMENSION(:) :: t1, t2 ! temporary variables aimed
! to store pres. layer thickness and height when transforming from MCGA
! to SP-CAM3 vertical grid

! Variables and constants used by cloud condensate estimative
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: e_s, q_s, cond
REAL(r8), PARAMETER :: e_s0 = 6.112, t0 = 273.16, rwv = rair/epsilo
!REAL(r8) :: tcel


! Roughly speaking, for code clarification purposes, we have three 
! 'classes' of variables. The first are variables used to get data stored 
! in file. The second are the variables inherited from SP-CAM's 
! tphysbc.F90, they are basically the variables passed as arguments to 
! crm subroutine. The third and last are local variables used to combine 
! the procedures of data readind and crm subroutine initialization.

!*************************************************************************
! Format specifiers, used for screen output
   11 FORMAT (A11,I4)
   12 FORMAT (A11,F11.3)
   13 FORMAT (A5,I4,A8,I4,A8,A12,A1)
   14 FORMAT (I3,F6.0,F7.2,F9.2,F8.2,X,2E12.4,2F8.2,2(X,E11.4),F3.6)
!  14 FORMAT (I3,F6.0,F7.2,F9.2,F8.2,F12.8,2F8.2,2F15.10)
   15 FORMAT (6(ES11.4,X)) !(5F15.10,F15.8)
   16 FORMAT (4(ES11.4,X)) !(4F15.10,4F7.4)

! Printing some information to tell the user what's happening
!WRITE(*,*) "Semi-prognostic test of SP-CAM superparameterization"
WRITE(*,*) "**********************"
WRITE(*,*) "SP-SAM Diagnostic Test"
WRITE(*,*) "**********************"
WRITE(*,*) 
WRITE(*,*) "CRM grid description"
WRITE(*,11) "crm_nx  = ", crm_nx
WRITE(*,11) "crm_ny  = ", crm_ny
WRITE(*,11) "crm_nz  = ", crm_nz
WRITE(*,12) "crm_dx  = ", crm_dx
WRITE(*,12) "crm_dy  = ", crm_dy 
WRITE(*,12) "crm_dt  = ", crm_dt
WRITE(*,*)  "YES3D?  = ", YES3D, "(1 - 3D CRM, 0 - 2D CRM)"
!WRITE(*,*)


! Reading time-varying input data stored in SEMIPRO_IN, step 1:
! Check number of vertical levels and time records
OPEN(31,FILE='SEMIPROG_IN',STATUS='OLD',FORM='formatted')
READ(31,*) nz, tstep, long, lati, topo, lsmk
nt = 0
DO WHILE(.TRUE.)
   READ(31,ERR=44,END=44,FMT=*) tlev, iotmp, iotmp, iotmp, iotmp, iotmp, &
           iotmp, iotmp, iotmp
   READ(31,ERR=44,END=44,FMT=*) iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, &
           iotmp, iotmp
   DO k=1,nz
      READ(31,*) iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, &
              iotmp
   END DO
   nt = nt + 1
END DO

44 CONTINUE

! Reading time-varying input data stored in SEMIPROG_IN, step 2:
! Effectivelly read the data
nt = nt + 1
REWIND(31)
READ(31,*)

ALLOCATE(timestamp(nt), pres(nz))
ALLOCATE(pslc(nt), usst(nt), vsst(nt), cssf(nt), clsf(nt), ocis(nt), &
        oces(nt), iswf(nt), roce(nt), olis(nt), oles(nt), role(nt), &
        swtc(nt),ocic(nt),lwtc(nt),lwbc(nt))
ALLOCATE(temp(nz,nt), umes(nz,nt), liqm(nz,nt), icem(nz,nt), uvel(nz,nt),&
        vvel(nz,nt), swrh(nz,nt), lwrh(nz,nt))

DO l=1,nt-1
   READ(31,*) timestamp(l), pslc(l), usst(l), vsst(l), cssf(l), clsf(l), &
           ocis(l), oces(l), iswf(l)
   READ(31,*) roce(l), olis(l), oles(l), role(l), swtc(l), ocic(l), &
           lwtc(l), lwbc(l)
   DO k=1,nz
      READ(31,*) pres(k), temp(k,l), umes(k,l), liqm(k,l), icem(k,l), &
              uvel(k,l), vvel(k,l), swrh(k,l), lwrh(k,l)
   END DO
END DO


WRITE(*,*) "Start time = ", timestamp(1)
WRITE(*,*) "Number of global time iterations = ", nt + 1
WRITE(*,*)


nt = nt - 1
plev = nz
pplev = nz + 1
dt_gl = tstep !3600._r8 !tstep
!dt_gl = 86400. !10800. !7200. !3600.
phis = topo * gravit
ocnfrac = 0. !lsmk


!*************************************************************************
! Pressure and pressure layer thickness (t1, pdel), constants
ALLOCATE(t1(plev), pmid(plev), pdel(plev))
t1(1) = pres(2) - pres(1)
DO k=2,plev-1
!  t1(k) = pres(k) - pres(k-1)
!  t1(k) = 0.5_r8 * (pres(k) - pres(k-1)) + 0.5_r8 ( pres(k+1) - pres(k) )
!  t1(k) = 0.5_r8 * ( pres(k) - pres(k-1) + pres(k+1) - pres(k) )
   t1(k) = 0.5_r8 * ( pres(k+1) - pres(k-1) )
END DO
t1(plev) = pres(plev) - pres(plev-1)
! Reversing vertical level counting
DO k=1,plev
   m = plev - k + 1
   pmid(k) = pres(m)
   pdel(k) = t1(m) * (-1)
END DO

!*************************************************************************
! Following a strongly simplified approach to estimate cloud liquid water
! and ice contents from Clausius-Clapeyron equation for saturation vapor
! pressure and the definition of saturation specific humidity. For details
! check chapter 2 of Rogers & Yau: A Short Course on Cloud Physics
ALLOCATE(e_s(nz,nt), q_s(nz,nt), cond(nz,nt))
DO l=1,nt
   DO k=1,nz
      e_s(k,l) = e_s0 * exp ( (latvap/rwv) * (1/t0 - 1/temp(k,l) ) )
      !q_s(k,l) = epsilo * e_s(k,l) / ( (pmid(k)-(1-epsilo)*e_s(k,l) ) )
      q_s(k,l) = epsilo * e_s(k,l) / pmid(k)
! The current approach assumes that the difference between the observed
! specific humidity and saturation specific humidity is the actual amount
! of saturated (liquid + ice) water
!     cond(k,l) = umes(k,l) - q_s(k,l)
      cond(k,l) = q_s(k,l) - umes(k,l)
      IF (cond(k,l).le.0.) cond(k,l) = 0.
   END DO
END DO
!PRINT*, latvap/1000._r8, rwv, epsilo

! This point we allocate the variables regarding the grid dimensions
ALLOCATE(t2(plev), zmid(plev), zint(pplev))
ALLOCATE(tl(plev), ql(plev), qcl(plev), qci(plev), &
        ul(plev), vl(plev), qrs(plev), qrl(plev), &
        ultend(plev), vltend(plev), qltend(plev), &
        qcltend(plev), qcitend(plev), sltend(plev), &
        cld(plev), cldr(plev), cldtop(plev), gicewp(plev), &
        gliqwp(plev), mc(plev), mcup(plev), mcdn(plev), &
        mcuup(plev), mcudn(plev), crm_qc(plev), crm_qi(plev), &
        crm_qs(plev), crm_qg(plev), crm_qr(plev), tkez(plev), &
        tkesgsz(plev), flux_u(plev), flux_v(plev), flux_qt(plev), &
        fluxsgs_qt(plev), flux_qp(plev), pflx(plev), qt_ls(plev), &
        qt_trans(plev), qp_trans(plev), qp_fall(plev), qp_evp(plev), &
        qp_src(plev), t_ls(plev) )

! SEMIPROG_OUT, the file which receives the semiprognostic test output
! ... declare here the file opening

OPEN(33,FILE='SEMIPROG_OUT',STATUS='unknown',ACTION='write',IOSTAT=io_out)

WRITE(33,*) timestamp(1), crm_nx, crm_ny, nz, nt, dt_gl, long, lati, topo, ocnfrac
WRITE(33,*) pmid
WRITE(33,*) pdel

!*************************************************************************
! The first time iteration uses zero mass flux profile in condensation
! subroutine

DO l=1,nt ! the main time loop
   WRITE(*,13) "Step ", l, " out of ", nt, " (time: ", timestamp(l), ")"

! Initializing vertical level heights
   tave = 0.5_r8  * ( temp(2,l) + temp(1,l) )
   t2(1) = 0.5_r8 * (rair*tave/gravit) * log( pres(1)/pres(2) )
! Pressure layer height (t2)
   DO k=2,plev-1
      tave = 0.5_r8 * ( temp(k+1,l) + temp(k-1,l) )
      t2(k) = t2(k-1) + (rair*tave/gravit) * log(pres(k-1)/pres(k))
!     write(*,*) k, tave, temp(k,l), pres(k), t2(k)
   END DO
   tave = 0.5_r8 * ( temp(plev,l) + temp(plev-1,l) )
   t2(plev) = t2(plev-1) + (rair*tave/gravit) * log(pres(plev-1)/pres(plev))
! Pressure layers interface height
   zint(pplev) = 0.0_r8
   DO k=2,pplev-1
      m = pplev - k + 1
      zint(k) = 0.5_r8 * (t2(m) + t2(m-1))
   END DO
   zint(1) = 2.0_r8 * t2(nz) - t2(plev)

! Initializing large-scale profiles
! SEMIPROG_IN profiles are vertically reversed in relation to CAM3, which 
! are originally read by CRM superparm.
   DO k=1,plev
      m = plev - k + 1
      tl(k)   = temp(m,l)
      ql(k)   = umes(m,l)
      qcl(k)  = 0._r8 !cond(m,l) !.1_r8 !liqm(m,l) !0.25 !cond(m,l) !liqm(m,l) !!+ .1_r8 ! add this to cause clouds and rain
      qci(k)  = 0._r8 !icem(m,l) !!+ .1_r8
      ul(k)   = uvel(m,l)
      vl(k)   = vvel(m,l)
      qrs(k)  = swrh(m,l)
      qrl(k)  = lwrh(m,l)
      zmid(k) = t2(m)
   END DO

! Surface pressure
   ps = pslc(l)

! L-S tendencies; according to CAM3's tphysbc.F90, they all starts null 
! every large-scale (global) time iteration
   DO k=1,plev
      ultend(k)  = 0.
      vltend(k)  = 0.
      qltend(k)  = 0.
      qcltend(k) = 0.
      qcitend(k) = 0.
      sltend(k)  = 0.
   END DO

! CRM grid variables; they are reversed in relation to CAM3's vert. grid
   DO k=1,plev
      m = plev - k + 1
      u_crm(:,:,k)   = ul(m)
      v_crm(:,:,k)   = vl(m)
      w_crm(:,:,k)   = 0. 
      t_crm(:,:,k)   = tl(m)
      q_crm(:,:,k)   = ql(m) + qcl(m) + qci(m) ! total water
      qn_crm(:,:,k)  = qcl(m) + qci(m) ! cloud water and ice
      qp_crm(:,:,k)  = 0. ! precipitation
      qc_crm(:,:,k)  = 0.
      qi_crm(:,:,k)  = 0.
      qpc_crm(:,:,k) = 0.
      qpi_crm(:,:,k) = 0.
      qrs_crm(:,:,k) = qrs(m) / (cpair * pdel(m))
      qrl_crm(:,:,k) = qrl(m) / (cpair * pdel(m))
   END DO
   prec_crm(:,:) = 0

   open(73, FILE='q_crm.txt', STATUS='unknown', ACTION='write', IOSTAT=io_out)
   write(73,*) q_crm
   close(73)

! Initializing radiation variables
   fsds_crm(:,:)    = ocis(l)
   fsns_crm(:,:)    = ocis(l) - oces(l)
   fsntoa_crm(:,:)  = iswf(l) - roce(l)
   fsutoa_crm(:,:)  = roce(l)
   flwds_crm(:,:)   = olis(l)
   flns_crm(:,:)    = oles(l)
   flut_crm(:,:)    = role(l)
   fsntoac_crm(:,:) = iswf(l) - swtc(l)
   fsdsc_crm(:,:)   = ocic(l)
   flutc_crm(:,:)   = lwtc(l)
   flnsc_crm(:,:)   = lwbc(l)

! Diagnostic total cloud cover
   precc  = 0.
   precl  = 0.
   precsc = 0.
   precsl = 0.
   cltot  = 0.
   clhgh  = 0.
   clmed  = 0.
   cllow  = 0.

   stat_buffer(:) = 0.

! According to SP-CAM3's tphysbc.F90, all those profiles are null when 
! passed to crm subroutine
   cld(:)        = 0.
   cldr(:)       = 0.
   cldtop(:)     = 0.
   gicewp(:)     = 0.
   gliqwp(:)     = 0.
   mc(:)         = 0.
   mcup(:)       = 0.
   mcdn(:)       = 0.
   mcuup(:)      = 0.
   mcudn(:)      = 0.
   crm_qc(:)     = 0.
   crm_qi(:)     = 0.
   crm_qs(:)     = 0.
   crm_qg(:)     = 0.
   crm_qr(:)     = 0.
   tkez(:)       = 0.
   tkesgsz(:)    = 0.
   flux_u(:)     = 0.
   flux_v(:)     = 0.
   flux_qt(:)    = 0.
   fluxsgs_qt(:) = 0.
   flux_qp(:)    = 0.
   pflx(:)       = 0.
   qt_ls(:)      = 0.
   qt_trans(:)   = 0.
   qp_trans(:)   = 0.
   qp_fall(:)    = 0.
   qp_evp(:)     = 0.
   qp_src(:)     = 0.
   t_ls(:)       = 0.

! Surface characteristics 
   wnd      = sqrt(ul(plev)**2 + vl(plev)**2)
   wnd_r4   = wnd
   tau00    = sqrt(usst(l)**2 + vsst(l)**2)
   bflxls   = cssf(l)/cpair + 0.61 * tl(plev) * clsf(l) / latvap
   taux_crm = 0.
   tauy_crm = 0.
   z0m      = 0.

! Checking input profiles
!  WRITE(*,*) "  1    3.   7.00 39645.50  239.11  0.00001540  -17.74   -3.49   0.0000000000  -0.0000448116"
   WRITE(*,*) " k pres    pdel   zmid      tl    ql         qcl           ul      vl    qrs         qrl"
   DO k=1,plev
      WRITE(*,14) k, pmid(k), pdel(k), zmid(k), tl(k), ql(k), &
                  qcl(k), ul(k), vl(k), qrs(k), qrl(k), q_s(k,l)
   END DO   

! The crm subrotine call
   CALL crm(lchnk, icol, &
      tl, ql, qcl, qci, ul, vl, &
      ps, pmid, pdel, phis, &
      zmid, zint, dt_gl, plev, &
      ultend, vltend, qltend, qcltend, qcitend, sltend, &
      u_crm, v_crm, w_crm, &
      t_crm, q_crm, qn_crm, qp_crm, &
      qc_crm, qi_crm, qpc_crm, qpi_crm, &
      prec_crm, qrs_crm, qrl_crm, &
      fsds_crm, fsns_crm, fsntoa_crm, fsutoa_crm, &
      flwds_crm, flns_crm, flut_crm, &
      fsntoac_crm, fsdsc_crm, flutc_crm, flnsc_crm, &
      t_rad, qv_rad, qc_rad, qi_rad, &
      precc, precl, precsc, precsl, &
      cltot, clhgh, clmed, cllow, &
      stat_buffer, cld, cldr, cldtop, &
      gicewp, gliqwp, &
      mc, mcup, mcdn, mcuup, mcudn, &
      crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, &
      tkez, tkesgsz, flux_u, flux_v, flux_qt, fluxsgs_qt, flux_qp, &
      pflx, qt_ls, qt_trans, qp_trans, qp_fall, &
      qp_evp, qp_src, t_ls, prectend, precstend, &
      ocnfrac, wnd_r4, tau00, bflxls, &
      taux_crm, tauy_crm, z0m, timing_factor)

! Writing some output to stdout
   WRITE(*,*)
   WRITE(*,*) "Check some output:"
   WRITE(*,*) "ultend      vltend      qltend      qcltend     qcitend     sltend"
   WRITE(*,15) ultend(plev), vltend(plev), qltend(plev), &
           qcltend(plev), qcitend(plev), sltend(plev)
   WRITE(*,*) "precc       precl       precsc      precsl"
   WRITE(*,16) precc, precl, precsc, precsl
   WRITE(*,*) "cltot       clhgh       clmed       cllow"
   WRITE(*,16) cltot, clhgh, clmed, cllow
   WRITE(*,*) 

! File output: SEMIPROG_OUT
   WRITE(33,*) precc, precl, precsc, precsl, cltot, clhgh, clmed, cllow, &
           taux_crm, tauy_crm, z0m, prectend, precstend
   WRITE(33,*) zmid
   WRITE(33,*) ultend
   WRITE(33,*) vltend
   WRITE(33,*) qltend
   WRITE(33,*) qcltend
   WRITE(33,*) qcitend
   WRITE(33,*) sltend
   WRITE(33,*) cld
   WRITE(33,*) cldr
   WRITE(33,*) cldtop
   WRITE(33,*) gicewp
   WRITE(33,*) gliqwp
   WRITE(33,*) mc
   WRITE(33,*) mcup
   WRITE(33,*) mcdn
   WRITE(33,*) mcuup
   WRITE(33,*) mcudn
   WRITE(33,*) crm_qc
   WRITE(33,*) crm_qi
   WRITE(33,*) crm_qs
   WRITE(33,*) crm_qg
   WRITE(33,*) crm_qr
   WRITE(33,*) tkez
   WRITE(33,*) tkesgsz
   WRITE(33,*) flux_u
   WRITE(33,*) flux_v
   WRITE(33,*) flux_qt
   WRITE(33,*) fluxsgs_qt
   WRITE(33,*) flux_qp
   WRITE(33,*) pflx
   WRITE(33,*) qt_ls
   WRITE(33,*) qt_trans
   WRITE(33,*) qp_trans
   WRITE(33,*) qp_fall
   WRITE(33,*) qp_evp
   WRITE(33,*) qp_src
   WRITE(33,*) t_ls

END DO   

END PROGRAM MCGA_CRM_SEMIPROG
