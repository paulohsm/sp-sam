 PROGRAM crm_semiprog

 USE shr_kind_mod, ONLY: r8 => shr_kind_r8 !, r4 => shr_kind_r4
 USE physconst,    ONLY: gravit, latvap, latice, cpair, rair
 USE crmdims,      ONLY: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy, crm_dt, YES3D

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
 REAL(r8) :: shf  ! Surface sensible heat flux (W/m2)
 REAL(r8) :: lhf  ! Surface latent heat flux (W/m2)
 REAL(r8) :: taux ! Surface stress zonal (N/m2)
 REAL(r8) :: tauy ! Surface stress merid (N/m2)

! AGCM parameters, read from file
 REAL(r8) :: dt_gl ! Global model's time step
 INTEGER :: plev, pplev   ! Number of AGCM vertical levels
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
 REAL(r8) :: fsds  ! Flux SW downwelling surface
 REAL(r8) :: fsns  ! Surface solar absorbed flux
 REAL(r8) :: fsnt  ! Net column abs solar flx at model top
 REAL(r8) :: fsut  ! Flux SW upwelling TOA
 REAL(r8) :: flwds ! Surface LW down flux
 REAL(r8) :: flns  ! Surface LW cooling (up-down) flux
 REAL(r8) :: flut  ! Outgoing LW flux at model top
 REAL(r8) :: fsntc ! Clear sky total column abs solar flux
 REAL(r8) :: fsdsc ! Clr. sky downward solar flux surface
 REAL(r8) :: flutc ! Clr. sky outgoing LW flx at model top
 REAL(r8) :: flnsc ! Clr. sky LW flux at srf (up-down)

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

! CRM input/output - originally provided by AGCM's cloud diagnostics package,
! but vanished before superparameterization CRM subroutine call. That is, 
! locally initialized ( all = 0. )
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
 REAL(r8) :: long, lati, topo, pslc, oces, iswf, swtc, tave
 INTEGER :: k, m!, i, j
 CHARACTER(len=64) :: FMT1, FMT2, FMT3
! Temporary variable for reading AGCM profiles from file
 REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tmp

! Screen output format specification
 FMT1 = "(A11,I4)"
 FMT2 = "(A11,F11.3)"
 FMT3 = "(2X,I3,999F11.3)"

 WRITE(*,*) "Semi-prognostic test of SP-CAM superparameterization"
 WRITE(*,*) 
 WRITE(*,*) "CRM grid description"
 WRITE(*,FMT1) "crm_nx  = ", crm_nx
 WRITE(*,FMT1) "crm_ny  = ", crm_ny
 WRITE(*,FMT1) "crm_nz  = ", crm_nz
 WRITE(*,FMT2) "crm_dx  = ", crm_dx
 WRITE(*,FMT2) "crm_dy  = ", crm_dy 
 WRITE(*,FMT2) "crm_dt  = ", crm_dt
 WRITE(*,*)    "YES3D?  = ", YES3D, "(1 - 3D CRM, 0 - 2D CRM)"

! Reading AGCM parameters, surface variables and profiles
 OPEN(11,FILE='SEMIPROG_IN',STATUS='old')
 READ(11,*) pplev, dt_gl, long, lati, topo, pslc, fsds, oces, iswf, fsut
 READ(11,*) flut, swtc, fsdsc, flnsc, ocnfrac, shf, lhf, taux, tauy

 WRITE(*,*)
 WRITE(*,*) "AGCM parameters (read from file)"
 WRITE(*,FMT1) "lchnk   = ", lchnk
 WRITE(*,FMT1) "icol    = ", icol
 WRITE(*,FMT1) "pplev   = ", pplev
 WRITE(*,FMT2) "dt_gl   = ", dt_gl
 WRITE(*,FMT2) "long    = ", long
 WRITE(*,FMT2) "lati    = ", lati

 WRITE(*,*)
 WRITE(*,*) "AGCM surface fields check"
 WRITE(*,FMT2) "topo    = ", topo
 WRITE(*,FMT2) "pslc    = ", pslc
 WRITE(*,FMT2) "fsds    = ", fsds
 WRITE(*,FMT2) "oces    = ", oces
 WRITE(*,FMT2) "iswf    = ", iswf
 WRITE(*,FMT2) "fsut    = ", fsut
 WRITE(*,FMT2) "flut    = ", flut
 WRITE(*,FMT2) "swtc    = ", swtc
 WRITE(*,FMT2) "fsdsc   = ", fsdsc
 WRITE(*,FMT2) "flnsc   = ", flnsc
 WRITE(*,FMT2) "ocnfrac = ", ocnfrac
 WRITE(*,FMT2) "shf     = ", shf
 WRITE(*,FMT2) "lhf     = ", lhf
 WRITE(*,FMT2) "taux    = ", taux
 WRITE(*,FMT2) "tauy    = ", tauy

 ALLOCATE(tmp(12,pplev))

 DO k=1,pplev
    READ(11,*) tmp(1,k), tmp(2,k), tmp(3,k), tmp(4,k), tmp(5,k), &
               tmp(6,k), tmp(7,k), tmp(8,k), tmp(9,k)
 END DO

 CLOSE(11)

! Deriving L-S variables from INPE AGCM file variables
 plev  = pplev-1
 ps    = pslc * 10
 phis  = topo * gravit
 fsns  = fsds - oces
 fsnt  = iswf - fsut
 fsntc = iswf - swtc
 tmp(10,1) = ps
 tmp(12,1) = 0
 DO k=2,pplev
    tmp(10,k) = ps * tmp(1,k)
    tmp(11,k) = tmp(10,k) - tmp(10,k-1)
    tave       = 0.5 * (tmp(2,k) + tmp(2,k-1))
    tmp(12,k) = tmp(12,k-1) + (rair*tave/gravit) * log(tmp(10,k-1)/tmp(10,k))
 END DO

 ALLOCATE(tl(plev),ql(plev),qcl(plev),qci(plev),ul(plev),vl(plev), &
    qrs(plev),qrl(plev),pmid(plev),pdel(plev),zmid(plev),zint(pplev))

 DO k=1,plev
    m = plev - k + 1
    tl(k)   = tmp(2,m+1)
    ql(k)   = tmp(3,m+1)
    qcl(k)  = tmp(4,m+1)
    qci(k)  = tmp(5,m+1)
    ul(k)   = tmp(6,m+1)
    vl(k)   = tmp(7,m+1)
    qrs(k)  = tmp(8,m+1)
    qrl(k)  = tmp(9,m+1)
    pmid(k) = tmp(10,m+1)
    pdel(k) = tmp(11,m+1) * (-1)
    zmid(k) = tmp(12,m+1)
 END DO

 WRITE(*,*)
 WRITE(*,*) "Checking large-scale profiles just after read from file"
 WRITE(*,*) "plev         tl         ql        qcl        qci         ul         vl       pmid       pdel       zmid"
 DO k=1,plev
    WRITE(*,FMT3) k, tl(k), ql(k), qcl(k), qci(k), ul(k), vl(k), pmid(k), pdel(k), zmid(k)
 END DO

!CAM3/SP-SAM interface height (zint)
 zint(pplev) = topo
 DO k=2,plev
    m = plev - k + 1
    zint(k) = 0.5 * (tmp(12,m+1) + tmp(12,m))
 END DO
 zint(1) = 2*tmp(12,pplev) - tmp(12,plev)

 WRITE(*,*)
 WRITE(*,*) "Interface height"
 DO k=1,pplev
    WRITE(*,'(I3,F11.3)') k, zint(k)
 END DO

 DEALLOCATE(tmp)

! Assuming null tendencies
 ALLOCATE(ultend(plev),vltend(plev),qltend(plev),qcltend(plev),qcitend(plev),sltend(plev))
 DO k=1,plev
    m = plev - k + 1
    ultend(m)  = 0.
    vltend(m)  = 0.
    qltend(m)  = 0.
    qcltend(m) = 0.
    qcitend(m) = 0.
    sltend(m)  = 0.
 END DO

 WRITE(*,*)
 WRITE(*,*) "Large-scale tendencies"
 WRITE(*,*) "plev     ultend     vltend     qltend    qcltend    qcitend     sltend"
 DO k=1,plev
    WRITE(*,FMT3) k, ultend(k), vltend(k), qltend(k), qcltend(k), qcitend(k), sltend(k)
 END DO

! Initializing other variables locally
 DO k = 1, crm_nz
    m = plev - k + 1
    u_crm(:,:,k)   = ul(m)
    v_crm(:,:,k)   = vl(m)
    w_crm(:,:,k)   = 0.
    t_crm(:,:,k)   = tl(m)
    q_crm(:,:,k)   = ql(m) + qcl(m) + qci(m)  ! total water
    qn_crm(:,:,k)  = qcl(k) + qci(m) ! cloud water/ice
    qp_crm(:,:,k)  = 0. ! precipitation
    qc_crm(:,:,k)  = 0.
    qi_crm(:,:,k)  = 0.
    qpc_crm(:,:,k) = 0.
    qpi_crm(:,:,k) = 0.

    qrs_crm(:,:,k) = qrs(m) / (cpair * pdel(m))
    qrl_crm(:,:,k) = qrl(m) / (cpair * pdel(m))
 END DO
 prec_crm(:,:) = 0.

! Initializing radiation variables
 fsds_crm(:,:)    = fsds
 fsns_crm(:,:)    = fsns
 fsntoa_crm(:,:)  = fsnt
 fsutoa_crm(:,:)  = fsut
 flwds_crm(:,:)   = flwds
 flns_crm(:,:)    = flns
 flut_crm(:,:)    = flut
 fsntoac_crm(:,:) = fsntc
 fsdsc_crm(:,:)   = fsdsc
 flutc_crm(:,:)   = flutc
 flnsc_crm(:,:)   = flnsc

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

 ALLOCATE(cld(plev),cldr(plev),cldtop(plev),gicewp(plev),gliqwp(plev),   &
    mc(plev),mcup(plev),mcdn(plev),mcuup(plev),mcudn(plev),crm_qc(plev), &
    crm_qi(plev),crm_qs(plev),crm_qg(plev),crm_qr(plev),tkez(plev),      &
    tkesgsz(plev),flux_u(plev),flux_v(plev),flux_qt(plev),               &
    fluxsgs_qt(plev),flux_qp(plev),pflx(plev),qt_ls(plev),qt_trans(plev),&
    qp_trans(plev),qp_fall(plev),qp_evp(plev),qp_src(plev),t_ls(plev))

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
 wnd      = sqrt(ul(1)**2 + vl(1)**2)
 tau00    = sqrt(taux**2 + tauy**2)
 bflxls   = shf/cpair + 0.61*tl(1)*lhf/latvap
 taux_crm = 0.
 tauy_crm = 0.
 z0m      = 0.

 call crm(lchnk, icol, & ! YES
          tl, ql, qcl, qci, ul, vl, & ! read from file (RFF)
          ps, pmid, pdel, phis, & ! RFF
          zmid, zint, dt_gl, plev, & ! RFF
          ultend, vltend, qltend, qcltend, qcitend, sltend, & ! RFF
          u_crm, v_crm, w_crm, & ! YES
          t_crm, q_crm, qn_crm, qp_crm, & ! YES
          qc_crm, qi_crm, qpc_crm, qpi_crm, & ! YES
          prec_crm, qrs_crm, qrl_crm, & ! YES
          fsds_crm, fsns_crm, fsntoa_crm, fsutoa_crm, & ! YES
          flwds_crm, flns_crm, flut_crm, & ! YES
          fsntoac_crm, fsdsc_crm, flutc_crm, flnsc_crm, & ! YES
          t_rad, qv_rad, qc_rad, qi_rad, & ! N/A
          precc, precl, precsc, precsl, & ! YES
          cltot, clhgh, clmed, cllow, & ! YES
          stat_buffer, cld, cldr, cldtop, & ! YES
          gicewp, gliqwp, & ! YES
          mc, mcup, mcdn, mcuup, mcudn, & ! YES
          crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, & ! YES
          tkez, tkesgsz, flux_u, flux_v, flux_qt, fluxsgs_qt, flux_qp, & ! YES
          pflx, qt_ls, qt_trans, qp_trans, qp_fall, & ! YES
          qp_evp, qp_src, t_ls, prectend, precstend, & ! YES
          ocnfrac, wnd, tau00, bflxls, & ! YES
          taux_crm, tauy_crm, z0m, timing_factor) ! YES

! 'tphysbc.F90' imediatelly follows computing cloud water and ice particle sizes
! using 'cldefr' routine which may be found in CAM3 code at
! cam1/models/atm/cam/src/physics/cam1/pkg_cldoptics.F90

! Output to screen some large-scale fields
 WRITE(*,*)
 WRITE(*,*) "Checking large-scale profiles after modifyed by superparameterization"
 WRITE(*,*) "plev         tl         ql        qcl        qci         ul         vl       pmid       pdel       zmid"
 DO k=1,plev
    WRITE(*,FMT3) k, tl(k), ql(k), qcl(k), qci(k), ul(k), vl(k), pmid(k), pdel(k), zmid(k)
 END DO

! Write output to file
 OPEN(13,FILE='SEMIPROG_OUT',STATUS='unknown',ACTION='write', &
    FORM='unformatted',ACCESS='direct',RECL=crm_nx*crm_ny*crm_nz*8)

 WRITE(13,rec=1) u_crm
 WRITE(13,rec=2) v_crm
 WRITE(13,rec=3) w_crm
 WRITE(13,rec=4) t_crm
 WRITE(13,rec=5) q_crm
 WRITE(13,rec=6) qn_crm
 WRITE(13,rec=7) qp_crm
 WRITE(13,rec=8) qc_crm
 WRITE(13,rec=9) qi_crm
 WRITE(13,rec=10) qpc_crm
 WRITE(13,rec=11) qpi_crm
 WRITE(13,rec=12) qrs_crm
 WRITE(13,rec=13) qrl_crm
 
 CLOSE(13)


 DEALLOCATE(tl,ql,qcl,qci,ul,vl,qrs,qrl,pmid,pdel,zmid,zint)

 DEALLOCATE(ultend,vltend,qltend,qcltend,qcitend,sltend)

 DEALLOCATE(cld,cldr,cldtop,gicewp,gliqwp,mc,mcup,mcdn,mcuup,mcudn,  &
    crm_qc,crm_qi,crm_qs,crm_qg,crm_qr,tkez,tkesgsz,flux_u,flux_v,   &
    flux_qt,fluxsgs_qt,flux_qp,pflx,qt_ls,qt_trans,qp_trans,qp_fall, &
    qp_evp,qp_src,t_ls)

 END PROGRAM crm_semiprog
