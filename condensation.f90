    SUBROUTINE CONDENSATION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,            &
                           PPABS, PZZ, PT, PRV, PRC, PRI, PMFLX, PCLDFR, LUSERI )
!   ############################################################################
!
!!
!!    PURPOSE
!!    -------
!!**  Routine to diagnose cloud fraction and liquid and ice condensate mixing ratios
!!     
!!    
!!**  METHOD
!!    ------
!!    Based on the large-scale fields of temperature, water vapor, and possibly
!!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed 
!!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
!!
!!    The total variance is parameterized as the sum of  stratiform/turbulent variance 
!!    and a convective variance.
!!    The turbulent variance is parameterized as a function of first-order moments, and
!!    the convective variance is modelled as a function of the convective mass flux (units kg/s m^2)
!!    as provided by a mass flux convection scheme.
!!
!!    Nota: if the host model does not use prognostic values for liquid and solid condensate
!!    or does not provide a convective mass flux, put all these values to zero.
!!    Also, it is supposed that vertical model levels are numbered from
!!    1 to KLEV, where 1 is the first model level above the surface
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!      INI_CST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST       : contains physical constants
!!
!!    REFERENCE
!!    ---------
!!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/06/2001
!!      modified    20/03/2002 : add convective Sigma_s and improve turbulent
!!                               length-scale in boundary-layer and near tropopause
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN)   :: KLON    ! horizontal dimension
INTEGER,                    INTENT(IN)   :: KLEV    ! vertical dimension
INTEGER,                    INTENT(IN)   :: KIDIA   ! value of the first point in x
                                                    ! default=1
INTEGER,                    INTENT(IN)   :: KFDIA   ! value of the last point in x
                                                    ! default=KLON
INTEGER,                    INTENT(IN)   :: KBDIA   ! vertical  computations start at
!                                                   ! KBDIA that is at least 1
INTEGER,                    INTENT(IN)   :: KTDIA   ! vertical computations can be
                                                    ! limited to KLEV + 1 - KTDIA
                                                    ! default=1
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PT     ! grid scale T  (K)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
LOGICAL                                   :: LUSERI ! logical switch to compute both
						    ! liquid and solid condensate (LUSERI=.TRUE.)
						    ! or only liquid condensate (LUSERI=.FALSE.)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PMFLX  ! convective mass flux (kg/(s m^2))
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   :: PCLDFR ! fractional cloudiness (between 0 and 1)
!
!						 
!*       0.2   Declarations of local variables :
!
INTEGER  :: JI, JK, JKT, JKP, JKM     ! loop index
REAL, DIMENSION(KLON,KLEV) :: ZTLK, ZRT       ! work arrays for T_l, r_t 
REAL, DIMENSION(KLON,KLEV) :: ZL              ! length-scale
INTEGER, DIMENSION(KLON)   :: ITPL    ! top levels of tropopause/highest inversion
REAL, DIMENSION(KLON)      :: ZTMIN   ! min Temp. related to ITPL
!
REAL :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
REAL :: ZLL, DZZ, ZZZ ! length scales
REAL :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
REAL :: ZSIG_CONV                                   ! convective part of Sig_s
!
!*       0.3  Definition of constants :  
!
!-------------------------------------------------------------------------------
!
REAL :: ZL0     = 600.        ! tropospheric length scale
                              ! changed to 600 m instead of 900 m to give a consistent
                              ! value (linear increase) in general 500 m deep oceanic
                              ! mixed layer - but could be put back to 900 m if wished
REAL :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
REAL :: ZCSIG_CONV = 0.30E-2  ! scaling factor for ZSIG_CONV as function of mass flux 
!
!-------------------------------------------------------------------------------
!
REAL :: XP00   ! reference pressure
REAL :: XPI    ! Pi
REAL ::  XG    ! gravity constant
REAL :: XMD    ! molecular weight of dry air
REAL :: XMV    ! molecular weight of water vapor
REAL :: XRD    ! gaz constant for dry air
REAL :: XRV    ! gaz constant for water vapor
REAL :: XCPD   ! specific heat of dry air
REAL :: XCPV   ! specific heat of water vapor
REAL :: XRHOLW ! density of liquid water
REAL :: XCL    ! specific heat of liquid water
REAL :: XCI    ! specific heat of ice
REAL :: XTT    ! triple point temperature
REAL :: XLVTT  ! latent heat of vaporisation at XTT
REAL :: XLSTT  ! latent heat of sublimation at XTT 
REAL :: XLMTT  ! latent heat of melting at XTT
REAL :: XESTT  ! saturation pressure at XTT
REAL :: XALPW  ! constants in saturation pressure over liquid water
REAL :: XBETAW 
REAL :: XGAMW 
REAL :: XALPI  ! constants in saturation pressure over ice
REAL :: XBETAI 
REAL :: XGAMI 

! Initialize thermodynamic constants in module MODD_CST

XP00   = 1.E5        ! reference pressure
XPI    = 3.141592654 ! Pi
 XG    = 9.80665     ! gravity constant
XMD    = 28.9644E-3  ! molecular weight of dry air
XMV    = 18.0153E-3  ! molecular weight of water vapor
XRD    = 287.05967   ! gaz constant for dry air
XRV    = 461.524993  ! gaz constant for water vapor
XCPD   = 1004.708845 ! specific heat of dry air
XCPV   = 1846.1      ! specific heat of water vapor
XRHOLW = 1000.       ! density of liquid water
XCL    = 4218.       ! specific heat of liquid water
XCI    = 2106.       ! specific heat of ice
XTT    = 273.16      ! triple point temperature
XLVTT  = 2.5008E6    ! latent heat of vaporisation at XTT
XLSTT  = 2.8345E6    ! latent heat of sublimation at XTT 
XLMTT  = 0.3337E6    ! latent heat of melting at XTT
XESTT  = 611.14      ! saturation pressure at XTT
XALPW  = 60.22416    ! constants in saturation pressure over liquid water
XBETAW = 6822.459384
XGAMW  = 5.13948
XALPI  = 32.62116    ! constants in saturation pressure over ice
XBETAI = 6295.421
XGAMI  = 0.56313


PCLDFR(:,:) = 0. ! Initialize values

JKT = KLEV+1-KTDIA
DO JK=KBDIA,JKT
DO JI=KIDIA,KFDIA
   ZTEMP  = PT(JI,JK)
    !latent heat of vaporisation/sublimation
   ZLV    = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
   ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

    !store temperature at saturation and total water mixing ratio
   ZRT(JI,JK)   = PRV(JI,JK) + PRC(JI,JK) + PRI(JI,JK)
   ZCPD         = XCPD + ZRT(JI,JK) * XCPV
   ZTLK(JI,JK)  = ZTEMP - ZLV*PRC(JI,JK)/ZCPD - ZLS*PRI(JI,JK)/ZCPD
END DO
END DO

!-------------------------------------------------------------------------------
! Determine tropopause/inversion  height from minimum temperature

ITPL(:)  = KBDIA+1
ZTMIN(:) = 400.
DO JK = KBDIA+1,JKT-1
   DO JI=KIDIA,KFDIA 
         IF ( PT(JI,JK) < ZTMIN(JI) ) THEN
              ZTMIN(JI) = PT(JI,JK)
              ITPL(JI) = JK
         END IF
   END DO
END DO

! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

ZL(:,KBDIA) = 20.
DO JK = KBDIA+1,JKT
DO JI=KIDIA,KFDIA 
      ! free troposphere
   ZL(JI,JK) = ZL0
   JKP = ITPL(JI)
   ZZZ =  PZZ(JI,JK) -  PZZ(JI,KBDIA)
      ! approximate length for boundary-layer : linear increase
   IF ( ZL0 > ZZZ )  ZL(JI,JK) = ZZZ
      ! gradual decrease of length-scale near and above tropopause/top inversion
   IF ( ZZZ > 0.9*(PZZ(JI,JKP)-PZZ(JI,KBDIA)) ) &
        ZL(JI,JK) = .6 * ZL(JI,JK-1) 
END DO
END DO
!-------------------------------------------------------------------------------


DO JK=KBDIA+1,JKT-1
   JKP=JK+1
   JKM=JK-1
DO JI=KIDIA,KFDIA
   ZTEMP  = PT(JI,JK)
    !latent heat of vaporisation/sublimation
   ZLV    = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
   ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

   ZCPD   = XCPD + ZRT(JI,JK) * XCPV
    !temperature at saturation
   ZTL    = ZTEMP - ZLV*PRC(JI,JK)/ZCPD - ZLS*PRI(JI,JK)/ZCPD
    !saturated water vapor mixing ratio over liquid water
   ZPV    = EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) )
   ZQSL   = XRD / XRV * ZPV / ( PPABS(JI,JK) - ZPV )

    !saturated water vapor mixing ratio over ice
   ZPIV   = EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) )
   ZQSI   = XRD / XRV * ZPIV / ( PPABS(JI,JK) - ZPIV )

    !interpolate between liquid and solid as function of temperature
    ! glaciation interval is specified here to 20 K
   ZFRAC = ( ZTL  - 250.16 ) / ( XTT - 250.16 )  ! liquid/solid fraction
   ZFRAC = MAX( 0., MIN(1., ZFRAC ) )
   ZFRAC = ZFRAC * ZFRAC
   IF(.NOT. LUSERI) ZFRAC=1.
   ZQSL = ( 1. - ZFRAC ) * ZQSI + ZFRAC * ZQSL
   ZLV  = ( 1. - ZFRAC ) * ZLS  + ZFRAC * ZLV

    !coefficients a and b
   ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )
   ZA   = 1. / ( 1. + ZLV/ZCPD * ZAH )
   ZB   = ZAH * ZA

    !parameterize Sigma_s with first_order closure

   DZZ    =  PZZ(JI,JKP)  - PZZ(JI,JKM)
   ZDRW   =  ZRT(JI,JKP)  - ZRT(JI,JKM)
   ZDTL   =  ZTLK(JI,JKP) - ZTLK(JI,JKM) + XG/ZCPD * DZZ
   ZLL    =  ZL(JI,JK)

   ZSIG_CONV = ZCSIG_CONV * PMFLX(JI,JK) / ZA ! standard deviation due to convection
   ZSIGMA =  SQRT( MAX( 1.E-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
             ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL + ZB*ZB*ZDTL*ZDTL  ) &
                                       + ZSIG_CONV * ZSIG_CONV ) )
    !zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
   ZSIGMA = MAX( ZSIGMA, 1.E-12 )
   
    !normalized saturation deficit
   ZSBAR = ZA * ( ZRT(JI,JK) - ZQSL )
   ZQ1   = ZSBAR / ZSIGMA 

    !cloud fraction
   PCLDFR(JI,JK) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

    !total condensate
   IF (ZQ1 > 0. .AND. ZQ1 <= 2. ) THEN
      ZCOND = EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
   ELSE IF (ZQ1 > 2.) THEN
      ZCOND = ZQ1
   ELSE
      ZCOND = EXP( 1.2*ZQ1-1 )
   END IF
   ZCOND = ZCOND * ZSIGMA

   if ( zcond<1.e-6) then
       zcond = 0.
       pcldfr(ji,jk) = 0.
   end if

   PRC(JI,JK) = ZFRAC * ZCOND ! liquid condensate
   IF (LUSERI) THEN
      PRI(JI,JK) = (1.-ZFRAC) * ZCOND   ! solid condensate
   END IF


    ! compute s'rl'/Sigs^2
    ! used in w'rl'= w's' * s'rl'/Sigs^2
!  PSIGRC(JI,JK) = PCLDFR(JI,JK)   ! Gaussian relation 

END DO
END DO
!
END SUBROUTINE CONDENSATION
