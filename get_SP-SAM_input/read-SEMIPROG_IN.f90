PROGRAM GET_SEMIPROG_IN
IMPLICIT NONE

INTEGER :: k, l, nz, nt, plev, pplev
REAL :: iotmp, tstep, longi, lati, topo, lsmk
CHARACTER(12) :: tlev
CHARACTER(12), ALLOCATABLE, DIMENSION(:) :: timestamp
REAL, ALLOCATABLE, DIMENSION(:) :: pslc, usst, vsst, cssf, clsf, ocis, oces, iswf, &
                                   roce, olis, oles, role, swtc, ocic, lwtc, lwbc, &
                                   pres
REAL, ALLOCATABLE, DIMENSION(:,:) :: temp, umes, liqm, icem, uvel, vvel, swrh, lwrh
! Calculate saturation specific humidity
REAL, ALLOCATABLE, DIMENSION(:,:) :: e_s, q_s, ew, qsat, cond
REAL, PARAMETER :: e_s0 = 6.11, latvap = 2.5 * 10**6, rgas = 461., t0 = 273., eps = 0.622, latice=677.5
REAL :: lath

!check the number of vertical levels and time records before effectively reading the data
OPEN(31,FILE='SEMIPROG_IN',STATUS='OLD',FORM='formatted')
READ(31,*) nz, tstep, longi, lati, topo, lsmk
nt = 0
DO WHILE(.TRUE.)
   READ(31,ERR=44,END=44,FMT=*) tlev, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp
   READ(31,ERR=44,END=44,FMT=*) iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp
   DO k=1,nz
      READ(31,*) iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp, iotmp
   END DO
   nt = nt + 1
END DO

44 CONTINUE

!read the data
nt = nt + 1
REWIND(31)
READ(31,*)

ALLOCATE(timestamp(nt), pres(nz))
ALLOCATE(pslc(nt),usst(nt),vsst(nt),cssf(nt),clsf(nt),ocis(nt),oces(nt),iswf(nt),&
         roce(nt),olis(nt),oles(nt),role(nt),swtc(nt),ocic(nt),lwtc(nt),lwbc(nt))
ALLOCATE(temp(nz,nt),umes(nz,nt),liqm(nz,nt),icem(nz,nt),uvel(nz,nt),vvel(nz,nt),swrh(nz,nt),lwrh(nz,nt))

DO l=1,nt-1
   READ(31,*) timestamp(l), pslc(l), usst(l), vsst(l), cssf(l), clsf(l), ocis(l), oces(l), iswf(l)
   READ(31,*) roce(l), olis(l), oles(l), role(l), swtc(l), ocic(l), lwtc(l), lwbc(l)
   DO k=1,nz
      READ(31,*) pres(k), temp(k,l), umes(k,l), liqm(k,l), icem(k,l), uvel(k,l), vvel(k,l), swrh(k,l), lwrh(k,l)
   END DO
END DO


! Calculate saturation specific humidity (q_s) [kg/kg]
!e_s0 = 6.11 ! millibar
!latvap = 2.5 * 10^3 ! J g-1
!rgas = 461 ! J kg-1 K-1
!t0 = 273 ! K
!eps = 0.622
!ew = 6.1121*(1.0007+3.46e-6*Pa).*exp((17.502*Ta)./(240.97+Ta))! % in mb
!qsat  = 0.62197*(ew./(Pa-0.378*ew)) !;                         % mb -> kg/kg
ALLOCATE(e_s(nz,nt), q_s(nz,nt), ew(nz,nt), qsat(nz,nt), cond(nz,nt))
DO l=1,nt-1
   WRITE(*,*) timestamp(l), "       TEMP,            UMES,            e_s,             q_s,             cond"
   DO k=1,nz
      lath = latvap
      IF (temp(k,l).lt.0.) lath = latice
      e_s(k,l) = e_s0 * exp( (lath/rgas) * ( 1/t0 - 1/temp(k,l) ) )
!     q_s(k,l) = eps*e_s(k,l) / ( pres(k) - (1-eps)*e_s(k,l) )
      q_s(k,l) = eps * e_s(k,l) / pres(k)
      cond(k,l) = q_s(k,l)-umes(k,l)
!     cond(k,l) = umes(k,l)-q_s(k,l)
      IF (cond(k,l).le.0.) cond(k,l) = 0.
      WRITE(*,*) pres(k), temp(k,l), umes(k,l), e_s(k,l), q_s(k,l), cond(k,l)*1000!, umes(k,l)-q_s(k,l)
   END DO
END DO


END PROGRAM GET_SEMIPROG_IN
