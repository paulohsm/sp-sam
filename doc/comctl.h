!----------------------------------------------------------------------- 
! 
! Purpose: Model control variables
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

      integer itsst             ! Sea surf. temp. update freq. (iters)
      integer nsrest            ! Restart flag
      integer iradsw            ! Iteration freq. for shortwave radiation
      integer iradlw            ! Iteration freq. for longwave radiation
      integer iradae            ! Iteration freq. for absorptivity/emissivity
      integer irad_always       ! Length of time SW/LW radiation will be run
                                ! continuously from the start of an initial run
      common /comctl/ itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae, irad_always

      integer nrefrq            ! Restart write freq.
      common /comctl/ nrefrq

      logical anncyc            ! true => do annual cycle (otherwise perpetual)
      logical nlend             ! true => end of run
      logical nlres             ! true => continuation run
      logical nlhst             ! true => regeneration run
      logical lbrnch            ! true => branch run
      common /comctl/ anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch

      logical aeres             ! true => read/write a/e data to/from restart file
      logical sstcyc            ! true => cycle sst dataset
      logical icecyc            ! true => cycle ice fraction dataset
      common /comctl/ aeres, sstcyc, icecyc

      logical ideal_phys        ! true => run "idealized" model configuration
      logical adiabatic         ! true => no physics
      logical flxave            ! true => send to coupler only on radiation time steps
      logical aqua_planet       ! Flag to run model in "aqua planet" mode
      common /comctl/ ideal_phys, adiabatic, flxave, aqua_planet

! f-v dynamics specific
! _ord = 1: first order upwind
! _ord = 2: 2nd order van Leer (Lin et al 1994)
! _ord = 3: standard PPM 
! _ord = 4: enhanced PPM (default)
      integer nsplit            ! Lagrangian time splits (Lin-Rood only)
      integer iord              ! scheme to be used in E-W direction
      integer jord              ! scheme to be used in N-S direction
      integer kord              ! scheme to be used for vertical mapping
      logical use_eta           ! Flag to use a's and b's set by dynamics/lr/set_eta.F90
      common /comctl/ nsplit, iord, jord, kord, use_eta

      logical doRamp_scon       ! true => turn on ramping for scon
      common /comctl/ doRamp_scon

      logical print_step_cost   ! true => print per-timestep cost info
      common /comctl/ print_step_cost

      logical doabsems          ! True => abs/emiss calculation this timestep
      logical dosw              ! True => shortwave calculation this timestep
      logical dolw              ! True => longwave calculation this timestep
      logical indirect          ! True => include indirect radiative effects of sulfate aerosols
      common /comctl/ doabsems, dosw, dolw, indirect

      integer som_conschk_frq   ! Energy conservation check frequency in SOM code
      integer ice_conschk_frq   ! Energy conservation check frequency in CSIM4 code
      common /comctl/ som_conschk_frq, ice_conschk_frq

      real(r8) divdampn         ! Number of days to invoke divergence damper
      real(r8) precc_thresh     ! Precipitation threshold for PRECCINT and PRECCFRQ
      real(r8) precl_thresh     ! Precipitation threshold for PRECLINT and PRECLFRQ
      common /comctl_r8/ divdampn, precc_thresh, precl_thresh
