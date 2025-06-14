!WRF:MODEL_LAYER:PHYSICS
!
! Current version is tagged as V3.5
! Reference is Morrison et al. (2005), JAS, Morrison et al. (2009), MWR

! changes with respect to V1.4

! V1.5
! 1) more pathways to allow hail to form (only affects IHAIL=1 option), from collisions of snow/cloud water
! 2) bug fix to PGAM calculation (multiplication instead of division by air density)

! V1.6
! 1) added parameter TMELT for all calculations involving melting point
! 2) replaced hard-wired gas constant for air with parameter value 'R'

! V1.7
! 1) modification to minimum mixing ratio in dry conditions, change from 10^-6 to 10^-8 kg/kg
!   to improve reflectivity at low mixing ratio amounts
! 2) bug fix to prevent possible division by zero error involving LAMI
! 3) change for liquid saturation vapor pressure, replace old formula with Flatau et al. 1992

! V2
! 1) bug fix to maximum-allowed particle fallspeeds (air density correction factor considered)
! 2) change to comments

! V2.1   
! 1) addition of rain drop breakup following Verlinde and Cotton (1993)
! 2) change to minimum allowed lambda (slope parameter) for rain
! 3) addition of accelerated melting of graupel/hail/snow due to collision with rain

! V3
! minor revisions by Andy Ackerman
! 1) replaced kinematic with dynamic viscosity 
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! Additional fixes
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NO SUBTRACTED SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SIN
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS 
! 7) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 8) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 9) BUG FIX TO IGRAUP SWITCH FOR NO GRAUPEL/HAIL

! hm bug fix 3/16/12

! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted

! hm, changes 3/4/13 for V3.3

! 1) removed second initialization of evpms (non-answer-changing)
! 2) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass 
!    collected by rain
! 3) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!    cm-3. This was done to address the problem of excessive and persistent
!    anvil cirrus produced by the scheme, and was found to greatly improve forecasts over
!    at convection-permitting scales over the central U.S. in summertime.

! hm, changes 7/25/13 for V3.4

! 1) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!    Hailong Wang (PNNL)
! 2) very minor fix to immersion freezing rate formulation (negligible impact)
! 3) clarifications to code comments
! 4) minor change to shedding of rain, remove limit so that the number of
!    collected drops can smaller than number of shed drops
! 5) change of specific heat of liquid water from 4218 to 4187 J/kg/K

! hm, changes 1/20/15 for version 3.5

! 1) minor bug fix to diagnostic steady-state supersaturation equation for droplet activation in
!    cloud interior for option IBASE=1, minus sign is missing from the pressure term in this equation,
!    note this only has ~10% effect relative to the DQSDT term so the impact of the fix is small.
!    Bug reported by Xiaowen Li (NASA Goddard).
! 2) minor bug fix to melting of snow and graupel, an extra factor of air density (RHO) was removed
!    from the calculation of PSMLT and PGMLT
! 3) redundant initialization of PSMLT (non answer-changing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE module_mp_GRAUPEL
!bloss   USE     module_wrf_error
!bloss      USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
!bloss      USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT

!  USE module_state_description
  ! parameters from SAM and options from wrapper routine.
   use params, only: lcond, lsub, cp, rgas, rv
   use grid, only: dx
   use micro_params
   use hoppel_transfer, only: mass_fraction, hoppel_aitken_accum_transfer
   
   IMPLICIT NONE

   REAL, PARAMETER :: PI = 3.1415926535897932384626434
   REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297

   PUBLIC  ::  MP_GRAUPEL
   PUBLIC  ::  POLYSVP

   PRIVATE :: GAMMA, DERF1
   PRIVATE :: PI, SQRTPI
   PUBLIC :: M2005MICRO_GRAUPEL !bloss

! SWITCHES FOR MICROPHYSICS SCHEME
! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

     INTEGER, PRIVATE ::  IACT

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

     INTEGER, PRIVATE ::  INUM

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
     REAL, PRIVATE ::      NDCNST

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

     INTEGER, PRIVATE ::  ILIQ

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS

     INTEGER, PRIVATE ::  INUC

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

     INTEGER, PRIVATE ::  IBASE

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

     INTEGER, PRIVATE ::  ISUB      

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

     INTEGER, PRIVATE ::  IGRAUP

! HM ADDED NEW OPTION FOR HAIL V1.3
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

     INTEGER, PRIVATE ::  IHAIL

! HM ADDED 8/1/08, v1.4
! SWITCH FOR WARM RAIN SCHEME
! IRAIN = 0, WARM RAIN (AUTO, ACC, SELF-COLL) FROM KHAIROUTIDNOV AND KOGAN (2000)
! IRAIN = 1, WARM RAIN (AUTO, ACC, SELF-COLL) FROM SEIFERT AND BEHENG (2001)

     INTEGER, PRIVATE ::  IRAIN      

! PB ADDED 4/13/09
! SWITCH TO TURN ON/OFF CLOUD LIQUID WATER SATURATION ADJUSTMENT
! WHEN USING TOTAL WATER FORMULATION IN SAM, THE SATURATION 
! ADJUSTMENT IS PERFORMED BEFORE CALLING M2005MICRO_GRAUPEL.
! THIS OPTION ALLOWS US TO AVOID PERFORMING IT IN M2005MICRO_GRAUPEL
! UNDER THE THEORY THAT THE OTHER MICROPHYSICAL PROCESSES WILL NOT
! DRIVE IT FAR FROM SATURATION.
! ISATADJ = 0, SATURATION ADJUSTMENT PEROFORMED IN M2005MICRO_GRAUPEL
! ISATADJ = 1, SATURATION ADJUSTMENT _NOT_ PEROFORMED IN M2005MICRO_GRAUPEL

     INTEGER, PRIVATE :: ISATADJ

! BRNR ADDED 8/31/11
! SWITCH TO TURN ON/OFF EVAPORATION TENDENCY ON NC. 
! IEVPNC = 0, EVAPORATION HAS NO TENDENCY ON NC
! IEVPNC = 1, NC IS REMOVED IN 1:1 RATIO WITH QN

     INTEGER, PRIVATE :: IEVPNC

! BRNR ADDED 9/8/11
! SWITCH TO TURN ON 3D AEROSOL TREATMENT
! IPRGAER = 0, USE FIXED AEROSOL NUMBER AND MEAN RADIUS
! IPRGAER = 1, USE PROGNOSTIC AEROSOL. PRESUMES IACT = 2 (dopredictNc = .true.)

     INTEGER, PRIVATE :: IPRGAER

! BRNR ADDED 9/16/11
! SWITCH TO TURN OFF PRECIP BY SHUTTING DOWN AUTOCONVERSION
! IPRECOFF = 0, AUTOCONVERSION ALLOWED
! IPRECOFF = 1, AUTOCONVERSION DISABLED

     INTEGER, PRIVATE :: IPRECOFF

! BRNR ADDED 10/31/11
! SWITCH TO TURN OFF SEDIMENTATION PROCESSES FOR CLOUD WATER
! ISEDOFF = 0, CLOUD WATER SEDIMENTATION IS ALLOWED
! ISEDOFF = 1, CLOUD WATER SEDIMENTATION IS DISABLED

     INTEGER, PRIVATE :: ISEDOFF

! CLOUD MICROPHYSICS CONSTANTS

     REAL, PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL, PRIVATE ::      R           ! GAS CONSTANT FOR AIR
!bloss     REAL, PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
!bloss     REAL, PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
     REAL, PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
     REAL, PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
     REAL, PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
     REAL, PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
     REAL, PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
     REAL, PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
     REAL, PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
     REAL, PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
     REAL, PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
     REAL, PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
     REAL, PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
     REAL, PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
     REAL, PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
     REAL, PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
     REAL, PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
     REAL, PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
     REAL, PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
! V1.6
     REAL, PRIVATE ::      TMELT     ! melting temp (K)
! hm, add for V2.1
     REAL, PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

! CCN SPECTRA FOR IACT = 1

     REAL, PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
     REAL, PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

! AEROSOL PARAMETERS FOR IACT = 2

     REAL, PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
     REAL, PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
     REAL, PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
     REAL, PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
     REAL, PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
     REAL, PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
     REAL, PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
     REAL, PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
     REAL, PRIVATE ::      BACT    ! ACTIVATION PARAMETER
     REAL, PRIVATE ::      BACT_coarse
     REAL, PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
     REAL, PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
     REAL, PRIVATE ::      RM3     ! ", fake coarse mode
     REAL, PRIVATE ::      Ncoarse ! Number concentration /mg coarse
     REAL, PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
     REAL, PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
     REAL, PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
     REAL, PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
     REAL, PRIVATE ::      SIG3    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 3
     REAL, PRIVATE ::      SG1
     REAL, PRIVATE ::      SG2
     REAL, PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F13
     REAL, PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL, PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL, PRIVATE ::      F23     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2     
     REAL, PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
     REAL, PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING


     REAL, PRIVATE ::      NA1     !Total number mode 1
     REAL, PRIVATE ::      QA1     !Total mass mode 1

! ASSUMED MUR WHEN SWITCH IMUR IS SET TO CONSTANT MUR
     REAL, PRIVATE ::      MURCNST

! AEROSOL PARAMETERS FOR IPRGAER = 1 !BRNR 9/8/11

     REAL, PRIVATE ::      MAER1   ! TOTAL AEROSOL MASS MIXING RATIO, MODE 1 (KG/KG) 
     REAL, PRIVATE ::      MAER2   ! TOTAL AEROSOL MASS MIXING RATIO, MODE 2 (KG/KG)

! CONSTANTS TO IMPROVE EFFICIENCY

     REAL, PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
     REAL, PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
     REAL, PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
     REAL, PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
     REAL, PRIVATE :: CONS41

! v1.4
     REAL, PRIVATE :: dnu(16)

!..Various radar related variables, from GT

!..Lookup table dimensions
      INTEGER, PARAMETER, PRIVATE:: nbins = 100
      INTEGER, PARAMETER, PRIVATE:: nbr = nbins
      INTEGER, PARAMETER, PRIVATE:: nbs = nbins
      INTEGER, PARAMETER, PRIVATE:: nbg = nbins
      REAL(8), DIMENSION(nbins+1):: ddx
      REAL(8), DIMENSION(nbr):: Dr, dtr
      REAL(8), DIMENSION(nbs):: Dds, dts
      REAL(8), DIMENSION(nbg):: Ddg, dtg
      REAL(8), PARAMETER, PRIVATE:: lamda_radar = 0.10         ! in meters
      REAL(8), PRIVATE:: K_w, PI5, lamda4
      COMPLEX*16, PRIVATE:: m_w_0, m_i_0
      REAL(8), DIMENSION(nbins+1), PRIVATE:: simpson
      REAL(8), DIMENSION(3), PARAMETER, PRIVATE:: basis =      &
                           (/1.d0/3.d0, 4.d0/3.d0, 1.d0/3.d0/)

      INTEGER, PARAMETER, PRIVATE:: slen = 20
      CHARACTER(len=slen), PRIVATE::                                    &
              mixingrulestring_s, matrixstring_s, inclusionstring_s,    &
              hoststring_s, hostmatrixstring_s, hostinclusionstring_s,  &
              mixingrulestring_g, matrixstring_g, inclusionstring_g,    &
              hoststring_g, hostmatrixstring_g, hostinclusionstring_g

      REAL, PARAMETER, PRIVATE:: D0r = 50.E-6
      REAL, PARAMETER, PRIVATE:: D0s = 100.E-6
      REAL, PARAMETER, PRIVATE:: D0g = 100.E-6
      CHARACTER*256:: mp_debug

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRAUPEL_INIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      integer n,i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
! SET PRIOR TO CODE COMPILATION

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

      INUM = 1 !bloss: use flag in prm file
      if(dopredictNc) then
         INUM = 0
      end if

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)

      NDCNST = Nc0 !bloss: use value from prm file (default=100.)

! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      if(dospecifyaerosol) then !bloss: specify using flag from prm file
         IACT = 2
      else
         IACT = 1
      end if

! IPRGAER = 0, USE FIXED AEROSOL MOMENTS FOR RAZZAK GHAN ACTIVATION
! IPRGAER = 1, USE PROGNOSTIC AEROSOL

      IPRGAER = 0 !brnr: turn on prognostic aerosol from prm file
      if(doprogaerosol) then 
         IPRGAER = 1
      end if

! IPRECOFF = 0, ALLOW AUTOCONVERSION
! IPRECOFF = 1, DISABLE AUTOCONVERSION

      IPRECOFF = 0 !brnr: turn off autoconversion from prm file
      if(doprecoff) then 
         IPRECOFF = 1
      end if 

! ISEDOFF = 0, ALLOW CLOUD DROPLET SEDIMENTATION
! ISEDOFF = 1, DISABLE CLOUD DROPLET SEDIMENTATION

      ISEDOFF = 0 !brnr: turn off sedimentation from prm file
      if(dosedoff) then
         ISEDOFF = 1
      end if


! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      if(docloudedgeactivation) then
         IBASE = 2
      else
         IBASE = 1
      end if

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      if(dosubgridw) then
         ISUB = 0
      else
         ISUB = 1      
      end if

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      if(doicemicro) then !bloss: specify using flag from prm file
         ILIQ = 0
      else
         ILIQ = 1
      end if

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

      if(doarcticicenucl) then !bloss: specify using flag from prm file
         INUC = 1
      else
         INUC = 0
      end if

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

      if(dograupel) then
         IGRAUP = 0
      else
         IGRAUP = 1
      end if

! HM ADDED 11/7/07, V1.3
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL

      if(dohail) then
         IHAIL = 1
      else
         IHAIL = 0
      end if
 
! HM ADDED 8/1/08, v1.4
! SWITCH FOR WARM RAIN SCHEME
! IRAIN = 0, WARM RAIN (AUTO, ACC, SELF-COLL) FROM KHAIROUTIDNOV AND KOGAN (2000)
! IRAIN = 1, WARM RAIN (AUTO, ACC, SELF-COLL) FROM SEIFERT AND BEHENG (2001)

      if(dosb_warm_rain) then
        IRAIN = 1
      else
        IRAIN = 0
      end if


! PB ADDED 4/13/09.  TURN OFF SATURATION ADJUSTMENT WITHIN M2005MICRO_GRAUPEL
! IN TOTAL WATER VERSION.  IT NOW TAKES PLACE BEFORE M2005MICRO_GRAUPEL IS CALLED.

      if(dototalwater) then
        ISATADJ = 1 ! total water version -- saturation adjustment outside this routine
      else
        ISATADJ = 0 ! separate vapor and cloud liquid version -- saturation adjustment here
      end if

! BRNR ADDED 8/31/11. OPTION TO TURN ON EVAPORATION TENDENCY FOR NC


      if(doevapnc) then
         IEVPNC = 1
      else
         IEVPNC = 0
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET PHYSICAL CONSTANTS

! FALLSPEED PARAMETERS (V=AD^B)
         AI = 700.
         AC = 3.E7
         AS = 11.72
         AR = 841.99667
         BI = 1.
         BC = 2.
         BS = 0.41
         BR = 0.8
! V1.3
         IF (IHAIL.EQ.0) THEN
	 AG = 19.3
	 BG = 0.37
         ELSE ! (MATSUN AND HUGGINS 1980)
         AG = 114.5 
         BG = 0.5
         END IF

! CONSTANTS AND PARAMETERS
         !bloss: use values from params module
         R = rgas
!bloss         R = 287.15
!bloss         RV = 465.5
!bloss         CP = 1005.
! V1.6
         TMELT = 273.15
! V1.6
         RHOSU = 85000./(R*TMELT)
         RHOW = rho_water !bloss 997.
         RHOI = rho_cloud_ice !bloss 500.
         RHOSN = rho_snow !bloss 100.
! V1.3
         IF (IHAIL.EQ.0) THEN
	 RHOG = 400.
         ELSE
         RHOG = 900.
         END IF
         AIMM = 0.66
         BIMM = 100.
         ECR = 1.
         DCS = 125.E-6
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3
	 MG0 = 1.6E-10
         F1S = 0.86
         F2S = 0.28
         F1R = 0.78
! V3 5/27/11
!         F2R = 0.32
! AA revision 4/1/11
         F2R = 0.308
         G = 9.806
         QSMALL = 1.E-14
         EII = 0.1
         ECI = 0.7
! HM, ADD FOR V3.2
! hm, 7/23/13
!         CPW = 4218.
         CPW = 4187.

! SIZE DISTRIBUTION PARAMETERS

         CI = RHOI*PI/6.
         DI = 3.
         CS = RHOSN*PI/6.
         DS = 3.
         CG = RHOG*PI/6.
         DG = 3.

! RADIUS OF CONTACT NUCLEI
         RIN = 0.1E-6

         MMULT = 4./3.*PI*RHOI*(5.E-6)**3

! SIZE LIMITS FOR LAMBDA

         LAMMAXI = 1./1.E-6
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./20.E-6
!         LAMMINR = 1./500.E-6
! V2.1
         LAMMINR = 1./2800.E-6

         LAMMAXS = 1./10.E-6
         LAMMINS = 1./2000.E-6
         LAMMAXG = 1./20.E-6
         LAMMING = 1./2000.E-6

! CCN SPECTRA FOR IACT = 1

! MARITIME
! MODIFIED FROM RASMUSSEN ET AL. 2002
! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

              K1 = ccnexpnt !bloss: specify using values from prm file
              C1 = ccnconst !bloss

!bloss              K1 = 0.4
!bloss              C1 = 120. 

! CONTINENTAL

!              K1 = 0.5
!              C1 = 1000. 

! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

         MW = 0.018
         OSM = 1.
         VI = 3.
         EPSM = 0.7
         RHOA = 1777.
         MAP = 0.132
         MA = 0.0284
         RR = 8.3187
         BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)
         BACT_coarse = 1.1

! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE
! MODE 1

         RM1 = rm_accum !bloss: specify using values from prm file
         SIG1 = sigma_accum
         SG1 = LOG(sigma_accum)
         NANEW1 = N_accum
!bloss         RM1 = 0.052E-6
!bloss         SIG1 = 2.04
!bloss         NANEW1 = 100.0E6
         F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
         F21 = 1.+0.25*LOG(SIG1)

! MODE 2

         RM2 = rm_aitken !bloss: specify using values from prm file
         SIG2 = sigma_aitken
         SG2 = LOG(sigma_aitken)
         NANEW2 = N_aitken
!bloss         RM2 = 1.3E-6
!bloss         SIG2 = 2.5
!bloss         NANEW2 = 1.E6
         F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
         F22 = 1.+0.25*LOG(SIG2)

         ! Use fixed coarse values
         RM3 = 1.3e-6
         SIG3 = 2.5
         F13 = 0.5*EXP(2.5*(LOG(SIG3))**2)
         F23 = 1.+0.25*LOG(SIG3)
         Ncoarse = 1.5 * 10./ (PI * (RM3*2*1.E6)**2  * EXP(2*(LOG(SIG3))**2))  ! 10 m/s  units /mg
! CONSTANTS FOR EFFICIENCY

         CONS1=GAMMA(1.+DS)*CS
         CONS2=GAMMA(1.+DG)*CG
         CONS3=GAMMA(4.+BS)/6.
         CONS4=GAMMA(4.+BR)/6.
         CONS5=GAMMA(1.+BS)
         CONS6=GAMMA(1.+BR)
         CONS7=GAMMA(4.+BG)/6.
         CONS8=GAMMA(1.+BG)
         CONS9=GAMMA(5./2.+BR/2.)
         CONS10=GAMMA(5./2.+BS/2.)
         CONS11=GAMMA(5./2.+BG/2.)
         CONS12=GAMMA(1.+DI)*CI
         CONS13=GAMMA(BS+3.)*PI/4.*ECI
         CONS14=GAMMA(BG+3.)*PI/4.*ECI
         CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
         CONS16=GAMMA(BI+3.)*PI/4.*ECI
         CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
         CONS18=RHOSN*RHOSN
         CONS19=RHOW*RHOW
         CONS20=20.*PI*PI*RHOW*BIMM
         CONS21=4./(DCS*RHOI)
         CONS22=PI*RHOI*DCS**3/6.
         CONS23=PI/4.*EII*GAMMA(BS+3.)
         CONS24=PI/4.*ECR*GAMMA(BR+3.)
         CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
         CONS26=PI/6.*RHOW
         CONS27=GAMMA(1.+BI)
         CONS28=GAMMA(4.+BI)/6.
         CONS29=4./3.*PI*RHOW*(25.E-6)**3
         CONS30=4./3.*PI*RHOW
         CONS31=PI*PI*ECR*RHOSN
         CONS32=PI/2.*ECR
         CONS33=PI*PI*ECR*RHOG
         CONS34=5./2.+BR/2.
         CONS35=5./2.+BS/2.
         CONS36=5./2.+BG/2.
         CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
         CONS38=PI*PI/3.*RHOW
         CONS39=PI*PI/36.*RHOW*BIMM
         CONS40=PI/6.*BIMM
         CONS41=PI*PI*ECR*RHOW

! v1.4
         dnu(1) = -0.557
         dnu(2) = -0.557
         dnu(3) = -0.430
         dnu(4) = -0.307
         dnu(5) = -0.186
         dnu(6) = -0.067
         dnu(7) = 0.050
         dnu(8) = 0.167
         dnu(9) = 0.282
         dnu(10) = 0.397
         dnu(11) = 0.512
         dnu(12) = 0.626
         dnu(13) = 0.739
         dnu(14) = 0.853
         dnu(15) = 0.966
         dnu(16) = 0.966

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables for radar reflecitivity calculations
!..Create bins of rain (from min diameter up to 5 mm).
      ddx(1) = D0r*1.0d0
      ddx(nbr+1) = 0.005d0
      do n = 2, nbr
         ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbr) &
                  *DLOG(ddx(nbr+1)/ddx(1)) +DLOG(ddx(1)))
      enddo
      do n = 1, nbr
         Dr(n) = DSQRT(ddx(n)*ddx(n+1))
         dtr(n) = ddx(n+1) - ddx(n)
      enddo

!..Create bins of snow (from min diameter up to 2 cm).
      Ddx(1) = D0s*1.0d0
      Ddx(nbs+1) = 0.02d0
      do n = 2, nbs
         Ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbs) &
                  *DLOG(Ddx(nbs+1)/Ddx(1)) +DLOG(Ddx(1)))
      enddo
      do n = 1, nbs
         Dds(n) = DSQRT(Ddx(n)*Ddx(n+1))
         dts(n) = Ddx(n+1) - Ddx(n)
      enddo

!..Create bins of graupel (from min diameter up to 5 cm).
      Ddx(1) = D0g*1.0d0
      Ddx(nbg+1) = 0.05d0
      do n = 2, nbg
         Ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbg) &
                  *DLOG(Ddx(nbg+1)/Ddx(1)) +DLOG(Ddx(1)))
      enddo   
      do n = 1, nbg
         Ddg(n) = DSQRT(Ddx(n)*Ddx(n+1))
         dtg(n) = Ddx(n+1) - Ddx(n)
      enddo

      do i = 1, 256
         mp_debug(i:i) = char(0)
      enddo

      call radar_init

END SUBROUTINE GRAUPEL_INIT

!interface copied from new thompson interface
!and added NC, NS, NR, and NG variables.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE M2005MICRO_GRAUPEL) 
! WHICH OPERATES ON 1D VERTICAL COLUMNS.
! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE

! ******IMPORTANT******
! THIS CODE ASSUMES THE DRIVER MODEL USES PROCESS-SPLITTING FOR SOLVING THE TIME-DEPENDENT EQS.
! THUS, MODEL VARIABLES ARE UPDATED WITH MICROPHYSICS TENDENCIES INSIDE OF THE MICROPHYSICS
! SCHEME. THESE UPDATED VARIABLES ARE PASSED BACK TO DRIVER MODEL. THIS IS WHY THERE
! ARE NO TENDENCIES PASSED BACK AND FORTH BETWEEN DRIVER AND THE INTERFACE SUBROUTINE

! AN EXCEPTION IS THE TURBULENT MIXING TENDENCIES FOR DROPLET AND CLOUD ICE NUMBER CONCENTRATIONS
! (NCTEND, NITEND BELOW). FOR APPLICATION IN MODELS OTHER THAN WRF, TURBULENT MIXING TENDENCIES
! CAN BE ADDED TO THE VARIABLES ELSEWHERE (IN DRIVER OR PBL ROUTINE), AND THEN DON'T
! NEED TO BE PASSED INTO THE SUBROUTINE HERE.....

! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MP_GRAUPEL(ITIMESTEP,                       &
                TH, QV, QC, QR, QI, QS, QG, NI, NC, NS, NR, NG, TKE, NCTEND, &
                NITEND,KZH,  &
                RHO, PII, P, DT_IN, DZ, HT, W,          &
                RAINNC, RAINNCV, SR                     &
               ,EFFCS,EFFIS                             & ! HM ADD 4/13/07 
               ,refl_10cm                           & ! GT
!bloss               ,grid_clock                          & ! GT
!bloss               ,grid_alarms                         & ! GT
               ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
               ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
               ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
                                            )
 
! QV - water vapor mixing ratio (kg/kg)
! QC - cloud water mixing ratio (kg/kg)
! QR - rain water mixing ratio (kg/kg)
! QI - cloud ice mixing ratio (kg/kg)
! QS - snow mixing ratio (kg/kg)
! QG - graupel mixing ratio (KG/KG)
! NI - cloud ice number concentration (1/kg)
! NC - Droplet Number concentration (1/kg)
! NS - Snow Number concentration (1/kg)
! NR - Rain Number concentration (1/kg)
! NG - Graupel number concentration (1/kg)
! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
! P - AIR PRESSURE (PA)
! W - VERTICAL AIR VELOCITY (M/S)
! TH - POTENTIAL TEMPERATURE (K)
! PII - exner function - used to convert potential temp to temp
! DZ - difference in height over interface (m)
! DT_IN - model time step (sec)
! ITIMESTEP - time step counter
! RAINNC - accumulated grid-scale precipitation (mm)
! RAINNCV - one time step grid scale precipitation (mm/time step)
! SR - one time step mass ratio of snow to total precip
! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
!................................
! GRID_CLOCK, GRID_ALARMS - parameters to limit radar reflectivity calculation only when needed
! otherwise radar reflectivity calculation every time step is too slow
! only needed for coupling with WRF, see code below for details

! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

! ADDITIONAL INPUT NEEDED BY MICRO
! ********NOTE: WVAR IS SHOULD BE USED IN DROPLET ACTIVATION
! FOR CASES WHEN UPDRAFT IS NOT RESOLVED, EITHER BECAUSE OF
! LOW MODEL RESOLUTION OR CLOUD TYPE

! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                                       ims, ime, jms, jme, kms, kme , &
                                       its, ite, jts, jte, kts, kte
! Temporary changed from INOUT to IN

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          qv, qc, qr, qi, qs, qg, ni, nc, ns, nr, TH, NG, effcs, effis

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                          pii, p, dz, rho, w, tke, nctend, nitend,kzh
   REAL, INTENT(IN):: dt_in
   INTEGER, INTENT(IN):: ITIMESTEP

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          RAINNC, RAINNCV, SR
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &  ! GT
                          refl_10cm

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

!bloss      TYPE (WRFU_Clock):: grid_clock                  ! GT
!bloss      TYPE (WRFU_Alarm), POINTER:: grid_alarms(:)     ! GT

   ! LOCAL VARIABLES

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)::                     &
                      effi, effs, effr, EFFG

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)::                     &
                      T, WVAR, EFFC

   REAL, DIMENSION(kts:kte) ::                                                                & 
                            QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, NC_TEND1D,           &
                            NI_TEND1D, NS_TEND1D, NR_TEND1D,                                  &
                            QC1D, QI1D, QR1D, NC1D,NI1D, NS1D, NR1D, QS1D,                    &
                            T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, RHO1D, W1D, WVAR1D,         &
                            EFFC1D, EFFI1D, EFFS1D, EFFR1D,DZ1D,   &
   ! HM ADD GRAUPEL
                            QG_TEND1D, NG_TEND1D, QG1D, NG1D, EFFG1D, &

! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
                            QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN, &

! HM add reflectivity, bloss add lamc, pgam
                            dbz,  LAMC1D, PGAM1D
                          
   REAL, DIMENSION(kts:kte,nmicro_proc) :: micro_proc_rates

   REAL PRECPRT1D, SNOWRT1D

   INTEGER I,K,J
   
   REAL DT
   LOGICAL:: dBZ_tstep ! GT

! set dbz logical based on grid_clock
!+---+
! only calculate reflectivity when it is needed for output
! in this instance, logical dbz_tstep is set to .true.
! *******NOTE: FOR COUPLING WITH DRIVER MODEL OTHER THAN WRF,
! THIS BLOCK OF CODE WILL NEED TO BE MODIFIED TO CORRECTLY
! SET WHEN REFLECTIVIITY CALCULATION IS MADE

      dBZ_tstep = .false.
!bloss      if ( Is_alarm_tstep(grid_clock, grid_alarms(HISTORY_ALARM)) ) then
!bloss         dBZ_tstep = .true.
!bloss      endif

   ! Initialize tendencies (all set to 0) and transfer
   ! array to local variables
   DT = DT_IN   
   do I=ITS,ITE
   do J=JTS,JTE
   DO K=KTS,KTE
       T(I,K,J)        = TH(i,k,j)*PII(i,k,j)

! wvar is the ST. DEV. OF sub-grid vertical velocity, used for calculating droplet 
! activation rates.
! WVAR BE DERIVED EITHER FROM PREDICTED TKE (AS IN MYJ PBL SCHEME),
! OR FROM EDDY DIFFUSION COEFFICIENT KZH (AS IN YSU PBL SCHEME),
! DEPENDING ON THE PARTICULAR pbl SCHEME DRIVER MODEL IS COUPLED WITH
! NOTE: IF MODEL HAS HIGH ENOUGH RESOLUTION TO RESOLVE UPDRAFTS, WVAR IS 
! PROBABLY NOT NEEDED 

! for MYJ pbl scheme:
!       WVAR(I,K,J)     = (0.667*tke(i,k,j))**0.5
! for YSU pbl scheme:
       WVAR(I,K,J) = KZH(I,K,J)/20.
       WVAR(I,K,J) = MAX(0.1,WVAR(I,K,J))
       WVAR(I,K,J) = MIN(4.,WVAR(I,K,J))

! add tendency from pbl to droplet and cloud ice concentration
! NEEDED FOR WRF TEMPORARILY!!!!
! OTHER DRIVER MODELS MAY ADD TURBULENT DIFFUSION TENDENCY FOR
! SCALARS SOMEWHERE ELSE IN THE MODEL (I.E, NOT IN THE MICROPHYSICS)
! IN THIS CASE THESE 2 LINES BELOW MAY BE REMOVED
       nc(i,k,j) = nc(i,k,j)+nctend(i,k,j)*dt
       ni(i,k,j) = ni(i,k,j)+nitend(i,k,j)*dt
   END DO
   END DO
   END DO

   do i=its,ite      ! i loop (east-west)
   do j=jts,jte      ! j loop (north-south)
!
   ! Transfer 3D arrays into 1D for microphysical calculations
!

! hm , initialize 1d tendency arrays to zero

      do k=kts,kte   ! k loop (vertical)

          QC_TEND1D(k)  = 0.
          QI_TEND1D(k)  = 0.
          QNI_TEND1D(k) = 0.
          QR_TEND1D(k)  = 0.
          NC_TEND1D(k)  = 0.
          NI_TEND1D(k)  = 0.
          NS_TEND1D(k)  = 0.
          NR_TEND1D(k)  = 0.
          T_TEND1D(k)   = 0.
          QV_TEND1D(k)  = 0.

          QC1D(k)       = QC(i,k,j)
          QI1D(k)       = QI(i,k,j)
          QS1D(k)       = QS(i,k,j)
          QR1D(k)       = QR(i,k,j)

          NC1D(k)       = NC(i,k,j)
          NI1D(k)       = NI(i,k,j)

          NS1D(k)       = NS(i,k,j)
          NR1D(k)       = NR(i,k,j)
! HM ADD GRAUPEL
          QG1D(K)       = QG(I,K,j)
          NG1D(K)       = NG(I,K,j)
          QG_TEND1D(K)  = 0.
          NG_TEND1D(K)  = 0.

          T1D(k)        = T(i,k,j)
          QV1D(k)       = QV(i,k,j)
          P1D(k)        = P(i,k,j)
          RHO1D(k)      = P1D(K)/(R*T1D(K))
          DZ1D(k)       = DZ(i,k,j)
          W1D(k)        = W(i,k,j)
          WVAR1D(k)     = WVAR(i,k,j)
      end do

      !bloss: add extra argument for rho for consistency with below subroutine.
      !         done by repeating p1z.
      !         diable routine to make sure it is not used.
      STOP 'in mp_graupel wrapper routine.  Only use m2005micro_graupel()'

!bloss: comment out the call to m2005_graupel since the interface has changed.

!!$      call m2005micro_graupel(QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, NC_TEND1D,            &
!!$       NI_TEND1D, NS_TEND1D, NR_TEND1D,                                                  &
!!$       QC1D, QI1D, QS1D, QR1D, NC1D,NI1D, NS1D, NR1D,                                    &
!!$       T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, RHO1D, DZ1D, W1D, WVAR1D,                   &
!!$       PRECPRT1D,SNOWRT1D,                                                               &
!!$       EFFC1D,EFFI1D,EFFS1D,EFFR1D,DT,                                                   &
!!$                                            IMS,IME, JMS,JME, KMS,KME,                   &
!!$                                            ITS,ITE, JTS,JTE, KTS,KTE,                   & ! HM ADD GRAUPEL
!!$                                    QG_TEND1D,NG_TEND1D,QG1D,NG1D,EFFG1D, &
!!$! ADD SEDIMENTATION TENDENCIES
!!$                                  QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, &
!!$                                  PGAM1D, LAMC1D,MICRO_PROC_RATES,.true.)

   !
   ! Transfer 1D arrays back into 3D arrays
   !
      do k=kts,kte

! hm, add tendencies to update global variables 
! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

          QC(i,k,j)        = QC1D(k)
          QI(i,k,j)        = QI1D(k)
          QS(i,k,j)        = QS1D(k)
          QR(i,k,j)        = QR1D(k)
          NC(i,k,j)        = NC1D(k)
          NI(i,k,j)        = NI1D(k)
          NS(i,k,j)        = NS1D(k)          
          NR(i,k,j)        = NR1D(k)
	  QG(I,K,j)        = QG1D(K)
          NG(I,K,j)        = NG1D(K)

          T(i,k,j)         = T1D(k)
          TH(I,K,J)        = T(i,k,j)/PII(i,k,j) ! CONVERT TEMP BACK TO POTENTIAL TEMP
          QV(i,k,j)        = QV1D(k)

          EFFC(i,k,j)      = EFFC1D(k)
          EFFI(i,k,j)      = EFFI1D(k)
          EFFS(i,k,j)      = EFFS1D(k)
          EFFR(i,k,j)      = EFFR1D(k)
	  EFFG(I,K,j)      = EFFG1D(K)

! EFFECTIVE RADIUS FOR RADIATION CODE
! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
! LIMITS ARE FROM THE CAM MODEL APPLIED BY ANDREW GETTELMAN
          EFFCS(I,K,J)     = MIN(EFFC(I,K,J),16.)
          EFFCS(I,K,J)     = MAX(EFFCS(I,K,J),4.)
          EFFIS(I,K,J)     = MIN(EFFI(I,K,J),130.)
          EFFIS(I,K,J)     = MAX(EFFIS(I,K,J),13.)

      end do

! hm modified so that m2005 precip variables correctly match wrf precip variables
      RAINNC(i,j) = RAINNC(I,J)+PRECPRT1D
      RAINNCV(i,j) = PRECPRT1D
      SR(i,j) = SNOWRT1D/(PRECPRT1D+1.E-12)

! add reflectivity calculations
! only calculate if logical parameter dbz_tstep = .true.

         if (dBZ_tstep) then
          call calc_refl10cm (qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ,    &
                      kts, kte, i, j, nr1d, ns1d, ng1d)
          do k = kts, kte
             refl_10cm(i,k,j) = dBZ(k)
          enddo
         endif

   end do
   end do   

END SUBROUTINE MP_GRAUPEL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE M2005MICRO_GRAUPEL(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,QAD3DTEN,QAD2_3DTEN, QAW3DTEN,QAR3DTEN, &
       NC3DTEN,NI3DTEN,NS3DTEN,NR3DTEN,NAD3DTEN,NAD2_3DTEN,QC3D,QI3D,QNI3D,QR3D,QAD3D,QAD2_3D, QAW3D,QAR3D, &
       NC3D,NI3D,NS3D,NR3D,NAD3D,NAD2_3D,NUC3D,NUR3D,T3DTEN,QV3DTEN,T3D,QV3D,PRES,RHO,DZQ,W3D,WVAR,PRECRT,SNOWRT, &
       EFFC,EFFI,EFFS,EFFR,DT,                                                   &
                                            IMS,IME, JMS,JME, KMS,KME,           &
                                            ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
                                            QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,QGSTEN,QRSTEN, &
                                            QADSTEN,NADSTEN,QAWSTEN,QARSTEN,QISTEN,QNISTEN,QCSTEN,& !brnr: add aerosol tendencies
                                            PGAM, LAMC, MICRO_PROC_RATES, & !bloss: add radiative outputs
                                            DO_ACCUMULATE_MICRO_PROC_RATES )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
! MORRISON ET AL. 2005 JAS; MORRISON AND PINTO 2005 JAS.
! ADDITIONAL CHANGE IS ADDITION OF GRAUPEL MICROPHYSICS.
! SCHEME IS DESCRIBED IN DETAIL BY MORRISON ET AL. (MONTHLY WEATHER REVIEW, IN PREP.)

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.

! CODE STRUCTURE: MAIN SUBROUTINE IS 'M2005MICRO_GRAUPEL'. ALSO INCLUDED IN THIS FILE IS
! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
! 'FUNCTION GAMMA'.

! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! DECLARATIONS

      IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
! DEFINE ARRAY SIZES

! INPUT NUMBER OF GRID CELLS

! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
      INTEGER, INTENT( IN)  :: IMS,IME, JMS,JME, KMS,KME,          &
                               ITS,ITE, JTS,JTE, KTS,KTE

      REAL, DIMENSION(KMS:KME) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)

!brnr: make tendencies for aerosol mass available
      REAL, DIMENSION(KMS:KME) ::  QAD3DTEN           ! DRY AEROSOL ACCUM MASS MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QAD2_3DTEN           ! DRY AEROSOL AITKEN MASS MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QAW3DTEN           ! WET AEROSOL MASS MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QAR3DTEN           ! RAIN AEROSOL MASS MIXING RATIO TENDENCY (KG/KG/S)

      REAL, DIMENSION(KMS:KME) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)

!brnr: make tendencies for aerosol number available
      REAL, DIMENSION(KMS:KME) ::  NAD3DTEN           ! DRY AEROSOL ACCUM NUMBER CONCENTRATION TENDENCY (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NAD2_3DTEN           ! DRY AEROSOL AITKEN NUMBER CONCENTRATION TENDENCY (1/KG/S)

      REAL, DIMENSION(KMS:KME) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)

!brnr: add vars for gamma exponents for use in scavenging
      REAL, DIMENSION(KMS:KME) :: NUC3D               ! CLOUD GAMMA EXPONENT
      REAL, DIMENSION(KMS:KME) :: NUR3D               ! RAIN GAMMA EXPNONET (currently 0)

!brnr: add mass for wet (rain and cloud) and dry aerosol     
      REAL, DIMENSION(KMS:KME) ::  QAD3D              ! DRY AEROSOL ACCUM MASS MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QAD2_3D              ! DRY AEROSOL AITKEN MASS MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QAW3D              ! WET AEROSOL MASS MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QAR3D              ! RAIN AEROSOL MASS MIXING RATIO (KG/KG)

      REAL, DIMENSION(KMS:KME) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)

!brnr: add number for wet and dry aerosol
      REAL, DIMENSION(KMS:KME) ::  NAD3D              ! DRY AEROSOL ACCUM NUMBER CONCENTRATION (1/KG)

      REAL, DIMENSION(KMS:KME) ::  NAD2_3D              ! DRY AEROSOL AITKEN NUMBER CONCENTRATION (1/KG)

      REAL, DIMENSION(KMS:KME) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
      REAL, DIMENSION(KMS:KME) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  T3D                ! TEMPERATURE (K)
      REAL, DIMENSION(KMS:KME) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
!bloss: make rho an input argument
      REAL, DIMENSION(KMS:KME), INTENT(IN) ::  RHO   ! AIR DENSITY 
      REAL, DIMENSION(KMS:KME) ::  DZQ                ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
      REAL, DIMENSION(KMS:KME) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
      REAL, DIMENSION(KMS:KME) ::  WVAR               ! SUB-GRID VERTICAL VELOCITY (M/S)

! HM ADDED GRAUPEL VARIABLES
      REAL, DIMENSION(KMS:KME) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

      REAL, DIMENSION(KMS:KME) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
     
!brnr: make sedimentation tendencies available for aerosol
      REAL, DIMENSION(KMS:KME) :: QADSTEN
      REAL, DIMENSION(KMS:KME) :: QAWSTEN
      REAL, DIMENSION(KMS:KME) :: QARSTEN
      REAL, DIMENSION(KMS:KME) :: NADSTEN
      REAL, DIMENSION(KMS:KME) :: NCSTEN
      REAL, DIMENSION(KMS:KME) :: NRSTEN

      REAL, DIMENSION(KMS:KME) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)      

! OUTPUT VARIABLES

        REAL PRECRT                ! TOTAL PRECIP PER TIME STEP (mm)
        REAL SNOWRT                ! SNOW PER TIME STEP (mm)

        REAL, DIMENSION(KMS:KME) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)

! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

        REAL DT         ! MODEL TIME STEP (SEC)

!.....................................................................................................
! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
! REST OF THE MODEL.

! SIZE PARAMETER VARIABLES

     REAL, DIMENSION(KMS:KME), INTENT(OUT) :: LAMC          ! SLOPE PARAMETER FOR DROPLETS (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMI          ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMS          ! SLOPE PARAMETER FOR SNOW (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMR          ! SLOPE PARAMETER FOR RAIN (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMG          ! SLOPE PARAMETER FOR GRAUPEL (M-1)
     REAL, DIMENSION(KMS:KME) :: CDIST1        ! PSD PARAMETER FOR DROPLETS
     REAL, DIMENSION(KMS:KME) :: N0I           ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: N0S           ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: N0RR          ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: N0G           ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
     REAL, DIMENSION(KMS:KME), INTENT(OUT) :: PGAM          ! SPECTRAL SHAPE PARAMETER FOR DROPLETS

     REAL, DIMENSION(KMS:KME,nmicro_proc), INTENT(INOUT) :: &
          MICRO_PROC_RATES ! ARRAY FOR ACCUMULATION OF MICROPHYSICS PROC RATES
     LOGICAL, INTENT(IN) :: DO_ACCUMULATE_MICRO_PROC_RATES ! FLAG

! MICROPHYSICAL PROCESSES

     REAL, DIMENSION(KMS:KME) ::  NSUBC     ! LOSS OF NC DURING EVAP
     REAL, DIMENSION(KMS:KME) ::  NSUBI     ! LOSS OF NI DURING SUB.
     REAL, DIMENSION(KMS:KME) ::  NSUBS     ! LOSS OF NS DURING SUB.
     REAL, DIMENSION(KMS:KME) ::  NSUBR     ! LOSS OF NR DURING EVAP
     REAL, DIMENSION(KMS:KME) ::  PRD       ! DEP CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  PRE       ! EVAP OF RAIN
     REAL, DIMENSION(KMS:KME) ::  PRDS      ! DEP SNOW
     REAL, DIMENSION(KMS:KME) ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     REAL, DIMENSION(KMS:KME) ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     REAL, DIMENSION(KMS:KME) ::  PRA       ! ACCRETION DROPLETS BY RAIN
     REAL, DIMENSION(KMS:KME) ::  PRC       ! AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KMS:KME) ::  PCC       ! COND/EVAP DROPLETS
     REAL, DIMENSION(KMS:KME) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KMS:KME) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KMS:KME) ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
     REAL, DIMENSION(KMS:KME) ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
     REAL, DIMENSION(KMS:KME) ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
     REAL, DIMENSION(KMS:KME) ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
     REAL, DIMENSION(KMS:KME) ::  NSAGG     ! SELF-COLLECTION OF SNOW
     REAL, DIMENSION(KMS:KME) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KMS:KME) ::  NPRC1     ! CHANGE NR AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KMS:KME) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  PRCI      ! CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
     REAL, DIMENSION(KMS:KME) ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
     REAL, DIMENSION(KMS:KME) ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
     REAL, DIMENSION(KMS:KME) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
     REAL, DIMENSION(KMS:KME) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
     REAL, DIMENSION(KMS:KME) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
     REAL, DIMENSION(KMS:KME) ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
     REAL, DIMENSION(KMS:KME) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
     REAL, DIMENSION(KMS:KME) ::  NSMLTS    ! CHANGE N MELTING SNOW
     REAL, DIMENSION(KMS:KME) ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN

! BRNR ADDED 2/15/12 EXTRA DIAGNOSTIC AEROSOL TENDENCIES 
     REAL, DIMENSION(KMS:KME) ::  NARG1      ! ARG predicted N1
     REAL, DIMENSION(KMS:KME) ::  NARG2      ! ARG predicted N2
     
     REAL, DIMENSION(KMS:KME) ::  NACTRATE   ! Activation rate of N
     REAL, DIMENSION(KMS:KME) ::  QACTRATE   ! Activation rate of N
     REAL, DIMENSION(KMS:KME) ::  NACTDIFF   !  Difference between ARG N activaed and OLD NC
     REAL, DIMENSION(KMS:KME) ::  NATRANSFER  ! TRANSFER FROM ACCUM TO AITKEN DRY TO COUNTER ACTIVATION
     REAL, DIMENSION(KMS:KME) ::  QATRANSFER
     REAL, DIMENSION(KMS:KME) ::  DC1       ! Critical diameter of activation
     REAL, DIMENSION(KMS:KME) ::  DC2       ! Critical diameter of activation
     REAL, DIMENSION(KMS:KME) ::  DG1       ! Modal diameter
     REAL, DIMENSION(KMS:KME) ::  DG2
     REAL, DIMENSION(KMS:KME) ::  ISACT     ! keep track of which levels are going through activation
     REAL, DIMENSION(KMS:KME) ::  SSPK      ! ARG Peak Supersaturation
     REAL, DIMENSION(KMS:KME) ::  QAPRA     ! CHANGE Q DUE TO AUTOCONVERSION
     REAL, DIMENSION(KMS:KME) ::  QAPRC     ! CHANGE Q DUE TO ACCRETION
     REAL, DIMENSION(KMS:KME) ::  QAPRE     ! CHANGE Q DUE TO EVAP OF RAIN
     REAL, DIMENSION(KMS:KME) ::  QASUBC    ! CHANGE Q DUE TO EVAP OF CLOUD

     REAL, DIMENSION(KMS:KME,30) :: TEN_ARRAY ! ARRAY TO HOLD ALL DIAGNOSTIC TENDENCIES

! BRNR ADDED 2/13/12
     REAL, DIMENSION(KMS:KME) ::  NCPOSLIM  ! POS CHANGE NC DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  NCNEGLIM  ! NEG CHANGE NC DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  NRPOSLIM  ! POS CHANGE NR DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  NRNEGLIM  ! NEG CHANGE NR DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  NADPOSLIM ! POS CHANGE NAD DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  NADNEGLIM ! NEG CHANGE NAD DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  QVPOSLIM  ! POS CHANGE QV DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  QCNEGLIM  ! NEG CHANGE QC DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  QRNEGLIM  ! NEG CHANGE QR DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  QAWNEGLIM ! NEG CHANGE QAW DUE TO LIMITERS
     REAL, DIMENSION(KMS:KME) ::  QARNEGLIM ! NEG CHANGE QAR DUE TO LIMITERS
! HM ADDED 12/13/06
     REAL, DIMENSION(KMS:KME) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
     REAL, DIMENSION(KMS:KME) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
     REAL, DIMENSION(KMS:KME) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
     REAL, DIMENSION(KMS:KME) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KMS:KME) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KMS:KME) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KMS:KME) ::  EPRD      ! SUBLIMATION CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  EPRDS     ! SUBLIMATION SNOW
! HM ADDED GRAUPEL PROCESSES
     REAL, DIMENSION(KMS:KME) ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PRDG    ! DEP OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  EPRDG    ! SUB OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     REAL, DIMENSION(KMS:KME) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
     REAL, DIMENSION(KMS:KME) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
     REAL, DIMENSION(KMS:KME) ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
! hm new process rate output
     REAL, DIMENSION(KMS:KME) ::  QHOMOC    ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
     REAL, DIMENSION(KMS:KME) ::  QHOMOR    ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF RAIN
     REAL, DIMENSION(KMS:KME) ::  NHOMOC    ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
     REAL, DIMENSION(KMS:KME) ::  NHOMOR    ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF RAIN
     REAL, DIMENSION(KMS:KME) ::  QMELTI    ! CHANGE Q DUE TO MELTING OF CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NMELTI    ! CHANGE N DUE TO MELTING OF CLOUD ICE


! TIME-VARYING ATMOSPHERIC PARAMETERS

     REAL, DIMENSION(KMS:KME) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
     REAL, DIMENSION(KMS:KME) ::   EVS   ! SATURATION VAPOR PRESSURE
     REAL, DIMENSION(KMS:KME) ::   EIS   ! ICE SATURATION VAPOR PRESSURE
     REAL, DIMENSION(KMS:KME) ::   QVS   ! SATURATION MIXING RATIO
     REAL, DIMENSION(KMS:KME) ::   QVI   ! ICE SATURATION MIXING RATIO
     REAL, DIMENSION(KMS:KME) ::   QVQVS ! SAUTRATION RATIO
     REAL, DIMENSION(KMS:KME) ::   QVQVSI! ICE SATURAION RATIO
     REAL, DIMENSION(KMS:KME) ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR
     REAL, DIMENSION(KMS:KME) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
     REAL, DIMENSION(KMS:KME) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
     REAL, DIMENSION(KMS:KME) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
     REAL, DIMENSION(KMS:KME) ::   MU    ! VISCOCITY OF AIR
     REAL, DIMENSION(KMS:KME) ::   SC    ! SCHMIDT NUMBER
     REAL, DIMENSION(KMS:KME) ::   XLF   ! LATENT HEAT OF FREEZING
!bloss     REAL, DIMENSION(KMS:KME) ::   RHO   ! AIR DENSITY
     REAL, DIMENSION(KMS:KME) ::   AB    ! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
     REAL, DIMENSION(KMS:KME) ::   ABI    ! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

! TIME-VARYING MICROPHYSICS PARAMETERS

     REAL, DIMENSION(KMS:KME) ::   DAP    ! DIFFUSIVITY OF AEROSOL
     REAL    NACNT                    ! NUMBER OF CONTACT IN
     REAL    FMULT                    ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
     REAL    COFFI                    ! ICE AUTOCONVERSION PARAMETER

! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

      REAL, DIMENSION(KMS:KME) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
      REAL UNI, UMI,UMR
      REAL, DIMENSION(KMS:KME) ::    FR, FI, FNI,FG,FNG
      REAL RGVM
      REAL, DIMENSION(KMS:KME) ::   FALOUTR,FALOUTI,FALOUTNI
      REAL FALTNDR,FALTNDI,FALTNDNI,RHO2
      REAL, DIMENSION(KMS:KME) ::   DUMQS,DUMFNS
      REAL UMS,UNS
      REAL, DIMENSION(KMS:KME) ::   FS,FNS, FALOUTS,FALOUTNS,FALOUTG,FALOUTNG
      REAL FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
      REAL, DIMENSION(KMS:KME) ::    DUMC,DUMFNC
      REAL UNC,UMC,UNG,UMG
      REAL, DIMENSION(KMS:KME) ::   FC,FALOUTC,FALOUTNC
      REAL FALTNDC,FALTNDNC
      REAL, DIMENSION(KMS:KME) ::   FNC,DUMFNR,FALOUTNR
      REAL FALTNDNR
      REAL, DIMENSION(KMS:KME) ::   FNR
      REAL, DIMENSION(KMS:KME) ::   DUMQAW,DUMQAR,FALOUTQAW,FALOUTQAR
      REAL FALTNDQAW,FALTNDQAR

! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

      REAL, DIMENSION(KMS:KME) ::    AIN,ARN,ASN,ACN,AGN

! EXTERNAL FUNCTION CALL RETURN VARIABLES

!      REAL GAMMA,      ! EULER GAMMA FUNCTION
!      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
!      REAL DERF1        ! ERROR FUNCTION

! DUMMY VARIABLES

     REAL DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS

! DEBUG VARIABLES
      
     INTEGER DBGIND

     REAL, DIMENSION(6,KMS:KME) :: DBGQAR3DTEN, DBGQAR3D
     LOGICAL :: QNEG

! PROGNOSTIC SUPERSATURATION

     REAL DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
     REAL DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
     REAL EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
     REAL EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
     REAL EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
     REAL EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

! NEW DROPLET ACTIVATION VARIABLES
     REAL TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
     REAL TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
     REAL TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
     REAL TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
     REAL TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
     REAL DUMACT,DUM3,DUMACTM, DUMACT1, DUMACT2, DUMACTM1, DUMACTM2, D3

! COUNTING/INDEX VARIABLES

! V1.3 DIFFERENT NSTEP FOR EACH SPECIES
     INTEGER K,N,NSTEPR,NSTEPI,NSTEPS,NSTEPC,NSTEPG ! ,I

! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN, 
!               = 1, HYDROMETEORS IN COLUMN

      INTEGER LTRUE

! DROPLET ACTIVATION/FREEZING AEROSOL


     REAL    CT      ! DROPLET ACTIVATION PARAMETER
     REAL    TEMP1   ! DUMMY TEMPERATURE
     REAL    SAT1    ! DUMMY SATURATION
     REAL    SIGVL   ! SURFACE TENSION LIQ/VAPOR
     REAL    KEL     ! KELVIN PARAMETER
     REAL    KC2     ! TOTAL ICE NUCLEATION RATE

       REAL CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

! MORE WORKING/DUMMY VARIABLES

     REAL DUMQI,DUMNI,DC0,DS0,DG0
     REAL DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

! EFFECTIVE VERTICAL VELOCITY  (M/S)
     REAL WEF

! WORKING PARAMETERS FOR ICE NUCLEATION

      REAL ANUC,BNUC

! WORKING PARAMETERS FOR AEROSOL ACTIVATION

      REAL AACT,GAMM,GG,PSI,ETA1,ETA2,ETA3,SM1,SM2,SM3,SMAX,UU1,UU2,ALPHA
      REAL NTRANS, QTRANS, TEMP_NAD, TEMP_QAD
      REAL NCNEW, QAWNEW, NADNEW
! DUMMY SIZE DISTRIBUTION PARAMETERS

        REAL DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN

        INTEGER IDROP

        INTEGER idx !bloss: index for process rate array

! v1.4
! new variables for seifert and beheng warm rain scheme
      REAL, DIMENSION(KMS:KME) :: nu
      integer dumii

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! SET LTRUE INITIALLY TO 0

         LTRUE = 0

! V.13 initialize effective radii to default values (from P. Blossey)
         effc(kts:kte) = 25.         
         effi(kts:kte) = 25.         
         effs(kts:kte) = 25.         
         effr(kts:kte) = 25.         
         effg(kts:kte) = 25.         

!bloss: Initialize pgam and lamc so that it is defined everywhere
         PGAM(KTS:KTE) = -10.
         LAMC(KTS:KTE) = -10.

! INIT DIAGNOSTIC PROFILES OF CLOUD AND RAIN SHAPE PARAMETER
            NUC3D(KTS:KTE) = 0.
            NUR3D(KTS:KTE) = 0.

! INIT LIMITER PROFILES TO ZERO
            NCPOSLIM(KTS:KTE) = 0.
            NCNEGLIM(KTS:KTE) = 0. 
            NRPOSLIM(KTS:KTE) = 0. 
            NRNEGLIM(KTS:KTE) = 0. 
            NADPOSLIM(KTS:KTE) = 0.
            NADNEGLIM(KTS:KTE) = 0.
            QVPOSLIM(KTS:KTE) = 0. 
            QCNEGLIM(KTS:KTE) = 0. 
            QRNEGLIM(KTS:KTE) = 0. 
            QAWNEGLIM(KTS:KTE) = 0.
            QARNEGLIM(KTS:KTE) = 0.

            ! INIT PROCESS PROFILES TO ZERO
            !bloss: Copied from Berner code -- there may be some duplicate process initialization here.
            NARG1(KTS:KTE) = 0. 
            NARG2(KTS:KTE) = 0. 
            NACTRATE(KTS:KTE) = 0.
            QACTRATE(KTS:KTE) = 0. 
            NACTDIFF(KTS:KTE) = 0. 
            NATRANSFER(KTS:KTE) = 0.
            QATRANSFER(KTS:KTE) = 0.
            DC1(KTS:KTE) = 0.
            DC2(KTS:KTE) = 0.
            DG1(KTS:KTE) = 0.
            DG2(KTS:KTE) = 0.
            ISACT(KTS:KTE) = 0.
            DC2(KTS:KTE) = 0.
            SSPK(KTS:KTE) = 0.
            QAPRA(KTS:KTE) = 0.
            QAPRC(KTS:KTE) = 0.
            QAPRE(KTS:KTE) = 0.
            QASUBC(KTS:KTE) = 0.
            PRC(KTS:KTE) = 0.
            NPRC(KTS:KTE) = 0.
            NPRC1(KTS:KTE) = 0.
            PRA(KTS:KTE) = 0.
            NPRA(KTS:KTE) = 0.
            NRAGG(KTS:KTE) = 0.
            PSMLT(KTS:KTE) = 0.
            NSMLTS(KTS:KTE) = 0.
            NSMLTR(KTS:KTE) = 0.
            EVPMS(KTS:KTE) = 0.
            PCC(KTS:KTE) = 0.
            PRE(KTS:KTE) = 0.
            NSUBC(KTS:KTE) = 0.
            NSUBR(KTS:KTE) = 0.
            PRACG(KTS:KTE) = 0.
            NPRACG(KTS:KTE) = 0.
            PSMLT(KTS:KTE) = 0.
            EVPMS(KTS:KTE) = 0.
            PGMLT(KTS:KTE) = 0.
            EVPMG(KTS:KTE) = 0.
            PRACS(KTS:KTE) = 0.
            NPRACS(KTS:KTE) = 0.
            NGMLTG(KTS:KTE) = 0.
            NGMLTR(KTS:KTE) = 0.

         if(do_accumulate_micro_proc_rates) then
           NSUBC    = 0. ! LOSS OF NC DURING EVAP
           NSUBI    = 0. ! LOSS OF NI DURING SUB.
           NSUBS    = 0. ! LOSS OF NS DURING SUB.
           NSUBR    = 0. ! LOSS OF NR DURING EVAP
           PRD      = 0. ! DEP CLOUD ICE
           PRE      = 0. ! EVAP OF RAIN
           PRDS     = 0. ! DEP SNOW
           NNUCCC   = 0. ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
           MNUCCC   = 0. ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
           PRA      = 0. ! ACCRETION DROPLETS BY RAIN
           PRC      = 0. ! AUTOCONVERSION DROPLETS
           PCC      = 0. ! COND/EVAP DROPLETS
           NNUCCD   = 0. ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
           MNUCCD   = 0. ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
           MNUCCR   = 0. ! CHANGE Q DUE TO CONTACT FREEZ RAIN
           NNUCCR   = 0. ! CHANGE N DUE TO CONTACT FREEZ RAIN
           NPRA     = 0. ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
           NRAGG    = 0. ! SELF-COLLECTION OF RAIN
           NSAGG    = 0. ! SELF-COLLECTION OF SNOW
           NPRC     = 0. ! CHANGE NC AUTOCONVERSION DROPLETS
           NPRC1     = 0. ! CHANGE NR AUTOCONVERSION DROPLETS
           PRAI     = 0. ! CHANGE Q ACCRETION CLOUD ICE
           PRCI     = 0. ! CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW
           PSACWS   = 0. ! CHANGE Q DROPLET ACCRETION BY SNOW
           NPSACWS  = 0. ! CHANGE N DROPLET ACCRETION BY SNOW
           PSACWI   = 0. ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
           NPSACWI  = 0. ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
           NPRCI    = 0. ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
           NPRAI    = 0. ! CHANGE N ACCRETION CLOUD ICE
           NMULTS   = 0. ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
           NMULTR   = 0. ! ICE MULT DUE TO RIMING RAIN BY SNOW
           QMULTS   = 0. ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
           QMULTR   = 0. ! CHANGE Q DUE TO ICE RAIN/SNOW
           PRACS    = 0. ! CHANGE Q RAIN-SNOW COLLECTION
           NPRACS   = 0. ! CHANGE N RAIN-SNOW COLLECTION
           PCCN     = 0. ! CHANGE Q DROPLET ACTIVATION
           PSMLT    = 0. ! CHANGE Q MELTING SNOW TO RAIN
           EVPMS    = 0. ! CHNAGE Q MELTING SNOW EVAPORATING
           NSMLTS   = 0. ! CHANGE N MELTING SNOW
           NSMLTR   = 0. ! CHANGE N MELTING SNOW TO RAIN
           PIACR    = 0. ! CHANGE QR, ICE-RAIN COLLECTION
           NIACR    = 0. ! CHANGE N, ICE-RAIN COLLECTION
           PRACI    = 0. ! CHANGE QI, ICE-RAIN COLLECTION
           PIACRS    = 0. ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
           NIACRS    = 0. ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
           PRACIS    = 0. ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
           EPRD     = 0. ! SUBLIMATION CLOUD ICE
           EPRDS    = 0. ! SUBLIMATION SNOW
           PRACG   = 0. ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
           PSACWG   = 0. ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
           PGSACW   = 0. ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
           PGRACS   = 0. ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
           PRDG   = 0. ! DEP OF GRAUPEL
           EPRDG   = 0. ! SUB OF GRAUPEL
           EVPMG   = 0. ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
           PGMLT   = 0. ! CHANGE Q MELTING OF GRAUPEL
           NPRACG   = 0. ! CHANGE N COLLECTION RAIN BY GRAUPEL
           NPSACWG   = 0. ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
           NSCNG   = 0. ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
           NGRACS   = 0. ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
           NGMLTG   = 0. ! CHANGE N MELTING GRAUPEL
           NGMLTR   = 0. ! CHANGE N MELTING GRAUPEL TO RAIN
           NSUBG   = 0. ! CHANGE N SUB/DEP OF GRAUPEL
           PSACR   = 0. ! CONVERSION DUE TO COLL OF SNOW BY RAIN
           NMULTG   = 0. ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
           NMULTRG   = 0. ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
           QMULTG   = 0. ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
           QMULTRG   = 0. ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
! hm new process rate output
           QHOMOC = 0.    ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
           QHOMOR = 0.   ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF RAIN
           NHOMOC = 0.    ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
           NHOMOR = 0.    ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF RAIN
           QMELTI = 0.    ! CHANGE Q DUE TO MELTING OF CLOUD ICE
           NMELTI = 0.    ! CHANGE N DUE TO MELTING OF CLOUD ICE

   end if

! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
         DO K = KTS,KTE
      
! LATENT HEAT OF VAPORATION

            XXLV(K) = lcond !bloss 3.1484E6-2370.*T3D(K)

! LATENT HEAT OF SUBLIMATION

            XXLS(K) = lsub !bloss 3.15E6-2370.*T3D(K)+0.3337E6

            CPM(K) = cp !bloss CP*(1.+0.887*QV3D(K))

! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = .622*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = .622*EIS(K)/(PRES(K)-EIS(K))

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! V1.3, change limit from 10^-7 to 10^-6
! V1.7 7/9/09 change limit from 10^-6 to 10^-8
! this improves reflectivity at low mixing ratios

             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QVPOSLIM(K) = QVPOSLIM(K)+QR3D(K)/DT !Limiter adds QR to QV
                  QRNEGLIM(K) = QRNEGLIM(K)-QR3D(K)/DT !Limiter zeros QR
                  NRNEGLIM(K) = NRNEGLIM(K)-NR3D(K)/DT !Limiter zeros NR
                  NADPOSLIM(K) = NADPOSLIM(K)+NR3D(K)/DT !Limiter adds NR to NAD
                  QARNEGLIM(K) = QARNEGLIM(K)-QAR3D(K)/DT ! Limiter zeros QAR
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
                  QAD3D(K)=QAD3D(K)+QAR3D(K) !move evaporated aerosol mass
                  QAR3D(K)=0.
                  NAD3D(K)=NAD3D(K)+NR3D(K) !move evaporated number
                  NR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QVPOSLIM(K) = QVPOSLIM(K)+QC3D(K)/DT !Limiter adds QC to QV
                  QCNEGLIM(K) = QCNEGLIM(K)-QC3D(K)/DT !Limiter zeros QC
                  NCNEGLIM(K) = NCNEGLIM(K)-NC3D(K)/DT !Limiter zeros NC
                  NADPOSLIM(K) = NADPOSLIM(K)+NC3D(K)/DT !Limiter adds NC to NAD
                  QAWNEGLIM(K) = QAWNEGLIM(K)-QAW3D(K)/DT !Limiter zeros QAW
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
                  QAD3D(K)=QAD3D(K)+QAW3D(K) !move evaporated aerosol mass
                  QAW3D(K)=0.
                  NAD3D(K)=NAD3D(K)+NC3D(K) !move evaporated number
                  NC3D(K)=0.
               END IF
             END IF

             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
             END IF

! AIR DENSITY

!bloss: now an input argument            RHO(K) = PRES(K)/(R*T3D(K))

! HEAT OF FUSION

            XLF(K) = XXLS(K)-XXLV(K)

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

            IF (QC3D(K).LT.QSMALL) THEN
               QCNEGLIM(K) = QCNEGLIM(K)-QC3D(K)/DT
               NCNEGLIM(K) = NCNEGLIM(K)-NC3D(K)/DT
               QAWNEGLIM(K) = QAWNEGLIM(K)-QAW3D(K)/DT
               QV3D(K) = MAX(1.e-10, QV3D(K) + QC3D(K)) !bloss: Conserve!!!
               NAD3D(K) = MAX(0., NAD3D(K) + NC3D(K)) !bloss: Conserve!!!
               QAD3D(K) = MAX(0., QAD3D(K) + QAW3D(K)) !bloss: Conserve!!!
               QC3D(K) = 0.
               NC3D(K) = 0.
               QAW3D(K) = 0.
            END IF
            IF (QR3D(K).LT.QSMALL) THEN
               QRNEGLIM(K) = QRNEGLIM(K)-QR3D(K)/DT
               NRNEGLIM(K) = NRNEGLIM(K)-NR3D(K)/DT
               QARNEGLIM(K) = QARNEGLIM(K)-QAR3D(K)/DT
               QV3D(K) = MAX(1.e-10, QV3D(K) + QR3D(K)) !bloss: Conserve!!!
               NAD3D(K) = MAX(0., NAD3D(K) + NR3D(K)) !bloss: Conserve!!!
               QAD3D(K) = MAX(0., QAD3D(K) + QAR3D(K)) !bloss: Conserve!!!
               QR3D(K) = 0.
               NR3D(K) = 0.
               QAR3D(K) = 0.
            END IF
            IF (QI3D(K).LT.QSMALL) THEN
               QI3D(K) = 0.
               NI3D(K) = 0.
            END IF
            IF (QNI3D(K).LT.QSMALL) THEN
               QNI3D(K) = 0.
               NS3D(K) = 0.
            END IF
            IF (QG3D(K).LT.QSMALL) THEN
               QG3D(K) = 0.
               NG3D(K) = 0.
            END IF
            
! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO
            
            QRSTEN(K) = 0.
            QISTEN(K) = 0.
            QNISTEN(K) = 0.
            QCSTEN(K) = 0.
            QGSTEN(K) = 0.
!brnr
            QADSTEN(K) = 0.
            QAWSTEN(K) = 0.
            QARSTEN(K) = 0.

! INITIALIZE SEDIMENTATION TENDENCIES FOR NUMBER CONCENTRATION
            NADSTEN(K) = 0. 
            NCSTEN(K) = 0.
            NRSTEN(K) = 0.

!..................................................................
! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

! DYNAMIC VISCOSITY OF AIR

            MU(K) = 1.496E-6*T3D(K)**1.5/(T3D(K)+120.)
            
! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

            DUM = (RHOSU/RHO(K))**0.54
! v3 5/27/11
!            AIN(K) = DUM*AI
! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction 
! hm fix 11/18/11
            AIN(K) = (RHOSU/RHO(K))**0.35*AI
            ARN(K) = DUM*AR
            ASN(K) = DUM*AS
!            ACN(K) = DUM*AC
! AA revision 4/1/11: temperature-dependent Stokes fall speed
            ACN(K) = G*RHOW/(18.*MU(K))
! HM ADD GRAUPEL 8/28/06
            AGN(K) = DUM*AG

! V1.7
! bug fix 7/10/09 
!hm 4/15/09 bug fix, initialize lami to prevent later division by zero
            LAMI(K)=0.

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
! FOR THIS LEVEL

            IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) THEN
                 IF (T3D(K).LT.TMELT.AND.QVQVSI(K).LT.0.999) GOTO 200
                 IF (T3D(K).GE.TMELT.AND.QVQVS(K).LT.0.999) GOTO 200
            END IF

! THERMAL CONDUCTIVITY FOR AIR

! v3 5/27/11
            KAP(K) = 1.414E3*MU(K)

! DIFFUSIVITY OF WATER VAPOR

            DV(K) = 8.794E-5*T3D(K)**1.81/PRES(K)

! SCHMIT NUMBER

! v3 5/27/11
            SC(K) = MU(K)/(RHO(K)*DV(K))

! PSYCHOMETIC CORRECTIONS

! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

            DUM = (RV*T3D(K)**2)
            
            DQSDT = XXLV(K)*QVS(K)/DUM
            DQSIDT =  XXLS(K)*QVI(K)/DUM
            
            ABI(K) = 1.+DQSIDT*XXLS(K)/CPM(K)
            AB(K) = 1.+DQSDT*XXLV(K)/CPM(K)

! 
!.....................................................................
!.....................................................................
! CASE FOR TEMPERATURE ABOVE FREEZING

            IF (T3D(K).GE.TMELT) THEN

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

               IF (INUM.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
                  NC3D(K)=NDCNST*1.E6/RHO(K)
               END IF

! GET SIZE DISTRIBUTION PARAMETERS

! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
               IF (QNI3D(K).LT.1.E-6) THEN
                  QR3D(K)=QR3D(K)+QNI3D(K)
                  NR3D(K)=NR3D(K)+NS3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
                  QNI3D(K) = 0.
                  NS3D(K) = 0.
               END IF
               IF (QG3D(K).LT.1.E-6) THEN
                  QR3D(K)=QR3D(K)+QG3D(K)
                  NR3D(K)=NR3D(K)+NG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
                  QG3D(K) = 0.
                  NG3D(K) = 0.
               END IF

               IF (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8) GOTO 300

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE
! BRNR 2/13/12 Add tracking of tendency due to this limiter

               NS3D(K) = MAX(0.,NS3D(K))
               
               DUM = NC3D(K)
               NC3D(K) = MAX(0.,NC3D(K))
               NCPOSLIM(K) = NCPOSLIM(K)+(NC3D(K)-DUM)/DT
               
               DUM = NR3D(K)
               NR3D(K) = MAX(0.,NR3D(K))
               NRPOSLIM(K) = NRPOSLIM(K)+(NR3D(K)-DUM)/DT
               
               NG3D(K) = MAX(0.,NG3D(K))

               DUM = NAD3D(K)
               NAD3D(K) = MAX(0.,NAD3D(K))
               NADPOSLIM(K) = NADPOSLIM(K)+(NAD3D(K)-DUM)/DT

!......................................................................
! RAIN

               IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
                 END IF
               END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

               IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
                  IF(dofix_pgam) THEN
                     pgam(k) = pgam_fixed
                  ELSE

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
                     PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
                     PGAM(K)=1./(PGAM(K)**2)-1.
                     PGAM(K)=MAX(PGAM(K),2.)
                     PGAM(K)=MIN(PGAM(K),10.)
                     
                  END IF
! v1.4
! interpolate
                  dumii=int(pgam(k))
                  nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
                       (pgam(k)-real(dumii))

! CALCULATE LAMC

                  LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                       (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

! Brnr modify the bound check routines to conserve Na+Nc
                  LAMMIN = (PGAM(K)+1.)/60.E-6
                  LAMMAX = (PGAM(K)+1.)/1.E-6

                  DUM1 = NC3D(K)
                  DUM2 = NAD3D(K)

                  IF (LAMC(K).LT.LAMMIN) THEN
                     LAMC(K) = LAMMIN
                     IF (IPRGAER.EQ.1) THEN
                        DUM = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                             LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                        !bloss(2018-02): prevent creation of number
                        !bloss(2019-03): prevent problems if total cloud+aerosol number is zero or negative
                        NC3D(K) = MAX(0.1*DUM, MIN(DUM, DUM1+DUM2) )                        
                        NAD3D(K) = MAX(0., DUM1+DUM2 -NC3D(K) )
                        ! RE-CALCULATE LAMC
                        LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                             (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
!bloss(2018-02)                        NAD3D(K) = MAX(0.,NAD3D(K)-(DUM-NC3D(K)))
!bloss(2018-02)                        NC3D(K) = DUM
                     ELSE
                        NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                             LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26              
                     END IF !IPRGAER.EQ.1
                  ELSE IF (LAMC(K).GT.LAMMAX) THEN
                     LAMC(K) = LAMMAX
                     IF (IPRGAER.EQ.1) THEN
            
                        DUM = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                             LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                        !bloss(2018-02): prevent creation of number
                        NC3D(K) = MIN(DUM, DUM1+DUM2)
                        NAD3D(K) = MAX(0., DUM1+DUM2 -NC3D(K) )
                        ! RE-CALCULATE LAMC
                        LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                             (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
!bloss(2018-02)                        NAD3D(K) = MAX(0.,NAD3D(K)-(DUM-NC3D(K)))
!bloss(2018-02)                        NC3D(K) = DUM
                     ELSE
                        NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                             LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                     END IF !IPRGAER.EQ.1
                  END IF !LAMC.LT.LAMMIN
                  
                  !BRNR Tracking for limiters
                  DUM = (NC3D(K)-DUM1)/DT
                  IF (DUM.LT.0.) THEN
                     NCNEGLIM(K) = NCNEGLIM(K) + DUM
                     NADPOSLIM(K) = NADPOSLIM(K) + (NAD3D(K)-DUM2)/DT
                  ELSE
                     NCPOSLIM(K) = NCPOSLIM(K) + DUM
                     NADNEGLIM(K) = NADNEGLIM(K) + (NAD3D(K)-DUM2)/DT
                  END IF
               
                END IF !QC3D.LT.QSMALL

!!$                IF(NAD3D(K).GT.NAD3D(1)) write(*,*) 'Line 1995: NAD3D(',k,') = ', NAD3D(k)

!......................................................................
! SNOW

               IF (QNI3D(K).GE.QSMALL) THEN
                  LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
                  N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

                  IF (LAMS(K).LT.LAMMINS) THEN
                     LAMS(K) = LAMMINS
                     N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

                     NS3D(K) = N0S(K)/LAMS(K)

                  ELSE IF (LAMS(K).GT.LAMMAXS) THEN

                     LAMS(K) = LAMMAXS
                     N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

                     NS3D(K) = N0S(K)/LAMS(K)
                  END IF
               END IF

!......................................................................
! GRAUPEL

               IF (QG3D(K).GE.QSMALL) THEN
                  LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
                  N0G(K) = NG3D(K)*LAMG(K)

! ADJUST VARS

                  IF (LAMG(K).LT.LAMMING) THEN
                     LAMG(K) = LAMMING
                     N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2
                     
                     NG3D(K) = N0G(K)/LAMG(K)
                     
                  ELSE IF (LAMG(K).GT.LAMMAXG) THEN

                     LAMG(K) = LAMMAXG
                     N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2
                     
                     NG3D(K) = N0G(K)/LAMG(K)
                  END IF
               END IF

!.....................................................................
! ZERO OUT PROCESS RATES

               NARG1(K) = 0. 
               NARG2(K) = 0. 
               NACTRATE(K) = 0.
               QACTRATE(K) = 0. 
               NACTDIFF(K) = 0. 
               NATRANSFER(K) = 0.
               QATRANSFER(K) = 0.
               DC1(K) = 0.
               DC2(K) = 0.
               DG1(K) = 0.
               DG2(K) = 0.
               ISACT(K) = 0.
               SSPK(K) = 0.
               QAPRA(K) = 0.
               QAPRC(K) = 0.
               QAPRE(K) = 0.
               QASUBC(K) = 0.
               PRC(K) = 0.
               NPRC(K) = 0.
               NPRC1(K) = 0.
               PRA(K) = 0.
               NPRA(K) = 0.
               NRAGG(K) = 0.
               PSMLT(K) = 0.
               NSMLTS(K) = 0.
               NSMLTR(K) = 0.
               EVPMS(K) = 0.
               PCC(K) = 0.
               PRE(K) = 0.
               NSUBC(K) = 0.
               NSUBR(K) = 0.
               PRACG(K) = 0.
               NPRACG(K) = 0.
               PSMLT(K) = 0.
               EVPMS(K) = 0.
               PGMLT(K) = 0.
               EVPMG(K) = 0.
               PRACS(K) = 0.
               NPRACS(K) = 0.
               NGMLTG(K) = 0.
               NGMLTR(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

               IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                  IF (IRAIN.EQ.0) THEN

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

!     Put a fudge fuctor for dependence on horizontal resolution - Marat
	   if(do_scale_dependence_of_autoconv) then
	      PRC(K) = PRC(K)*min(1.,100./dx)**2
           end if

! note: nprc1 is change in Nr,
! nprc is change in Nc

                     NPRC1(K) = PRC(K)/CONS29
                     NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

                  ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! replace with seifert and beheng

                     dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
                     dum1 = 600.*dum**0.68*(1.-dum**0.68)**3
                     
                     prc(k) = 9.44e9/(20.*2.6e-7)* &
                          (nu(k)+2.)*(nu(k)+4.)/(nu(k)+1.)**2* &
                          (rho(k)*qc3d(k)/1000.)**4/(rho(k)*nc3d(k)/1.e6)**2* &
                          (1.+dum1/(1.-dum)**2)*1000./rho(k)

                     nprc(k) = prc(k)*2./2.6e-7*1000.
                     nprc1(k) = 0.5*nprc(k)
                     
                   END IF
                 END IF
                 
                 IF (IPRECOFF.EQ.1) THEN
                   NPRC(K) = 0.
                   NPRC1(K) = 0.
                   PRC(K) = 0.
                 END IF
!.......................................................................
! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING

                 IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

                   UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
                   UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
                   dum=(rhosu/rho(k))**0.54
                   UMS=MIN(UMS,1.2*dum)
                   UNS=MIN(UNS,1.2*dum)
                   UMR=MIN(UMR,9.1*dum)
                   UNR=MIN(UNR,9.1*dum)

! hm fix, 3/4/13
! for above freezing conditions to get accelerated melting of snow,
! we need collection of rain by snow (following Lin et al. 1983)
!            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
!                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
!                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
!                  (5./(LAMS(K)**3*LAMR(K))+                    &
!                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
!                  0.5/(LAMS(K)*LAMR(K)**3)))

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

! v3 5/27/11 npracs no longer used
!            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
!                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
!                (1./(LAMR(K)**3*LAMS(K))+                      &
!                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
!                 1./(LAMR(K)*LAMS(K)**3))

                 END IF

! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
! ASSUME SHED DROPS ARE 1 MM IN SIZE

                 IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

                   UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
                   UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
                   dum=(rhosu/rho(k))**0.54
                   UMG=MIN(UMG,20.*dum)
                   UNG=MIN(UNG,20.*dum)
                   UMR=MIN(UMR,9.1*dum)
                   UNR=MIN(UNR,9.1*dum)
                   
! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

! ASSUME 1 MM DROPS ARE SHED, GET NUMBER CONC (KG-1) SHED PER SEC

                   DUM = PRACG(K)/5.2E-7

! GET NUMBER CONC OF RAIN DROPS COLLECTED

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! hm 7/15/13, remove limit so that the number of collected drops can smaller than
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
            NPRACG(K)=NPRACG(K)-DUM

                 END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

                 IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

                   IF (IRAIN.EQ.0) THEN
                     
                     DUM=(QC3D(K)*QR3D(K))
                     PRA(K) = 67.*(DUM)**1.15
                     NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))
                     
                   ELSE IF (IRAIN.EQ.1) THEN
              
! v1.4
! seifert and beheng (2001) formulation

                     dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
                     dum1 = (dum/(dum+5.e-4))**4
                     pra(k) = 5.78e3*rho(k)/1000.*qc3d(k)*qr3d(k)*dum1
                     npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
                          (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)

                   END IF
                 END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

! v1.4, replace with seifert and beheng (2001)

                 IF (QR3D(K).GE.1.E-8) THEN
! v1.4
! seifert and beheng
! include breakup, V2.1
                   dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
                     dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
                   end if
            nragg(k) = -5.78*dum*qr3d(k)*nr3d(k)*rho(k)

                 END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF RAIN

                 IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                          F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
                 ELSE
                   EPSR = 0.
                 END IF

! NO CONDENSATION ONTO RAIN, ONLY EVAP

               IF (QV3D(K).LT.QVS(K)) THEN
                  PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
                  PRE(K) = MIN(PRE(K),0.)
               ELSE
                  PRE(K) = 0.
               END IF

!.......................................................................
! MELTING OF SNOW

! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

               IF (QNI3D(K).GE.1.E-8) THEN

! v3 5/27/11 bug fix
!             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
                 DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACS(K)

! hm fix 1/20/15
!             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(TMELT-T3D(K))/       &
!                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
!                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS10/                   &
!                   (LAMS(K)**CONS35))+DUM
                 PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(TMELT-T3D(K))/       &
                    XLF(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                      F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                      SC(K)**(1./3.)*CONS10/                   &
                      (LAMS(K)**CONS35))+DUM

! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

                 IF (QVQVS(K).LT.1.) THEN
                   EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                        (F1S/(LAMS(K)*LAMS(K))+                       &
                        F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                        SC(K)**(1./3.)*CONS10/                   &
                        (LAMS(K)**CONS35))
! bug fix V1.4
                   EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)    
                   EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
                   PSMLT(K) = PSMLT(K)-EVPMS(K)
                 END IF
               END IF

!.......................................................................
! MELTING OF GRAUPEL

! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

               IF (QG3D(K).GE.1.E-8) THEN

! v3 5/27/11 bug fix
!             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
                 DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACG(K)

! hm fix 1/20/15
!             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(TMELT-T3D(K))/ 		 &
!                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
!                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS11/                   &
!                   (LAMG(K)**CONS36))+DUM
             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(TMELT-T3D(K))/ 		 &
                    XLF(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                      F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                      SC(K)**(1./3.)*CONS11/                   &
                      (LAMG(K)**CONS36))+DUM

! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

                 IF (QVQVS(K).LT.1.) THEN
                   EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                        (F1S/(LAMG(K)*LAMG(K))+                               &
                        F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                        SC(K)**(1./3.)*CONS11/                   &
                        (LAMG(K)**CONS36))
! bug fix V1.4
                   EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
                   EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
                   PGMLT(K) = PGMLT(K)-EVPMG(K)
                 END IF
               END IF

! HM, V2.1
! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

               PRACG(K) = 0.
               PRACS(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
! CALCULATION

! CONSERVATION OF QC

               DUM = (PRC(K)+PRA(K))*DT

               IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN

                 RATIO = QC3D(K)/DUM
                 
                 PRC(K) = PRC(K)*RATIO
                 PRA(K) = PRA(K)*RATIO
                  
! RECALC NUMBER TENDENCIES BASED ON NEW PRC AND PRA

! CLOUD TENDENCIES
                 IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                   IF (IRAIN.EQ.0) THEN

! note: nprc1 is change in Nr,
! nprc is change in Nc

                     NPRC1(K) = PRC(K)/CONS29
                     NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))
                     
                     NPRC1(K) = MIN(NPRC1(K),NC3D(K)/DT)

                     IF (NPRC1(K).GT.NPRC(K)) NPRC1(K) = NPRC(K)
                     
                   ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! replace with seifert and beheng

                     nprc(k) = prc(k)*2./2.6e-7*1000.
                     nprc1(k) = 0.5*nprc(k)
                     
                   END IF
                 END IF
                 
                 IF (IPRECOFF.EQ.1) THEN
                   NPRC(K) = 0.
                   NPRC1(K) = 0.
                   PRC(K) = 0.
                 END IF

! RAIN TENDENCIES
                 IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

                   IF (IRAIN.EQ.0) THEN
                     
                     NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))
                     
                   ELSE IF (IRAIN.EQ.1) THEN
                     
! v1.4
! seifert and beheng (2001) formulation

                     npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
                          (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)
                     
                   END IF
                 END IF
                 
               END IF


               !IF(QC3D(K).GT.0.) THEN
               !   DUM = MIN(QAW3D(K)*((PRC(K)+PRA(K))/QC3D(K)),QAW3D(K)/DT)
               !   QAW3DTEN(K) = QAW3DTEN(K)-DUM
               !   QAR3DTEN(K) = QAR3DTEN(K)+DUM
               !END IF


! CONSERVATION OF SNOW

               DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

               IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

! NO SOURCE TERMS FOR SNOW AT T > FREEZING
                  RATIO = QNI3D(K)/DUM

                  PSMLT(K) = PSMLT(K)*RATIO
                  EVPMS(K) = EVPMS(K)*RATIO
                  PRACS(K) = PRACS(K)*RATIO

               END IF

! CONSERVATION OF GRAUPEL

               DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT
               
               IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
                  RATIO = QG3D(K)/DUM
                  
                  PGMLT(K) = PGMLT(K)*RATIO
                  EVPMG(K) = EVPMG(K)*RATIO
                  PRACG(K) = PRACG(K)*RATIO

               END IF

! CONSERVATION OF QR
! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

               DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT
               
               IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN
                  
                  RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                       (-PRE(K))
                  PRE(K) = PRE(K)*RATIO
        
               END IF

!....................................

               QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

               T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)

               QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
               QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
               QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
               QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
               
               IF((IPRGAER.EQ.1).AND.(IPRECOFF.EQ.0)) THEN
                  IF(QC3D(K).GE.QSMALL) THEN
                     QAPRA(K) = QAW3D(K)*PRA(K)/QC3D(K)
                     QAPRC(K) = QAW3D(K)*PRC(K)/QC3D(K)
                     DUM = MIN((PRA(K)+PRC(K))*DT/QC3D(K),1.)
                  ELSE
                     DUM = 0.
                     QAPRA(K) = 0.
                     QAPRC(K) = 0.
                  END IF
                  QAW3DTEN(K) = QAW3DTEN(K)-QAPRA(K)-QAPRC(K)
                  
                  IF(QR3D(K).GE.QSMALL) THEN
                     DUM1 = MAX(PRE(K)*DT/QR3D(K),-1.) !Revisit tendency due to evap
                     QAPRE(K) = QAR3D(K)*DUM1/DT
                  ELSE
                     QAPRE(K) = 0.
                     DUM1 = 0.
                  END IF

                  QAR3DTEN(K) = QAR3DTEN(K)+QAR3D(K)*(DUM1/DT)+QAW3D(K)*(DUM/DT)
                  QAD3DTEN(K) = QAD3DTEN(K)+QAR3D(K)*(-DUM1/DT)
               END IF

!add declarations for QAPRA, QAPRC, QAPRE 

!DBG var updates
!               DBGIND = 3
!               DBGQAR3DTEN(DBGIND,K) = QAR3DTEN(K)
!               DBGQAR3D(DBGIND,K) = DBGQAR3D(2,K) + DBGQAR3DTEN(DBGIND,K)*DT
!
!               DBGQAW3DTEN(DBGIND,K) = QAW3DTEN(K)
!               DBGQAW3D(DBGIND,K) = DBGQAW3D(1,K) + DBGQAW3DTEN(DBGIND,K)*DT
!               DBGQAD3DTEN(DBGIND,K) = QAD3DTEN(K)
!               DBGQAD3D(DBGIND,K) = DBGQAD3D(1,K) + DBGQAD3DTEN(DBGIND,K)*DT
!
!               IF((MINVAL(DBGQAR3D).LT.0)) THEN
!                   QNEG = .TRUE.
!                  print*,'error, mass is negative!!! DBGIND = 3'
!               END IF


! v3 5/27/11
!      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
!      NG3DTEN(K) = NG3DTEN(K)
               NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
               NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))

               IF (PRE(K).LT.0.) THEN
                 DUM = PRE(K)*DT/QR3D(K)
                 DUM = MAX(-1.,DUM)
                 NSUBR(K) = DUM*NR3D(K)/DT
               END IF

! V1.3 move code below to before saturation adjustment
               IF (EVPMS(K)+PSMLT(K).LT.0.) THEN
                 DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
                 DUM = MAX(-1.,DUM)
                 NSMLTS(K) = DUM*NS3D(K)/DT
               END IF
               IF (PSMLT(K).LT.0.) THEN
                 DUM = PSMLT(K)*DT/QNI3D(K)
                 DUM = MAX(-1.0,DUM)
                 NSMLTR(K) = DUM*NS3D(K)/DT
               END IF
               IF (EVPMG(K)+PGMLT(K).LT.0.) THEN
                 DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
                 DUM = MAX(-1.,DUM)
                 NGMLTG(K) = DUM*NG3D(K)/DT
               END IF
               IF (PGMLT(K).LT.0.) THEN
                 DUM = PGMLT(K)*DT/QG3D(K)
                 DUM = MAX(-1.0,DUM)
                 NGMLTR(K) = DUM*NG3D(K)/DT
               END IF

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

               NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
               NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
               NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))
               NAD3DTEN(K) = NAD3DTEN(K)-NSUBR(K) !brnr dry aerosol tendency due to evaporation or rain

300            CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               IF(ISATADJ.EQ.0) THEN !PB 4/13/09

! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = 0.622*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

                  DUMS = DUMQV-DUMQSS
                  PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
                  IF (PCC(K)*DT+DUMQC.LT.0.) THEN
                     PCC(K) = -DUMQC/DT
                  END IF

                  QV3DTEN(K) = QV3DTEN(K)-PCC(K)
                  T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
                  QC3DTEN(K) = QC3DTEN(K)+PCC(K)

               END IF

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS

!brnr modify below to also calculate aerosol mass in activated CCN
!following Abdul-Razzak et al Pt. 2

!bloss: only do activation if droplet number is predicted
!bloss      IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) THEN

           
! SET UP AEROSOL MASS AND NUMBER FIELDS IF USING PROGNOSTIC OPTION !brnr
            
!               DBGNCTND(1) = DBGNCTND(1) + NC3DTEN(K)
!               DBGNADTND(1) = DBGNADTND(1) + NAD3DTEN(K)
!               DBGNC(1) = DBGNC(1) + NC3D(K)
!               DBGNAD(1) = DBGNAD(1) + NAD3D(K)

               IF (IPRGAER.EQ.1) THEN 
                  NANEW1 = (NAD3D(K) + NC3D(K))*RHO(K) !Accum
                  NANEW2 = NAD2_3D(K) * RHO(K) ! Aitken 
                  MAER1 = (QAD3D(K) + QAW3D(K))*RHO(K)
                  MAER2 = QAD2_3D(K) * RHO(K)

                  IF (NANEW1.GT.0.) THEN
                     RM1 = (((0.75*MAER1/(PI*RHOA*NANEW1))*EXP(-1.*9./2*LOG(SIG1)**2)))**(1/3.) ! calculate geometric mean radius (mode 1)
                  ELSE
                     RM1 = 0.
                     MAER1 = 0.
                  END IF
                  IF (NANEW2.GT.0.) THEN
                     RM2 = (((0.75*MAER2/(PI*RHOA*NANEW2))*EXP(-1.*9./2*LOG(SIG2)**2)))**(1/3.) ! calculate geometric mean radius (mode 2)
                  ELSE
                     RM2 = 0.
                     MAER2 = 0.
                  END IF
                  DG1(K) = RM1*2.
                  DG2(K) = RM2*2.
               END IF
               
               IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL.AND.INUM.EQ.0) THEN
!                IF (QC3D(K)+QC3DTEN(K)*DT.GE.1.E-6.AND.INUM.EQ.0) THEN
! EFFECTIVE VERTICAL VELOCITY (M/S)

                  IF (ISUB.EQ.0) THEN
! ADD SUB-GRID VERTICAL VELOCITY
                     DUM = W3D(K)+WVAR(K)

! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
                     DUM = MAX(DUM,0.10)

                  ELSE IF (ISUB.EQ.1) THEN 
                    DUM=W3D(K)
         if(do_scale_dependence_of_activation) then
            ! Marat:take into account dependence of W on hor gridsize
            DUM=W3D(K)*max(1.,dx/1000.)
         end if
                  END IF

! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
                 IF (DUM.GE.0.001) THEN

                     IF (IBASE.EQ.1) THEN

! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

                        IDROP=0

! V1.3 USE CURRENT VALUE OF QC FOR IDROP
                        IF (QC3D(K).LE.0.05E-3/RHO(K)) THEN
                           IDROP=1
                        END IF
                        IF (K.EQ.1) THEN
                           IDROP=1
                        ELSE IF (K.GE.2) THEN
                           IF (QC3D(K).GT.0.05E-3/RHO(K).AND. &
                                QC3D(K-1).LE.0.05E-3/RHO(K-1)) THEN
                              IDROP=1
                          END IF
                        END IF

                        IF (IDROP.EQ.1) THEN
! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

                           IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

                              DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
                              DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
                              DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
                              DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                              DUM2 = (DUM2-NC3D(K))/DT
                              DUM2 = MAX(0.,DUM2)
                              NARG1(K) = DUM2 !brnr track activation tendency for NC
                              NC3DTEN(K) = NC3DTEN(K)+DUM2

                           ELSE IF (IACT.EQ.2) THEN
! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

                              IF ((NANEW1.GT.0.AND.RM1.GT.0).OR.(NANEW2.GT.0.AND.RM2.GT.0)) THEN
                              
                                 SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                                 AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
                                 ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
                                 GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))
                                 
                                 GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                      (T3D(K)*RR)-1.))
                                 
                                 PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT
                                 
                                 IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                                    ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
                                    SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                                    DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
                                 ELSE
                                    DUM1 = 0.
                                 END IF

                                 IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                                    ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)
                                    SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                                    DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)
                                 ELSE 
                                    DUM2 = 0.
                                 END IF
 
                                D3 = 0.                                
                                 IF (dofixedcoarsemode) THEN                                   
                                    ETA3 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*Ncoarse*RHO(K)*1.e6)  ! convert Ncoarse from /mg to /m^3
                                    SM3 = 2./BACT_coarse**0.5*(AACT/(3.*RM3))**1.5
                                    D3 = 1./SM3**2*(F13*(PSI/ETA3)**1.5+F23*(SM3**2/(ETA3+3.*PSI))**0.75)
                                 END IF
                                    
                                 SMAX = 1./(DUM1+DUM2+D3)**0.5

                                 IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                                    UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                                    DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                                 ELSE
                                    UU1 = 0.
                                 END IF
                            
                                 IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                                    UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                                    DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                                 ELSE
                                    UU2 = 0.
                                 END IF

                                 ISACT(K) = 1.
                                 SSPK(K) = SMAX ! save statistic
                                 DC1(K) = DG1(K) * (SM1/SMAX)**(2./3.) !Critical diameters
                                 DC2(K) = DG2(K) * (SM2/SMAX)**(2./3.)

                                 NA1 = NANEW1/RHO(K)  ! units of /kg  total number of mode 1
                                 QA1 = MAER1/RHO(K)   ! units of kg/kg total mass of mode 1
                                 
                                 NARG1(K) = DUM1/RHO(K) ! number of prescribed activated in /kg
                                 NARG2(K) = DUM2/RHO(K)
                                 
                                 NARG1(K) = MIN(NARG1(K), NA1) ! is this really necessary?
                                 NARG2(K) = MIN(NARG2(K), NANEW2/RHO(K))

                                 NACTDIFF(K) = (NARG1(K) + NARG2(K) - NC3D(K)) ! saved stat

                                 NCNEW = NC3D(K) + MAX(0., NACTDIFF(K)) ! only increase nc 
                                 
				 NCNEW = MIN(NA1, NCNEW)  ! cap activation could make this optional

                                 NTRANS = 0.
                                 QTRANS = 0.
                                 IF (doacttransfer) THEN
                                    call hoppel_aitken_accum_transfer(NAD2_3D(K), NA1, QAD2_3D(K), QA1, &
                                         SG2, SG1, DC2(K), DC1(K), RHOA, RHOA, NTRANS, QTRANS)
                                 END IF   

                                 ! transfer to accum from aitken
                                 NA1 = NA1 + NTRANS
                                 QA1 = QA1 + QTRANS

                                 ! create tendencies of NAD, NC, QAD, QAW,  based on new values
                                 ! disregarding old tendencies?
                                 NADNEW = NA1 - NCNEW
                                 NAD3DTEN(K) = NAD3DTEN(K) + (NADNEW - NAD3D(K))/DT
                                 NACTRATE(K) = (NCNEW - NC3D(K))/DT
                                 NC3DTEN(K) = NC3DTEN(K) + NACTRATE(K)                                 

                                 ! activated mass diagnosed
                                 QAWNEW = QA1 * (1. - mass_fraction(NADNEW/NA1, sigma_accum))
                                 QAW3DTEN(K) = QAW3DTEN(K) + (QAWNEW - QAW3D(K))/DT
                                 QAD3DTEN(K) = QAD3DTEN(K) + (QA1 - QAWNEW - QAD3D(K))/DT ! unactivated is what is leftover
                                 
                                 NAD2_3DTEN(K) = NAD2_3DTEN(K) - NTRANS/DT
                                 QAD2_3DTEN(K) = QAD2_3DTEN(K) - QTRANS/DT

                                 NATRANSFER(K) = NTRANS/DT  ! save stat
                                 QATRANSFER(K) = QTRANS/DT  ! save stat                             
                                 
                              END IF ! (NANEW1.GT.0).OR.(NANEW2.GT.0)
                           END IF  ! IACT
                       
!.............................................................................
                        ELSE IF (IDROP.EQ.0) THEN
! ACTIVATE IN CLOUD INTERIOR
! FIND EQUILIBRIUM SUPERSATURATION

                           TAUC=1./(2.*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.)/LAMC(K))
                           IF (EPSR.GT.1.E-8) THEN
                              TAUR=1./EPSR
                           ELSE
                              TAUR=1.E8
                           END IF

! hm fix 1/20/15
!           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=(-QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
                           DUM3=DUM3*TAUC*TAUR/(TAUC+TAUR)
                           
                           IF (DUM3/QVS(K).GE.1.E-6) THEN
                              IF (IACT.EQ.1) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

                                 DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
                                 DUMACT = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))

! USE POWER LAW CCN SPECTRA

! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
                                 DUM3=DUM3/QVS(K)*100.

                                 DUM2=C1*DUM3**K1
! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
                                 DUM2=MIN(DUM2,DUMACT)
                                 DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
                                 DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                                 DUM2 = (DUM2-NC3D(K))/DT
                                 DUM2 = MAX(0.,DUM2)
                                 NARG1(K) = DUM2 !brnr track activation number tendency
                                 NC3DTEN(K) = NC3DTEN(K)+DUM2

                              ELSE IF (IACT.EQ.2) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

                                 IF ((NANEW1.GT.0.AND.RM1.GT.0).OR.(NANEW2.GT.0.AND.RM2.GT.0)) THEN
                                 
                                    SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                                    AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
                                    ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
                                    GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))
                                    
                                    GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                         (T3D(K)*RR)-1.))
                                    
                                    PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT
                                    
                                    IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                                       ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
                                       SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                                       DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
                                    ELSE
                                       DUM1 = 0.
                                       SM1 = 1.
                                    END IF

                                    IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                                       ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)
                                       SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                                       DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)
                                    ELSE 
                                       DUM2 = 0.
                                       SM2 = 1.
                                    END IF

                                    D3 = 0.                                
                                    IF (dofixedcoarsemode) THEN                                   
                                       ETA3 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*Ncoarse*RHO(K)*1.e6)  ! convert Ncoarse from /mg to /m^3
                                       SM3 = 2./BACT_coarse**0.5*(AACT/(3.*RM3))**1.5
                                       D3 = 1./SM3**2*(F13*(PSI/ETA3)**1.5+F23*(SM3**2/(ETA3+3.*PSI))**0.75)
                                    END IF
                                    SMAX = 1./(DUM1+DUM2+D3)**0.5

                                    
                                    UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                                    UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                                    DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                                    DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                                    
                                    DUM1 = DUM1/RHO(K)
                                    DUM2 = DUM2/RHO(K)    !CONVERT TO KG-1

                                    SSPK(K) = SMAX
                                    
! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL
                                    DUMACT1 = MIN(DUM1, NANEW1/RHO(K))
                                    DUMACT2 = MIN(DUM2, NANEW2/RHO(K))
!                                   DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! AEROSOL MASS FOR EQUILIBRIUM SUPER SATURATION USING PROGNOSTIC AEROSOL
                              

! USE LOGNORMAL AEROSOL
                                    SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                                    AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)

! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
                                    SMAX = DUM3/QVS(K)

                                    
                                    IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                                       SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                                       UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                                       DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                                    ELSE
                                       DUM1 = 0.
                                       UU1 = 0.
                                    END IF
                                 
                                    IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                                       SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                                       UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                                       DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                                    ELSE
                                       DUM2 = 0.
                                       UU2 = 0.
                                    END IF
                                    
                                    DUM1 = DUM1/RHO(K)
                                    DUM2 = DUM2/RHO(K)  !CONVERT TO KG-1
                                 
                                    
! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL
! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
                                    DUM1 = MIN(DUM1, NANEW1/RHO(K), DUMACT1)
                                    DUM2 = MIN(DUM2,NANEW2/RHO(K), DUMACT2)

                                    ISACT(K) = 1.
                                    SSPK(K) = MIN(SSPK(K), SMAX)
                                    DC1(K) = DG1(K) * (SM1/SSPK(K))**(2./3.) !Critical diameters
                                    DC2(K) = DG2(K) * (SM2/SSPK(K))**(2./3.)

                                    NA1 = NANEW1/RHO(K)
                                    QA1 = MAER1/RHO(K)

                                    NARG1(K) = DUM1
                                    NARG2(K) = DUM2

                                    NARG1(K) = MIN(NARG1(K), NA1) ! is this really necessary?
                                    NARG2(K) = MIN(NARG2(K), NANEW2/RHO(K))
                                    
                                    NACTDIFF(K) = (NARG1(K) + NARG2(K) - NC3D(K)) ! saved stat
                                    
                                    NCNEW = NC3D(K) + MAX(0., NACTDIFF(K)) ! only increase nc 
                                 
                                    NCNEW = MIN(NA1, NCNEW)  ! cap activation could make this optional

                                    NTRANS = 0.
                                    QTRANS = 0.
                                    IF (doacttransfer) THEN
                                       call hoppel_aitken_accum_transfer(NAD2_3D(K), NA1, QAD2_3D(K), QA1, &
                                            SG2, SG1, DC2(K), DC1(K), RHOA, RHOA, NTRANS, QTRANS)
                                    END IF

                                    ! transfer to accum from aitken
                                    NA1 = NA1 + NTRANS
                                    QA1 = QA1 + QTRANS

                                    ! create tendencies of NAD, NC, QAD, QAW,  based on new values

                                    NADNEW = NA1 - NCNEW
                                    NACTRATE(K) = (NCNEW - NC3D(K))/DT

                                    ! activated mass diagnosed
                                    QAWNEW = QA1 * (1. - mass_fraction(NADNEW/NA1, sigma_accum))
                                    QACTRATE(K) = (QAWNEW - QAW3D(K))/DT
                                    QACTRATE(K) = min(QACTRATE(K), (QAD3D(K) + QATRANSFER(K)*DT)/DT)

                                    NATRANSFER(K) = NTRANS/DT  
                                    QATRANSFER(K) = QTRANS/DT              

                                    NAD2_3DTEN(K) = NAD2_3DTEN(K) - NATRANSFER(K)
                                    QAD2_3DTEN(K) = QAD2_3DTEN(K) - QATRANSFER(K)

                                    NAD3DTEN(K) = NAD3DTEN(K) + NATRANSFER(K) - NACTRATE(K)
                                    QAD3DTEN(K) = QAD3DTEN(K) + QATRANSFER(K) - QACTRATE(K)

                                    NC3DTEN(K) = NC3DTEN(K) + NACTRATE(K)
                                    QAW3DTEN(K) = QAW3DTEN(K) + QACTRATE(K)
                                 
                                 END IF ! (NANEW1.GT.0).OR.(NANEW2.GT.0)
                              END IF ! IACT
                           END IF ! DUM3/QVS > 1.E-6
                        END IF  ! IDROP = 1

!.......................................................................
                     ELSE IF (IBASE.EQ.2) THEN

                        IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

                           DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
                           DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
                           DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
                           DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                           DUM2 = (DUM2-NC3D(K))/DT
                           DUM2 = MAX(0.,DUM2)
                           NARG1(K) = DUM2 !brnr track activated number tendency
                           NC3DTEN(K) = NC3DTEN(K)+DUM2
                           
                        ELSE IF (IACT.EQ.2) THEN
                           
                           IF ((NANEW1.GT.0.AND.RM1.GT.0).OR.(NANEW2.GT.0.AND.RM2.GT.0)) THEN
                              
                              SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                              AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
                              ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
                              GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))
                              
                              GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                   (T3D(K)*RR)-1.))
                              
                              PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT
                              
                              IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                                 ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
                                 SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                                 DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
                              ELSE                                 
                                 DUM1 = 0.
                              END IF
                              
                              IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                                 ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)
                                 SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                                 DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)
                              ELSE 
                                 DUM2 = 0.
                              END IF

                              D3 = 0.                                
                              IF (dofixedcoarsemode) THEN                                   
                                 ETA3 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*Ncoarse*RHO(K)*1.e6)  ! convert Ncoarse from /mg to /m^3
                                 SM3 = 2./BACT_coarse**0.5*(AACT/(3.*RM3))**1.5
                                 D3 = 1./SM3**2*(F13*(PSI/ETA3)**1.5+F23*(SM3**2/(ETA3+3.*PSI))**0.75)
                              END IF
                              SMAX = 1./(DUM1+DUM2+D3)**0.5
                              
                              IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                                 UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                                 DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                              ELSE
                                 UU1 = 0.
                              END IF
                              
                              IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                                 UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                                 DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                              ELSE
                                 UU2 = 0.
                              END IF

                              ISACT(K) = 1.
                              SSPK(K) = SMAX
                              DC1(K) = DG1(K) * (SM1/SMAX)**(2./3.) !Critical diameters
                              DC2(K) = DG2(K) * (SM1/SMAX)**(2./3.)

                              NA1 = NANEW1/RHO(K)  ! units of /kg  total number of mode 1
                              QA1 = MAER1/RHO(K)   ! units of kg/kg total mass of mode 1

                              NARG1(K) = DUM1/RHO(K) ! number of prescribed activated in /kg
                              NARG2(K) = DUM2/RHO(K)

                              NARG1(K) = MIN(NARG1(K), NA1) ! is this really necessary?
                              NARG2(K) = MIN(NARG2(K), NANEW2/RHO(K))

                              NACTDIFF(K) = (NARG1(K) + NARG2(K) - NC3D(K)) ! saved stat

                              NCNEW = NC3D(K) + MAX(0., NACTDIFF(K)) ! only increase nc 

                              NCNEW = MIN(NA1, NCNEW)  ! cap activation could make this optional

                              NTRANS = 0.
                              QTRANS = 0.
                              IF (doacttransfer) THEN
                                 call hoppel_aitken_accum_transfer(NAD2_3D(K), NA1, QAD2_3D(K), QA1, &
                                      SG2, SG1, DC2(K), DC1(K), RHOA, RHOA, NTRANS, QTRANS)
                              END IF

                              ! transfer to accum from aitken
                              NA1 = NA1 + NTRANS
                              QA1 = QA1 + QTRANS

                              NADNEW = NA1 - NCNEW
                              NAD3DTEN(K) = NAD3DTEN(K) + (NADNEW - NAD3D(K))/DT
                              NACTRATE(K) = (NCNEW - NC3D(K))/DT
                              NC3DTEN(K) = NC3DTEN(K) + NACTRATE(K)                                 

                              QAWNEW = QA1 * (1. - mass_fraction(NADNEW/NA1, sigma_accum))
                              QAW3DTEN(K) = QAW3DTEN(K) + (QAWNEW - QAW3D(K))/DT
                              QAD3DTEN(K) = QAD3DTEN(K) + (QA1 - QAWNEW - QAD3D(K))/DT ! unactivated is what is leftover
                                 
                              NAD2_3DTEN(K) = NAD2_3DTEN(K) - NTRANS/DT
                              QAD2_3DTEN(K) = QAD2_3DTEN(K) - QTRANS/DT
                              NATRANSFER(K) = NTRANS/DT  ! save stat
                              QATRANSFER(K) = QTRANS/DT  ! save stat            

                           END IF ! (NANEW1.GT.0).OR.(NANEW2.GT.0)
                        END IF  ! IACT
                     END IF  ! IBASE
                  END IF  ! W > 0.001
               END IF  ! QC3D > QSMALL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION
!brnr: make the Nc tendency due to evaporation switchable from prm file
    
            IF (PCC(K).LT.0.) THEN
               IF (IEVPNC.EQ.1) THEN
                  DUM = PCC(K)*DT/QC3D(K)
                  !DUM = MAX(-1.,DUM)
                  NSUBC(K) = DUM*NC3D(K)/DT
                  NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
                  IF (IPRGAER.EQ.1) THEN
                     NAD3DTEN(K) = NAD3DTEN(K)-NSUBC(K) !brnr dry aerosol tendency is opposite of cloud droplet tendency
                     DUM1 = DUM*QAW3D(K)/DT
                     IF((QAW3D(K)+(QAW3DTEN(K)+DUM1)*DT).LT.0) THEN
                        QAD3DTEN(K) = QAD3DTEN(K)+QAW3DTEN(K)+QAW3D(K)/DT                       
                        QASUBC(K) = -QAW3DTEN(K)-QAW3D(K)/DT
                        QAW3DTEN(K) = -QAW3D(K)/DT 
                     ELSE
                        QAW3DTEN(K) = QAW3DTEN(K)+DUM1
                        QAD3DTEN(K) = QAD3DTEN(K)-DUM1
                        QASUBC(K) = DUM1
                     END IF
                  END IF
               ELSEIF ((QC3D(K)+QC3DTEN(K)*DT).LT.QSMALL) THEN !brnr if evaporation in a gridbox will essentially evaporate all 
                  NSUBC(K) = -NC3DTEN(K)-NC3D(K)/DT
                  IF (IPRGAER.EQ.1) THEN
                     NAD3DTEN(K) = NAD3DTEN(K)-NSUBC(K) !     and nucleated aerosol mass to dry aerosol mass
                     QASUBC(K) = -QAW3DTEN(K)-QAW3D(K)/DT                    
                     QAW3DTEN(K) = -QAW3D(K)/DT
                     QAD3DTEN(K) = QAD3DTEN(K)-QASUBC(K)
                  END IF
                  NC3DTEN(K) = -NC3D(K)/DT   !     cloud, transfer all cloud number back to dry aerosol number
               END IF
            END IF

           
!            DBGNCTND(3) = DBGNCTND(3) + NC3DTEN(K)
!            DBGNADTND(3) = DBGNADTND(3) + NAD3DTEN(K)
!            DBGNC(3) = DBGNC(3) + NC3D(K)
!            DBGNAD(3) = DBGNAD(3) + NAD3D(K)

!DBG var updates
!            DBGIND = 5
!            DBGQAR3DTEN(DBGIND,K) = QAR3DTEN(K)
!            DBGQAR3D(DBGIND,K) = DBGQAR3D(2,K) + DBGQAR3DTEN(DBGIND,K)*DT
!            
!            DBGQAW3DTEN(DBGIND,K) = QAW3DTEN(K)
!            DBGQAW3D(DBGIND,K) = DBGQAW3D(1,K) + DBGQAW3DTEN(DBGIND,K)*DT
!            
!            DBGQAD3DTEN(DBGIND,K) = QAD3DTEN(K)
!            DBGQAD3D(DBGIND,K) = DBGQAD3D(1,K) + DBGQAD3DTEN(DBGIND,K)*DT
!
!            IF((MINVAL(DBGQAR3D).LT.0)) THEN
!               QNEG = .TRUE.
!               print*,'error, mass is negative!!! DBGIND = 4'
!            END IF

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

!.....................................................................
!.....................................................................
          ELSE  ! TEMPERATURE < 273.15

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

            IF (INUM.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
              NC3D(K)=NDCNST*1.E6/RHO(K)
            END IF

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

            NI3D(K) = MAX(0.,NI3D(K))
            NS3D(K) = MAX(0.,NS3D(K))
            NC3D(K) = MAX(0.,NC3D(K))
            NR3D(K) = MAX(0.,NR3D(K))
            NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

            IF (QI3D(K).GE.QSMALL) THEN
               LAMI(K) = (CONS12*                 &
                    NI3D(K)/QI3D(K))**(1./DI)
               N0I(K) = NI3D(K)*LAMI(K)

! CHECK FOR SLOPE

! ADJUST VARS

               IF (LAMI(K).LT.LAMMINI) THEN
                  
                  LAMI(K) = LAMMINI
                  
                  N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12
                  
                  NI3D(K) = N0I(K)/LAMI(K)
               ELSE IF (LAMI(K).GT.LAMMAXI) THEN
                  LAMI(K) = LAMMAXI
                  N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12
                  
                  NI3D(K) = N0I(K)/LAMI(K)
               END IF
            END IF

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF

      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

            IF (QC3D(K).GE.QSMALL) THEN
               
         !bloss: option for fixing pgam
               if(dofix_pgam) then
                  pgam(k) = pgam_fixed
               else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
                  PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
                  PGAM(K)=1./(PGAM(K)**2)-1.
                  PGAM(K)=MAX(PGAM(K),2.)
                  PGAM(K)=MIN(PGAM(K),10.)

               end if
! v1.4
! interpolate
               dumii=int(pgam(k))
               nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
                    (pgam(k)-real(dumii))

! CALCULATE LAMC

               LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                    (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON
               
               LAMMIN = (PGAM(K)+1.)/60.E-6
               LAMMAX = (PGAM(K)+1.)/1.E-6

! Brnr modify the bound check routines to conserve Na+Nc
               DUM1 = NC3D(K)
               DUM2 = NAD3D(K)

               IF (LAMC(K).LT.LAMMIN) THEN
                  LAMC(K) = LAMMIN
                  IF (IPRGAER.EQ.1) THEN
                     DUM = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                          LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                     !DUM = MIN(DUM,(NC3D(K)+NAD3D(K)))
                        !bloss(2018-02): prevent creation of number
                        !bloss(2019-03): prevent problems if total cloud+aerosol number is zero or negative
                        NC3D(K) = MAX(0.1*DUM, MIN(DUM, DUM1+DUM2) )
                        NAD3D(K) = MAX(0., DUM1+DUM2 -NC3D(K) )
                        ! RE-CALCULATE LAMC
                        LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                             (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
!bloss(2018-02)                        NAD3D(K) = MAX(0.,NAD3D(K)-(DUM-NC3D(K)))
!bloss(2018-02)                        NC3D(K) = DUM
                  ELSE
                     NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                          LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                  END IF!IPRGAER.EQ.1
               ELSE IF (LAMC(K).GT.LAMMAX) THEN
                  LAMC(K) = LAMMAX
                  IF (IPRGAER.EQ.1) THEN
                     DUM = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                          LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                     !DUM = MIN(DUM,(NC3D(K)+NAD3D(K)))
                        !bloss(2018-02): prevent creation of number
                        NC3D(K) = MIN(DUM, DUM1+DUM2)
                        NAD3D(K) = MAX(0., DUM1+DUM2 -NC3D(K) )
                        ! RE-CALCULATE LAMC
                        LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                             (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
!bloss(2018-02)                        NAD3D(K) = MAX(0.,NAD3D(K)-(DUM-NC3D(K)))
!bloss(2018-02)                        NC3D(K) = DUM
                  ELSE
                     NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                          LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                     END IF!IPRGAER.EQ.1
               END IF

! TO CALCULATE DROPLET FREEZING

               CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)

            END IF

!!$                IF(NAD3D(K).GT.NAD3D(1)) write(*,*) 'Line 3281: NAD3D(',k,') = ', NAD3D(k)

!......................................................................
! SNOW

            IF (QNI3D(K).GE.QSMALL) THEN
               LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
               N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

               IF (LAMS(K).LT.LAMMINS) THEN
                  LAMS(K) = LAMMINS
                  N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

                  NS3D(K) = N0S(K)/LAMS(K)
                  
               ELSE IF (LAMS(K).GT.LAMMAXS) THEN
                  
                  LAMS(K) = LAMMAXS
                  N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1
                  
                  NS3D(K) = N0S(K)/LAMS(K)
               END IF
            END IF
            
!......................................................................
! GRAUPEL

            IF (QG3D(K).GE.QSMALL) THEN
               LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
               N0G(K) = NG3D(K)*LAMG(K)
               
! CHECK FOR SLOPE

! ADJUST VARS

               IF (LAMG(K).LT.LAMMING) THEN
                  LAMG(K) = LAMMING
                  N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

                  NG3D(K) = N0G(K)/LAMG(K)
                  
               ELSE IF (LAMG(K).GT.LAMMAXG) THEN
                  
                  LAMG(K) = LAMMAXG
                  N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2
                  
                  NG3D(K) = N0G(K)/LAMG(K)
               END IF
            END IF

!.....................................................................
! ZERO OUT PROCESS RATES

            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PGSACW(K) = 0.
            PGRACS(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES
! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
!.......................................................................
! FREEZING OF CLOUD DROPLETS
! ONLY ALLOWED BELOW -4 C
            IF (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15) THEN

! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

! MEYERS CURVE

               NACNT = EXP(-2.80+0.262*(TMELT-T3D(K)))*1000.

! COOPER CURVE
!        NACNT =  5.*EXP(0.304*(TMELT-T3D(K)))

! FLECTHER
!     NACNT = 0.01*EXP(0.6*(TMELT-T3D(K)))

! CONTACT FREEZING

! MEAN FREE PATH

               DUM = 7.37*T3D(K)/(288.*10.*PRES(K))/100.

! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
! BASED ON BROWNIAN DIFFUSION

               DAP(K) = CONS37*T3D(K)*(1.+DUM/RIN)/MU(K)
               
               MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                    LOG(GAMMA(PGAM(K)+5.))-4.*LOG(LAMC(K)))
               NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
                    GAMMA(PGAM(K)+2.)/                         &
                    LAMC(K)

! IMMERSION FREEZING (BIGG 1953)

!           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
!                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
!                   EXP(AIMM*(273.15-T3D(K)))

!           NNUCCC(K) = NNUCCC(K)+                                  &
!            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
!                *EXP(AIMM*(273.15-T3D(K)))

! hm 7/15/13 fix for consistency w/ original formula
               MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                    EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                   (EXP(AIMM*(273.15-T3D(K)))-1.)

               NNUCCC(K) = NNUCCC(K)+                                  &
                    CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                *(EXP(AIMM*(273.15-T3D(K)))-1.)

! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
! MIXING RATIO, SINCE CONSERVATION NOT APPLIED TO NUMBER CONC

               NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

            END IF
        
!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

            IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

               IF (IRAIN.EQ.0) THEN

!     Put a fudge fuctor for dependence on horizontal resolution - Marat
                PRC(K)=min(1.,100./dx)**2*1350.*QC3D(K)**2.47*  &
!                PRC(K)=1350.*QC3D(K)**2.47*  &
                       (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

                  NPRC1(K) = PRC(K)/CONS29
                  NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

! hm bug fix 3/20/12
                  NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))
                  
               ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! replace with seifert and beheng

                  dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
                  dum1 = 600.*dum**0.68*(1.-dum**0.68)**3
                  
                  prc(k) = 9.44e9/(20.*2.6e-7)* &
                       (nu(k)+2.)*(nu(k)+4.)/(nu(k)+1.)**2* &
                       (rho(k)*qc3d(k)/1000.)**4/(rho(k)*nc3d(k)/1.e6)**2* &
                       (1.+dum1/(1.-dum)**2)*1000./rho(k)

                  nprc(k) = prc(k)*2./2.6e-7*1000.
                  nprc1(k) = 0.5*nprc(k)

               END IF
            END IF
            
            IF (IPRECOFF.EQ.1) THEN
               NPRC(K) = 0.
               NPRC1(K) = 0.
               PRC(K) = 0.
            END IF
!.......................................................................
            ! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME


! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

            IF (QNI3D(K).GE.1.E-8) THEN
               NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
                    ((2.+BS)/3.)*QNI3D(K)**((2.+BS)/3.)*                  &
                    (NS3D(K)*RHO(K))**((4.-BS)/3.)/                       &
                    (RHO(K))
            END IF

!.......................................................................
! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
! HERE USE CONTINUOUS COLLECTION EQUATION WITH
! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

! SNOW

            IF (QNI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN
               
               PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                    N0S(K)/                        &
                    LAMS(K)**(BS+3.)
               NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                    N0S(K)/                        &
                    LAMS(K)**(BS+3.)

            END IF

!............................................................................
! COLLECTION OF CLOUD WATER BY GRAUPEL

            IF (QG3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN
               
               PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                    N0G(K)/                        &
                    LAMG(K)**(BG+3.)
               NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                    N0G(K)/                        &
                    LAMG(K)**(BG+3.)
	    END IF

!.......................................................................
! HM, ADD 12/13/06
! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
! BEFORE RIMING CAN OCCUR
! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
! TO HALLET-MOSSOP SPLINTERING

            IF (QI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
! FROM THOMPSON ET AL. 2004, MWR

               IF (1./LAMI(K).GE.100.E-6) THEN

                  PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                       N0I(K)/                        &
                       LAMI(K)**(BI+3.)
                  NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                       N0I(K)/                        &
                       LAMI(K)**(BI+3.)
               END IF
            END IF

!.......................................................................
! ACCRETION OF RAIN WATER BY SNOW
! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                (1./(LAMR(K)**3*LAMS(K))+                      &
                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
                 1./(LAMR(K)*LAMS(K)**3))

! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

! V1.3
! ASSUME COLLECTION OF SNOW BY RAIN PRODUCES GRAUPEL NOT HAIL

! V1.5
!            IF (IHAIL.EQ.0) THEN
            IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
            PSACR(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                 N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                  (5./(LAMS(K)**3*LAMR(K))+                    &
                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                  0.5/(LAMS(K)*LAMR(K)**3)))            
            END IF
!            END IF

         END IF

!.......................................................................

! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990, 
! USED BY REISNER ET AL 1998
         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

	    END IF

!.......................................................................
! RIME-SPLINTERING - SNOW
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

!v1.4
         IF (QNI3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWS(K).GT.0..OR.PRACS(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

               IF (PSACWS(K).GT.0.) THEN
                  NMULTS(K) = 35.E4*PSACWS(K)*FMULT*1000.
                  QMULTS(K) = NMULTS(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                  QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                  PSACWS(K) = PSACWS(K)-QMULTS(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACS(K).GT.0.) THEN
                   NMULTR(K) = 35.E4*PRACS(K)*FMULT*1000.
                   QMULTR(K) = NMULTR(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                   QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                   PRACS(K) = PRACS(K)-QMULTR(K)

               END IF

            END IF
         END IF
         END IF
         END IF

!.......................................................................
! RIME-SPLINTERING - GRAUPEL 
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

! V1.3
! ONLY CALCULATE FOR GRAUPEL NOT HAIL
! V1.5
!         IF (IHAIL.EQ.0) THEN
! v1.4
         IF (QG3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWG(K).GT.0..OR.PRACG(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

               IF (PSACWG(K).GT.0.) THEN
                  NMULTG(K) = 35.E4*PSACWG(K)*FMULT*1000.
                  QMULTG(K) = NMULTG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                  QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                  PSACWG(K) = PSACWG(K)-QMULTG(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACG(K).GT.0.) THEN
                   NMULTRG(K) = 35.E4*PRACG(K)*FMULT*1000.
                   QMULTRG(K) = NMULTRG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                   QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                   PRACG(K) = PRACG(K)-QMULTRG(K)

               END IF
               
            END IF
         END IF
      END IF
   END IF
!         END IF

!........................................................................
! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL
! ASSUME CONVERTED SNOW FORMS GRAUPEL NOT HAIL
! HAIL ASSUMED TO ONLY FORM BY FREEZING OF RAIN
! OR COLLISIONS OF RAIN WITH CLOUD ICE

! V1.3
! V1.5
!           IF (IHAIL.EQ.0) THEN
   IF (PSACWS(K).GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
      IF (QNI3D(K).GE.0.1E-3.AND.QC3D(K).GE.0.5E-3) THEN

! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
         PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
              ASN(K)*ASN(K)/ &
              (RHO(K)*LAMS(K)**(2.*BS+2.))) 

! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
	     DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.) 

! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
	     NSCNG(K) = DUM/MG0*RHO(K)
! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
             NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

! PORTION OF RIMING LEFT FOR SNOW
             PSACWS(K) = PSACWS(K) - PGSACW(K)
             END IF
	   END IF

! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

	   IF (PRACS(K).GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
	      DUM = CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3 &    
                   /(CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3+ &  
                   CONS19*(4./LAMR(K))**3*(4./LAMR(K))**3)
              DUM=MIN(DUM,1.)
              DUM=MAX(DUM,0.)
	      PGRACS(K) = (1.-DUM)*PRACS(K)
            NGRACS(K) = (1.-DUM)*NPRACS(K)
! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
            NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
            NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

! AMOUNT LEFT FOR SNOW PRODUCTION
            PRACS(K) = PRACS(K) - PGRACS(K)
            NPRACS(K) = NPRACS(K) - NGRACS(K)
! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
            PSACR(K)=PSACR(K)*(1.-DUM)
            END IF
	   END IF
!           END IF

!.......................................................................
! FREEZING OF RAIN DROPS
! FREEZING ALLOWED BELOW -4 C

         IF (T3D(K).LT.269.15.AND.QR3D(K).GE.QSMALL) THEN

! IMMERSION FREEZING (BIGG 1953)
!            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3 &
!                 /LAMR(K)**3

!            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3

! hm fix 7/15/13 for consistency w/ original formula
            MNUCCR(K) = CONS20*NR3D(K)*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3 &
                 /LAMR(K)**3

            NNUCCR(K) = PI*NR3D(K)*BIMM*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3

! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
            NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

         END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

      IF (IRAIN.EQ.0) THEN

         DUM=(QC3D(K)*QR3D(K))
         PRA(K) = 67.*(DUM)**1.15
         NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))
         
      ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! seifert and beheng (2001) formulation

         dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
         dum1 = (dum/(dum+5.e-4))**4
         pra(k) = 5.78e3*rho(k)/1000.*qc3d(k)*qr3d(k)*dum1
         npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
              (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)
         
      END IF
   END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

! v1.4 replace with seifert and beheng (2001)

         IF (QR3D(K).GE.1.E-8) THEN
! v1.4
! seifert and beheng
! include breakup, V2.1
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
            nragg(k) = -5.78*dum*qr3d(k)*nr3d(k)*rho(k)

         END IF

!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

         IF (QI3D(K).GE.1.E-8 .AND.QVQVSI(K).GE.1.) THEN

!           COFFI = 2./LAMI(K)
!           IF (COFFI.GE.DCS) THEN
              NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
              PRCI(K) = CONS22*NPRCI(K)
              NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

!           END IF
         END IF

!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

         IF (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) THEN
            PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                     LAMS(K)**(BS+3.)
            NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                  RHO(K)*N0S(K)/                                 &
                  LAMS(K)**(BS+3.)
            NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
         END IF

!.......................................................................
! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
! FOLLOWS REISNER ET AL. 1998
! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

   IF (QR3D(K).GE.1.E-8.AND.QI3D(K).GE.1.E-8.AND.T3D(K).LE.TMELT) THEN

! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
! OTHERWISE ADD TO SNOW

            IF (QR3D(K).GE.0.1E-3) THEN
            NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
            NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
            ELSE 
            NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
            NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
            END IF
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

   IF (INUC.EQ.0) THEN

! FREEZING OF AEROSOL ONLY ALLOWED BELOW -5 C
! AND ABOVE DELIQUESCENCE THRESHOLD OF 80%
! AND ABOVE ICE SATURATION

! add threshold according to Greg Thomspon

      if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
           QVQVSI(K).ge.1.08) then

! hm, modify dec. 5, 2006, replace with cooper curve
         kc2 = 0.005*exp(0.304*(TMELT-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
         kc2 = min(kc2,500.e3)
         kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1

         IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
            NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
            MNUCCD(K) = NNUCCD(K)*MI0
         END IF

      END IF
      
   ELSE IF (INUC.EQ.1) THEN

      IF (T3D(K).LT.TMELT.AND.QVQVSI(K).GT.1.) THEN

             KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
             IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
                NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
                MNUCCD(K) = NNUCCD(K)*MI0
             END IF
          END IF
          
       END IF
       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

101    CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

! NO VENTILATION FOR CLOUD ICE

       IF (QI3D(K).GE.QSMALL) THEN

          EPSI = 2.*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

       ELSE
          EPSI = 0.
       END IF
       
       IF (QNI3D(K).GE.QSMALL) THEN
          EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
               (F1S/(LAMS(K)*LAMS(K))+                       &
               F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
               SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
       ELSE
          EPSS = 0.
       END IF
       
       IF (QG3D(K).GE.QSMALL) THEN
          EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
               (F1S/(LAMG(K)*LAMG(K))+                               &
               F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
               SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))


      ELSE
      EPSG = 0.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
               F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
        EPSR = 0.
      END IF

! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

       IF (QI3D(K).GE.QSMALL) THEN              
          DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
          PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
       ELSE
          DUM=0.
       END IF
! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
       IF (QNI3D(K).GE.QSMALL) THEN
          PRDS(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
               EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
! OTHERWISE ADD TO CLOUD ICE
       ELSE
          PRD(K) = PRD(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
       END IF
! VAPOR DPEOSITION ON GRAUPEL
       PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)
       
! NO CONDENSATION ONTO RAIN, ONLY EVAP
       
       IF (QV3D(K).LT.QVS(K)) THEN
          PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
          PRE(K) = MIN(PRE(K),0.)
       ELSE
          PRE(K) = 0.
       END IF

! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
! FORMULA FROM REISNER 2 SCHEME

       DUM = (QV3D(K)-QVI(K))/DT
       
       FUDGEF = 0.9999
       SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)
       
       IF( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
            (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
          MNUCCD(K) = FUDGEF*MNUCCD(K)*DUM/SUM_DEP
          PRD(K) = FUDGEF*PRD(K)*DUM/SUM_DEP
          PRDS(K) = FUDGEF*PRDS(K)*DUM/SUM_DEP
          PRDG(K) = FUDGEF*PRDG(K)*DUM/SUM_DEP
       ENDIF
       
! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

       IF (PRD(K).LT.0.) THEN
          EPRD(K)=PRD(K)
          PRD(K)=0.
       END IF
       IF (PRDS(K).LT.0.) THEN
          EPRDS(K)=PRDS(K)
          PRDS(K)=0.
       END IF
       IF (PRDG(K).LT.0.) THEN
          EPRDG(K)=PRDG(K)
          PRDG(K)=0.
       END IF

!.......................................................................
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CONSERVATION OF WATER
! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.
! THIS SECTION IS SEPARATED INTO TWO PARTS, IF T < 0 C, T > 0 C
! DUE TO DIFFERENT PROCESSES THAT ACT DEPENDING ON FREEZING/ABOVE FREEZING

! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

! NOTE: CONSERVATION NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
! STEP

!****SENSITIVITY - NO ICE

       IF (ILIQ.EQ.1) THEN
          MNUCCC(K)=0.
          NNUCCC(K)=0.
          MNUCCR(K)=0.
          NNUCCR(K)=0.
          MNUCCD(K)=0.
          NNUCCD(K)=0.
       END IF

! ****SENSITIVITY - NO GRAUPEL
       IF (IGRAUP.EQ.1) THEN
          PRACG(K) = 0.
          PSACR(K) = 0.
          PSACWG(K) = 0.
          PGSACW(K) = 0.
          PGRACS(K) = 0.
          PRDG(K) = 0.
          EPRDG(K) = 0.
          EVPMG(K) = 0.
          PGMLT(K) = 0.
          NPRACG(K) = 0.
          NPSACWG(K) = 0.
          NSCNG(K) = 0.
          NGRACS(K) = 0.
          NSUBG(K) = 0.
          NGMLTG(K) = 0.
          NGMLTR(K) = 0.
! v3 5/27/11
          PIACRS(K)=PIACRS(K)+PIACR(K)
          PIACR(K) = 0.
! fix 070713
	    PRACIS(K)=PRACIS(K)+PRACI(K)
	    PRACI(K) = 0.
	    PSACWS(K)=PSACWS(K)+PGSACW(K)
	    PGSACW(K) = 0.
	    PRACS(K)=PRACS(K)+PGRACS(K)
	    PGRACS(K) = 0.
       END IF

! CONSERVATION OF QC

       DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT
       
       IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN
          RATIO = QC3D(K)/DUM

          PRC(K) = PRC(K)*RATIO
          PRA(K) = PRA(K)*RATIO
          MNUCCC(K) = MNUCCC(K)*RATIO
          PSACWS(K) = PSACWS(K)*RATIO
          PSACWI(K) = PSACWI(K)*RATIO
          QMULTS(K) = QMULTS(K)*RATIO
          QMULTG(K) = QMULTG(K)*RATIO
          PSACWG(K) = PSACWG(K)*RATIO
          PGSACW(K) = PGSACW(K)*RATIO
       END IF
 
! CONSERVATION OF QI

       DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
            -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT
       
       IF (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) THEN
          
          RATIO = (QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
               MNUCCD(K)+PSACWI(K))/ &
               (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))

          PRCI(K) = PRCI(K)*RATIO
          PRAI(K) = PRAI(K)*RATIO
          PRACI(K) = PRACI(K)*RATIO
          PRACIS(K) = PRACIS(K)*RATIO
          EPRD(K) = EPRD(K)*RATIO
          
       END IF

! CONSERVATION OF QR

       DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
            PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT
       
       IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN
          
          RATIO = (QR3D(K)/DT+PRC(K)+PRA(K))/ &
               (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))

          PRE(K) = PRE(K)*RATIO
          PRACS(K) = PRACS(K)*RATIO
          QMULTR(K) = QMULTR(K)*RATIO
          QMULTRG(K) = QMULTRG(K)*RATIO
          MNUCCR(K) = MNUCCR(K)*RATIO
          PIACR(K) = PIACR(K)*RATIO
          PIACRS(K) = PIACRS(K)*RATIO
          PGRACS(K) = PGRACS(K)*RATIO
          PRACG(K) = PRACG(K)*RATIO

       END IF
       
! CONSERVATION OF QNI
! CONSERVATION FOR GRAUPEL SCHEME

       IF (IGRAUP.EQ.0) THEN
          
          DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT
          
          IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN
             
             RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K))
             
             EPRDS(K) = EPRDS(K)*RATIO
             PSACR(K) = PSACR(K)*RATIO
             
          END IF

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       ELSE IF (IGRAUP.EQ.1) THEN
          
          DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT
          
          IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN
             
             RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K))
             
             EPRDS(K) = EPRDS(K)*RATIO
             PSACR(K) = PSACR(K)*RATIO
             
          END IF
          
       END IF
       
! CONSERVATION OF QG

       DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT
       
       IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN
          
          RATIO = (QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
               PIACR(K)+PRACI(K))/(-EPRDG(K))
          
          EPRDG(K) = EPRDG(K)*RATIO

       END IF

! TENDENCIES

       QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))
       
! v3 5/27/11 bug fix
       T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
            *XXLV(K)+(PRD(K)+PRDS(K)+                            &
            MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
            (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
            QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
            +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*XLF(K))/CPM(K)
       
       QC3DTEN(K) = QC3DTEN(K)+                                      &
            (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
            PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
       QI3DTEN(K) = QI3DTEN(K)+                                      &
            (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
            PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
       QR3DTEN(K) = QR3DTEN(K)+                                      &
            (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
            -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))
       
       IF (IGRAUP.EQ.0) THEN
          
          QNI3DTEN(K) = QNI3DTEN(K)+                                    &
               (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
          NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
          QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
               PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
          NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       ELSE IF (IGRAUP.EQ.1) THEN

          QNI3DTEN(K) = QNI3DTEN(K)+                                    &
               (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
          NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))
          
       END IF
       
       !bloss(2018-02): Patch in warm cloud aerosol-in-rain treatment here
       !  Should work here, even if we are neglecting aerosol in ice.
       IF((IPRGAER.EQ.1).AND.(IPRECOFF.EQ.0)) THEN
         IF(QC3D(K).GE.QSMALL) THEN
           QAPRA(K) = QAW3D(K)*PRA(K)/QC3D(K)
           QAPRC(K) = QAW3D(K)*PRC(K)/QC3D(K)
           DUM = MIN((PRA(K)+PRC(K))*DT/QC3D(K),1.)
         ELSE
           DUM = 0.
           QAPRA(K) = 0.
           QAPRC(K) = 0.
         END IF
         QAW3DTEN(K) = QAW3DTEN(K)-QAPRA(K)-QAPRC(K)

         IF(QR3D(K).GE.QSMALL) THEN
           DUM1 = MAX(PRE(K)*DT/QR3D(K),-1.) !Revisit tendency due to evap
           QAPRE(K) = QAR3D(K)*DUM1/DT
         ELSE
           QAPRE(K) = 0.
           DUM1 = 0.
         END IF

         QAR3DTEN(K) = QAR3DTEN(K)+QAR3D(K)*(DUM1/DT)+QAW3D(K)*(DUM/DT)
         QAD3DTEN(K) = QAD3DTEN(K)+QAR3D(K)*(-DUM1/DT)
       END IF

       NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
            -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))
       
       NI3DTEN(K) = NI3DTEN(K)+                                      &
            (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
            NNUCCD(K)-NIACR(K)-NIACRS(K))
       
       NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
            +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))

! V1.3 move code below to before saturation adjustment
       IF (EPRD(K).LT.0.) THEN
          DUM = EPRD(K)*DT/QI3D(K)
          DUM = MAX(-1.,DUM)
          NSUBI(K) = DUM*NI3D(K)/DT
       END IF
       IF (EPRDS(K).LT.0.) THEN
          DUM = EPRDS(K)*DT/QNI3D(K)
          DUM = MAX(-1.,DUM)
          NSUBS(K) = DUM*NS3D(K)/DT
       END IF
       IF (PRE(K).LT.0.) THEN
          DUM = PRE(K)*DT/QR3D(K)
          DUM = MAX(-1.,DUM)
          NSUBR(K) = DUM*NR3D(K)/DT
       END IF
       IF (EPRDG(K).LT.0.) THEN
          DUM = EPRDG(K)*DT/QG3D(K)
          DUM = MAX(-1.,DUM)
          NSUBG(K) = DUM*NG3D(K)/DT
       END IF

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

       NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
       NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
       NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
       NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF(ISATADJ.EQ.0) THEN !PB 4/13/09

! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = 0.622*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

          DUMS = DUMQV-DUMQSS
          PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
          IF (PCC(K)*DT+DUMQC.LT.0.) THEN
             PCC(K) = -DUMQC/DT
          END IF
          
          QV3DTEN(K) = QV3DTEN(K)-PCC(K)
          T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
          QC3DTEN(K) = QC3DTEN(K)+PCC(K)
          
       END IF

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS

!bloss: only do activation if droplet number is predicted
!bloss      IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) THEN
       IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL.AND.INUM.EQ.0) THEN
          
! EFFECTIVE VERTICAL VELOCITY (M/S)

          IF (ISUB.EQ.0) THEN
! ADD SUB-GRID VERTICAL VELOCITY
             DUM = W3D(K)+WVAR(K)
             
! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
             DUM = MAX(DUM,0.10)

          ELSE IF (ISUB.EQ.1) THEN
!         DUM=W3D(K)
         DUM=W3D(K)*max(1.,dx/1000.)  ! Marat:take into account dependence of W on hor gridsize
          END IF

! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
          IF (DUM.GE.0.001) THEN
             
             IF (IBASE.EQ.1) THEN

! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)
                
                IDROP=0
                
! V1.3 USE CURRENT VALUE OF QC FOR IDROP
                IF (QC3D(K).LE.0.05E-3/RHO(K)) THEN
                   IDROP=1
                END IF
                IF (K.EQ.1) THEN
                   IDROP=1
                ELSE IF (K.GE.2) THEN
                   IF (QC3D(K).GT.0.05E-3/RHO(K).AND. &
                        QC3D(K-1).LE.0.05E-3/RHO(K-1)) THEN
                      IDROP=1
                   END IF
                END IF

                IF (IDROP.EQ.1) THEN
! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

                   IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

                      DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
                      DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
                      DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
                      DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                      DUM2 = (DUM2-NC3D(K))/DT
                      DUM2 = MAX(0.,DUM2)
                      NC3DTEN(K) = NC3DTEN(K)+DUM2
                      
                   ELSE IF (IACT.EQ.2) THEN

! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)
                      
                      IF ((NANEW1.GT.0.AND.RM1.GT.0).OR.(NANEW2.GT.0.AND.RM2.GT.0)) THEN
                              
                         SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                         AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
                         ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
                         GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))
                         
                         GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                              (T3D(K)*RR)-1.))
                         
                         PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT
                         
                         IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                            ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
                            SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                            DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
                         ELSE
                            DUM1 = 0.
                         END IF
                         
                         IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                            ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)
                            SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                            DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)
                         ELSE 
                            DUM2 = 0.
                         END IF

                         D3 = 0.                                
                         IF (dofixedcoarsemode) THEN                                   
                            ETA3 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*Ncoarse*RHO(K)*1.e6)  ! convert Ncoarse from /mg to /m^3
                            SM3 = 2./BACT_coarse**0.5*(AACT/(3.*RM3))**1.5
                            D3 = 1./SM3**2*(F13*(PSI/ETA3)**1.5+F23*(SM3**2/(ETA3+3.*PSI))**0.75)
                         END IF
                         
                         SMAX = 1./(DUM1+DUM2+D3)**0.5
                         
                         IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                            UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                            DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                         ELSE
                            UU1 = 0.
                         END IF
                         
                         IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                            UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                            DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                         ELSE
                            UU2 = 0.
                         END IF

                         ISACT(K) = 1.
                         SSPK(K) = SMAX
                         DC1(K) = DG1(K) * (SM1/SMAX)**(2./3.) !Critical diameters
                         DC2(K) = DG2(K) * (SM1/SMAX)**(2./3.)
                         
                         NA1 = NANEW1/RHO(K)  ! units of /kg  total number of mode 1
                         QA1 = MAER1/RHO(K)   ! units of kg/kg total mass of mode 1

                         NARG1(K) = DUM1/RHO(K) ! number of prescribed activated in /kg
                         NARG2(K) = DUM2/RHO(K)

                         NARG1(K) = MIN(NARG1(K), NA1) ! is this really necessary?
                         NARG2(K) = MIN(NARG2(K), NANEW2/RHO(K))

                         NACTDIFF(K) = (NARG1(K) + NARG2(K) - NC3D(K)) ! saved stat

                         NCNEW = NC3D(K) + MAX(0., NACTDIFF(K)) ! only increase nc 

                         NCNEW = MIN(NA1, NCNEW)  ! cap activation could make this optional

                         NTRANS = 0.
                         QTRANS = 0.
                         IF (doacttransfer) THEN
                            call hoppel_aitken_accum_transfer(NAD2_3D(K), NA1, QAD2_3D(K), QA1, &
                                 SG2, SG1, DC2(K), DC1(K), RHOA, RHOA, NTRANS, QTRANS)
                         END IF

                         ! transfer to accum from aitken
                         NA1 = NA1 + NTRANS
                         QA1 = QA1 + QTRANS

                         NADNEW = NA1 - NCNEW
                         NACTRATE(K) = (NCNEW - NC3D(K))/DT

                         ! activated mass diagnosed
                         QAWNEW = QA1 * (1. - mass_fraction(NADNEW/NA1, sigma_accum))
                         QACTRATE(K) = (QAWNEW - QAW3D(K))/DT
                         QACTRATE(K) = min(QACTRATE(K), (QAD3D(K) + QATRANSFER(K)*DT)/DT)

                         NATRANSFER(K) = NTRANS/DT  
                         QATRANSFER(K) = QTRANS/DT              

                         NAD2_3DTEN(K) = NAD2_3DTEN(K) - NATRANSFER(K)
                         QAD2_3DTEN(K) = QAD2_3DTEN(K) - QATRANSFER(K)

                         NAD3DTEN(K) = NAD3DTEN(K) + NATRANSFER(K) - NACTRATE(K)
                         QAD3DTEN(K) = QAD3DTEN(K) + QATRANSFER(K) - QACTRATE(K)

                         NC3DTEN(K) = NC3DTEN(K) + NACTRATE(K)
                         QAW3DTEN(K) = QAW3DTEN(K) + QACTRATE(K)
                                 
                     
                      END IF ! (NANEW1.GT.0).OR.(NANEW2.GT.0)           
                   END IF  ! IACT

!.............................................................................
                ELSE IF (IDROP.EQ.0) THEN
! ACTIVATE IN CLOUD INTERIOR
! FIND EQUILIBRIUM SUPERSATURATION

                   TAUC=1./(2.*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.)/LAMC(K))
                   IF (EPSR.GT.1.E-8) THEN
                      TAUR=1./EPSR
                   ELSE
                      TAUR=1.E8
                   END IF
                   IF (EPSI.GT.1.E-8) THEN
                      TAUI=1./EPSI
                   ELSE
                      TAUI=1.E8
                   END IF
                   IF (EPSS.GT.1.E-8) THEN
                      TAUS=1./EPSS
                   ELSE
                      TAUS=1.E8
                   END IF
                   IF (EPSG.GT.1.E-8) THEN
                      TAUG=1./EPSG
                   ELSE
                      TAUG=1.E8
                   END IF

! EQUILIBRIUM SS INCLUDING BERGERON EFFECT

! hm fix 1/20/15
!           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=(-QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
                   DUM3=(DUM3*TAUC*TAUR*TAUI*TAUS*TAUG- &
                        (QVS(K)-QVI(K))*(TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS))/ &
                        (TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS+ &
                        TAUR*TAUI*TAUS*TAUG+TAUC*TAUI*TAUS*TAUG)
                   
                   IF (DUM3/QVS(K).GE.1.E-6) THEN
                      IF (IACT.EQ.1) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

                         DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
                         DUMACT = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
                         
! USE POWER LAW CCN SPECTRA

! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
                         DUM3=DUM3/QVS(K)*100.

                         DUM2=C1*DUM3**K1
! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
                         DUM2=MIN(DUM2,DUMACT)
                         DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
                         DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                         DUM2 = (DUM2-NC3D(K))/DT
                         DUM2 = MAX(0.,DUM2)
                         NC3DTEN(K) = NC3DTEN(K)+DUM2
                         
                      ELSE IF (IACT.EQ.2) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

                         IF ((NANEW1.GT.0.AND.RM1.GT.0).OR.(NANEW2.GT.0.AND.RM2.GT.0)) THEN
                                 
                            SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                            AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
                            ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
                            GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))
                            
                            GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                 (T3D(K)*RR)-1.))
                            
                            PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT
                            
                            IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                               ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
                               SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                               DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
                            ELSE
                               DUM1 = 0.
                               SM1 = 1.
                            END IF
                            
                            IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                               ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)
                               SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                               DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)
                            ELSE 
                               DUM2 = 0.
                               SM2 = 1.
                            END IF

                            D3 = 0.                                
                            IF (dofixedcoarsemode) THEN                                   
                               ETA3 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*Ncoarse*RHO(K)*1.e6)  ! convert Ncoarse from /mg to /m^3
                               SM3 = 2./BACT_coarse**0.5*(AACT/(3.*RM3))**1.5
                               D3 = 1./SM3**2*(F13*(PSI/ETA3)**1.5+F23*(SM3**2/(ETA3+3.*PSI))**0.75)
                            END IF

                            SMAX = 1./(DUM1+DUM2+D3)**0.5
                            
                            UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                            UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                            DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                            DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                            
                            DUM1 = DUM1/RHO(K)
                            DUM2 = DUM2/RHO(K) !CONVERT TO KG-1

                            SSPK(K) = SMAX
! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL
                            DUMACT1 = MIN(DUM1, NANEW1/RHO(K))
                            DUMACT2 = MIN(DUM2, NANEW2/RHO(K))

! AEROSOL MASS FOR EQUILIBRIUM SUPER SATURATION USING PROGNOSTIC AEROSOL

! USE LOGNORMAL AEROSOL
                            SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                            AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)

! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
                            SMAX = DUM3/QVS(K)

                            
                            IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                               SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                               UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                               DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                            ELSE
                               DUM1 = 0.
                               UU1 = 0.
                            END IF
                            
                            IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                               SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                               UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                               DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                            ELSE
                               DUM2 = 0.
                               UU2 = 0.
                            END IF
                            
                            DUM1 = DUM1/RHO(K)
                            DUM2 = DUM2/RHO(K) !CONVERT TO KG-1
                           
                            
! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL AND
! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
                            DUM1 = MIN(DUM1, NANEW1/RHO(K), DUMACT1)
                            DUM2 = MIN(DUM2, NANEW2/RHO(K), DUMACT2)

                            ISACT(K) = 1.
                            SSPK(K) = MIN(SSPK(K), SMAX)
                            DC1(K) = DG1(K) * (SM1/SSPK(K))**(2./3.) !Critical diameters
                            DC2(K) = DG2(K) * (SM2/SSPK(K))**(2./3.)

                            NA1 = NANEW1/RHO(K)
                            QA1 = MAER1/RHO(K)

                            NARG1(K) = DUM1
                            NARG2(K) = DUM2

                            NARG1(K) = MIN(NARG1(K), NA1) ! is this really necessary?
                            NARG2(K) = MIN(NARG2(K), NANEW2/RHO(K))

                            NACTDIFF(K) = (NARG1(K) + NARG2(K) - NC3D(K)) ! saved stat

                            NCNEW = NC3D(K) + MAX(0., NACTDIFF(K)) ! only increase nc 

                            NCNEW = MIN(NA1, NCNEW)  ! cap activation could make this optional

                            NTRANS = 0.
                            QTRANS = 0.
                            IF (doacttransfer) THEN
                               call hoppel_aitken_accum_transfer(NAD2_3D(K), NA1, QAD2_3D(K), QA1, &
                                    SG2, SG1, DC2(K), DC1(K), RHOA, RHOA, NTRANS, QTRANS)
                            END IF

                            ! transfer to accum from aitken
                            NA1 = NA1 + NTRANS
                            QA1 = QA1 + QTRANS

                            NADNEW = NA1 - NCNEW
                            NACTRATE(K) = (NCNEW - NC3D(K))/DT

                            QAWNEW = QA1 * (1. - mass_fraction(NADNEW/NA1, sigma_accum))
                            QACTRATE(K) = (QAWNEW - QAW3D(K))/DT
                            QACTRATE(K) = min(QACTRATE(K), (QAD3D(K) + QATRANSFER(K)*DT)/DT)

                            NATRANSFER(K) = NTRANS/DT  
                            QATRANSFER(K) = QTRANS/DT              

                            NAD2_3DTEN(K) = NAD2_3DTEN(K) - NATRANSFER(K)
                            QAD2_3DTEN(K) = QAD2_3DTEN(K) - QATRANSFER(K)

                            NAD3DTEN(K) = NAD3DTEN(K) + NATRANSFER(K) - NACTRATE(K)
                            QAD3DTEN(K) = QAD3DTEN(K) + QATRANSFER(K) - QACTRATE(K)

                            NC3DTEN(K) = NC3DTEN(K) + NACTRATE(K)
                            QAW3DTEN(K) = QAW3DTEN(K) + QACTRATE(K)                          
                     

                                    

                         END IF ! (NANEW1.GT.0).OR.(NANEW2.GT.0)
                      END IF ! IACT
                   END IF ! DUM3/QVS > 1.E-6
                END IF  ! IDROP = 1
                
!.......................................................................
             ELSE IF (IBASE.EQ.2) THEN
                
                IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

                   DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
                   DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
                   DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
                   DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                   DUM2 = (DUM2-NC3D(K))/DT
                   DUM2 = MAX(0.,DUM2)
                   NC3DTEN(K) = NC3DTEN(K)+DUM2

                ELSE IF (IACT.EQ.2) THEN

! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)
                           
                   IF ((NANEW1.GT.0.AND.RM1.GT.0).OR.(NANEW2.GT.0.AND.RM2.GT.0)) THEN
                           
                      SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
                      AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
                      ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
                      GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))
                      
                      GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                           (T3D(K)*RR)-1.))
                      
                      PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT
                      
                      IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                         ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
                         SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
                         DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
                      ELSE                         
                         DUM1 = 0.
                      END IF
                      
                      IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                         ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)
                         SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5
                         DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)
                      ELSE 
                         DUM2 = 0.
                      END IF

                      D3 = 0.                                
                      IF (dofixedcoarsemode) THEN                                   
                         ETA3 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*Ncoarse*RHO(K)*1.e6)  ! convert Ncoarse from /mg to /m^3
                         SM3 = 2./BACT_coarse**0.5*(AACT/(3.*RM3))**1.5
                         D3 = 1./SM3**2*(F13*(PSI/ETA3)**1.5+F23*(SM3**2/(ETA3+3.*PSI))**0.75)
                      END IF
                      
                      SMAX = 1./(DUM1+DUM2+D3)**0.5
                      
                      IF (NANEW1.GT.0.AND.RM1.GT.0) THEN
                         UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
                         DUM1 = NANEW1/2.*(1.-DERF1(UU1))
                      ELSE
                         UU1 = 0.
                      END IF
                      
                      IF (NANEW2.GT.0.AND.RM2.GT.0.AND.doaitkenactivate) THEN
                         UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
                         DUM2 = NANEW2/2.*(1.-DERF1(UU2))
                      ELSE
                         UU2 = 0.
                      END IF

                      ISACT(K) = 1.
                      SSPK(K) = SMAX ! save statistic
                      DC1(K) = DG1(K) * (SM1/SMAX)**(2./3.) !Critical diameters
                      DC2(K) = DG2(K) * (SM2/SMAX)**(2./3.)

                      NA1 = NANEW1/RHO(K)  ! units of /kg  total number of mode 1
                      QA1 = MAER1/RHO(K)   ! units of kg/kg total mass of mode 1

                      NARG1(K) = DUM1/RHO(K) ! number of prescribed activated in /kg
                      NARG2(K) = DUM2/RHO(K)

                      NARG1(K) = MIN(NARG1(K), NA1) ! is this really necessary?
                      NARG2(K) = MIN(NARG2(K), NANEW2/RHO(K))

                      NACTDIFF(K) = (NARG1(K) + NARG2(K) - NC3D(K)) ! saved stat

                      NCNEW = NC3D(K) + MAX(0., NACTDIFF(K)) ! only increase nc 

                      NCNEW = MIN(NA1, NCNEW)  ! cap activation could make this optional

                      NTRANS = 0.
                      QTRANS = 0.
                      IF (doacttransfer) THEN
                         call hoppel_aitken_accum_transfer(NAD2_3D(K), NA1, QAD2_3D(K), QA1, &
                              SG2, SG1, DC2(K), DC1(K), RHOA, RHOA, NTRANS, QTRANS)
                      END IF

                      ! transfer to accum from aitken
                      NA1 = NA1 + NTRANS
                      QA1 = QA1 + QTRANS

                      ! activated mass diagnosed
                      QAWNEW = QA1 * (1. - mass_fraction(NADNEW/NA1, sigma_accum))
                      QACTRATE(K) = (QAWNEW - QAW3D(K))/DT
                      QACTRATE(K) = min(QACTRATE(K), (QAD3D(K) + QATRANSFER(K)*DT)/DT)

                      NATRANSFER(K) = NTRANS/DT  
                      QATRANSFER(K) = QTRANS/DT              

                      NAD2_3DTEN(K) = NAD2_3DTEN(K) - NATRANSFER(K)
                      QAD2_3DTEN(K) = QAD2_3DTEN(K) - QATRANSFER(K)

                      NAD3DTEN(K) = NAD3DTEN(K) + NATRANSFER(K) - NACTRATE(K)
                      QAD3DTEN(K) = QAD3DTEN(K) + QATRANSFER(K) - QACTRATE(K)

                      NC3DTEN(K) = NC3DTEN(K) + NACTRATE(K)
                      QAW3DTEN(K) = QAW3DTEN(K) + QACTRATE(K)




                   END IF ! (NANEW1.GT.0).OR.(NANEW2.GT.0)           
                END IF  ! IACT
             END IF  ! IBASE
          END IF  ! W > 0.001
        END IF  ! QC3D > QSMALL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

        IF (PCC(K).LT.0) THEN  
          IF (IEVPNC.EQ.1) THEN  !brnr modified to be turned on with flag in prm file 
            DUM = PCC(K)*DT/QC3D(K)
            DUM = MAX(-1.,DUM)
            NSUBC(K) = DUM*NC3D(K)/DT
            NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
           IF (IPRGAER.EQ.1) THEN
             NAD3DTEN(K) = NAD3DTEN(K)-NSUBC(K) !brnr dry aerosol tendency is opposite of cloud droplet tendency
             DUM1 = DUM*QAW3D(K)/DT
             QAW3DTEN(K) = QAW3DTEN(K)+DUM1 !brnr use Ivanova et al 2008 assumption 
             QAD3DTEN(K) = QAD3DTEN(K)-DUM1
           END IF
         ELSEIF ((QC3D(K)+PCC(K)*DT).LT.QSMALL) THEN !brnr if evaporation in a gridbox will essentially evaporate all 
           NC3DTEN(K) = NC3DTEN(K)-NC3D(K)/DT   !     cloud, transfer all cloud number back to dry aerosol number
           IF (IPRGAER.EQ.1) THEN
             NAD3DTEN(K) = NAD3DTEN(K)+NC3D(K)/DT !     and nucleated aerosol mass to dry aerosol mass
             QAW3DTEN(K) = -QAW3D(K)/DT
             QAD3DTEN(K) = QAD3DTEN(K)+QAW3D(K)/DT !WORKING HERE
           END IF
         END IF
       END IF
       
       IF (EPRD(K).LT.0.) THEN
         DUM = EPRD(K)*DT/QI3D(K)
         DUM = MAX(-1.,DUM)
         NSUBI(K) = DUM*NI3D(K)/DT
       END IF
       IF (EPRDS(K).LT.0.) THEN
         DUM = EPRDS(K)*DT/QNI3D(K)
         DUM = MAX(-1.,DUM)
         NSUBS(K) = DUM*NS3D(K)/DT
       END IF
       IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
         DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
       END IF
       IF (EPRDG(K).LT.0.) THEN
         DUM = EPRDG(K)*DT/QG3D(K)
         DUM = MAX(-1.,DUM)
         NSUBG(K) = DUM*NG3D(K)/DT
       END IF
       
!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
         NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
         NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
         NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
         NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)

       END IF !!!!!! TEMPERATURE

! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
         LTRUE = 1

 200     CONTINUE

       END DO

! V1.3 move precip initialization to here
! INITIALIZE PRECIP AND SNOW RATES
       PRECRT = 0.
       SNOWRT = 0.

! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

       IF (LTRUE.EQ.0) GOTO 400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!.......................................................................
! CALCULATE SEDIMENATION
! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
! STABILITY, I.E. COURANT# < 1

!.......................................................................

 NSTEPR=1
 NSTEPI=1
 NSTEPS=1
 NSTEPC=1
 NSTEPG=1
! v3 5/27/11
 DO K = KTE,KTS,-1
    
    DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
    DUMQS(K)  = QNI3D(K)+QNI3DTEN(K)*DT
    DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
    DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
    DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
    DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
    DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
    DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
    DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
    DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT

    DUMQAW(K) = QAW3D(K)+QAW3DTEN(K)*DT
    DUMQAR(K) = QAR3D(K)+QAR3DTEN(K)*DT
         
!ADDSTUFF3
    
    
! SWITCH FOR CONSTANT DROPLET NUMBER
    IF (INUM.EQ.1) THEN
       DUMFNC(K) = NC3D(K)
    END IF
    
! GET DUMMY LAMDA

! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
    DUMFNI(K) = MAX(0.,DUMFNI(K))
    DUMFNS(K) = MAX(0.,DUMFNS(K))
    DUMFNC(K) = MAX(0.,DUMFNC(K))
    DUMFNR(K) = MAX(0.,DUMFNR(K))
    DUMFNG(K) = MAX(0.,DUMFNG(K))
    
!......................................................................
! CLOUD ICE

    IF (DUMI(K).GE.QSMALL) THEN
       DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
       DLAMI=MAX(DLAMI,LAMMINI)
       DLAMI=MIN(DLAMI,LAMMAXI)
    END IF
!......................................................................
! RAIN

      IF (DUMR(K).GE.QSMALL) THEN
        DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1./3.)
        DLAMR=MAX(DLAMR,LAMMINR)
        DLAMR=MIN(DLAMR,LAMMAXR)
      END IF
!......................................................................
! CLOUD DROPLETS

    IF (DUMC(K).GE.QSMALL) THEN
       !bloss: option for fixing pgam
       if(dofix_pgam) then
          pgam(k) = pgam_fixed
       else
          
!         DUM = PRES(K)/(R*T3D(K))
! V1.5
          PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
          PGAM(K)=1./(PGAM(K)**2)-1.
          PGAM(K)=MAX(PGAM(K),2.)
          PGAM(K)=MIN(PGAM(K),10.)
          
       end if

       DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.)/(DUMC(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
       LAMMIN = (PGAM(K)+1.)/60.E-6
       LAMMAX = (PGAM(K)+1.)/1.E-6
       DLAMC=MAX(DLAMC,LAMMIN)
       DLAMC=MIN(DLAMC,LAMMAX)
    END IF
!......................................................................
! SNOW
    
    IF (DUMQS(K).GE.QSMALL) THEN
       DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1./DS)
       DLAMS=MAX(DLAMS,LAMMINS)
       DLAMS=MIN(DLAMS,LAMMAXS)
    END IF
!......................................................................
! GRAUPEL
    
    IF (DUMG(K).GE.QSMALL) THEN
       DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1./DG)
       DLAMG=MAX(DLAMG,LAMMING)
       DLAMG=MIN(DLAMG,LAMMAXG)
    END IF

!......................................................................
! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

! CLOUD WATER

      IF (DUMC(K).GE.QSMALL) THEN
      UNC =  ACN(K)*GAMMA(1.+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.))
      UMC = ACN(K)*GAMMA(4.+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.))
      ELSE
      UMC = 0.
      UNC = 0.
      END IF

      IF (DUMI(K).GE.QSMALL) THEN
      UNI =  AIN(K)*CONS27/DLAMI**BI
      UMI = AIN(K)*CONS28/(DLAMI**BI)
      ELSE
      UMI = 0.
      UNI = 0.
      END IF

      IF (DUMR(K).GE.QSMALL) THEN
      UNR = ARN(K)*CONS6/DLAMR**BR
      UMR = ARN(K)*CONS4/(DLAMR**BR)
      ELSE
        UMR = 0.
        UNR = 0.
      END IF

      IF (DUMQS(K).GE.QSMALL) THEN
        UMS = ASN(K)*CONS3/(DLAMS**BS)
        UNS = ASN(K)*CONS5/DLAMS**BS
      ELSE
        UMS = 0.
        UNS = 0.
      END IF

      IF (DUMG(K).GE.QSMALL) THEN
        UMG = AGN(K)*CONS7/(DLAMG**BG)
        UNG = AGN(K)*CONS8/DLAMG**BG
      ELSE
        UMG = 0.
        UNG = 0.
      END IF

! SET REALISTIC LIMITS ON FALLSPEED

! bug fix, 10/08/09
      dum=(rhosu/rho(k))**0.54
      UMS=MIN(UMS,1.2*dum)
      UNS=MIN(UNS,1.2*dum)
! v3 5/27/11
! fix for correction by AA 4/6/11
      UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)
      UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)
      UMR=MIN(UMR,9.1*dum)
      UNR=MIN(UNR,9.1*dum)
      UMG=MIN(UMG,20.*dum)
      UNG=MIN(UNG,20.*dum)
    
      FR(K) = UMR
      FI(K) = UMI
      FNI(K) = UNI
      FS(K) = UMS
      FNS(K) = UNS
      FNR(K) = UNR
      FC(K) = UMC
      FNC(K) = UNC
      FG(K) = UMG
      FNG(K) = UNG

! v3 5/27/11 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

      IF (K.LE.KTE-1) THEN
        IF (FR(K).LT.1.E-10) THEN
          FR(K)=FR(K+1)
        END IF
        IF (FI(K).LT.1.E-10) THEN
          FI(K)=FI(K+1)
        END IF
        IF (FNI(K).LT.1.E-10) THEN
          FNI(K)=FNI(K+1)
        END IF
        IF (FS(K).LT.1.E-10) THEN
          FS(K)=FS(K+1)
        END IF
        IF (FNS(K).LT.1.E-10) THEN
          FNS(K)=FNS(K+1)
        END IF
        IF (FNR(K).LT.1.E-10) THEN
          FNR(K)=FNR(K+1)
        END IF
        IF (FC(K).LT.1.E-10) THEN
          FC(K)=FC(K+1)
        END IF
        IF (FNC(K).LT.1.E-10) THEN
          FNC(K)=FNC(K+1)
        END IF
        IF (FG(K).LT.1.E-10) THEN
          FG(K)=FG(K+1)
        END IF
        IF (FNG(K).LT.1.E-10) THEN
          FNG(K)=FNG(K+1)
        END IF
      END IF ! K LE KTE-1
    
      RGVM = MAX(FR(K),FNR(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPR = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEPR)
    
      RGVM = MAX(FI(K),FNI(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPI = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEPI)

      RGVM = MAX(FS(K),FNS(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPS = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEPS)

      RGVM = MAX(FC(K),FNC(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPC = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEPC)

      RGVM = MAX(FG(K),FNG(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPG = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEPG)

! MULTIPLY VARIABLES BY RHO
      DUMR(k) = DUMR(k)*RHO(K)
      DUMI(k) = DUMI(k)*RHO(K)
      DUMFNI(k) = DUMFNI(K)*RHO(K)
      DUMQS(k) = DUMQS(K)*RHO(K)
      DUMFNS(k) = DUMFNS(K)*RHO(K)
      DUMFNR(k) = DUMFNR(K)*RHO(K)
      DUMC(k) = DUMC(K)*RHO(K)
      DUMFNC(k) = DUMFNC(K)*RHO(K)
      DUMG(k) = DUMG(K)*RHO(K)
      DUMFNG(k) = DUMFNG(K)*RHO(K)
      
      DUMQAW(K) = DUMQAW(K)*RHO(K)
      DUMQAR(K) = DUMQAR(K)*RHO(K)
      
    END DO

! V1.3, change so that sub-stepping is done
! individually for each species


! RAIN

 DO N = 1,NSTEPR
    
    DO K = KTS,KTE
       FALOUTR(K) = FR(K)*DUMR(K)
       FALOUTNR(K) = FNR(K)*DUMFNR(K)
       FALOUTQAR(K) = FR(K)*DUMQAR(K) !brnr fall out of aerosol in rain
    END DO

! TOP OF MODEL
    K = KTE
    FALTNDR = FALOUTR(K)/DZQ(k)
    FALTNDNR = FALOUTNR(K)/DZQ(k)
    FALTNDQAR = FALOUTQAR(K)/DZQ(k)
    
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
    QRSTEN(K) = QRSTEN(K)-FALTNDR/NSTEPR/RHO(k)
    NRSTEN(K) = NRSTEN(K)-FALTNDNR/NSTEPR/RHO(k)
    NR3DTEN(K) = NR3DTEN(K)-FALTNDNR/NSTEPR/RHO(k)
    DUMR(K) = DUMR(K)-FALTNDR*DT/NSTEPR
    DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT/NSTEPR
    DUMQAR(K) = DUMQAR(K)-FALTNDQAR*DT/NSTEPR
    
    DO K = KTE-1,KTS,-1
       FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)
       FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)
       QRSTEN(K) = QRSTEN(K)+FALTNDR/NSTEPR/RHO(k)
       NRSTEN(K) = NRSTEN(K)+FALTNDNR/NSTEPR/RHO(k)
       NR3DTEN(K) = NR3DTEN(K)+FALTNDNR/NSTEPR/RHO(k)
       DUMR(K) = DUMR(K)+FALTNDR*DT/NSTEPR
       DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT/NSTEPR
       IF((IPRGAER.EQ.1).AND.(IPRECOFF.EQ.0)) THEN
          FALTNDQAR = (FALOUTQAR(K+1)-FALOUTQAR(K))/DZQ(K)
          QARSTEN(K) = QARSTEN(K)+FALTNDQAR/NSTEPR/RHO(k)
          DUMQAR(K) = DUMQAR(K)+FALTNDQAR*DT/NSTEPR
       END IF
    END DO
    PRECRT = PRECRT+(FALOUTR(KTS))  &
         *DT/NSTEPR
 END DO
 
! CLOUD ICE

 DO N = 1,NSTEPI
    
    DO K = KTS,KTE
       FALOUTI(K) = FI(K)*DUMI(K)
       FALOUTNI(K) = FNI(K)*DUMFNI(K)
    END DO

! TOP OF MODEL
    K = KTE
    FALTNDI = FALOUTI(K)/DZQ(k)
    FALTNDNI = FALOUTNI(K)/DZQ(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
    QISTEN(K) = QISTEN(K)-FALTNDI/NSTEPI/RHO(k)
    NI3DTEN(K) = NI3DTEN(K)-FALTNDNI/NSTEPI/RHO(k)
    DUMI(K) = DUMI(K)-FALTNDI*DT/NSTEPI
    DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT/NSTEPI
    DO K = KTE-1,KTS,-1
       FALTNDI = (FALOUTI(K+1)-FALOUTI(K))/DZQ(K)
       FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))/DZQ(K)
       QISTEN(K) = QISTEN(K)+FALTNDI/NSTEPI/RHO(k)
       NI3DTEN(K) = NI3DTEN(K)+FALTNDNI/NSTEPI/RHO(k)
       DUMI(K) = DUMI(K)+FALTNDI*DT/NSTEPI
       DUMFNI(K) = DUMFNI(K)+FALTNDNI*DT/NSTEPI
    END DO
    PRECRT = PRECRT+(FALOUTI(KTS))  &
         *DT/NSTEPI
    SNOWRT = SNOWRT+(FALOUTI(KTS))*DT/NSTEPI
 END DO

! CLOUD DROPLETS

 DO N = 1,NSTEPC
    
    DO K = KTS,KTE
       FALOUTC(K) = FC(K)*DUMC(K)
       FALOUTNC(K) = FNC(K)*DUMFNC(K)
       FALOUTQAW(K) = FC(K)*DUMQAW(K) ! sediment activated aerosol mass at mass weighted droplet fall rate
    END DO
! TOP OF MODEL

    K = KTE
    FALTNDC = FALOUTC(K)/DZQ(k)
    FALTNDNC = FALOUTNC(K)/DZQ(k)
    FALTNDQAW = FALOUTQAW(K)/DZQ(K) ! Tendency due to top of model flux

! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
    QCSTEN(K) = QCSTEN(K)-FALTNDC/NSTEPC/RHO(k)
    NCSTEN(K) = NCSTEN(K)-FALTNDNC/NSTEPC/RHO(k)
    NC3DTEN(K) = NC3DTEN(K)-FALTNDNC/NSTEPC/RHO(k)
    DUMC(K) = DUMC(K)-FALTNDC*DT/NSTEPC
    DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT/NSTEPC
    IF (IPRGAER.EQ.1) THEN
       QAWSTEN(K) = QAWSTEN(K)-FALTNDQAW/NSTEPC/RHO(K)
       DUMQAW(K) = DUMQAW(K)-FALTNDQAW*DT/NSTEPC
    END IF
    
    DO K = KTE-1,KTS,-1
       FALTNDC = (FALOUTC(K+1)-FALOUTC(K))/DZQ(K)
       FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))/DZQ(K)
       QCSTEN(K) = QCSTEN(K)+FALTNDC/NSTEPC/RHO(k)
       NCSTEN(K) = NCSTEN(K)+FALTNDNC/NSTEPC/RHO(k)
       NC3DTEN(K) = NC3DTEN(K)+FALTNDNC/NSTEPC/RHO(k)
       DUMC(K) = DUMC(K)+FALTNDC*DT/NSTEPC
       DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT/NSTEPC
       IF (IPRGAER.EQ.1) THEN
          FALTNDQAW = (FALOUTQAW(K+1)-FALOUTQAW(K))/DZQ(K)                  
          QAWSTEN(K) = QAWSTEN(K)+FALTNDQAW/NSTEPC/RHO(K) ! sediment 
          DUMQAW(K) = DUMQAW(K)+FALTNDQAW*DT/NSTEPC
       END IF
    END DO
    PRECRT = PRECRT+(FALOUTC(KTS))  &
         *DT/NSTEPC
 END DO
 
! SNOW

 DO N = 1,NSTEPS
    
    DO K = KTS,KTE
       FALOUTS(K) = FS(K)*DUMQS(K)
       FALOUTNS(K) = FNS(K)*DUMFNS(K)
    END DO

! TOP OF MODEL
    K = KTE
    FALTNDS = FALOUTS(K)/DZQ(k)
    FALTNDNS = FALOUTNS(K)/DZQ(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
    QNISTEN(K) = QNISTEN(K)-FALTNDS/NSTEPS/RHO(k)
    NS3DTEN(K) = NS3DTEN(K)-FALTNDNS/NSTEPS/RHO(k)
    DUMQS(K) = DUMQS(K)-FALTNDS*DT/NSTEPS
    DUMFNS(K) = DUMFNS(K)-FALTNDNS*DT/NSTEPS
    DO K = KTE-1,KTS,-1
       FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
       FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
       QNISTEN(K) = QNISTEN(K)+FALTNDS/NSTEPS/RHO(k)
       NS3DTEN(K) = NS3DTEN(K)+FALTNDNS/NSTEPS/RHO(k)
       DUMQS(K) = DUMQS(K)+FALTNDS*DT/NSTEPS
       DUMFNS(K) = DUMFNS(K)+FALTNDNS*DT/NSTEPS
    END DO
    PRECRT = PRECRT+(FALOUTS(KTS))  &
         *DT/NSTEPS
    SNOWRT = SNOWRT+(FALOUTS(KTS))*DT/NSTEPS
 END DO

! GRAUPEL

 DO N = 1,NSTEPG

    DO K = KTS,KTE
       FALOUTG(K) = FG(K)*DUMG(K)
       FALOUTNG(K) = FNG(K)*DUMFNG(K)
    END DO

! TOP OF MODEL
    K = KTE
    FALTNDG = FALOUTG(K)/DZQ(k)
    FALTNDNG = FALOUTNG(K)/DZQ(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
    QGSTEN(K) = QGSTEN(K)-FALTNDG/NSTEPG/RHO(k)
    NG3DTEN(K) = NG3DTEN(K)-FALTNDNG/NSTEPG/RHO(k)
    DUMG(K) = DUMG(K)-FALTNDG*DT/NSTEPG
    DUMFNG(K) = DUMFNG(K)-FALTNDNG*DT/NSTEPG
    DO K = KTE-1,KTS,-1
       FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
       FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)
       QGSTEN(K) = QGSTEN(K)+FALTNDG/NSTEPG/RHO(k)
       NG3DTEN(K) = NG3DTEN(K)+FALTNDNG/NSTEPG/RHO(k)
       DUMG(K) = DUMG(K)+FALTNDG*DT/NSTEPG
       DUMFNG(K) = DUMFNG(K)+FALTNDNG*DT/NSTEPG
    END DO
    PRECRT = PRECRT+(FALOUTG(KTS))  &
         *DT/NSTEPG
      SNOWRT = SNOWRT+(FALOUTG(KTS))*DT/NSTEPG
   END DO

   DO K=KTS,KTE

! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

      QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
      QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
      QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
      QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
      QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)
      IF (IPRGAER.EQ.1) THEN
         QAW3DTEN(K)=QAW3DTEN(K)+QAWSTEN(K)
         QAR3DTEN(K)=QAR3DTEN(K)+QARSTEN(K)
      END IF
      
!DBG var updates
!            DBGIND = 6
!            DBGQAR3DTEN(DBGIND,K) = QAR3DTEN(K)
!            DBGQAR3D(DBGIND,K) = DBGQAR3D(2,K) + DBGQAR3DTEN(DBGIND,K)*DT
!            
!            DBGQAW3DTEN(DBGIND,K) = QAW3DTEN(K)
!            DBGQAW3D(DBGIND,K) = DBGQAW3D(1,K) + DBGQAW3DTEN(DBGIND,K)*DT
!            
!            DBGQAD3DTEN(DBGIND,K) = QAD3DTEN(K)
!            DBGQAD3D(DBGIND,K) = DBGQAD3D(1,K) + DBGQAD3DTEN(DBGIND,K)*DT
!
!
!            IF((MINVAL(DBGQAR3D).LT.0)) THEN
!               QNEG = .TRUE.
!               print*,'error, mass is negative!!! DBGIND = 5'
!            END IF


! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

! V1.7
!hm 7/9/09 bug fix
!        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
      IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.TMELT.AND.LAMI(K).GE.1.E-10) THEN
         
         IF (1./LAMI(K).GE.2.*DCS) THEN
            QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
            NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
            QI3DTEN(K) = -QI3D(K)/DT
            NI3DTEN(K) = -NI3D(K)/DT
         END IF
      END IF
      
! hm add tendencies here, then call sizeparameter
! to ensure consisitency between mixing ratio and number concentration
      
      QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
      QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
      QNI3D(k)       = QNI3D(k)+QNI3DTEN(k)*DT
      QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
      NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
      NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
      NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
      NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT
      
      IF (IGRAUP.EQ.0) THEN
         QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
         NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
      END IF
      
      IF (IPRGAER.EQ.1) THEN
         NAD3D(K)    = NAD3D(K)+NAD3DTEN(K)*DT
         NAD2_3D(K)  = NAD2_3D(K) + NAD2_3DTEN(K)*DT
         QAD3D(K)    = QAD3D(K)+QAD3DTEN(K)*DT
         QAD2_3D(K)  = QAD2_3D(k) + QAD2_3DTEN(K)*DT
         QAW3D(K)    = QAW3D(K)+QAW3DTEN(K)*DT
         QAR3D(K)    = QAR3D(K)+QAR3DTEN(K)*DT
      END IF

!!$      IF(NAD3D(K).GT.NAD3D(1)) write(*,*) 'Line 5410: NAD3D(',k,') = ', NAD3D(k)

!IF((MINVAL(DBGQAR3D).LT.0)) THEN
!   QNEG = .TRUE.
!  print*,'error, mass is negative!!!, DBGIND = 6'
!END IF
      
!ADDSTUFF5
!      DBGNCTND(6) = DBGNCTND(6)+NC3DTEN(K) 
!      DBGNADTND(6) = DBGNADTND(6)+NAD3DTEN(K)
!      DBGNC(6) = DBGNC(6)+NC3D(K)
!      DBGNAD(6) =DBGNAD(6)+NAD3D(K)

!NEED TO COMPENSATE AEROSOL BUDGET FOR THE EFFECT OF LIMITERS BELOW

! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
      T3D(K)         = T3D(K)+T3DTEN(k)*DT
      QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT
      
! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
      EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
      EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

      
! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

      IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)
      
      QVS(K) = .622*EVS(K)/(PRES(K)-EVS(K))
      QVI(K) = .622*EIS(K)/(PRES(K)-EIS(K))
      
      QVQVS(K) = QV3D(K)/QVS(K)
      QVQVSI(K) = QV3D(K)/QVI(K)
      
! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER

! V1.3, change limit from 10^-7 to 10^-6
! V1.7 7/9/09 change limit from 10^-6 to 10^-8

      IF (QVQVS(K).LT.0.9) THEN
         IF (QR3D(K).LT.1.E-8) THEN
            QVPOSLIM(K) = QVPOSLIM(K)+QR3D(K)/DT !Limiter adds QR to QV
            QRNEGLIM(K) = QRNEGLIM(K)-QR3D(K)/DT !Limiter zeros QR
            QV3D(K)=QV3D(K)+QR3D(K)
            T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
            QR3D(K)=0.
         END IF
         IF (QC3D(K).LT.1.E-8) THEN
            QVPOSLIM(K) = QVPOSLIM(K)+QC3D(K)/DT !Limiter adds QC to QV
            QCNEGLIM(K) = QCNEGLIM(K)-QC3D(K)/DT !Limiter zeros QC
            QV3D(K)=QV3D(K)+QC3D(K)
            T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
            QC3D(K)=0.
         END IF
      END IF
      
      IF (QVQVSI(K).LT.0.9) THEN
         IF (QI3D(K).LT.1.E-8) THEN
            QV3D(K)=QV3D(K)+QI3D(K)
            T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
            QI3D(K)=0.
         END IF
         IF (QNI3D(K).LT.1.E-8) THEN
            QV3D(K)=QV3D(K)+QNI3D(K)
            T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
            QNI3D(K)=0.
         END IF
         IF (QG3D(K).LT.1.E-8) THEN
            QV3D(K)=QV3D(K)+QG3D(K)
            T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
            QG3D(K)=0.
         END IF
      END IF

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

      IF (QC3D(K).LT.QSMALL) THEN
         QCNEGLIM(K) = QCNEGLIM(K)-QC3D(K)/DT
         QAWNEGLIM(K) = QAWNEGLIM(K)-QAW3D(K)/DT
         NCNEGLIM(K) = NCNEGLIM(K)-NC3D(K)/DT
         QC3D(K) = 0.
         QAW3D(K) = 0.
         NC3D(K) = 0.
      END IF
      IF (QR3D(K).LT.QSMALL) THEN
         QRNEGLIM(K) = QRNEGLIM(K)-QR3D(K)/DT
         QARNEGLIM(K) = QARNEGLIM(K)-QAR3D(K)/DT
         NRNEGLIM(K) = NRNEGLIM(K)-NR3D(K)/DT
         QR3D(K) = 0.
         QAR3D(K) = 0.
         NR3D(K) = 0.
      END IF
      IF (QI3D(K).LT.QSMALL) THEN
         QI3D(K) = 0.
         NI3D(K) = 0.
      END IF
      IF (QNI3D(K).LT.QSMALL) THEN
         QNI3D(K) = 0.
         NS3D(K) = 0.
      END IF
      IF (QG3D(K).LT.QSMALL) THEN
         QG3D(K) = 0.
         NG3D(K) = 0.
      END IF

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

      IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
           .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) GOTO 500
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE INSTANTANEOUS PROCESSES

! ADD MELTING OF CLOUD ICE TO FORM RAIN

      IF (QI3D(K).GE.QSMALL.AND.T3D(K).GE.TMELT) THEN
         QR3D(K) = QR3D(K)+QI3D(K)
         T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
! hm new process rate output
           QMELTI(K)=QI3D(K)/DT
           NMELTI(K)=NI3D(K)/DT
         QI3D(K) = 0.
         NR3D(K) = NR3D(K)+NI3D(K)
         NI3D(K) = 0.
      END IF

! ****SENSITIVITY - NO ICE
      IF (ILIQ.EQ.1) GOTO 778
      
! HOMOGENEOUS FREEZING OF CLOUD WATER

      IF (T3D(K).LE.233.15.AND.QC3D(K).GE.QSMALL) THEN
         QI3D(K)=QI3D(K)+QC3D(K)
         T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
! hm new process rate output
           QHOMOC(K)=QC3D(K)/DT
           NHOMOC(K)=NC3D(K)/DT
         QC3D(K)=0.
         NI3D(K)=NI3D(K)+NC3D(K)
         NC3D(K)=0.
      END IF

! HOMOGENEOUS FREEZING OF RAIN
      
      IF (IGRAUP.EQ.0) THEN
         
         IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
            QG3D(K) = QG3D(K)+QR3D(K)
            T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
! hm new process rate output
           QHOMOR(K)=QR3D(K)/DT
           NHOMOR(K)=NR3D(K)/DT
            QR3D(K) = 0.
            NG3D(K) = NG3D(K)+ NR3D(K)
            NR3D(K) = 0.
         END IF
         
      ELSE IF (IGRAUP.EQ.1) THEN
         
         IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
            QNI3D(K) = QNI3D(K)+QR3D(K)
            T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
! hm new process rate output
           QHOMOR(K)=QR3D(K)/DT
           NHOMOR(K)=NR3D(K)/DT
            QR3D(K) = 0.
            NS3D(K) = NS3D(K)+NR3D(K)
            NR3D(K) = 0.
         END IF
         
      END IF
      
778   CONTINUE
      
! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      
      DUM = NC3D(K)
      NC3D(K) = MAX(0.,NC3D(K))
      NCPOSLIM(K) = NCPOSLIM(K)+(NC3D(K)-DUM)/DT
      
      DUM = NR3D(K)
      NR3D(K) = MAX(0.,NR3D(K))
      NRPOSLIM(K) = NRPOSLIM(K)+(NR3D(K)-DUM)/DT

      NG3D(K) = MAX(0.,NG3D(K))
      
      IF(IPRGAER.EQ.1) THEN
         DUM = NAD3D(K)
         NAD3D(K) = MAX(0.,NAD3D(K))
         NADPOSLIM(K) = NADPOSLIM(K)+(NAD3D(K)-DUM)/DT
      END IF
      !......................................................................
! CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)

! CHECK FOR SLOPE

! ADJUST VARS

         IF (LAMI(K).LT.LAMMINI) THEN
            
            LAMI(K) = LAMMINI
            
            N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12
            
            NI3D(K) = N0I(K)/LAMI(K)
         ELSE IF (LAMI(K).GT.LAMMAXI) THEN
            LAMI(K) = LAMMAXI
            N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12
            
            NI3D(K) = N0I(K)/LAMI(K)
         END IF
      END IF

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
        DUM2 = NR3D(K)
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN
        DUM2 = NR3D(K)

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)
      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF

      DUM1 = (NR3D(K)-DUM2)/DT
      IF (DUM1.LT.0.) THEN
        NRNEGLIM(K) = NRNEGLIM(K) + DUM1
        NADPOSLIM(K) = NADPOSLIM(K) - DUM1
        NAD3D(K) = NAD3D(K) - DUM1*DT
      ELSE
        NRPOSLIM(K) = NRPOSLIM(K) + DUM1 
        DUM1 = MAX(-NAD3D(K)/DT,-DUM1)
        NADNEGLIM(K) = NADNEGLIM(K) + DUM1
        NAD3D(K) = NAD3D(K) + DUM1*DT
      END IF
    END IF

!!$      IF(NAD3D(K).GT.NAD3D(1)) write(*,*) 'Line 5761: NAD3D(',k,') = ', NAD3D(k)

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         IF(dofix_pgam) THEN
            pgam(k) = pgam_fixed
         ELSE

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
            PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
            PGAM(K)=1./(PGAM(K)**2)-1.
            PGAM(K)=MAX(PGAM(K),2.)
            PGAM(K)=MIN(PGAM(K),10.)
            
         END IF
         
! CALCULATE LAMC
         
         LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
              (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

! Brnr modify the bound check routines to conserve Na+Nc
                  
         LAMMIN = (PGAM(K)+1.)/60.E-6
         LAMMAX = (PGAM(K)+1.)/1.E-6
         
         DUM1 = NC3D(K)
         DUM2 = NAD3D(K)
         
         IF (LAMC(K).LT.LAMMIN) THEN
            LAMC(K) = LAMMIN
            IF (IPRGAER.EQ.1) THEN
               DUM = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                    LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                        !bloss(2018-02): prevent creation of number
                        !bloss(2019-03): prevent problems if total cloud+aerosol number is zero or negative
                        NC3D(K) = MAX(0.1*DUM, MIN(DUM, DUM1+DUM2) )                        
                        NAD3D(K) = MAX(0., DUM1+DUM2 -NC3D(K) )
                        ! RE-CALCULATE LAMC
                        LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                             (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
!bloss(2018-02)                        NAD3D(K) = MAX(0.,NAD3D(K)-(DUM-NC3D(K)))
!bloss(2018-02)                        NC3D(K) = DUM
            ELSE
               NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                    LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26              
            END IF !IPRGAER.EQ.1
         ELSE IF (LAMC(K).GT.LAMMAX) THEN
            LAMC(K) = LAMMAX
            IF (IPRGAER.EQ.1) THEN
               DUM = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                    LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                        !bloss(2018-02): prevent creation of number
                        NC3D(K) = MIN(DUM, DUM1+DUM2)    
                        NAD3D(K) = MAX(0., DUM1+DUM2 -NC3D(K) )
                        ! RE-CALCULATE LAMC
                        LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                             (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
!bloss(2018-02)                        NAD3D(K) = MAX(0.,NAD3D(K)-(DUM-NC3D(K)))
!bloss(2018-02)                        NC3D(K) = DUM
            ELSE
               NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                    LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
            END IF !IPRGAER.EQ.1
         END IF !LAMC.LT.LAMMIN
         
         !BRNR Tracking for limiters
         DUM = (NC3D(K)-DUM1)/DT
         IF (DUM.LT.0.) THEN
            NCNEGLIM(K) = NCNEGLIM(K) + DUM
            NADPOSLIM(K) = NADPOSLIM(K) + (NAD3D(K)-DUM2)/DT
         ELSE
            NCPOSLIM(K) = NCPOSLIM(K) + DUM
            NADNEGLIM(K) = NADNEGLIM(K) + (NAD3D(K)-DUM2)/DT
         END IF
  
      END IF !QC3D.LT.QSMALL

!!$      IF(NAD3D(K).GT.NAD3D(1)) write(*,*) 'Line 5758: NAD3D(',k,') = ', NAD3D(k)

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
         LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
         
! CHECK FOR SLOPE

! ADJUST VARS
         
         IF (LAMS(K).LT.LAMMINS) THEN
            LAMS(K) = LAMMINS
            N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1
            
            NS3D(K) = N0S(K)/LAMS(K)
            
         ELSE IF (LAMS(K).GT.LAMMAXS) THEN
            
            LAMS(K) = LAMMAXS
            N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1
            NS3D(K) = N0S(K)/LAMS(K)
         END IF

      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
         LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
         
! CHECK FOR SLOPE

! ADJUST VARS

         IF (LAMG(K).LT.LAMMING) THEN
            LAMG(K) = LAMMING
            N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2
            
            NG3D(K) = N0G(K)/LAMG(K)
            
         ELSE IF (LAMG(K).GT.LAMMAXG) THEN
            
            LAMG(K) = LAMMAXG
            N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2
            
            NG3D(K) = N0G(K)/LAMG(K)
         END IF
         
      END IF
      
500   CONTINUE

! CALCULATE EFFECTIVE RADIUS

      IF (QI3D(K).GE.QSMALL) THEN
         EFFI(K) = 3./LAMI(K)/2.*1.E6
      ELSE
         EFFI(K) = 25.
      END IF
      
      IF (QNI3D(K).GE.QSMALL) THEN
         EFFS(K) = 3./LAMS(K)/2.*1.E6
      ELSE
         EFFS(K) = 25.
      END IF
      
      IF (QR3D(K).GE.QSMALL) THEN
         EFFR(K) = 3./LAMR(K)/2.*1.E6
      ELSE
         EFFR(K) = 25.
      END IF

      IF (QC3D(K).GE.QSMALL) THEN
         EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
              GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
      ELSE
         EFFC(K) = 25.
      END IF
      
      IF (QG3D(K).GE.QSMALL) THEN
         EFFG(K) = 3./LAMG(K)/2.*1.E6
      ELSE
         EFFG(K) = 25.
      END IF

! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
!          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
! HM, 3/4/13, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
! OF EXCESSIVE AND PERSISTENT ANVIL
! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
          NI3D(K) = MIN(NI3D(K),0.3E6/RHO(K))

! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
      IF (INUM.EQ.0.AND.IACT.EQ.2) THEN
         DUM = NC3D(K)
         NC3D(K) = MAX(0.,NC3D(K)) !MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K)))
         NCPOSLIM(K) = NCPOSLIM(K)+(NC3D(K)-DUM)/DT
      END IF
      IF (IPRGAER.EQ.1) THEN
         DUM = NAD3D(K)
         NAD3D(K) = MAX(0.,NAD3D(K))
         NADPOSLIM(K) = NADPOSLIM(K)+(NAD3D(K)-DUM)/DT
      END IF
! SWITCH FOR CONSTANT DROPLET NUMBER
      IF (INUM.EQ.1) THEN
! CHANGE NDCNST FROM CM-3 TO KG-1
         NC3D(K) = NDCNST*1.E6/RHO(K)
      END IF
      
!bloss: preserve Berner output of gamma exponent (even though it duplicates pgam)
      NUC3D(K) = PGAM(K)

   END DO !!! K LOOP
   
400 CONTINUE

      if(do_accumulate_micro_proc_rates) then
        !bloss: Accumulate microphysics process rates
        idx = 0
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PCC       ! COND/EVAP DROPLETS
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PCCN      ! CHANGE Q DROPLET ACTIVATION
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSUBC     ! LOSS OF NC DURING EVAP
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRE       ! EVAP OF RAIN
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSUBR     ! LOSS OF NR DURING EVAP
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRA       ! ACCRETION DROPLETS BY RAIN
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRC       ! AUTOCONVERSION DROPLETS
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
        idx = idx + 1
        micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NRAGG     ! SELF-COLLECTION OF RAIN
        if (IPRGAER.eq.1) then
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NARG1     ! change in Nd from activation mode 1 (not rate)
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NARG2     ! change in Nd from activation mode 2 (not rate)
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NACTRATE     ! activation - rate of transfer of aerosol mass from unactivated to activated 
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QACTRATE     ! activation - rate of transfer of aerosol mass from unactivated to activated 
          idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NACTDIFF     ! activation - difference between ARG suggested NC and existing NC (not rate)
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NATRANSFER ! transfer rate from Aitken to accum
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QATRANSFER ! transfer rate from Aitken to accum
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + ISACT !  activation is run on particular point
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + DC1*ISACT !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + DC2*ISACT !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + DG1*ISACT !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + DG2*ISACT !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + SSPK*ISACT !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + DG1 !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:, idx) = micro_proc_rates(:,idx) + DG2 !  Critical diameter of activation
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QAPRA     ! autoconversion - xfer of aerosol mass from qaw to qar
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QAPRC     ! accretion -xfer of aerosol mass from qaw to qar.
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QAPRE     ! source of qad due to rain evaporatoin
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QASUBC     ! source of qad due to cloud evaporation
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NCPOSLIM  
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NCNEGLIM  
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NRPOSLIM  
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NRNEGLIM  
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NADPOSLIM  
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NADNEGLIM
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QVPOSLIM  
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QCNEGLIM    
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QRNEGLIM    
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QAWNEGLIM    
           idx = idx + 1
           micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QARNEGLIM 
       end if   
        if(iliq.eq.0) then
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSUBI     ! LOSS OF NI DURING SUB.
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSUBS     ! LOSS OF NS DURING SUB.
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRD       ! DEP CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRDS      ! DEP SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSAGG     ! SELF-COLLECTION OF SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRAI      ! CHANGE Q ACCRETION CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRCI      ! CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRAI     ! CHANGE N ACCRETION CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSMLTS    ! CHANGE N MELTING SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NIACR     ! CHANGE N, ICE-RAIN COLLECTION
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + EPRD      ! SUBLIMATION CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + EPRDS     ! SUBLIMATION SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PRDG    ! DEP OF GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + EPRDG    ! SUB OF GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PGMLT    ! CHANGE Q MELTING OF GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NGMLTG    ! CHANGE N MELTING GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
          ! hm new process rate output
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QHOMOC    ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QHOMOR    ! CHANGE Q DUE TO HOMOGENEOUS FREEZING OF RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NHOMOC    ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NHOMOR    ! CHANGE N DUE TO HOMOGENEOUS FREEZING OF RAIN
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + QMELTI    ! CHANGE Q DUE TO MELTING OF CLOUD ICE
          idx = idx + 1
          micro_proc_rates(:,idx) = micro_proc_rates(:,idx) + NMELTI    ! CHANGE N DUE TO MELTING OF CLOUD ICE
        end if
      end if

! ALL DONE !!!!!!!!!!!
   RETURN
 END  SUBROUTINE M2005MICRO_GRAUPEL
     
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL FUNCTION POLYSVP (T,TYPE)

!-------------------------------------------

!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

      IMPLICIT NONE

      REAL DUM
      REAL T
      INTEGER TYPE

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8 

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE.EQ.1) THEN

!         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
!          LOG10(6.1071))*100.


      dt = max(-80.,t-273.16)
      polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
      polysvp = polysvp*100.

      END IF

! LIQUID

      IF (TYPE.EQ.0) THEN

       dt = max(-80.,t-273.16)
       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp = polysvp*100.

!         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
!             5.02808*LOG10(373.16/T)-									&
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
!             LOG10(1013.246))*100.

         END IF


      END FUNCTION POLYSVP

!------------------------------------------------------------------------------

      REAL FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      INTEGER I,N
      LOGICAL PARITY
      REAL                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,				   &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA


      REAL FUNCTION DERF1(X)
      IMPLICIT NONE
      REAL X
      REAL, DIMENSION(0 : 64) :: A, B
      REAL W,T,Y
      INTEGER K,I
      DATA A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,	                          &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,							  &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0, 							  &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,							  &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      DATA (B(I), I = 0, 12) /                                  &
         -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
         -0.00000000640788827665E0, -0.00000001667820132100E0,  &
         -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
          0.00001612722157047886E0, -0.00025616361025506629E0, 	&
          0.00015380842432375365E0,  0.00815533022524927908E0, 	&
         -0.01402283663896319337E0, -0.19746892495383021487E0,  &
          0.71511720328842845913E0 /
      DATA (B(I), I = 13, 25) /                                 &
         -0.00000000001951073787E0, -0.00000000032302692214E0,  &
          0.00000000522461866919E0,  0.00000000342940918551E0, 	&
         -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
          0.00002687044575042908E0, -0.00011843240273775776E0, 	&
         -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
          0.00909530922354827295E0, -0.20160072778491013140E0, 	&
          0.51169696718727644908E0 /
      DATA (B(I), I = 26, 38) /                                 &
         0.00000000003147682272E0, -0.00000000048465972408E0,   &
         0.00000000063675740242E0,  0.00000003377623323271E0, 	&
        -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
         0.00001947204525295057E0,  0.00002854147231653228E0, 	&
        -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
         0.02328095035422810727E0, -0.16725021123116877197E0, 	&
         0.32490054966649436974E0 /
      DATA (B(I), I = 39, 51) /                                 &
         0.00000000002319363370E0, -0.00000000006303206648E0,   &
        -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
         0.00000011371857327578E0, -0.00000211211337219663E0, 	&
         0.00000368797328322935E0,  0.00009823686253424796E0, 	&
        -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
         0.02585434424202960464E0, -0.11637092784486193258E0, 	&
         0.18267336775296612024E0 /
      DATA (B(I), I = 52, 64) /                                 &
        -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
        -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
         0.00000018006992266137E0, -0.00000078441223763969E0, 	&
        -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
        -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
         0.02064129023876022970E0, -0.06905562880005864105E0,   &
         0.09084526782065478489E0 /
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
              A(K + 11)) * T + A(K + 12)) * W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF (X .LT. 0) Y = -Y
      DERF1 = Y
      END FUNCTION DERF1
!+---+-----------------------------------------------------------------+
!
      subroutine radar_init

      IMPLICIT NONE
      INTEGER:: n
      PI5 = PI*PI*PI*PI*PI
      lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
      m_w_0 = m_complex_water_ray (lamda_radar, 0.0d0)
      m_i_0 = m_complex_ice_maetzler (lamda_radar, 0.0d0)
      K_w = (ABS( (m_w_0*m_w_0 - 1.0) /(m_w_0*m_w_0 + 2.0) ))**2

      do n = 1, nbins+1
         simpson(n) = 0.0d0
      enddo
      do n = 1, nbins-1, 2
         simpson(n) = simpson(n) + basis(1)
         simpson(n+1) = simpson(n+1) + basis(2)
         simpson(n+2) = simpson(n+2) + basis(3)
      enddo

      do n = 1, slen
         mixingrulestring_s(n:n) = char(0)
         matrixstring_s(n:n) = char(0)
         inclusionstring_s(n:n) = char(0)
         hoststring_s(n:n) = char(0)
         hostmatrixstring_s(n:n) = char(0)
         hostinclusionstring_s(n:n) = char(0)
         mixingrulestring_g(n:n) = char(0)
         matrixstring_g(n:n) = char(0)
         inclusionstring_g(n:n) = char(0)
         hoststring_g(n:n) = char(0)
         hostmatrixstring_g(n:n) = char(0)
         hostinclusionstring_g(n:n) = char(0)
      enddo

      mixingrulestring_s = 'maxwellgarnett'
      hoststring_s = 'air'
      matrixstring_s = 'water'
      inclusionstring_s = 'spheroidal'
      hostmatrixstring_s = 'icewater'
      hostinclusionstring_s = 'spheroidal'

      mixingrulestring_g = 'maxwellgarnett'
      hoststring_g = 'air'
      matrixstring_g = 'water'
      inclusionstring_g = 'spheroidal'
      hostmatrixstring_g = 'icewater'
      hostinclusionstring_g = 'spheroidal'

      end subroutine radar_init
!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION m_complex_water_ray(lambda,T)

!      Complex refractive Index of Water as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
!      after Ray (1972)

      IMPLICIT NONE
      REAL(8), INTENT(IN):: T,lambda
      REAL(8):: epsinf,epss,epsr,epsi
      REAL(8):: alpha,lambdas,sigma,nenner
      COMPLEX*16, PARAMETER:: i = (0d0,1d0)

      epsinf  = 5.27137d0 + 0.02164740d0 * T - 0.00131198d0 * T*T
      epss    = 78.54d+0 * (1.0 - 4.579d-3 * (T - 25.0)                 &
              + 1.190d-5 * (T - 25.0)*(T - 25.0)                        &
              - 2.800d-8 * (T - 25.0)*(T - 25.0)*(T - 25.0))
      alpha   = -16.8129d0/(T+273.16) + 0.0609265d0
      lambdas = 0.00033836d0 * exp(2513.98d0/(T+273.16)) * 1e-2

      nenner = 1.d0+2.d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PI*0.5) &
             + (lambdas/lambda)**(2d0-2d0*alpha)
      epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)   &
           * sin(alpha*PI*0.5)+1d0)) / nenner
      epsi = ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)            &
           * cos(alpha*PI*0.5)+0d0)) / nenner                           &
           + lambda*1.25664/1.88496
      
      m_complex_water_ray = SQRT(CMPLX(epsr,-epsi))
      
      END FUNCTION m_complex_water_ray

!+---+-----------------------------------------------------------------+
      
      COMPLEX*16 FUNCTION m_complex_ice_maetzler(lambda,T)
      
!      complex refractive index of ice as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.0001,30] m; T in [-250.0,0.0] C
!      Original comment from the Matlab-routine of Prof. Maetzler:
!      Function for calculating the relative permittivity of pure ice in
!      the microwave region, according to C. Maetzler, "Microwave
!      properties of ice and snow", in B. Schmitt et al. (eds.) Solar
!      System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
!      Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
!      TK = temperature (K), range 20 to 273.15
!      f = frequency in GHz, range 0.01 to 3000
         
      IMPLICIT NONE
      REAL(8), INTENT(IN):: T,lambda
      REAL(8):: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa

      c = 2.99d8
      TK = T + 273.16
      f = c / lambda * 1d-9

      B1 = 0.0207
      B2 = 1.16d-11
      b = 335.0d0
      deltabeta = EXP(-10.02 + 0.0364*(TK-273.16))
      betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1)**2) ) + B2*f*f
      beta = betam + deltabeta
      theta = 300. / TK - 1.
      alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
      m_complex_ice_maetzler = 3.1884 + 9.1e-4*(TK-273.16)
      m_complex_ice_maetzler = m_complex_ice_maetzler                   &
                             + CMPLX(0.0d0, (alfa/f + beta*f)) 
      m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))
      
      END FUNCTION m_complex_ice_maetzler
!+---+-----------------------------------------------------------------+

      subroutine rayleigh_soak_wetgraupel (x_g, a_geo, b_geo, fmelt,    &
                     meltratio_outside, m_w, m_i, lambda, C_back,       &
                     mixingrule,matrix,inclusion,                       &
                     host,hostmatrix,hostinclusion)

      IMPLICIT NONE

      REAL(8), INTENT(in):: x_g, a_geo, b_geo, fmelt, lambda,  &
                                     meltratio_outside
      REAL(8), INTENT(out):: C_back
      COMPLEX*16, INTENT(in):: m_w, m_i
      CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion,     &
                                     host, hostmatrix, hostinclusion

      COMPLEX*16:: m_core, m_air
      REAL(8):: D_large, D_g, rhog, x_w, xw_a, fm, fmgrenz,    &
                         volg, vg, volair, volice, volwater,            &
                         meltratio_outside_grenz, mra
      INTEGER:: error
      real :: rho_i, rho_w

      rho_i = 900.
      rho_w = 1000.


!     refractive index of air:
      m_air = (1.0d0,0.0d0)

!     Limiting the degree of melting --- for safety: 
      fm = DMAX1(DMIN1(fmelt, 1.0d0), 0.0d0)
!     Limiting the ratio of (melting on outside)/(melting on inside):
      mra = DMAX1(DMIN1(meltratio_outside, 1.0d0), 0.0d0)

!    ! The relative portion of meltwater melting at outside should increase
!    ! from the given input value (between 0 and 1)
!    ! to 1 as the degree of melting approaches 1,
!    ! so that the melting particle "converges" to a water drop.
!    ! Simplest assumption is linear:
      mra = mra + (1.0d0-mra)*fm

      x_w = x_g * fm

      D_g = a_geo * x_g**b_geo

      if (D_g .ge. 1d-12) then

       vg = PI/6. * D_g**3
       rhog = DMAX1(DMIN1(x_g / vg, DBLE(rho_i)), 10.0d0)
       vg = x_g / rhog
      
       meltratio_outside_grenz = 1.0d0 - rhog / rho_w

       if (mra .le. meltratio_outside_grenz) then
        !..In this case, it cannot happen that, during melting, all the
        !.. air inclusions within the ice particle get filled with
        !.. meltwater. This only happens at the end of all melting.
        volg = vg * (1.0d0 - mra * fm)
 
       else
        !..In this case, at some melting degree fm, all the air
        !.. inclusions get filled with meltwater.
        fmgrenz=(rho_i-rhog)/(mra*rho_i-rhog+rho_i*rhog/rho_w)

        if (fm .le. fmgrenz) then
         !.. not all air pockets are filled:
         volg = (1.0 - mra * fm) * vg
        else
         !..all air pockets are filled with meltwater, now the
         !.. entire ice sceleton melts homogeneously:
         volg = (x_g - x_w) / rho_i + x_w / rho_w
        endif

       endif

       D_large  = (6.0 / PI * volg) ** (1./3.)
       volice = (x_g - x_w) / (volg * rho_i)
       volwater = x_w / (rho_w * volg)
       volair = 1.0 - volice - volwater
      
       !..complex index of refraction for the ice-air-water mixture
       !.. of the particle:
       m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,      &
                         volwater, mixingrule, host, matrix, inclusion, &
                         hostmatrix, hostinclusion, error)
       if (error .ne. 0) then
        C_back = 0.0d0
        return
       endif

       !..Rayleigh-backscattering coefficient of melting particle: 
       C_back = (ABS((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2           &
                * PI5 * D_large**6 / lamda4

      else
       C_back = 0.0d0
      endif

      end subroutine rayleigh_soak_wetgraupel
!+---+-----------------------------------------------------------------+

      complex*16 function get_m_mix_nested (m_a, m_i, m_w, volair,      &
                     volice, volwater, mixingrule, host, matrix,        &
                     inclusion, hostmatrix, hostinclusion, cumulerror)

      IMPLICIT NONE

      REAL(8), INTENT(in):: volice, volair, volwater
      COMPLEX*16, INTENT(in):: m_a, m_i, m_w
      CHARACTER(len=*), INTENT(in):: mixingrule, host, matrix,          &
                     inclusion, hostmatrix, hostinclusion
      INTEGER, INTENT(out):: cumulerror

      REAL(8):: vol1, vol2
      COMPLEX*16:: mtmp
      INTEGER:: error

      !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
      !.. air, ice, or water

      cumulerror = 0
      get_m_mix_nested = CMPLX(1.0d0,0.0d0)

      if (host .eq. 'air') then

       if (matrix .eq. 'air') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, 0.0d0, vol1, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error
          
        if (hostmatrix .eq. 'air') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'icewater') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
         CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'ice') then

       if (matrix .eq. 'ice') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0d0, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'ice') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airwater') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         'air', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
         CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'water') then

       if (matrix .eq. 'water') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0d0,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'water') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,            &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airice') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,            &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
                           hostmatrix
         CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'none') then

       get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                     &
                       volair, volice, volwater, mixingrule,            &
                       matrix, inclusion, error)
       cumulerror = cumulerror + error
        
      else
       write(mp_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
       CALL wrf_debug(150, mp_debug)
       cumulerror = cumulerror + 1
      endif

      IF (cumulerror .ne. 0) THEN
       write(mp_debug,*) 'GET_M_MIX_NESTED: error encountered'
       CALL wrf_debug(150, mp_debug)
       get_m_mix_nested = CMPLX(1.0d0,0.0d0)    
      endif

      end function get_m_mix_nested

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION get_m_mix (m_a, m_i, m_w, volair, volice,     &
                     volwater, mixingrule, matrix, inclusion, error)

      IMPLICIT NONE

      REAL(8), INTENT(in):: volice, volair, volwater
      COMPLEX*16, INTENT(in):: m_a, m_i, m_w
      CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion
      INTEGER, INTENT(out):: error

      error = 0
      get_m_mix = CMPLX(1.0d0,0.0d0)

      if (mixingrule .eq. 'maxwellgarnett') then
       if (matrix .eq. 'ice') then
        get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,  &
                           m_i, m_a, m_w, inclusion, error)
       elseif (matrix .eq. 'water') then
        get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,  &
                           m_w, m_a, m_i, inclusion, error)
       elseif (matrix .eq. 'air') then
        get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,  &
                           m_a, m_w, m_i, inclusion, error)
       else
        write(mp_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
        CALL wrf_debug(150, mp_debug)
        error = 1
       endif

      else
       write(mp_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
       CALL wrf_debug(150, mp_debug)
       error = 2
      endif

      if (error .ne. 0) then
       write(mp_debug,*) 'GET_M_MIX: error encountered'
       CALL wrf_debug(150, mp_debug)
      endif

      END FUNCTION get_m_mix

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION m_complex_maxwellgarnett(vol1, vol2, vol3,    &
                     m1, m2, m3, inclusion, error)

      IMPLICIT NONE

      COMPLEX*16 :: m1, m2, m3
      REAL(8) :: vol1, vol2, vol3
      CHARACTER(len=*) :: inclusion

      COMPLEX*16 :: beta2, beta3, m1t, m2t, m3t
      INTEGER, INTENT(out) :: error

      error = 0

      if (DABS(vol1+vol2+vol3-1.0d0) .gt. 1d-6) then
       write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',       &
              'partial volume fractions is not 1...ERROR'
       CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=CMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m1t = m1**2
      m2t = m2**2
      m3t = m3**2

      if (inclusion .eq. 'spherical') then
       beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
       beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
      elseif (inclusion .eq. 'spheroidal') then
       beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
       beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
      else
       write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ',                  &
                         'unknown inclusion: ', inclusion
       CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=DCMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m_complex_maxwellgarnett = &
       SQRT(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
       (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

      END FUNCTION m_complex_maxwellgarnett

!+---+-----------------------------------------------------------------+
!..Compute radar reflectivity assuming 10 cm wavelength radar and using
!.. Rayleigh approximation.  Only complication is melted snow/graupel
!.. which we treat as water-coated ice spheres and use Uli Blahak's
!.. library of routines.  The meltwater fraction is simply the amount
!.. of frozen species remaining from what initially existed at the
!.. melting level interface.
!+---+-----------------------------------------------------------------+
      subroutine calc_refl10cm (qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ,  &
                          kts, kte, ii, jj, nr1d, ns1d, ng1d)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                qv1d, qr1d, qs1d, qg1d, t1d, p1d, nr1d, ns1d, ng1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs, rg,rnr,rns,rng

      REAL(8), DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g,ilams,n0_s

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

      REAL(8):: lamg
      REAL(8):: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0
      LOGICAL:: melti
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

!..Single melting snow/graupel particle 70% meltwater on external sfc
      REAL(8), PARAMETER:: melt_outside_s = 0.7d0
      REAL(8), PARAMETER:: melt_outside_g = 0.7d0

      REAL(8):: cback, x, eta, f_d

! hm added parameter
      REAL R1,t_0,dumlams,dumlamr,dumlamg,dumn0s,dumn0r,dumn0g,ocms,obms,ocmg,obmg

      integer n

      R1 = 1.E-12
      t_0 = 273.15

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         if (qr1d(k) .gt. R1) then
            rr(k) = qr1d(k)*rho(k)
            L_qr(k) = .true.
         else
            rr(k) = R1
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R1) then
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R1) then
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif

! hm add number concentration
         if (nr1d(k) .gt. R1) then
            rnr(k) = nr1d(k)*rho(k)
         else
            rnr(k) = R1
         endif
         if (ns1d(k) .gt. R1) then
            rns(k) = ns1d(k)*rho(k)
         else
            rns(k) = R1
         endif
         if (ng1d(k) .gt. R1) then
            rng(k) = ng1d(k)*rho(k)
         else
            rng(k) = R1
         endif

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte

! compute moments for snow

! calculate slope and intercept parameter

      dumLAMS = (CONS1*rns(K)/rs(K))**(1./DS)
      dumN0S = rns(K)*dumLAMS/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      IF (dumLAMS.LT.LAMMINS) THEN
      dumLAMS = LAMMINS
      dumN0S = dumLAMS**4*rs(K)/CONS1
      ELSE IF (dumLAMS.GT.LAMMAXS) THEN
      dumLAMS = LAMMAXS
      dumN0S = dumLAMS**4*rs(k)/CONS1
      end if

      ilams(k)=1./dumlams
      n0_s(k)=dumn0s

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

      do k = kte, kts, -1


! calculate slope and intercept parameter

      dumLAMg = (CONS2*rng(K)/rg(K))**(1./Dg)
      dumN0g = rng(K)*dumLAMg/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      IF (dumLAMg.LT.LAMMINg) THEN
      dumLAMg = LAMMINg
      dumN0g = dumLAMg**4*rg(K)/CONS2
      ELSE IF (dumLAMg.GT.LAMMAXg) THEN
      dumLAMg = LAMMAXg
      dumN0g = dumLAMg**4*rg(k)/CONS2
      end if

      ilamg(k)=1./dumlamg
      n0_g(k)=dumn0g

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept & slope values for rain.
!+---+-----------------------------------------------------------------+

      do k = kte, kts, -1

! calculate slope and intercept parameter

      dumLAMr = (PI*RHOW*rnr(K)/rr(K))**(1./3.)
      dumN0r = rnr(K)*dumLAMr/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      IF (dumLAMr.LT.LAMMINr) THEN
      dumLAMr = LAMMINr
      dumN0r = dumLAMr**4*rr(K)/(PI*RHOW)
      ELSE IF (dumLAMr.GT.LAMMAXr) THEN
      dumLAMr = LAMMAXr
      dumN0r = dumLAMr**4*rr(k)/(PI*RHOW)
      end if

      ilamr(k)=1./dumlamr
      n0_r(k)=dumn0r

      enddo

      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt. T_0) .and. (rr(k).gt. 0.001e-3) &
                   .and. ((rs(k+1)+rg(k+1)).gt. 0.01e-3) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*720.*ilamr(k)**7

         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (pi*rhosn/6./900.)*(pi*rhosn/6./900.) &
                                    * N0_s(k)*720.*ilams(k)**7
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (pi*rhog/6./900.)* (pi*rhog/6./900.)        &
                                    * N0_g(k)*720.*ilamg(k)**7
      enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.2) then
       do k = k_0-1, 1, -1

!..Reflectivity contributed by melting snow
          fmelt_s = DMIN1(1.0d0-rs(k)/rs(k_0), 1.0d0)
          if (fmelt_s.gt.0.01d0 .and. fmelt_s.lt.0.99d0 .and.           &
                         rs(k).gt.R1) then
           eta = 0.d0
           obms = 1./ds
           ocms = (1./(pi*rhosn/6.))**obms
           do n = 1, nbs
              x = pi*rhosn/6. * Dds(n)**3
              call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)* DEXP(-Dds(n)/ilams(k))
              eta = eta + f_d * CBACK * simpson(n) * dts(n)

           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif


!..Reflectivity contributed by melting graupel

          fmelt_g = DMIN1(1.0d0-rg(k)/rg(k_0), 1.0d0)
          if (fmelt_g.gt.0.01d0 .and. fmelt_g.lt.0.99d0 .and.           &
                         rg(k).gt.R1) then
           eta = 0.d0
           lamg = 1./ilamg(k)
           obmg = 1./dg
           ocmg = (1./(pi*rhog/6.))**obmg
           do n = 1, nbg
              x = pi*rhog/6. * Ddg(n)**3
              call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)* DEXP(-lamg*Ddg(n))
              eta = eta + f_d * CBACK * simpson(n) * dtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo


      end subroutine calc_refl10cm

 END MODULE module_mp_GRAUPEL
