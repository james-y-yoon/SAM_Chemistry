module micro_params
  !bloss: Holds microphysical parameter settings separately, so that they can be used
  !   inside both the microphysics module and the module(s) from WRF.

  logical :: doicemicro = .true. ! Turn on ice microphysical processes

  logical :: doaerosols = .false. ! Use Thompson-Eidhammer water- and ice-friendly aerosols

  logical :: doisotopes = .false. ! Enable water isotopologues

  !bloss: re-arrange rain-snow-snow-graupel-cloud ice collision processes
  !  so that each process has a single source and
  !  single destination.  Helpful in simplifying
  !  mass conservation and in defining transfer rates
  !  for water tracers or isotopologues.
  logical :: doRearrangeRainIceCollection = .false. 

  ! If using MSE chunk budgets, output budgets for water vapor and its isotopologues.
  logical :: do_vapor_chunk_budgets = .false. 

  real :: Nc0 = 100. ! initial/specified cloud droplet number concentration (#/cm3).

  !bloss(2021): Aerosol-aware Thompson-Eidhammer microphysics
  !  From p. 3638-3639 in Thompson & Eidhammer (2014), their activation scheme
  !  assumes that the underlying aerosol size distribution is lognormal with a 
  !  mean diameter of 80 nanometers and a geometric standard deviation of 1.8.
  !  The assumed hygroscopicity seems to be 0.4.
  !  Initial profile: nwfa(z) = CN_aloft + (CN_surf-CN_aloft)*exp(-z(k)/CN_DecayHeight)
  real :: CN_surf = 250.e6  ! Water-friendly aerosol number in #/kg at surface
  real :: CN_aloft = 50.e6  ! Water-friendly aerosol number in #/kg aloft
  real :: CN_DecayHeight = 800. ! DecayHeight in meters

  !  Profile: nifa(z) = IN_aloft + (IN_surf-IN_aloft)*exp(-z(k)/CN_DecayHeight)
  real :: IN_surf = 2.e6    ! Ice-friendly aerosol number in #/kg at surface
  real :: IN_aloft = 0.5e6  ! Ice-friendly aerosol number in #/kg aloft
  real :: IN_DecayHeight = 800. ! DecayHeight in meters

  ! option to allow the gamma exponent for rain, graupel and cloud ice to be specified.
  !   Note that mu=0 (the default) corresponds to an exponential distribution.
  real :: fixed_mu_r = 0.
  real :: fixed_mu_i = 0.
  real :: fixed_mu_g = 0.

  ! Fix the exponent in the gamma droplet size distribution for cloud liquid to a single value.
  logical :: dofix_mu_c = .false. 
  real :: fixed_mu_c = 10.3 ! fixed value from Geoffroy et al (2010).
  ! Full citation: Geoffroy, O., Brenguier, J.-L., and Burnet, F.:
  !   Parametric representation of the cloud droplet spectra for LES warm
  !   bulk microphysical schemes, Atmos. Chem. Phys., 10, 4835-4848,
  !   doi:10.5194/acp-10-4835-2010, 2010.

  !..Densities of rain, snow, graupel, and cloud ice. --> Needed in radiation
        REAL, PARAMETER :: rho_water = 1000.0
        REAL, PARAMETER :: rho_snow = 100.0
        REAL, PARAMETER :: rho_graupel = 500.0
        REAL, PARAMETER :: rho_cloud_ice = 890.0

  ! Fix rain number generation from melting snow (backported from WRF V3.9)      
  logical :: BrownEtAl2017_pnr_sml_fix = .true. 

  !bloss(2018-02): Enable choice of snow moment
  !parameterizations between Field et al (2005), the default,
  !and Field et al (2007).
  logical :: doFieldEtAl2007Snow = .false.
  
  ! Field et al (2007) has two snow size distributions: tropical and mid-latitude.
  !   The size distribution is used in the computation of sedimentation.
  logical :: TropicalSnow = .true. ! if false, use mid-latitude size distribution


  character*80 :: lookup_table_location = './RUNDATA/'

  !bloss(2021-03): Output individual process rates.
  !  Break them into groups:
  !     - mass and number tendencies
  !     - warm cloud tendencies (doicemicro=.false.)
  !     - cold cloud tendnecies (doicemicro=.true.)
  !     - AerosolAware tendnecies (doaerosols=.true.)
  logical :: do_output_process_rates = .false.
  integer ::  nproc_rates_mass = 1 , nproc_rates_number = 1 ! default value is one

  integer, parameter :: nmicro_process_rates_warm_mass = 4
  integer, parameter :: nmicro_process_rates_warm_number = 3
  integer, parameter :: nmicro_process_rates_warm_aerosol_number = 4

  integer, parameter :: nmicro_process_rates_cold_mass = 27
  integer, parameter :: nmicro_process_rates_cold_number = 14
  integer, parameter :: nmicro_process_rates_cold_aerosol_mass = 1 !pri_iha
  integer, parameter :: nmicro_process_rates_cold_aerosol_number = 8

end module micro_params
