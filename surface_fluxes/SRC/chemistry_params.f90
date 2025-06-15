module chemistry_params

   use chem_isoprene_Parameters, only: NVAR, NSPEC, NFIX
   use chem_isoprene_Monitor, only: SPC_NAMES

   implicit none

   logical :: do_iepox_droplet_chem=.false., do_iepox_aero_chem=.false., &
        hi_org=.true.

   real :: pHdrop=5.
   real :: pHaero=4.

   real :: deposition_rate = 0.01
   
   logical :: do_OH_diurnal = .true.
   real :: OH_night = 1.e5
   real :: OH_day_peak = 5.e6
   logical :: do_NO2_photolysis = .true.
   logical :: do_surface_Isoprene_diurnal = .true.

   logical, dimension(NVAR), public :: flag_gchemvar_out3D  ! which chem array to  output

   ! Define namelist variables
   character(len=15), dimension(NSPEC) :: gas_init_name  ! array of desired names
                                       !            for nonzero init
   real*8, dimension(NSPEC) :: gas_init_value      ! array of init values
                                               ! for corresponding gas_init_name
   character(len=15), dimension(NVAR) :: gas_out3D_name  ! array of desired gas names

   real :: p0 = 1013.25    ! pressure of 1atm in hPa
   real :: rhol = 1000.    ! kg/m^3


   
end module chemistry_params
