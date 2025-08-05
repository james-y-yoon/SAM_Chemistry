module chemistry
  
  use grid, only: nx, ny, nzm, nz, &
       dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
     z, zi, pres, adz, dz, dx, nx_gl, &
     time, dt, nstep, ncycle, nstat, nstatfrq, nrestart, day, &
     rank, dompi, masterproc, nsubdomains, save3Dbin, &
     case, caseid

  use vars, only: rho, t0, dtn, qcl, tabs0, pres, qv0
  use chem_isoprene_Parameters, only: NVAR, NSPEC, NFIX, ind_Isoprene, ind_IEPOX, ind_NO, ind_ISOPOOH, ind_OH, ind_ISOP1Nit, ind_ISOPDiNit, NREACT  ! NSPEC=NVAR+NFIX
  use chem_isoprene_Monitor, only: SPC_NAMES
  use chem_isoprene_Function, only: Fun
  use chem_isoprene_rate_functions, only: compute_rate_constants
  use chemistry_params, only: flag_gchemvar_out3D,  &
       gas_init_name, gas_init_value, gas_out3D_name, p0, rhol, &
       do_iepox_droplet_chem, do_iepox_aero_chem, hi_org, pHdrop, pHaero, &
       deposition_rate, do_OH_diurnal, OH_night, OH_day_peak, do_O3_photolysis, &
       do_NO2_photolysis, do_H2O2_photolysis, do_HCHO_photolysis, do_NO3_photolysis, &
       do_ISOPOOH_photolysis, do_surface_Isoprene_diurnal, do_only_tropospheric_chemistry

  use chem_aqueous, only: naqchem_fields, molwt, iIEPOX, iTETROL, iIEPOX_SO4,&
       iIP1NIT, iIPDINIT, aq_species_names, aq_gasprod_species_names, &
       flag_aqchemvar_out3D, flag_aqchemgasvar_out3D, &
       iepox_aqueous_tendencies, isop1nit_aqueous_tendencies
      
  use chem_aerosol, only:   iTETROLr, iIEPOX_SO4r,iepox_aero_transfer_rate, narchem_fields, molwt_ar, flag_archemvar_out3D, ar_species_names
  use chemistry_aux, only: get_OH
  
  use micro_params, only: MW_air, avgd, rho_aerosol, sigma_accum
  use microphysics, only: micro_field, iqcl, iqad, inad, incl, dBZ_cloudradar
  
  implicit none

  integer ngchem_fields, ngchem_fixed, ngchem_spec   ! equal to NVAR, NFIX, NSPEC  respectively

  integer nchem_fields_3Dsave 
  
  logical :: isallocatedCHEM = .false.

  real, allocatable, dimension(:,:,:,:) :: gchem_field  ! in ppv air
  real, allocatable, dimension(:,:,:,:) :: aqchem_field  ! in kg/kg
  real, allocatable, dimension(:,:,:,:) :: aqchem_gasprod_field  ! in kg/kg
  real, allocatable, dimension(:,:,:,:) :: archem_field  ! in kg/kg
  real, allocatable, dimension(:,:) :: gchem_profile_fixed ! in ppv 
  real, allocatable, dimension(:) :: M_profile  !  air density in molec/cm3
  real, allocatable, dimension(:,:) :: rate_const  ! array of gas reaction rate constants

  real, allocatable, dimension(:,:) :: gchem_horiz_mean_tend ! in ppv/s
  real, allocatable, dimension(:,:) :: aqchem_horiz_mean_tend ! in kg/kg/s
  real, allocatable, dimension(:,:) :: aqchem_gasprod_horiz_mean_tend ! in kg/kg/s
  real, allocatable, dimension(:,:) :: archem_horiz_mean_tend ! in kg/kg/s
  real, allocatable, dimension(:) :: g_depos_horiz_mean_tend_IEPOX ! in ppv/s
  real, allocatable, dimension(:) :: g_depos_horiz_mean_tend_ISOPOOH ! in ppv/s
  
  real, allocatable, dimension(:,:) :: & ! statistical arrays
       gchwle, &  ! resolved vertical flux
       gchadv, &  ! tendency due to vertical advection
       gchdiff, & ! tend. vertical diffusion
       gchwsb  ! SGS vertical flux

  real, allocatable, dimension(:,:) :: & ! statistical arrays
       aqchwle, &  ! resolved vertical flux
       aqchadv, &  ! tendency due to vertical advection
       aqchdiff, & ! tend. vertical diffusion
       aqchwsb  ! SGS vertical flux

  real, allocatable, dimension(:,:) :: & ! statistical arrays
       archwle, &  ! resolved vertical flux
       archadv, &  ! tendency due to vertical advection
       archdiff, & ! tend. vertical diffusion
       archwsb  ! SGS vertical flux

  real, allocatable, dimension(:,:,:) :: fluxbch, fluxtch  ! surface/top fluxes

  real gas_output_scale  ! convert all gas chem output to ppb
  real  gas_input_scale  ! convert all gas chem input to parts per unit air

  real, allocatable, dimension(:) :: Haq, NO3aq, SO4aq, HSO4aq      ! constant (for now) aqueous concentrations
  real, allocatable, dimension(:) :: Haero, SO4aero, HSO4aero      ! constant (for now) aerosol concentrations 
  real  actHaero   ! H activity in aerosols

   real OrgMF  ! Organic mass fraction of aerosol
   real FracTETROL  ! Fraction of IEPOXg to convert to TETROL on aerosol
   real FracIEPOX_SO4     ! Fraction of IEPOXg to convert to SO4 on aerosol
 
   CONTAINS

  subroutine chem_setparm()
    implicit none
   
    integer ierr, ios, ios_missing_namelist, place_holder
    ngchem_fields = NVAR    ! number of advected che fields
    ngchem_fixed = NFIX     ! number of fixed chem profiles
    ngchem_spec = NSPEC  ! = NVAR + NFIX
    
    gas_output_scale = 1.e9  ! convert all gas chem output to ppbv
    gas_input_scale = 1./gas_output_scale  ! convert all input from ppb to ppunit
    NAMELIST /CHEMISTRY/   do_iepox_droplet_chem, do_iepox_aero_chem, hi_org, pHdrop, pHaero, deposition_rate, do_OH_diurnal, OH_night, OH_day_peak, do_surface_Isoprene_diurnal, gas_init_name, gas_init_value, gas_out3D_name
        
    ! read in namelist
    NAMELIST /BNCUIODSBJCB/ place_holder
    open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
    read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
    rewind(55) !note that one must rewind before searching for new namelists
    read (55,CHEMISTRY,IOSTAT=ios)

    if (ios.ne.0) then
       if(ios.ne.ios_missing_namelist) then
           write(*,*) '****** ERROR: bad specification in CHEMISTRY namelist'
           rewind(55)
           read (55,CHEMISTRY) ! this should give a useful error message
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '****** No CHEMISTRY namelist in prm file *********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

    ! output namelist for documentation  
   if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.nml', form='formatted', position='append')    
      write (unit=55,nml=CHEMISTRY,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if
    
    ! allocate advection fields
   if(.not.isallocatedCHEM) then
       ! allocate isoprene gas chemistry fields and related variables
       allocate(gchem_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,NVAR))
       allocate(gchem_profile_fixed(nzm, ngchem_fixed), M_profile(nzm))
       allocate(gchem_horiz_mean_tend(nzm, NVAR))
       allocate(rate_const(nzm, NREACT))
       allocate(g_depos_horiz_mean_tend_IEPOX(nzm))
       allocate(g_depos_horiz_mean_tend_ISOPOOH(nzm))
       
       ! allocate aqueous IEPOX fields
       allocate(aqchem_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,naqchem_fields))
       allocate(aqchem_horiz_mean_tend(nzm, naqchem_fields))
       allocate(Haq(nzm), NO3aq(nzm), SO4aq(nzm), HSO4aq(nzm))
       allocate(Haero(nzm),SO4aero(nzm), HSO4aero(nzm))
       ! allocate gas product IEPOX fields
       allocate(aqchem_gasprod_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,naqchem_fields))
       allocate(aqchem_gasprod_horiz_mean_tend(nzm, naqchem_fields))

       ! allocate aerosol product fields
       allocate(archem_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,narchem_fields))
       allocate(archem_horiz_mean_tend(nzm, narchem_fields))
       
       allocate(gchwle(nz, ngchem_fields),gchadv(nz,ngchem_fields), &
            gchdiff(nz,ngchem_fields),gchwsb(nz,ngchem_fields))
       allocate(aqchwle(nz, naqchem_fields),aqchadv(nz,naqchem_fields), &
            aqchdiff(nz,naqchem_fields),aqchwsb(nz,naqchem_fields))
       ! Should have one here for gas products of aq chem FIX
       allocate(archwle(nz, narchem_fields),archadv(nz,narchem_fields), &
            archdiff(nz,narchem_fields),archwsb(nz,narchem_fields))

       allocate(fluxbch(nx, ny, ngchem_fields), fluxtch(nx, ny, ngchem_fields))
    end if   

    ! for now set these aqueous concentrations as constants, could be
    ! put into namelist input in future

    !Haq = 1.e-5 ! M H+ (pye et al)
    Haq = 10**(-pHdrop)
    NO3aq = 1.e-7 ! M NO3-
    SO4aq = 1.e-7 ! M SO4-2
    HSO4aq = 1.e-10 ! M HSO4-

    ! aerosol concentrations
    !Haero = 0.000038 ! M of H+ in aerosols based on RH = 0.8
    Haero = 10**(-pHaero)
    actHaero = 1.39  !  based on RH=0.8
    SO4aero = 3.8e-5  ! M of nucleophile SO4-2 
    HSO4aero = 6.3e-5  ! M of general acid HSO4-

    if (hi_org) then
       OrgMF = 0.85
    else
       OrgMF = 0.2
    end if
    
    FracTETROL = 0.85  ! portion of IEPOXg to convert to TETROL on aerosol 
    FracIEPOX_SO4 = 0.15     ! portion of IEPOXg to convert to SO4 on aerosol
    
end subroutine chem_setparm
  
subroutine chem_init()
  ! called at start of run or restart
  implicit none
  integer i,j,k
  integer v_selected ! index of namelist input variables
  integer v ! index of kpp variable
  logical match

  rate_const = 0.
  M_profile = 0.001 * RHO * avgd/ MW_air
  ! initialize gas fields
  gchem_profile_fixed = 0. 
  
  do v_selected = 1,ngchem_spec
     match = .false.
     do v=ngchem_fields+1, ngchem_spec  ! search only over fixed-variable names
        if(gas_init_name(v_selected)==trim(SPC_NAMES(v))) then
           match=.true.
           exit
        end if
     end do   
     if(match) then
        gchem_profile_fixed(:, v-ngchem_fields) = gas_init_value(v_selected)* &
             gas_input_scale
        if (masterproc) then
           write(*,*) 'SET FIXED CHEM PROFILE: ', gas_init_name(v_selected),                 gas_init_value(v_selected)
        end if   

     end if
    
  end do   

  if(nrestart.eq.0) then
     ! initialize gas chem fields
     gchem_field = 0.
     ! compute conversion profile 
     do v_selected = 1,ngchem_spec
        match = .false.   
        do v = 1,ngchem_fields   ! search only over variable names  
           if(gas_init_name(v_selected)==trim(SPC_NAMES(v))) then
              match=.true.
              exit
           end if
        end do  
        if(match) then
           do i = 1, nx
              do j = 1, ny
                 gchem_field(i,j, :, v) = gas_init_value(v_selected)*gas_input_scale
              end do
           end do   
        end if
        ! if (v_selected.eq.1.and.masterproc) write(*,*) 'Initialized ISOP102 to ', gchem_field(1,1,0,1)
      
     end do
     ! initialize aqchem and archem fields as zero
     aqchem_field = 0.
     aqchem_gasprod_field = 0.
     archem_field = 0.
  end if  ! restart       
  
  ! set flags for 3d output based on namelist input 
  flag_gchemvar_out3D(:)=.false.
  nchem_fields_3Dsave = 0.
  do v = 1,ngchem_fields
     if(any(gas_out3D_name==trim(SPC_NAMES(v)))) then
        flag_gchemvar_out3D(v) = .true.
        nchem_fields_3Dsave = nchem_fields_3Dsave + 1
        if (masterproc) write(*,*) &
             'Chem 3d output field added: ', trim(SPC_NAMES(v))
     end if
  end do

  ! for now save all aqeous variables and aerosol product fields
  flag_aqchemvar_out3D(:) = .true.
  flag_aqchemgasvar_out3D(:) = .true.
  flag_archemvar_out3D(:) = .true.
  nchem_fields_3Dsave = nchem_fields_3Dsave + naqchem_fields*2 + narchem_fields
  
  ! initialize some statistics profiles, not output yet
  gchwle = 0.
  gchadv = 0.
  gchdiff = 0.
  gchwsb = 0.

  aqchwle = 0.
  aqchadv = 0.
  aqchdiff = 0.
  aqchwsb = 0.

  archwle = 0.
  archadv = 0.
  archdiff = 0.
  archwsb = 0.

  ! top and bottom fluxes of fields
  fluxbch = 0.
  fluxtch = 0.

end subroutine chem_init  

subroutine chem_hbuf_init(namelist,deflist,unitlist,status,average_type,count,chemcount)

  
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*), average_type(*), count,chemcount, n

character*16 name
character*16 tend_name
character*80 longname
character*10 units

chemcount = 0  

do n = 1,ngchem_fields
    name = trim(SPC_NAMES(n))
    longname = trim(SPC_NAMES(n))
    units = 'ppbv'
    call add_to_namelist(count,chemcount,name,longname,units,0)

    longname = trim(SPC_NAMES(n)) // ' tendency due to gas reaction'
    units = 'ppbv/s'
    tend_name = trim(SPC_NAMES(n))//'+'
    call add_to_namelist(count,chemcount,tend_name,longname,units,0)
end do

do n = 1,ngchem_fixed
   name = trim(SPC_NAMES(n+ngchem_fields))
   longname = name // ' Fixed Species'
   units = 'ppbv'
   call add_to_namelist(count,chemcount,name,longname,units,0)
end do

do n = 1,naqchem_fields
    name = trim(aq_species_names(n))
    longname = trim(aq_species_names(n))
    units = 'kg/kg'
    call add_to_namelist(count,chemcount,name,longname,units,0)

    longname = trim(aq_species_names(n)) // ' tendency due to aqueous reaction'
    units = 'kg/kg/s'
    call add_to_namelist(count,chemcount,trim(aq_species_names(n))//'+',longname,units,0)
end do

do n = 1,naqchem_fields
    name = trim(aq_gasprod_species_names(n))
    longname = trim(aq_gasprod_species_names(n))
    units = 'ppbv'
    call add_to_namelist(count,chemcount,name,longname,units,0)

    longname = trim(aq_gasprod_species_names(n)) // 'gas tendency due to aqueous reaction'
    units = 'ppbv/s'
    call add_to_namelist(count,chemcount,trim(aq_gasprod_species_names(n))//'+',longname,units,0)
end do
    
do n = 1,narchem_fields
    name = trim(ar_species_names(n))
    longname = trim(ar_species_names(n))
    units = 'kg/kg'
    call add_to_namelist(count,chemcount,name,longname,units,0)

    longname = trim(ar_species_names(n)) // ' tendency due to aerosol chemistry'
    units = 'kg/kg/s'
    call add_to_namelist(count,chemcount,trim(ar_species_names(n))//'+',longname,units,0)
end do

    name = 'IPOOHd+'
    longname = 'IPOOH gas tendency due to dry deposition'
    units = 'ppbv/s'
    call add_to_namelist(count, chemcount, name, longname, units, 0)

    name = 'IEPOXd+'
    longname = 'IEPOX gas tendence due to dry deposition'
    units = 'ppbv/s'
    call add_to_namelist(count, chemcount, name, longname, units, 0)


if(masterproc) then
   write(*,*) 'Added ', chemcount, ' arrays to statistics for gas and aqueous chemical species'
end if

end subroutine chem_hbuf_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START CHEMISTRY FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
subroutine chem_flux()  
    use grid, only : day, time, dt, nstep
    use vars, only : precsfc, tabs, interactive_soil_wetness    ! precsfc (x,y), tabs (x,y,z) (JY), and interactive_soil_wetness (only used in chem_flux) (JY)
    use rad, only : swDownSurface, lwDownSurface, insolation_TOA         ! eventually converted to PAR

    implicit none

    !! For non-MEGAN ISOP diurnal cycle
    real :: t_solar_peak = 0.167  ! days
    real :: pi = 3.1415927
    real :: frac_of_peak

    ! Average fluxes
    real :: Iso_avg_flux = 1.e-9   ! converted from (ppb m/s) to  m/s 
    real :: conversion_between_isoprene_and_nox = 0.0072    ! to scale average isoprene flux

    ! Conversion factors
    real :: J_to_mol = 4.6                                           ! Approximate unit conversion between W m-2 and umol photons m-2 s-1
    real :: frac_PAR = 0.5                                           ! Fraction of downward shortwave that is PAR

    ! Allocatable arrays for emissions
    real, allocatable, dimension(:,:) :: soil_NOx ! soil NOx array (x,y)
    real, allocatable, dimension(:,:) :: ppfd ! to calculate PAR-dependence of isoprene emissions
    real, allocatable, dimension(:,:) :: radiation_activity_factor ! to calculate PAR-dependence of isoprene emissions
    real, allocatable, dimension(:,:) :: temperature_activity_factor ! to calculate PAR-dependence of isoprene emissions

    ! Other variables
    real :: LDF_i = 1                                                ! Set for isoprene

    ! Boolean controls
    logical :: do_megan_isoprene = .true.
    logical :: do_bdsnp_no = .true.
    logical :: do_lightning = .true.

    logical :: do_surface_Isoprene_diurnal = .false.

    if ( do_lightning ) then
        call lightning_decaria_ctg()
        ! call lightning_decaria_ic()
    endif

    !!!!! END LIGHTNING ATTEMPT !!!!!!

    ! Default values for emissions
    fluxbch = 0.
    fluxtch = 0.
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Isoprene Flux Module !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( do_megan_isoprene ) then 

        allocate(radiation_activity_factor(nx, ny), temperature_activity_factor(nx, ny), ppfd(nx, ny))
        radiation_activity_factor = 0.
        temperature_activity_factor = 0.
        ppfd = 0.

        !! MEGAN Isoprene Radiation Module !!
        ! Convert radiation from W m-2 to umol m-2 s-1 (https://search.r-project.org/CRAN/refmans/bigleaf/html/Rg.to.PPFD.html)
        if ( ALLOCATED(swDownSurface) ) then
            ! print*, "shortwave=", SUM(swDownSurface) / SIZE(swDownSurface)
            ! print*, "longwave=", SUM(lwDownSurface) / SIZE(lwDownSurface)
            ! print*, "clearsky shortwave=", SUM(insolation_TOA) / SIZE(insolation_TOA)
            ppfd(:,:) = swDownSurface(:,:) * J_to_mol * frac_PAR
            radiation_activity_factor(:,:) = calculate_megan_BVOC_radiation(ppfd(:,:), LDF_i)
            ! print*, "radiation_activity_factor=", SUM(radiation_activity_factor) / SIZE(radiation_activity_factor)
        endif

        !! MEGAN Isoprene Temperature Module !!
        ! Calculate temperature activity factor for isoprene
        temperature_activity_factor = calculate_megan_BVOC_temperature(tabs(:,:,1), LDF_i)
        ! print*, "temperature_activity_factor=", SUM(temperature_activity_factor) / SIZE(temperature_activity_factor)

        fluxbch(:,:,ind_Isoprene) =  Iso_avg_flux * temperature_activity_factor * radiation_activity_factor     ! Calculate isoprene fluxes using MEGAN
        deallocate(temperature_activity_factor, ppfd, radiation_activity_factor)

    elseif ( do_surface_Isoprene_diurnal ) then
        ! print*, 'day=', day, 'time=', time 
        frac_of_peak = MAX(0., cos((day - t_solar_peak)*2*pi))
        ! print *, 'frac_of_peak = ', frac_of_peak
        fluxbch(:,:,ind_Isoprene) =  frac_of_peak * Iso_avg_flux * (pi/2.)
    
    else
        fluxbch(:,:,ind_Isoprene) = Iso_avg_flux
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Soil NO Flux Module !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( do_bdsnp_no ) then
        allocate(soil_NOx(nx, ny))
        soil_NOx = 0.

        call calculate_bdsnp_NO(soil_NOx)
        fluxbch(:, :, ind_NO) = ( Iso_avg_flux * conversion_between_isoprene_and_nox ) * soil_NOx(:, :)
        deallocate(soil_NOx)

        ! print*, "nox=", SUM(fluxbch(:, :, ind_NO)) / SIZE(fluxbch(:, :, ind_NO))
    endif  
end subroutine chem_flux

subroutine lightning_decaria_ctg()
    ! Parametrizes (as close as possible!) to DeCaria et al. (2000), A cloud-scale model study of lightning-generated NOx...

    use grid, only : dt, nstep, nx, ny, nz, pres, dz                        ! pres is in mbar!
    use vars, only : tabs ! tabs (x,y,z) contains absolute temperature in K
    use params, only : land, ocean
    use microphysics, only: dBZ_cloudradar
    implicit none

    ! General parameters for both cloud-to-ground (CTG) and intracloud (IC) lightning
    real :: time_between_flashes = 180                                      ! 3 minutes
    integer :: lightning_time_step                                          ! How many time steps to skip before calculating lightning (function of dt)
    character(len = 15) :: progress_string = "CTG LIGHTNING!"
    real :: pi = 3.1415927
    real :: universal_gas_constant = 8.314                                  ! gas constant, in J mol-1 K-1
    real :: moist_adiabatic_lapse_rate = 5e-3                               ! an estimate of the moist adiabatic lapse rate, in C / m

    ! Cloud-to-ground lightning
    real :: cloud_to_ground_isotherm = -15                                  ! isoterm of the Gaussian mean for CTG lightning, in Celsius
    integer :: vertical_index_for_isotherm                                  ! index of the CTG mean isotherm
    real :: mu
    real :: std_dev                                                         ! in meters
    real :: x_area_of_storm                                                 ! area the storm takes up
    real :: temperature_of_closest_altitude

    ! Production rates
    real :: no_production_per_flash = 460                                   ! CTG NO production per flash, in mol/flash (460 is CTG from Price & Rind)

    real, allocatable, dimension(:) :: vertical_function_profile            ! f(x) for CTG in DeCaria
    real :: integral_of_vertical_function                                   ! Used to get LC from N_tot

    real, allocatable, dimension(:,:,:) :: change_in_mixing_ratio                   ! delta q_NO(z) in DeCaria

    logical :: price_and_rind = .false.                                       ! horizontal area using cloud top height
    logical :: decaria_reflectivity = .true.                                       ! horizontal area using cloud top height

    real :: cloud_top_height_threshold = 1300                                ! in meters
    real, allocatable, dimension(:, :) :: cloud_top_heights                  ! for use in Price & Rind parametrization
    real, allocatable, dimension(:, :) :: cloud_top_temps                    ! for use in Price & Rind parametrization
    real, allocatable, dimension(:, :) :: price_and_rind_flash_rates         ! for use in Price & Rind parametrization
    real, allocatable, dimension(:) :: number_of_20dbz_per_altitude

    integer :: i, j, k                                                          ! Counter variables

    lightning_time_step = time_between_flashes / dt

    ! If the desired time between flashes has passed, then do lightning
    if ( mod(nstep, lightning_time_step) == 0 ) then
        allocate(vertical_function_profile(nz), change_in_mixing_ratio(nx, ny, nz))
        change_in_mixing_ratio = 0.

        print*, "**************** Time for ", progress_string

        do i = 1, nx
            do j = 1, ny 

                ! Gives the closest to the isotherm, but there is no guarantee that this is close enough
                vertical_index_for_isotherm = minloc(abs((tabs(i,j,:) - 273.15) - cloud_to_ground_isotherm), dim = 1) ! dim=1 is required to return an integer, not an 1-element array

                mu = z(vertical_index_for_isotherm)
                temperature_of_closest_altitude = tabs(i, j, vertical_index_for_isotherm) - 273.15 

                ! If the closest temperature is farther than 1.5 deg C away from the isotherm, extrapolate using the MALR
                if ( abs( temperature_of_closest_altitude - cloud_to_ground_isotherm ) > 1.5 ) then
                    mu = mu - ( cloud_to_ground_isotherm - temperature_of_closest_altitude ) / moist_adiabatic_lapse_rate
                endif
                
                std_dev = mu / 3.
                vertical_function_profile = 1 / (sqrt(2 * pi) * std_dev) * exp(-1 * (z - mu)**2 / (2 * std_dev**2))
                
                ! Integration for denominator
                integral_of_vertical_function = 0.
                do k = 1, nz
                    integral_of_vertical_function = integral_of_vertical_function + ( vertical_function_profile(k) * pres(k) * 100 ) * dz           ! 100 is to convert from mbar to Pa
                enddo
                
                if ( integral_of_vertical_function > 0 ) then
                        change_in_mixing_ratio(i,j,:) = ( no_production_per_flash / integral_of_vertical_function ) * universal_gas_constant * vertical_function_profile * tabs(i, j, :)
                endif
                ! print*, "NO", SUM(gchem_field(i, j, :, ind_NO)) / SIZE(gchem_field(i, j, :, ind_NO))
            enddo
        enddo

        if ( price_and_rind ) then
            allocate(cloud_top_heights(nx, ny), cloud_top_temps(nx, ny), price_and_rind_flash_rates(nx, ny))
            
            cloud_top_heights = 0.
            cloud_top_temps = 0.

            call calculate_cloud_top_height(cloud_top_heights, cloud_top_temps)
            print*, "cloud_top_height=", MAXVAL(cloud_top_heights)
            print*, "cloud top temp=", MINVAL(cloud_top_temps)

            x_area_of_storm = 0.
            price_and_rind_flash_rates = 0.

            do i = 1, nx
                do j = 1, ny 
                    if ( cloud_top_heights(i, j) <= cloud_top_height_threshold ) then
                        cloud_top_heights(i, j) = 0.
                    else
                        x_area_of_storm = x_area_of_storm + 1

                        if ( land ) then
                                price_and_rind_flash_rates(i,j) = 3.44e-5 * ( cloud_top_heights(i,j) / 1000 )**4.9
                        elseif ( ocean ) then
                                price_and_rind_flash_rates(i,j) = 6.4e-4 * ( cloud_top_heights(i,j) / 1000 )**1.73
                        endif
                    endif
                end do
            end do

            do i = 1, nx
                do j = 1, ny 
                    if ( cloud_top_heights(i, j) < 1 ) then 
                        change_in_mixing_ratio(i,j,:)= 0.
                    else
                        do k = 1, nz
                            change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( x_area_of_storm * dx * dx ) ! eventually replace dx with dy
                            change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) * price_and_rind_flash_rates(i,j) / sum( price_and_rind_flash_rates(:,:) )
                        enddo
                    endif
                 enddo
             enddo

             print*, "maximum lightning***************", MAXVAL(change_in_mixing_ratio)

            gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
         endif



         if ( decaria_reflectivity ) then
            allocate(number_of_20dbz_per_altitude(nzm))
            number_of_20dbz_per_altitude = 0.

            do k = 1,nzm 
               do i = 1,nx
                  do j = 1,ny
                     if ( dBZ_cloudradar(i, j, k) >= 20 ) then
                        number_of_20dbz_per_altitude(k) = number_of_20dbz_per_altitude(k) + 1
                     else
                        change_in_mixing_ratio(i,j,k) = 0.
                     endif
                  enddo
               enddo
            enddo

            do i = 1, nx
                do j = 1, ny 
                     do k = 1, nzm
                        if ( number_of_20dbz_per_altitude(k) > 0. ) then 
                           change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( number_of_20dbz_per_altitude(k) * dx * 20000. * 100 ) ! eventually replace dx with dy
                        endif
                     enddo
                 enddo
             enddo

             print*, "maximum lightning***************", MAXVAL(change_in_mixing_ratio)

            gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
        endif
    endif
end subroutine lightning_decaria_ctg

subroutine lightning_decaria_ic()
    ! Parametrizes (as close as possible!) to DeCaria et al. (2000), A cloud-scale model study of lightning-generated NOx...

    use grid, only : dt, nstep, nx, ny, nz, pres, dz                        ! pres is in mbar!
    use vars, only : tabs, qcl, qpl, qci, qpi                                                   ! tabs (x,y,z) contains absolute temperature in K
    use params, only : land, ocean
    implicit none

    ! General parameters for both cloud-to-ground (CTG) and intracloud (IC) lightning
    real :: time_between_flashes = 180                                      ! 3 minutes
    integer :: lightning_time_step                                          ! How many time steps to skip before calculating lightning (function of dt)

    character(len = 15) :: progress_string = "IC LIGHTNING!"
    real :: pi = 3.1415927
    real :: universal_gas_constant = 8.314                                  ! gas constant, in J mol-1 K-1
    real :: moist_adiabatic_lapse_rate = 5e-3                               ! an estimate of the moist adiabatic lapse rate, in C / m

    ! Cloud-to-ground lightning
    real :: ic_isotherm_bottom = -15                                  ! isoterm of the Gaussian mean for CTG lightning, in Celsius
    integer :: vertical_index_for_isotherm_bottom                                  ! index of the CTG mean isotherm

    real :: ic_isotherm_top = -45                                  ! isoterm of the Gaussian mean for CTG lightning, in Celsius
    integer :: vertical_index_for_isotherm_top                                  ! index of the CTG mean isotherm

    real :: mu_bottom
    real :: mu_top

    real :: std_dev_bottom                                                         ! in meters
    real :: std_dev_top                                                         ! in meters

    real, allocatable, dimension(:) :: x_area_of_storm                                                 ! area the storm takes up
    real :: temperature_of_closest_altitude

    ! Production rates
    real :: no_production_per_flash = 460                                   ! CTG NO production per flash, in mol/flash (460 is CTG from Price & Rind)

    real, allocatable, dimension(:) :: vertical_function_profile            ! f(x) for CTG in DeCaria
    real :: integral_of_vertical_function                                   ! Used to get LC from N_tot

    real, allocatable, dimension(:,:,:) :: change_in_mixing_ratio                   ! delta q_NO(z) in DeCaria
    logical :: decaria_ic = .true.                                       ! horizontal area using cloud top height

    integer :: i, j, k                                                          ! Counter variables

    lightning_time_step = time_between_flashes / dt

    ! If the desired time between flashes has passed, then do lightning
    if ( mod(nstep, lightning_time_step) == 0 ) then
        allocate(vertical_function_profile(nz), change_in_mixing_ratio(nx, ny, nz))
        change_in_mixing_ratio = 0.

        print*, "**************** Time for ", progress_string

        do i = 1, nx
            do j = 1, ny 

                ! Gives the closest to the isotherm, but there is no guarantee that this is close enough
                vertical_index_for_isotherm_bottom = minloc(abs((tabs(i,j,:) - 273.15) - ic_isotherm_bottom), dim = 1) ! dim=1 is required to return an integer, not an 1-element array

                mu_bottom = z(vertical_index_for_isotherm_bottom)
                temperature_of_closest_altitude = tabs(i, j, vertical_index_for_isotherm_bottom) - 273.15 

                ! If the closest temperature is farther than 1.5 deg C away from the isotherm, extrapolate using the MALR
                if ( abs( temperature_of_closest_altitude - ic_isotherm_bottom ) > 1.5 ) then
                    mu_bottom = mu_bottom - ( ic_isotherm_bottom - temperature_of_closest_altitude ) / moist_adiabatic_lapse_rate
                endif
                std_dev_bottom = mu_bottom / 3.


                ! Gives the closest to the isotherm, but there is no guarantee that this is close enough
                vertical_index_for_isotherm_top = minloc(abs((tabs(i,j,:) - 273.15) - ic_isotherm_top), dim = 1) ! dim=1 is required to return an integer, not an 1-element array

                mu_top = z(vertical_index_for_isotherm_top)
                temperature_of_closest_altitude = tabs(i, j, vertical_index_for_isotherm_top) - 273.15 

                ! If the closest temperature is farther than 1.5 deg C away from the isotherm, extrapolate using the MALR
                if ( abs( temperature_of_closest_altitude - ic_isotherm_top ) > 1.5 ) then
                    mu_top = mu_top - ( ic_isotherm_top - temperature_of_closest_altitude ) / moist_adiabatic_lapse_rate
                endif
                std_dev_top = std_dev_bottom / 3.

                vertical_function_profile = ( (0.8 / (sqrt(2 * pi) * std_dev_bottom) * exp(-1 * (z - mu_bottom)**2 / (2 * std_dev_bottom**2))) + (1 / (sqrt(2 * pi) * std_dev_top) * exp(-1 * (z - mu_top)**2 / (2 * std_dev_top**2))) )
                
                ! Integration for denominator
                integral_of_vertical_function = 0.
                do k = 1, nz
                    integral_of_vertical_function = integral_of_vertical_function + ( vertical_function_profile(k) * pres(k) * 100 ) * dz           ! 100 is to convert from mbar to Pa
                enddo
                
                if ( integral_of_vertical_function > 0 ) then
                        change_in_mixing_ratio(i,j,:) = ( no_production_per_flash / integral_of_vertical_function ) * universal_gas_constant * vertical_function_profile * tabs(i, j, :)
                endif
            enddo
        enddo

         if ( decaria_ic ) then
            allocate(x_area_of_storm(nzm))
            x_area_of_storm = 0.

            print*, "check this is nonzero=", MAXVAL(qcl)
            print*, "similarly=", MAXVAL(qpl)

            do k = 1,nzm 
               do i = 1,nx
                  do j = 1,ny
                     if ( (qcl(i, j, k) + qpl(i, j, k) + qci(i, j, k) + qpi(i, j, k)) * 1000. > 0.01 ) then
                        x_area_of_storm(k) = x_area_of_storm(k) + 1
                     else
                        change_in_mixing_ratio(i,j,k) = 0.
                     endif
                  enddo
               enddo
            enddo

            do i = 1, nx
                do j = 1, ny 
                     do k = 1, nzm
                        if ( x_area_of_storm(k) > 0. ) then 
                           change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( x_area_of_storm(k) * dx * 20000. ) ! eventually replace dx with dy
                        endif
                     enddo
                 enddo
             enddo

             print*, "maximum IC lightning***************", MAXVAL(change_in_mixing_ratio)

            gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
        endif
    endif
end subroutine lightning_decaria_ic

subroutine calculate_cloud_top_height(cloud_top_heights, cloud_top_temps)
    ! Required because diagnose is run after chem_flux

    use vars, only : qcl, qci, rho, tabs
    use grid, only : adz, dz, nx, ny, nzm, z

    integer :: i, j, k                                                          ! Counter variables
    real :: tmp_lwp                                                             ! Temporary variable
    real, allocatable, dimension(:,:) :: cloud_top_heights
    real, allocatable, dimension(:,:) :: cloud_top_temps

    cloud_top_heights = 0.
    cloud_top_temps = 0.

    do j = 1,ny
       do i = 1,nx
          tmp_lwp = 0.
          do k = nzm,1,-1

             tmp_lwp = tmp_lwp + (qcl(i,j,k) + qci(i,j,k)) * rho(k) * dz * adz(k)
             
             if (tmp_lwp.gt.0.01) then
                cloud_top_heights(i,j) = z(k)
                cloud_top_temps(i,j) = tabs(i,j,k)
                exit
             end if
          end do
       end do
    end do
end subroutine calculate_cloud_top_height


subroutine calculate_bdsnp_NO(soil_NOx)
    ! Calculates soil NOx
    ! See Hudman et al. (2012) and Wang et al. (2021) for parametrization
    use vars, only : precsfc, tabs, interactive_soil_wetness ! precsfc (x,y), tabs (x,y,z) (JY), and interactive_soil_wetness (only used in chem_flux) (JY)

    implicit none

    ! BDSNP Parameters
    real :: a_bdsnp = 1.65
    real :: b_bdsnp = 3.3
    real :: temperature_in_celsius

    real, allocatable, dimension(:,:) :: soil_NOx                 ! temperature activity factor, light-dependent
    real, dimension(nx, ny) :: surface_temperature                 ! temperature at bottommost layer

    ! For interactive soil moisture
    real :: decay_in_soil_moisture = 0.00003                ! From MERRA2 Regressions
    real :: increase_in_soil_moisture = 8                   ! From MERRA2 Regressions

    ! Conversion factors in soil moisture calculations
    integer :: scale_between_kg_m2_s1_and_mm_day = 86400    ! MERRA2 precipitation is in kg m-2 s-1; precipitation in SAM is mm/day (?)
    real :: conversion_between_hour_and_second = 3600       ! MERRA2 is on an hourly time grid versus seconds for SAM

    integer :: i, j                                                  ! Counter variables

    ! Loop through the horizontal axes
    do i = 1, nx
        do j = 1, ny

            ! If there is precipitation, increase soil moisture
            if ( precsfc(i, j) > 0 ) then
                interactive_soil_wetness(i, j) = interactive_soil_wetness(i, j) + (increase_in_soil_moisture * dt / scale_between_kg_m2_s1_and_mm_day * precsfc(i, j))
            
            ! If there is no precipitation, decay the soil moisture
            else
                interactive_soil_wetness(i, j) = interactive_soil_wetness(i, j) - (decay_in_soil_moisture * dt / conversion_between_hour_and_second)
            endif

            !! Now, calculate soil NOx !!
            temperature_in_celsius = tabs(i, j, 1) - 273.15

            if ( temperature_in_celsius < 20 ) then
                soil_NOx(i, j) = exp( 0.103 * temperature_in_celsius ) * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
            elseif ( ( temperature_in_celsius >= 20 ) .and. ( temperature_in_celsius <= 40 ) ) then
                soil_NOx(i, j) = ( -0.009 * temperature_in_celsius**3 + 0.837 * temperature_in_celsius**2 - 22.52 * temperature_in_celsius + 196.149 ) * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
            else
                soil_NOx(i, j) = 58.269 * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
            endif
        end do
    end do
end subroutine calculate_bdsnp_NO

function calculate_megan_BVOC_temperature(surface_temperature, LDF_i) result(gamma_T)
    ! Returns the value of the MEGAN scaling factor (gamma_T) for temperature
    ! See Guenther et al. (2012) for parametrization
    implicit none

    ! parameters
    real, dimension(nx, ny) :: surface_temperature                   ! temperature at bottommost layer
    real :: LDF_i                                                    ! light-dependent fraction, = 1 for isoprene

    real, allocatable, dimension(:,:) :: gamma_T_ldf                 ! temperature activity factor, light-dependent
    real, allocatable, dimension(:,:) :: gamma_T_lif                 ! temperature activity factor, light-independent
    real, allocatable, dimension(:,:) :: gamma_T                     ! final temperature activity factor, to return

    ! MEGAN - light-independent fraction
    real :: beta_i = 0.13                                            ! empirically determined coefficient for each VOC, tuned to ISOP
    real :: T_s = 297                                                ! standard temperature conditions for leaf temperature [K]

    ! MEGAN - light-dependent fraction
    real :: C_eo = 2                                                 ! Changes with species (currently set at isoprene)
    real :: C_t1 = 95                                                ! Changes with species (currently set at isoprene)
    real :: C_t2 = 230                                               ! Empirical coefficient
    real :: T_24 = 300                                               ! Average leaf temperature of past 24 hours
    real :: T_240 = 300                                              ! Average leaf temperature of past 240 hours

    ! intermediate calculations
    real, allocatable, dimension(:,:) :: T_opt
    real, allocatable, dimension(:,:) :: E_opt
    real, allocatable, dimension(:,:) :: x

    allocate(gamma_T_ldf(nx, ny), gamma_T_lif(nx, ny), gamma_T(nx, ny))
    allocate(T_opt(nx, ny), E_opt(nx, ny), x(nx, ny))

    gamma_T_lif = exp(beta_i * (surface_temperature(:,:) - T_s))                                ! Light-independent -- similar to monoterpene flux

    T_opt = 313 + (0.6 * (T_240 - T_s))                                                         ! Optimal temperature
    E_opt = C_eo * exp(0.05 * (T_24 - T_s)) * exp(0.05 * (T_240 - T_s))
    x = ((1 / T_opt) - (1 / surface_temperature)) / 0.00831

    gamma_T_ldf = E_opt * (C_t2 * exp(C_t1 * x) / (C_t2 - C_t1 * (1 - exp(C_t2 * x))))          ! light dependent emission activity factor
    gamma_T = (1 - LDF_i) * gamma_T_lif + (LDF_i * gamma_T_ldf)                                 ! MEGAN emission activity factor, accounts for light dependent and light independent factors    
end function calculate_megan_BVOC_temperature

function calculate_megan_BVOC_radiation(ppfd, LDF_i) result(gamma_p)
    ! Returns the value of the MEGAN scaling factor (gamma_P) for temperature
    ! See Guenther et al. (2012) for parametrization
    implicit none

    ! constants
    real :: p_s = 200               ! Standard conditions for PPFD, equal to 200 umol m-2 s-1 for sunlit leaves, 50 for shaded leaves
    real :: p_24 = 600              ! average PPFD of past 24 hours
    real :: p_240 = 600             ! average PPFD of past 240 hours

    ! intermediate calculations
    real :: c_p 
    real :: alpha
    real, allocatable, dimension(:,:) :: gamma_p_ldf
    real, allocatable, dimension(:,:) :: gamma_p

    ! parameters
    real :: ppfd(:,:)                ! photosynthetic photon flux density, in umol m-2 s-1
    real :: LDF_i                    ! light-dependent fraction, = 1 for isoprene

    allocate(gamma_p_ldf(nx, ny))
    allocate(gamma_p(nx, ny))

    alpha = 0.004 - ( 0.0005 * log(p_240) )
    c_p = 0.0468 * exp(0.0005 * ( p_24 - p_s )) * p_240**(0.6)

    gamma_p_ldf = c_p * ((alpha * ppfd(:,:)) / (1 + (alpha**2  * ppfd(:,:)**2))**0.5)       ! light dependent emission activityfactor
    gamma_p = (1 - LDF_i) + ( LDF_i * gamma_p_ldf )                                         ! MEGAN emission activity factor, accounts for light dependent and light independent factors    
end function calculate_megan_BVOC_radiation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END CHEMISTRY FLUXES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine do_photolysis
   use grid, only : time, z, dz, nz, dtn
   ! use chemistry_aux, only : get_NO2_photolysis_rate
   use chem_Isoprene_Parameters, only : ind_NO2, ind_O3, ind_OH, ind_H2O2, ind_HCHO, ind_HO2, ind_CO, ind_NO, ind_NO3, ind_ISOPOOH
   use rad, only : insolation_TOA, swDownSurface

   implicit none

   real, allocatable, dimension(:,:,:) :: photolysis_rate
   real, allocatable, dimension(:,:,:) :: change_in_mixing_ratio

   integer :: i, j, k, n

   allocate(change_in_mixing_ratio(nz, ny, nzm))
   change_in_mixing_ratio = 0.
   photolysis_rate = 0.

   if (do_NO2_photolysis) then
      photolysis_rate = calculate_photolysis_rate(1.2E-2)
      change_in_mixing_ratio = gchem_field(:,:,:,ind_NO2) * ( photolysis_rate * dtn )
      gchem_field(:,:,:,ind_NO2) = gchem_field(:,:,:,ind_NO2) - change_in_mixing_ratio
      gchem_field(:,:,:,ind_NO) = gchem_field(:,:,:,ind_NO) + change_in_mixing_ratio
      gchem_field(:,:,:,ind_O3) = gchem_field(:,:,:,ind_O3) + change_in_mixing_ratio
   end if 

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! SECOND INSTANCE", gchem_field(i,j,k,n)
                  print*, "PHOTOLYSIS RATE=", photolysis_rate(i,j,k)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

   if (do_NO3_photolysis) then
      photolysis_rate = calculate_photolysis_rate(2.8E-4)
      change_in_mixing_ratio = gchem_field(:,:,:,ind_NO3) * ( photolysis_rate * dtn )
      gchem_field(:,:,:,ind_NO3) = gchem_field(:,:,:,ind_NO3) - change_in_mixing_ratio
      gchem_field(:,:,:,ind_NO2) = gchem_field(:,:,:,ind_NO2) + change_in_mixing_ratio
   end if 

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! THIRD INSTANCE", gchem_field(i,j,k,n)
                  print*, "PHOTOLYSIS RATE=", photolysis_rate(i,j,k)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

   if (do_ISOPOOH_photolysis) then
      photolysis_rate = calculate_photolysis_rate(9.9E-6)
      change_in_mixing_ratio = gchem_field(:,:,:,ind_ISOPOOH) * ( photolysis_rate * dtn )
      gchem_field(:,:,:,ind_ISOPOOH) = gchem_field(:,:,:,ind_ISOPOOH) - change_in_mixing_ratio
      gchem_field(:,:,:,ind_OH) = gchem_field(:,:,:,ind_OH) + 2 * change_in_mixing_ratio
      gchem_field(:,:,:,ind_HCHO) = gchem_field(:,:,:,ind_HCHO) + change_in_mixing_ratio
   end if 

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! FOURTH INSTANCE", gchem_field(i,j,k,n)
                  print*, "PHOTOLYSIS RATE=", photolysis_rate(i,j,k)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

   if (do_HCHO_photolysis) then
      photolysis_rate = calculate_photolysis_rate(9.8E-5)
      change_in_mixing_ratio = gchem_field(:,:,:,ind_HCHO) * ( photolysis_rate * dtn )
      gchem_field(:,:,:,ind_HCHO) = gchem_field(:,:,:,ind_HCHO) - change_in_mixing_ratio
      gchem_field(:,:,:,ind_HO2) = gchem_field(:,:,:,ind_HO2) + 2 * change_in_mixing_ratio
      gchem_field(:,:,:,ind_CO) = gchem_field(:,:,:,ind_CO) + change_in_mixing_ratio
   end if 

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! FIFTH INSTANCE", gchem_field(i,j,k,n)
                  print*, "PHOTOLYSIS RATE=", photolysis_rate(i,j,k)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

   if (do_H2O2_photolysis) then
      photolysis_rate = calculate_photolysis_rate(1.18E-6)
      change_in_mixing_ratio = gchem_field(:,:,:,ind_H2O2) * ( photolysis_rate * dtn )
      gchem_field(:,:,:,ind_H2O2) = gchem_field(:,:,:,ind_H2O2) - change_in_mixing_ratio
      gchem_field(:,:,:,ind_OH) = gchem_field(:,:,:,ind_OH) + 2 * change_in_mixing_ratio
   end if 

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! SIXTH INSTANCE", gchem_field(i,j,k,n)
                  print*, "PHOTOLYSIS RATE=", photolysis_rate(i,j,k)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

   if (do_O3_photolysis) then
      photolysis_rate = calculate_photolysis_rate_O3()
      change_in_mixing_ratio = gchem_field(:,:,:,ind_O3) * ( photolysis_rate * dtn )
      gchem_field(:,:,:,ind_O3) = gchem_field(:,:,:,ind_O3) - change_in_mixing_ratio
      gchem_field(:,:,:,ind_OH) = gchem_field(:,:,:,ind_OH) + 2 * change_in_mixing_ratio
   end if 

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! SEVENTH INSTANCE", gchem_field(i,j,k,n)
                  print*, "PHOTOLYSIS RATE=", photolysis_rate(i,j,k)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine do_photolysis

function calculate_radiation_flux() result(radiation_flux)
   ! Calculates the radiation at a certain layer using simple linear interpolation
   use rad, only : insolation_TOA, swDownSurface
   use grid, only : z, dz, nz, dtn

   implicit none

   ! Parameters
   real, allocatable, dimension(:,:,:) :: radiation_flux
   integer :: i, j, k 
   real :: slope = 0.
   
   allocate(radiation_flux(nx, ny, nz))
   radiation_flux = 0.
      
   if ( allocated(swDownSurface) .and. allocated(insolation_TOA) ) then 
      do i = 1, nx
         do j = 1, ny
            radiation_flux(i,j,1) = swDownSurface(i,j)                                             ! Bottommost layer is shortwave at surface
            radiation_flux(i,j,nz) = insolation_TOA(i,j)                                           ! Topmost layer is shortwave at ToA
            
            slope = ( radiation_flux(i,j,nz) - radiation_flux(i,j,1) ) / ( z(nz) - z(1) )          ! Slope to interpolate between the two

            do k = 1, nz
               radiation_flux(i,j,k) = radiation_flux(i,j,1) + ( slope * (z(k) - z(1)) )
            enddo
         enddo
      enddo
   endif
end function calculate_radiation_flux

function calculate_photolysis_rate(maximum) result(photolysis_rate)
   ! Calculates the photolysis loss frequency for a species
   use grid, only : z, dz, nz, nzm, dtn

   implicit none

   ! Parameters
   real :: maximum
   real, allocatable, dimension(:,:,:) :: photolysis_rate
   real, allocatable, dimension(:,:,:) :: radiation_flux

   integer :: i, j, k 
   real :: maximum_radiation_flux_TOA = 1000.
   
   allocate(radiation_flux(nx, ny, nz))
   allocate(photolysis_rate(nx, ny, nz))

   radiation_flux = 0.
   photolysis_rate = 0.

   radiation_flux = calculate_radiation_flux()
   
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz
            photolysis_rate(i,j,k) = MIN( (radiation_flux(i,j,k) / maximum_radiation_flux_TOA) * maximum, maximum )
               
            if ( photolysis_rate(i,j,k) .lt. 0 ) then
               photolysis_rate(i,j,k) = 0.
            endif

         enddo
      enddo
   enddo
end function calculate_photolysis_rate

function calculate_photolysis_rate_O3() result(photolysis_rate)
   ! Calculates the photolysis loss frequency for a species
   use vars, only : qv0
   implicit none

   ! Parameters
   real, allocatable, dimension(:,:,:) :: photolysis_rate
   real, allocatable, dimension(:,:,:) :: radiation_flux
   real, allocatable, dimension(:) :: maximum_photolysis_frequency

   integer :: i, j, k 
   real :: maximum_radiation_flux_TOA = 1000.

   real :: k1 = 2.1e-10       ! corresponds to O(1D) + H2O  --> 2 OH
   real :: k2 = 3.3e-11       ! corresponds to O(1D) + N2  --> O(3P) + N2
   
   allocate(radiation_flux(nx, ny, nz))
   allocate(photolysis_rate(nx, ny, nz))
   allocate(maximum_photolysis_frequency(nz))

   radiation_flux = 0.
   photolysis_rate = 0.

   radiation_flux = calculate_radiation_flux()
   
   do i = 1, nx
      do j = 1, ny

         maximum_photolysis_frequency = 6.52e-5 * (k1 * qv0 / 1000. * 28.97/18.02 *  M_profile) / (k2 * M_profile)
         
         do k = 1, nz
            photolysis_rate(i,j,k) = MIN( (radiation_flux(i,j,k) / maximum_radiation_flux_TOA) * maximum_photolysis_frequency(k), maximum_photolysis_frequency(k) )
               
            if ( photolysis_rate(i,j,k) .lt. 0 ) then
               photolysis_rate(i,j,k) = 0.
            endif
         enddo
      enddo
   enddo
end function calculate_photolysis_rate_O3

subroutine find_tropopause(tropopause_index)
   use grid, only : nx, ny, z
   use vars, only : tabs

   implicit none
   integer :: i,j,k
   
   real, allocatable, dimension(:, :) :: tropopause_index
   real, allocatable, dimension(:, :) :: tropopause_temp

   allocate(tropopause_temp(nx, ny))

   do i = 1,nx
      do j = 1,ny
         tropopause_temp(i, j) = minval(tabs(i,j,:), dim = 1)
         tropopause_index(i, j) = minloc (tabs(i,j,:), dim = 1)
      enddo
   enddo

   ! print*, "Tropopause Temperature=", SUM(tropopause_temp) / SIZE(tropopause_temp)
   ! print*, "Altitude Temperature=", z(MAXVAL(tropopause_index))
end subroutine find_tropopause

subroutine chem_finalize()
  ! deallocate
  implicit none
  integer :: ierr

  if(isallocatedCHEM) then
     deallocate(gchem_field, STAT=ierr)
     deallocate(aqchem_field, STAT=ierr)
     deallocate(aqchem_gasprod_field, STAT=ierr)
     deallocate(archem_field, STAT=ierr)
     deallocate(gchem_profile_fixed, gchwle, gchadv, gchdiff, gchwsb,M_profile,STAT=ierr)
     deallocate(aqchwle, aqchadv, aqchdiff, aqchwsb,STAT=ierr)
     deallocate(archwle, archadv, archdiff, archwsb,STAT=ierr)
     deallocate(gchem_horiz_mean_tend)
     deallocate(aqchem_horiz_mean_tend)
     deallocate(aqchem_gasprod_horiz_mean_tend)
     deallocate(archem_horiz_mean_tend)
     deallocate(g_depos_horiz_mean_tend_IEPOX)
     deallocate(g_depos_horiz_mean_tend_ISOPOOH)
     deallocate(Haq, NO3aq, SO4aq, HSO4aq)
     deallocate(Haero, SO4aero, HSO4aero)
     deallocate(rate_const)
     deallocate(fluxbch, fluxtch)

     if(ierr.ne.0) then
        write(*,*) 'Failed to deallocated chem arrays on proc ', rank
     end if
  end if
  
end subroutine chem_finalize

subroutine chem_proc()
  implicit none
  integer :: i,j,k,n, ispecies
  real, dimension(nzm,NVAR) :: var_profile
  real, dimension(nzm,NVAR) :: gas_column_tend_profile ! in molecules/cm3/s
  real, dimension(NVAR) :: adjusted_tendency
  real, dimension(naqchem_fields) :: aq_adjusted_tendency
  real :: aq_adj_tendency
  real :: min_aerosol_radius = 1.e-13
  real :: pi = 3.1415927
  real, dimension(nzm,NFIX) :: fixed_profile

  real, dimension(nzm, naqchem_fields) :: aq_tend  ! mean tendency profiles of aq species
  real, dimension(nzm, naqchem_fields) :: aq_gasprod_tend ! mean tendency profiles of gaseous products
  real, dimension(nzm, narchem_fields) :: ar_tend ! mean tendency profiles of aerosol chemistry products

  real, dimension(nzm, naqchem_fields) :: aqchem_conc(nzm, naqchem_fields)
  real, dimension(nzm, naqchem_fields) :: aqgas_conc(nzm, naqchem_fields)
  real, dimension(nzm, naqchem_fields) :: archem_conc(nzm, narchem_fields)

  real, dimension(nzm) :: rho_tot_aerosol(nzm) ! mean aerosol density (org+inorg component)
  real, dimension(nzm) :: aero_transfer_rate(nzm) !  fractional rate of transfer of IEPOXg to aerosol surface (/s)
  real, dimension(nzm) :: aero_radius(nzm)   ! interstitial accumulation mode area weighted radius
  real, dimension(nzm) :: rho_org_aerosol(nzm)
  real, dimension(nzm) :: num_conc(nzm)
  real, dimension(nzm) :: Rdrop  ! cloud droplet radius in m
  real, dimension(nzm) :: qcloud ! temporary cloud water array
  real, dimension(nzm) :: water_vol_frac ! temporary water volume array
  real, dimension(nzm) :: pressure_atm ! temporary pressure array
  real, dimension(nzm) :: dummy_water  ! debug variable

  real :: OH_conc ! OH gas concentration in molecules/cm3

  logical :: override_gamma = .true.
  logical :: do_debug_output

  real :: IEPOX_transfer_ppv
  
  real, allocatable, dimension(:, :) :: tropopause_index
  allocate(tropopause_index(nx, ny))

  call find_tropopause(tropopause_index)

  if ( do_only_tropospheric_chemistry ) then
   do i = 1, nx
      do j = 1,ny
         do k = 1,nzm
            if ( k .gt. tropopause_index(i, j) + 5 ) then
               gchem_field(i, j, k, :) = 0.
            endif
         enddo
      enddo
   enddo
  endif

  ! compute gas chem tendencies and apply them to gas chem fields
  call compute_rate_constants(nzm, t0, M_profile, qv0, rate_const)

   do k = 1, nzm
      do n = 1, NREACT
         if ( isnan(rate_const(k,n))) then
            print*, "THERE IS A NAN! THERE IS A NAN! RATE CONSTANT!!!!!!!!!", rate_const(k,n)
            print*, "i=", i 
            print*, "j=", j 
            print*, "k=", k 
            print*, "n=", n 
            STOP
         endif
      enddo
   enddo

  !do n = 1, NREACT
  !   write(*,*) 'rate_const ', n, rate_const(1,n)
  !end do
  gchem_horiz_mean_tend(:,:) = 0.
  g_depos_horiz_mean_tend_IEPOX(:) = 0.
  g_depos_horiz_mean_tend_ISOPOOH(:) = 0.
  
  !if (do_OH_diurnal) then
  !   OH_conc = get_OH(day)     
     ! set OH level  maybe set at molecules/cm3 and derive gchem (ppv) 
  !end if
  
  do j = 1,ny
     do i = 1,nx
        do k = 1,nzm
           var_profile(k,:) = gchem_field(i,j,k,:)*M_profile(k)
           fixed_profile(k,:) = gchem_profile_fixed(k,:)*M_profile(k)
        end do
        ! special treatment for OH
        ! if (do_OH_diurnal) then
        !   do k = 1,nzm
        !      fixed_profile(k,ind_OH) = OH_conc
        !      gchem_profile_fixed(k,ind_OH) = OH_conc/M_profile(k)
        !   end do
        ! end if
        gas_column_tend_profile(:,:) = 0.

        do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! BEFORE FUN!!!!!!!!!!", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  print*, "M_profile=", M_profile(k)
                  STOP
               endif
            enddo
         enddo

        call Fun(nzm, var_profile, fixed_profile, rate_const, gas_column_tend_profile)

        do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! AFTER FUN!!!!!!!!!!", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      
        do k = 1,nzm
            if ( k .gt. tropopause_index(i, j) + 5 ) then
               gas_column_tend_profile(k,:) = 0.
           endif
           adjusted_tendency = gas_column_tend_profile(k,:)

           do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! IN THE LOOP FIRST!!!!!!!!!!", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  print*, "adjusted_tendency=", adjusted_tendency(n)
                  STOP
               endif
            enddo

           where( gchem_field(i,j,k,:) + dtn * gas_column_tend_profile(k,:)/M_profile(k) < 0.)
              adjusted_tendency = -M_profile(k) * gchem_field(i,j,k,:)/dtn
           end where   

           do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! IN THE LOOP 1!!!!!!!!!!", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo

           gchem_field(i,j,k,:) = gchem_field(i,j,k,:) + dtn*adjusted_tendency/M_profile(k)

           do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! IN THE LOOP 2!!!!!!!!!!", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  print*, "adjusted_tendency=", adjusted_tendency(n)
                  print*, "original_tendency=", gas_column_tend_profile(k,n)
                  print*, "M_profile=", M_profile(k)
                  STOP
               endif
            enddo

           gchem_horiz_mean_tend(k,:) = gchem_horiz_mean_tend(k,:) + adjusted_tendency/M_profile(k)
        end do   
     end do
  end do   
  gchem_horiz_mean_tend = gchem_horiz_mean_tend/(nx*ny)  

  do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! AFTER GAS PHASE CHEMISTRY!!!!!!!!!", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
               endif
            enddo
         enddo
      enddo
   enddo

   do i = 1,nx
        do j = 1,ny
                do k = 1,nzm
                        do n = 1, ngchem_fields
                                if ( isnan(gchem_field(i,j,k,n)) ) then
                                        STOP
                                endif
                        enddo
                enddo
         enddo

  enddo

  ! compute aqueous chemistry tendencies and apply them to aqueous fields
  aqchem_horiz_mean_tend(:,:) = 0.
  aqchem_gasprod_horiz_mean_tend(:,:) = 0.

  ! copy IEPOX into IEPOXg before calling aqueous 
  aqchem_gasprod_field(:,:,:, iIEPOX) = gchem_field(:,:,:, ind_IEPOX)
  aqchem_gasprod_field(:,:,:, iIP1NIT) = gchem_field(:,:,:, ind_ISOP1Nit)
  aqchem_gasprod_field(:,:,:, iIPDINIT) = gchem_field(:,:,:, ind_ISOPDiNit)

  pressure_atm = pres/p0
  do j = 1,ny
     do i = 1,nx

        qcloud = micro_field(i,j,:,iqcl)
        water_vol_frac = qcloud*(rho/rhol)
        do_debug_output = (j.eq.1.and.i.eq.1)
        do_debug_output = .false.
        aq_tend(:,:) = 0.
        aq_gasprod_tend(:,:) = 0.

        ! compute radius - for now assume monodisperse
        do k = 1,nzm
           num_conc(k) = max(1., micro_field(i,j,k,incl))        
              
           Rdrop(k) = ((rho(k)/rhol) * (qcloud(k)/num_conc(k)) * (0.75/pi))**(1./3.)
           if (Rdrop(k).lt.1.e-8) then
               Rdrop(k) = min_aerosol_radius ! avoid division by zero in aqueous subroutines
           endif
        end do
        if (do_debug_output) then 
           do k = 1,nzm
              write(*,*) 'k, Rdrop, num_conc, qcl', k, Rdrop(k), num_conc(k), micro_field(i,j,k,iqcl)

           end do
        endif
 
        !Rdrop = 30.e-6
        !dummy_water = 8.e-7
        ! convert aqeous inputs kg/kg to M  (mol/L)

        ! if QC is zero but AQ field is nonzero, need to do something else
        do ispecies = 1,naqchem_fields
           
           aqchem_conc(:,ispecies) = aqchem_field(i,j,:,ispecies) * rhol/molwt(ispecies)
         !            aqchem_conc(:, ispecies) = aqchem_conc(:, ispecies)/dummy_water
           do k = 1,nzm
              if (qcloud(k).gt.0.00005) then
                aqchem_conc(k, ispecies) = aqchem_conc(k, ispecies)/qcloud(k)
              else
                aqchem_conc(k, ispecies) = 0. ! could convert this straggling aq to gas?
                 
              endif
           end do 
           aqgas_conc(:,ispecies) = aqchem_gasprod_field(i,j,:,ispecies)*pres(:)/p0
        end do
        
        ! override IEPOXg with IEPOX from Isoprene model
        !  aqgas_conc(:, iIEPOX) = gchem_field(i,j,:, ind_IEPOX) ! this was missing pres(:)/p0
        ! now not needed since we fill aqchem_gasprod_field(:,:,:,iIEPOX) with gchem(ind_IEPOX) every step
        if (do_iepox_droplet_chem) then 
            call iepox_aqueous_tendencies(nzm, tabs0, pressure_atm, Rdrop, &
                 water_vol_frac, Haq, NO3aq, SO4aq, HSO4aq, &
                 aqgas_conc(:,:), aqchem_conc(:,:), &  ! input conc fields
                 aq_gasprod_tend(:,:), aq_tend(:,:), do_debug_output)   ! output tend fields

            call isop1nit_aqueous_tendencies(nzm, tabs0, pressure_atm, Rdrop, &
               water_vol_frac, &
                 aqgas_conc(:,:), aqchem_conc(:,:), &  ! input conc fields
                 aq_gasprod_tend(:,:), aq_tend(:,:), do_debug_output)   ! output tend fields
        end if    
!         call iepox_aqueous_tendencies(nzm, tabs0, pres/p0, Rdrop, &
!             dummy_water, Haq, NO3aq, SO4aq, HSO4aq, &
!             aqgas_conc(:,:), aqchem_conc(:,:), &  ! input conc fields
!             aq_gasprod_tend(:,:), aq_tend(:,:), do_debug_output)   ! output tend fields

        do ispecies = 1,naqchem_fields
           do k = 1,nzm
              if (qcloud(k).le.0.00005) then
                 aq_tend(k, ispecies) = 0.
                 aq_gasprod_tend(k, ispecies) = 0.
              end if
           end do              
        end do

!        if (do_debug_output) then
!           write (*,*) 'preconv:  k, iepoxg, iepoxg_tend, iepoxa, iepoxa_tend'
!           do k=1,3
!              write(*,*) k, aqgas_conc(k, iIEPOX), aq_gasprod_tend(k, iIEPOX), aqchem_conc(k,iIEPOX), aq_tend(k, iIEPOX)
!           end do
!        end if
        
        ! convert aqueous output tendencies back to model dims        
        do ispecies = 1,naqchem_fields
           aq_tend(:, ispecies) = aq_tend(:, ispecies) * qcloud(:) * molwt(ispecies)/rhol   ! to kg/kg
           !aq_tend(:, ispecies) = aq_tend(:, ispecies) * dummy_water * molwt(ispecies)/rhol   ! to kg/kg
           aq_gasprod_tend(:, ispecies) = aq_gasprod_tend(:,ispecies) * (p0/pres(:)) ! to ppv       
        end do
   
!        if (do_debug_output) then
!           write (*,*) 'postconv:  k,  iepoxg, iepoxg_tend, iepoxa, iepoxa_tend'
!           do k=1,3
!              write(*,*) k, aqgas_conc(k, iIEPOX), aq_gasprod_tend(k, iIEPOX), aqchem_conc(k,iIEPOX), aq_tend(k, iIEPOX)
!           end do
!        end if
        
        do k = 1,nzm
           aq_adjusted_tendency = aq_tend(k,:)
           where (aqchem_field(i,j,k,:) + dtn*aq_tend(k,:) < 0.)
              aq_adjusted_tendency = -aqchem_field(i,j,k,:)/dtn
           end where   
           aqchem_field(i,j,k,:) = aqchem_field(i,j,k,:) + dtn*aq_adjusted_tendency
           aqchem_horiz_mean_tend(k,:) = aqchem_horiz_mean_tend(k,:) + aq_adjusted_tendency

           aq_adjusted_tendency = aq_gasprod_tend(k,:)
           where (aqchem_gasprod_field(i,j,k,:) + dtn*aq_gasprod_tend(k,:) < 0.)
              aq_adjusted_tendency = -aqchem_gasprod_field(i,j,k,:)/dtn
           end where   
           aqchem_gasprod_field(i,j,k,:) = aqchem_gasprod_field(i,j,k,:) + dtn*aq_adjusted_tendency
           aqchem_gasprod_horiz_mean_tend(k,:) = aqchem_gasprod_horiz_mean_tend(k,:) + aq_adjusted_tendency
!           aq_adj_tendency = aq_gasprod_tend(k,iIEPOX)
!           if ((gchem_field(i,j,k,ind_IEPOX) + dtn*aq_gasprod_tend(k,iIEPOX)).lt.0.) then
!              aq_adj_tendency = -gchem_field(i,j,k,ind_IEPOX)/dtn
!           end if
           if (do_iepox_droplet_chem) then
              !              gchem_field(i,j,k, ind_IEPOX) = gchem_field(i,j,k, ind_IEPOX) + dtn*aq_adj_tendency
              gchem_field(i,j,k, ind_IEPOX) = gchem_field(i,j,k, ind_IEPOX) + dtn*aq_adjusted_tendency(iIEPOX)
              gchem_field(i,j,k, ind_ISOP1Nit) = gchem_field(i,j,k, ind_ISOP1Nit) + dtn*aq_adjusted_tendency(iIP1NIT)
              gchem_field(i,j,k, ind_ISOPDiNit) = gchem_field(i,j,k, ind_ISOPDiNit) + dtn*aq_adjusted_tendency(iIPDINIT)
              
           end if
        end do
     end do
  end do

  do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! BEFORE DEPOSITION INSTANCE", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
               endif
            enddo
         enddo
      enddo
   enddo

   do i = 1,nx
        do j = 1,ny
                do k = 1,nzm
                        do n = 1, ngchem_fields
                                if ( isnan(gchem_field(i,j,k,n)) ) then
                                        STOP
                                endif
                        enddo
                enddo
         enddo

  enddo

  ! Apply deposition to selected variables at lowest grid level only
  do j = 1,ny
     do i = 1,nx
        gchem_field(i,j,1,ind_ISOPOOH) =  (1.-deposition_rate*dtn) * gchem_field(i,j,1,ind_ISOPOOH)
        gchem_field(i,j,1,ind_IEPOX) = (1.-deposition_rate*dtn) * gchem_field(i,j,1,ind_IEPOX)
        g_depos_horiz_mean_tend_ISOPOOH(1) = g_depos_horiz_mean_tend_ISOPOOH(1) - deposition_rate * gchem_field(i,j,1,ind_ISOPOOH)
        g_depos_horiz_mean_tend_IEPOX(1) = g_depos_horiz_mean_tend_IEPOX(1) - deposition_rate * gchem_field(i,j,1,ind_IEPOX)
     end do
  end do  

  do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! AFTER DEPOSITION INSTANCE", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo
  ! Do photolysis
  call do_photolysis

  do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "THERE IS A NAN! THERE IS A NAN! AFTER PHOTOLYSIS INSTANCE", gchem_field(i,j,k,n)
                  print*, "i=", i 
                  print*, "j=", j 
                  print*, "k=", k 
                  print*, "n=", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo

  archem_horiz_mean_tend(:,:) = 0.

  if (do_iepox_aero_chem) then
     do j = 1,ny
        do i = 1,nx
           do_debug_output = (j.eq.1.and.i.eq.1)
           do_debug_output = .false.
           ar_tend(:,:) = 0.
           aero_radius = min_aerosol_radius
           do k = 1,nzm
              if (micro_field(i,j,k,inad).gt.1) then
                 aero_radius(k) = 0.5 * ((1/rho_aerosol) * (micro_field(i,j,k,iqad)/micro_field(i,j,k,inad)) * (.75/pi)) **(1./3.)* &
                      EXP(3*LOG(sigma_accum)**2)
              end if
              if (do_debug_output) then
                 write(*,*) 'k, qad, nad, sigma, radius=' , k, micro_field(i,j,k,iqad), micro_field(i,j,k,inad), sigma_accum, aero_radius(k)
              end if
           end do           

           where(aero_radius.lt.min_aerosol_radius)
              aero_radius = min_aerosol_radius
           end where
           
           call iepox_aero_transfer_rate(nzm, tabs0, pressure_atm, aero_radius, Haero, &
                actHaero, SO4aero, HSO4aero, OrgMF, override_gamma, aero_transfer_rate, rho_org_aerosol, do_debug_output)
      
           if (do_debug_output) then
              do k = 1,nzm
                 write(*,*) 'k, radius, aero_transfer_rate= ', k, aero_radius(k), aero_transfer_rate(k)
              end do
           end if
           

           do k = 1,nzm
              ! Multiply aero_transfer_rate by aerosol number concentration
              aero_transfer_rate(k) = aero_transfer_rate(k) * micro_field(i,j,k,inad)
              ! limit IEPOXg loss
              if (aero_transfer_rate(k) * dtn.ge.1.) then
                 aero_transfer_rate(k) = 1./dtn
              end if
           end do   
            

           do k = 1,nzm
              if (aero_radius(k).ge.1.e-12) then  ! avoid division by zero for 0 size aerosol
                 ! apply IEPOXg loss
                 IEPOX_transfer_ppv = dtn * aero_transfer_rate(k)*gchem_field(i,j,k,ind_IEPOX)
                 gchem_field(i,j,k, ind_IEPOX) = gchem_field(i,j,k, ind_IEPOX) - IEPOX_transfer_ppv       ! *(p0/pres(k)) ! convert to ppv
                 ! distribute aerosol mass gain (converting to kg aerosol/kg air)
                 archem_field(i,j,k, iTETROLr) = archem_field(i,j,k, iTETROLr) + &
                   FracTETROL*IEPOX_transfer_ppv * molwt_ar(iTETROLr)/28.96  !  WHAT IS AIR MW CONSTANT CALLED - REPLACE HERE and next line
                 archem_field(i,j,k, iIEPOX_SO4r) = archem_field(i,j,k,iIEPOX_SO4r) + &
                   FracIEPOX_SO4 * IEPOX_transfer_ppv * molwt(iIEPOX_SO4r)/28.96 

                 archem_horiz_mean_tend(k,iTETROLr) =  archem_horiz_mean_tend(k,iTETROLr) + FracTETROL * IEPOX_transfer_ppv * molwt_ar(iTETROLr)/28.96
                 archem_horiz_mean_tend(k,iIEPOX_SO4r) =  archem_horiz_mean_tend(k,iIEPOX_SO4r) + FracIEPOX_SO4 * IEPOX_transfer_ppv * molwt_ar(iIEPOX_SO4r)/28.96
              
              end if                 
          
           end do   
        
        end do
     end do
  end if   

  ! normalize tendencies
  aqchem_horiz_mean_tend = aqchem_horiz_mean_tend/(nx*ny)
  aqchem_gasprod_horiz_mean_tend = aqchem_gasprod_horiz_mean_tend/(nx*ny)
  archem_horiz_mean_tend = archem_horiz_mean_tend/(nx*ny)

  g_depos_horiz_mean_tend_ISOPOOH = g_depos_horiz_mean_tend_ISOPOOH/(nx*ny)
  g_depos_horiz_mean_tend_IEPOX = g_depos_horiz_mean_tend_IEPOX/(nx*ny)
end subroutine chem_proc

subroutine chem_print()
  implicit none
end subroutine chem_print  

subroutine chem_statistics()
  use hbuffer, only: hbuf_put
  implicit none
  ! average fields in space for .stat file

  real, dimension(nzm) :: tr0, tendency
  real, dimension(nzm) :: zeros

  real factor_xy
  integer i,j,k,m, n, ii, jj, nn, ncond

  character*16 name
  character*16 tend_name
  factor_xy = 1./float(nx*ny)

  zeros(:) = 0.
  do n = 1,ngchem_fields
    ! compute horizontal mean of all gas chem fields  
    do k = 1,nzm
       tr0(k) = SUM(gchem_field(1:nx,1:ny,k,n))
    end do

    !if(n.eq.1) write(*,*) 'IP1O2 mean, (1,1,1) = ', tr0(1), gchem_field(1,1,1,1)*1.e9
    
    call hbuf_put(trim(SPC_NAMES(n)), tr0, gas_output_scale*factor_xy) ! factor is 1/(nx * ny)
    call hbuf_put(trim(SPC_NAMES(n))//'+', gchem_horiz_mean_tend(:, n), gas_output_scale)
  end do  

  if ( ngchem_fixed .ge. 1)  then
        do n = 1,ngchem_fixed
                name = trim(SPC_NAMES(n+ngchem_fields))
                call hbuf_put(name, gchem_profile_fixed(:, n), gas_output_scale)
        end do
  endif

  do n = 1,naqchem_fields
    ! compute horizontal mean of all aq chem fields  
    do k = 1,nzm
       tr0(k) = SUM(aqchem_field(1:nx,1:ny,k,n))
    end do

    call hbuf_put(trim(aq_species_names(n)), tr0, factor_xy) ! factor is 1/(nx * ny)
    call hbuf_put(trim(aq_species_names(n))//'+', aqchem_horiz_mean_tend(:, n), 1.)
 end do
 
 do n = 1,naqchem_fields
    do k = 1,nzm
       tr0(k) = SUM(aqchem_gasprod_field(1:nx,1:ny,k,n))
    end do
    
    call hbuf_put(trim(aq_gasprod_species_names(n)), tr0, factor_xy*gas_output_scale) ! factor is 1/(nx * ny)
    call hbuf_put(trim(aq_gasprod_species_names(n))//'+', aqchem_gasprod_horiz_mean_tend(:, n), gas_output_scale)   
 end do  

 do n = 1,narchem_fields
    do k = 1,nzm
        tr0(k) = SUM(archem_field(1:nx,1:ny,k,n))
     end do
     
     call hbuf_put(trim(ar_species_names(n)), tr0, factor_xy) ! factor is 1/(nx * ny)
     call hbuf_put(trim(ar_species_names(n))//'+', archem_horiz_mean_tend(:, n), 1.)
  end do

  call hbuf_put('IPOOHd+', g_depos_horiz_mean_tend_IEPOX, gas_output_scale)
  call hbuf_put('IEPOXd+', g_depos_horiz_mean_tend_IEPOX, gas_output_scale)
end subroutine chem_statistics  

subroutine chem_write_fields3D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1
  character *80 long_name
  character *16 name
  character *10 units
  integer :: i, j, k, f 
  real(4), dimension(nx,ny,nzm) :: tmp

  do f = 1,ngchem_fields
     if(flag_gchemvar_out3D(f)) then  
        nfields1=nfields1+1
        do k=1,nzm
           do j=1,ny
              do i=1,nx
                 tmp(i,j,k)=gchem_field(i,j,k,f)*gas_output_scale
              end do
           end do
        end do
        name=TRIM(SPC_NAMES(f))
        long_name=TRIM(SPC_NAMES(f))
        units='ppbv'
        call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains)
     end if
  end do

  do f = 1,naqchem_fields
     if(flag_aqchemvar_out3D(f)) then  
        nfields1=nfields1+1
        tmp = aqchem_field(1:nx,1:ny,:,f)
        name=TRIM(aq_species_names(f))
        long_name=TRIM(aq_species_names(f))
        units='kg/kg'
        call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains)
     end if
  end do

  do f = 1,naqchem_fields
     if(flag_aqchemgasvar_out3D(f)) then  
        nfields1=nfields1+1
        tmp = aqchem_gasprod_field(1:nx, 1:ny,:,f) ! account for ghost cells - specify 1:nx, 1,ny
        name=TRIM(aq_gasprod_species_names(f))
        long_name=TRIM(aq_gasprod_species_names(f))
        units='kg/kg'
        call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains)
     end if
  end do

  do f = 1,narchem_fields
     if(flag_archemvar_out3D(f)) then  
        nfields1=nfields1+1
        tmp = archem_field(1:nx,1:ny,:,f)
        name=TRIM(ar_species_names(f))
        long_name=TRIM(ar_species_names(f))
        units='kg/kg'
        call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains)
     end if
  end do
  
end subroutine chem_write_fields3D

end module chemistry
