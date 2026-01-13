module chemistry

  use grid, only: nx, ny, nzm, nz, &
       dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
     z, zi, pres, adz, dz, dx, nx_gl, &
     time, dt, nstep, ncycle, nstat, nstatfrq, nrestart, day, &
     rank, dompi, masterproc, nsubdomains, save3Dbin, &
     case, caseid

  use rad, only : swDown, swUp
  use vars, only: rho, dtn, qcl, tabs0, pres, qv0, tabs, qpl
  use cloudchem_Parameters, only: NVAR, NSPEC, NFIX, ind_IEPOX, ind_ISOP1Nit, ind_ISOPOOH, ind_ISOPDiNit, ind_OH, NREACT ! , ind_ISOP1Nit, ind_ISOPDiNit,    ! NSPEC=NVAR+NFIX

  !! From KPP-generated files
  use cloudchem_Monitor, only: SPC_NAMES
  use cloudchem_Function, only: Fun
  use cloudchem_Global, only: RCONST, TEMP, C, ATOL, RTOL, SUN, C_M, C_H2O
  use cloudchem_Rates, only: Update_RCONST, Update_PHOTO
  use cloudchem_Integrator, only: Integrate

  use chemistry_params, only: p0, rhol, do_only_tropospheric_chemistry, do_transport_loss, &
                           do_OH_diurnal, OH_night, OH_day_peak, do_rainout, do_washout, do_convective_scavenging, wet_deposition_time_step, wet_deposition_species, k0_constants, cr_constants, dry_deposition_species, dry_deposition_velocities, &
                           do_NO2_photolysis, do_O3_photolysis, do_H2O2_photolysis, do_HCHO_photolysis, do_NO3_photolysis, do_ISOPOOH_photolysis, &
                           do_megan_isoprene, do_surface_Isoprene_diurnal, do_bdsnp_no, &
                           do_IC_lightning, do_CTG_lightning, IC_decaria, CTG_decaria_reflectivity, CTG_price_and_rind, &
                           do_iepox_droplet_chem, do_iepox_aero_chem, hi_org, pHdrop, pHaero, & 
                           gas_init_name, gas_init_value, gas_out3D_name, flag_gchemvar_out3D

  use chem_aqueous, only: naqchem_fields, molwt, iIEPOX, iTETROL, iIEPOX_SO4,&
       iIP1NIT, iIPDINIT, aq_species_names, aq_gasprod_species_names, &
       flag_aqchemvar_out3D, flag_aqchemgasvar_out3D, &
       iepox_aqueous_tendencies, isop1nit_aqueous_tendencies
      
  use chem_aerosol, only:   iTETROLr, iIEPOX_SO4r,iepox_aero_transfer_rate, narchem_fields, molwt_ar, flag_archemvar_out3D, ar_species_names
  use chemistry_aux, only: get_OH
  
  use micro_params, only: MW_air, avgd, rho_aerosol, sigma_accum
  use microphysics, only: micro_field, iqcl, iqad, inad, incl

  use emissions, only : surface_emission_flux_driver, lightning_decaria_ic, lightning_decaria_ctg
  use photolysis, only : photolysis_driver
  use deposition, only : dry_deposition_driver
  use transport_loss, only : transport_removal_driver
  use wet_deposition, only : wet_deposition_driver

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
  real, allocatable, dimension(:) :: soil_NO_emission_flux 
  real, allocatable, dimension(:) :: isop_emission_flux


  ! For wet deposition
  real, allocatable, dimension(:,:,:) :: previous_qcl  !  air density in molec/cm3
  real, allocatable, dimension(:,:,:) :: change_in_qcl  !  air density in molec/cm3

  real, allocatable, dimension(:,:,:) :: previous_qpl  ! array of gas reaction rate constants
  real, allocatable, dimension(:,:,:) :: change_in_qpl  ! array of gas reaction rate constants
  
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
   
NAMELIST /CHEMISTRY/   do_only_tropospheric_chemistry, do_transport_loss, &
                           do_OH_diurnal, OH_night, OH_day_peak, do_rainout, do_washout, do_convective_scavenging, wet_deposition_species, wet_deposition_time_step, k0_constants, cr_constants, &
                           dry_deposition_species, dry_deposition_velocities, &
                           do_NO2_photolysis, do_O3_photolysis, do_H2O2_photolysis, do_HCHO_photolysis, do_NO3_photolysis, do_ISOPOOH_photolysis, &
                           do_megan_isoprene, do_surface_Isoprene_diurnal, do_bdsnp_no, &
                           do_IC_lightning, do_CTG_lightning, IC_decaria, CTG_decaria_reflectivity, CTG_price_and_rind, &
                           do_iepox_droplet_chem, do_iepox_aero_chem, hi_org, pHdrop, pHaero, & 
                           gas_init_name, gas_init_value, gas_out3D_name        
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
       allocate(soil_NO_emission_flux(nzm))
       allocate(isop_emission_flux(nzm))
       
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

       ! Wet deposition
      allocate(previous_qcl(nx, ny, nzm), previous_qpl(nx, ny, nzm))
      allocate(change_in_qcl(nx, ny, nzm), change_in_qpl(nx, ny, nzm))
    end if   

    previous_qcl(:,:,:) = -9999
    previous_qpl(:,:,:) = -9999

    change_in_qcl(:,:,:) = 0.
    change_in_qpl(:,:,:) = 0.

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
  M_profile = 0.001 * RHO * avgd / MW_air
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

  soil_NO_emission_flux(:) = 0.
  isop_emission_flux(:) = 0.

  ! top and bottom fluxes of fields
  fluxbch = 0.
  fluxtch = 0.

end subroutine chem_init 

subroutine chem_proc()
  implicit none

  integer :: i,j,k,n, ispecies
  real, dimension(nzm, NVAR) :: var_profile

  real, dimension(nzm, NVAR) :: gas_column_tend_profile ! in molecules/cm3/s
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

  ! For KPP !
   INTEGER                :: IERR                      ! KPP success or failure flag
   INTEGER                :: ICNTRL (20)
   INTEGER                :: ISTATUS(20)
   REAL(8)                   :: RCNTRL (20)
   REAL(8)                   :: RSTATE (20)
  
  real, allocatable, dimension(:, :) :: tropopause_index
  allocate(tropopause_index(nx, ny))
  
  IERR      = 0                        ! KPP success or failure flag
  ISTATUS   = 0.0                   ! Rosenbrock output
  ICNTRL    = 0                        ! Rosenbrock input (integer)
  RCNTRL    = 0.0                   ! Rosenbrock input (real)
  RSTATE    = 0.0                   ! Rosenbrock output
  
  ICNTRL(1) = 1
  ICNTRL(2) = 1
  ICNTRL(4) = 40
  ICNTRL(15) = -1
  ICNTRL(16) = 1

  ATOL = 1.0e-2
  RTOL = 0.5e-2

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

  gchem_horiz_mean_tend(:,:) = 0.
  g_depos_horiz_mean_tend_IEPOX(:) = 0.
  g_depos_horiz_mean_tend_ISOPOOH(:) = 0.

  do j = 1, ny
     do i = 1, nx
         gas_column_tend_profile(:,:) = 0.

        do k = 1, nzm 
            C_M = M_profile(k)                 ! Air density [molec/cm3]
            C_H2O = C_M * qv0(k) / 1000. * 28.97 / 18.02
            ! gchem_field(i, j, k, ind_H2O) = qv0(k) / 1000. * 28.97 / 18.02         ! Update water vapor mixing ratio (in mol mol-1 w.r.t. dry air)
            
            var_profile(k,:) = gchem_field(i,j,k,:) * M_profile(k)               ! Convert from mole fraction to [molec cm-3]
            fixed_profile(k,:) = gchem_profile_fixed(k,:) * M_profile(k)         ! Convert from mole fraction to [molec cm-3]

            TEMP = tabs(i, j, k)

            ! Set concentrations in KPP to the current box's concentration [molec cm-3]
            C(1:NVAR) = var_profile(k, :)
            ! C(NVAR+1:NSPEC) = fixed_profile(k, :)

            call Update_RCONST()
            if ( allocated(swDown) ) then
               SUN = MIN( ( swDown(i,k) + swUp(i, k) ) / 1500., 1.)
               call Update_PHOTO()
            endif

            rate_const(k, :) = RCONST                                            ! rate_const is the local version of RCONST
            call Fun(var_profile(k, :), fixed_profile(k, :), rate_const(k, :), gas_column_tend_profile(k, :))     ! Calculate derivatives
            CALL Integrate(0.0, dtn, ICNTRL, RCNTRL, ISTATUS, RSTATE, IERR)       ! Use RCONST and Func to integrate the concentrations a timestep

            IF ( IERR < 0 ) THEN
               print*, "Oops, IERR < 0!"
            ENDIF

            adjusted_tendency = ( C(1:NVAR)  / M_profile(k) ) - gchem_field(i,j,k,:)     ! Calculate the change in the species from one timestep to another
            gchem_field(i,j,k,:) = C(1:NVAR) / M_profile(k)
            gchem_horiz_mean_tend(k,:) = gchem_horiz_mean_tend(k,:) + adjusted_tendency ! / M_profile(k)
        enddo
        
        call check_nan
     end do
  end do   
   
  gchem_horiz_mean_tend = gchem_horiz_mean_tend / (nx * ny)  

  ! Do photolysis
  ! call photolysis_driver(gchem_field, M_profile)

  ! Do dry deposition
  call dry_deposition_driver(gchem_field, g_depos_horiz_mean_tend_ISOPOOH, g_depos_horiz_mean_tend_IEPOX)

  ! Transport removal to simulate large-scale advection
  if ( do_transport_loss ) then
   call transport_removal_driver(gchem_field)
  endif

  ! Do wet deposition
  if ( do_convective_scavenging .or. do_rainout  .or. do_washout ) then
   if ( mod(nstep, wet_deposition_time_step) .eq. 0 ) then
      if ( ( all(previous_qcl .eq. -9999) ) .or. ( all(previous_qpl .eq. -9999) ) ) then
         previous_qcl = qcl
         previous_qpl = qpl
      else
         change_in_qcl = qcl - previous_qcl
         change_in_qpl = qpl - previous_qpl

         call wet_deposition_driver(gchem_field, M_profile, change_in_qcl, change_in_qpl)

         previous_qcl = qcl
         previous_qpl = qpl
      endif
   endif
 endif

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

        qcloud = micro_field(i,j,:,iqcl) ! qcl(i,j,:)
        water_vol_frac = qcloud*(rho/rhol)
        do_debug_output = (j.eq.1.and.i.eq.1)
        do_debug_output = .false.
        aq_tend(:,:) = 0.
        aq_gasprod_tend(:,:) = 0.

        ! compute radius - for now assume monodisperse
        do k = 1,nzm
           num_conc(k) = max(1., micro_field(i,j,k,incl))        ! Talk to Peter
              
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
            aq_adj_tendency = aq_gasprod_tend(k,iIEPOX)
            if ((gchem_field(i,j,k,ind_IEPOX) + dtn*aq_gasprod_tend(k,iIEPOX)).lt.0.) then
               aq_adj_tendency = -gchem_field(i,j,k,ind_IEPOX)/dtn
            end if
           if (do_iepox_droplet_chem) then
              gchem_field(i,j,k, ind_IEPOX) = gchem_field(i,j,k, ind_IEPOX) + dtn*aq_adj_tendency
              ! gchem_field(i,j,k, ind_IEPOXD) = gchem_field(i,j,k, ind_IEPOXD) + dtn*aq_adjusted_tendency(iIEPOX)
              gchem_field(i,j,k, ind_ISOP1Nit) = gchem_field(i,j,k, ind_ISOP1Nit) + dtn*aq_adjusted_tendency(iIP1NIT)
              gchem_field(i,j,k, ind_ISOPDiNit) = gchem_field(i,j,k, ind_ISOPDiNit) + dtn*aq_adjusted_tendency(iIPDINIT)
              
           end if
        end do
     end do
  end do

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
end subroutine find_tropopause

subroutine check_nan
   implicit none

   real :: i, j, k, n

   do i = 1,nx
      do j = 1,ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( isnan(gchem_field(i,j,k,n)) ) then
                  print*, "******* There is a NaN in gchem_field. Check if your system is unstable! ****** ", gchem_field(i,j,k,n)
                  print*, "i = ", i 
                  print*, "j = ", j 
                  print*, "k = ", k 
                  print*, "n = ", n 
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo
end subroutine check_nan

subroutine chem_hbuf_init(namelist, deflist, unitlist, status, average_type, count, chemcount)
   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*), average_type(*), count, chemcount, n

   character*16 name
   character*16 tend_name
   character*80 longname
   character*10 units

   chemcount = 0

   do n = 1, ngchem_fields
      name = trim( SPC_NAMES(n) )
      longname = trim( SPC_NAMES(n) )
      units = 'ppbv'
      call add_to_namelist(count, chemcount, name, longname, units, 0)

      tend_name = trim(SPC_NAMES(n)) // '+'
      longname = trim(SPC_NAMES(n)) // ' tendency due to gas reaction'
      units = 'ppbv/s'
      call add_to_namelist(count, chemcount, tend_name, longname, units, 0)
   end do

   do n = 1, ngchem_fixed
      name = trim(SPC_NAMES(n + ngchem_fields))
      longname = name // ' Fixed Species'
      units = 'ppbv'
      call add_to_namelist(count, chemcount, name, longname, units, 0)
   end do

   if ( masterproc ) then
      write(*,*) 'Added ', chemcount, ' arrays to statistics for gaseous chemical species!'
   end if

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

   name = 'SOILNO_E'
   longname = 'Soil NOx Emission Flux'
   units = 'ppbv m/s'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   name = 'ISOP_E'
   longname = 'Isoprene Emission Flux'
   units = 'ppbv m/s'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   if ( masterproc ) then
      write(*,*) 'Added ', chemcount, ' arrays to statistics for gas and aqueous chemical species!'
   end if

end subroutine chem_hbuf_init

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

! chem_flux is a wrapper function to call emissions module
subroutine chem_flux()
   soil_NO_emission_flux(:) = 0.
   isop_emission_flux(:) = 0.

   call surface_emission_flux_driver(fluxbch, M_profile, isop_emission_flux, soil_NO_emission_flux)

   if ( do_CTG_lightning ) then
      call lightning_decaria_ctg(gchem_field)
   endif

   if ( do_IC_lightning ) then 
      call lightning_decaria_ic(gchem_field)
   endif

end subroutine chem_flux


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

  print*, "SOIL NO BEFORE WRITING: ", SUM(soil_NO_emission_flux)
  call hbuf_put('SOILNO_E', soil_NO_emission_flux, gas_output_scale * factor_xy)
  call hbuf_put('ISOP_E', isop_emission_flux, gas_output_scale * factor_xy)

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
