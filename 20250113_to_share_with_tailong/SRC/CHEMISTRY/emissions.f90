module emissions

    use grid, only : nx, ny, nz, nzm, dx, dy, dz, day, time, dt, pres, nstep               ! pres is in mbar!
    use cloudchem_Parameters, only : ind_NO, ind_ISOP
    use chemistry_params, only : do_megan_isoprene, do_surface_Isoprene_diurnal, do_bdsnp_no, CTG_decaria_reflectivity, CTG_price_and_rind, IC_decaria

    ! Common constants
    real :: pi = 3.1415927
    real :: universal_gas_constant = 8.314                                  ! gas constant, in J mol-1 K-1

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! SURFACE EMISSIONS !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! do_surface_emission_flux() provides a source of surface emissions to the bottommost model box
    ! you must send by reference the fluxbch (flux) (x,y) array to populate
    ! fluxbch is a variable in the chemistry module
    subroutine surface_emission_flux_driver(fluxbch, M_profile, isop_emission_flux, soil_NO_emission_flux)
        use vars, only : tabs
        use rad, only : swDownSurface, lwDownSurface, insolation_TOA                    ! eventually converted to PAR

        implicit none

        ! Average fluxes
        real :: iso_avg_flux = 1.575e-10 * 0.6  ! in kg ISOP m-2 s-1
        real :: conversion_between_isoprene_and_nox = 0.1765 * 0.8 * 0.15 !  0.0154454 !0.00077227 ! 0.0072    ! To scale average isoprene flux [ratio of kg m-2 s-1 to kg m-2 s-1]

        real :: isop_molar_mass = 68.12   ! in g/mol
        real :: no_molar_mass = 30.01     ! in g/mol
        real :: avogadro_number = 6.022e23       ! molecules/mol

        ! Conversion factors
        real :: J_to_mol = 4.6                                           ! Approximate unit conversion between W m-2 and umol photons m-2 s-1
        real :: frac_PAR = 0.5                                           ! Fraction of downward shortwave that is PAR

        ! Allocatable arrays for emissions
        real, allocatable, dimension(:,:,:) :: fluxbch                      ! parameter containing the surface fluxes
        real, allocatable, dimension(:) :: isop_emission_flux, soil_NO_emission_flux  ! parameter containing the surface fluxes

        real, allocatable, dimension(:,:) :: soil_NOx_activity_factor ! soil NOx array (x,y)
        real, allocatable, dimension(:,:) :: ppfd ! to calculate PAR-dependence of isoprene emissions
        real, allocatable, dimension(:,:) :: radiation_activity_factor ! to calculate PAR-dependence of isoprene emissions
        real, allocatable, dimension(:,:) :: temperature_activity_factor ! to calculate PAR-dependence of isoprene emissions
        real, allocatable, dimension(:)   :: M_profile                      ! parameter containing the surface fluxes

        ! Other variables
        real :: LDF_i = 1                                                ! Set for isoprene

        !! For non-MEGAN ISOP diurnal cycle
        real :: t_solar_peak = 0.167                                                ! days
        real :: frac_of_peak

        ! Default values for emissions
        fluxbch = 0.
    
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
                ppfd(:,:) = swDownSurface(:,:) * J_to_mol * frac_PAR
                radiation_activity_factor(:,:) = calculate_megan_BVOC_radiation(ppfd(:,:), LDF_i)
            endif

            !! MEGAN Isoprene Temperature Module !!
            ! Calculate temperature activity factor for isoprene
            temperature_activity_factor = calculate_megan_BVOC_temperature(tabs(:,:,1), LDF_i)

            fluxbch(:,:,ind_ISOP) =  ( iso_avg_flux * 1000 / isop_molar_mass / M_profile(1) * avogadro_number / 100**3 ) * temperature_activity_factor * radiation_activity_factor     ! Calculate isoprene fluxes using MEGAN
            isop_emission_flux(1) = SUM(fluxbch(:,:,ind_ISOP))
            deallocate(temperature_activity_factor, ppfd, radiation_activity_factor)

        elseif ( do_surface_Isoprene_diurnal ) then
            frac_of_peak = MAX(0., cos((day - t_solar_peak)*2*pi))
            fluxbch(:,:,ind_ISOP) =  frac_of_peak * Iso_avg_flux * (pi/2.)
        
        else
            fluxbch(:,:,ind_ISOP) = iso_avg_flux
        endif
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Soil NO Flux Module !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        if ( do_bdsnp_no ) then
            allocate(soil_NOx_activity_factor(nx, ny))
            soil_NOx_activity_factor = 0.                       ! Initialize to zero

            call calculate_bdsnp_NO(soil_NOx_activity_factor)
            fluxbch(:, :, ind_NO) = ( iso_avg_flux * conversion_between_isoprene_and_nox * 1000 / no_molar_mass / M_profile(1) * avogadro_number / 100**3 ) * soil_NOx_activity_factor(:, :)
            soil_NO_emission_flux(1) = SUM(fluxbch(:, :, ind_NO))
            
            print*, "SOIL NO, ACTUAL: ", SUM(fluxbch(:, :, ind_NO)) / SIZE(fluxbch(:, :, ind_NO))
            print*, "ACTIVITY FACTOR:  ", SUM(soil_NOx_activity_factor) / SIZE(soil_NOx_activity_factor)
            deallocate(soil_NOx_activity_factor)
        endif

    end subroutine surface_emission_flux_driver

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


    subroutine calculate_bdsnp_NO(soil_NOx_activity_factor)
        ! Calculates soil NOx
        ! See Hudman et al. (2012) and Wang et al. (2021) for parametrization
        use vars, only : precsfc, tabs, interactive_soil_wetness ! precsfc (x,y), tabs (x,y,z) (JY), and interactive_soil_wetness (only used in chem_flux) (JY)

        implicit none

        ! BDSNP Parameters
        real :: a_bdsnp = 1.65
        real :: b_bdsnp = 3.3
        real :: temperature_in_celsius

        real, allocatable, dimension(:,:) :: soil_NOx_activity_factor                 ! temperature activity factor, light-dependent
        real, dimension(nx, ny) :: surface_temperature                 ! temperature at bottommost layer

        ! For interactive soil moisture
        real :: decay_in_soil_moisture = 0.00007                    ! From MERRA2 Regressions
        real :: increase_in_soil_moisture = 5.6                   ! From MERRA2 Regressions

        ! Conversion factors in soil moisture calculations
        integer :: scale_between_kg_m2_s1_and_mm_day = 86400    ! MERRA2 precipitation is in kg m-2 s-1; precipitation in SAM is mm/day (?)
        real :: conversion_between_hour_and_second = 3600       ! MERRA2 is on an hourly time grid versus seconds for SAM

        integer :: i, j                                                  ! Counter variables

        ! Loop through the horizontal axes
        do i = 1, nx
            do j = 1, ny

                ! If there is precipitation, increase soil moisture
                ! Precsfc is outputted as mm/day, but in the model it is kg m-2 s-1!
                if ( ( precsfc(i, j) > 0 ) .and. ( precsfc(i, j) < 0.01 )  ) then
                    interactive_soil_wetness(i, j) = interactive_soil_wetness(i, j) + (increase_in_soil_moisture * dt * precsfc(i,j) / conversion_between_hour_and_second)
                
                ! If there is no precipitation, decay the soil moisture
                else
                    interactive_soil_wetness(i, j) = interactive_soil_wetness(i, j) - (decay_in_soil_moisture * dt / conversion_between_hour_and_second)
                endif

                if ( interactive_soil_wetness(i, j) .gt. 1 ) then
                    interactive_soil_wetness(i, j) = 1
                elseif ( interactive_soil_wetness(i, j) .lt. 0 ) then
                    interactive_soil_wetness(i, j) = 0
                endif

                !! Now, calculate soil NOx !!
                temperature_in_celsius = tabs(i, j, 1) - 273.15

                if ( temperature_in_celsius < 20 ) then
                    soil_NOx_activity_factor(i, j) = exp( 0.103 * temperature_in_celsius ) * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
                
                elseif ( ( temperature_in_celsius >= 20 ) .and. ( temperature_in_celsius <= 40 ) ) then
                    soil_NOx_activity_factor(i, j) = ( -0.009 * temperature_in_celsius**3 + 0.837 * temperature_in_celsius**2 - 22.52 * temperature_in_celsius + 196.149 ) * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
                
                else
                    soil_NOx_activity_factor(i, j) = 58.269 * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
                endif

            end do
        end do
        print*, "MEAN PREC: ", SUM(precsfc) / SIZE(precsfc)
        print*, "MAX PREC: ", MAXVAL(precsfc)
        print*, "SOIL WETNESS: ", SUM(interactive_soil_wetness) / SIZE(interactive_soil_wetness)

    end subroutine calculate_bdsnp_NO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! LIGHTNING EMISSIONS !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lightning_decaria_ctg(gchem_field)
        ! Parametrizes (as close as possible!) to DeCaria et al. (2000), 
        ! A cloud-scale model study of lightning-generated NOx...
        
        use grid, only : z
        use vars, only : tabs                       ! tabs (x,y,z) contains absolute temperature in K
        use params, only : land, ocean
        use microphysics, only: dBZ_cloudradar
        
        implicit none

        real, allocatable, dimension(:,:,:, :) :: gchem_field                      ! parameter containing the surface fluxes

        ! General parameters for both cloud-to-ground (CTG) and intracloud (IC) lightning
        real :: time_between_flashes = 180                                      ! 3 minutes
        integer :: lightning_time_step                                          ! How many time steps to skip before calculating lightning (function of dt)
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

        real :: cloud_top_height_threshold = 5000.                               ! in meters
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

                enddo
            enddo

            if ( CTG_price_and_rind ) then
                allocate(cloud_top_heights(nx, ny), cloud_top_temps(nx, ny), price_and_rind_flash_rates(nx, ny))
                
                cloud_top_heights = 0.
                cloud_top_temps = 0.

                call calculate_cloud_top_height(cloud_top_heights, cloud_top_temps)

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

                print*, "*************** Maximum CTG lightning = ", MAXVAL(change_in_mixing_ratio)

                gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
            endif

            if ( CTG_decaria_reflectivity ) then
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
                                change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( number_of_20dbz_per_altitude(k) * dx * 20000. * 500. ) ! eventually replace dx with dy; ADDED *2 11/19
                            endif
                        enddo
                    enddo
                enddo

                print*, "*************** Maximum CTG lightning = ", MAXVAL(change_in_mixing_ratio)

                gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
            endif
        endif
    end subroutine lightning_decaria_ctg

    subroutine lightning_decaria_ic(gchem_field)
        ! Parametrizes (as close as possible!) to DeCaria et al. (2000), A cloud-scale model study of lightning-generated NOx...
        use grid, only : z
        use vars, only : tabs, qcl, qpl, qci, qpi                                                   ! tabs (x,y,z) contains absolute temperature in K
        use params, only : land, ocean
        use microphysics, only: dBZ_cloudradar

        implicit none

        ! General parameters for both cloud-to-ground (CTG) and intracloud (IC) lightning
        real :: time_between_flashes = 180                                      ! 3 minutes
        integer :: lightning_time_step                                          ! How many time steps to skip before calculating lightning (function of dt)

        real, allocatable, dimension(:,:,:, :) :: gchem_field                      ! parameter containing the surface fluxes
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

        integer :: i, j, k                                                          ! Counter variables

        lightning_time_step = time_between_flashes / dt

        ! If the desired time between flashes has passed, then do lightning
        if ( mod(nstep, lightning_time_step) == 0 ) then
            allocate(vertical_function_profile(nz), change_in_mixing_ratio(nx, ny, nz))
            change_in_mixing_ratio = 0.

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

            if ( IC_decaria ) then
                allocate(x_area_of_storm(nzm))
                x_area_of_storm = 0.

                do k = 1,nzm 
                do i = 1,nx
                    do j = 1,ny
                        if ( ( (qcl(i, j, k) + qpl(i, j, k) + qci(i, j, k) + qpi(i, j, k)) * 1000. > 0.01 ) .and. ( maxval(dBZ_cloudradar) >= 20 ) ) then
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
                            change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( x_area_of_storm(k) * dx * 20000. * 1000. ) ! eventually replace dx with dy, ADDED *10 11/19
                            endif
                        enddo
                    enddo
                enddo

                print*, "*************** Maximum IC lightning", MAXVAL(change_in_mixing_ratio)

                gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
            endif
        endif
    end subroutine lightning_decaria_ic

    subroutine calculate_cloud_top_height(cloud_top_heights, cloud_top_temps)
        ! Copy of the function that calculates cloud top height in diagnose
        ! This function is required because diagnose is run after emissions/chemistry
        
        use grid, only : adz, dz, z
        use vars, only : qcl, qci, rho, tabs
        
        implicit none

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

end module emissions
