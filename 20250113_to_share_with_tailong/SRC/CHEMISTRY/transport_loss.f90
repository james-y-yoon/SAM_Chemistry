module transport_loss

    use grid, only : nx, ny, nzm, dx, dy, dz
    use vars, only : u, v, dtn

    CONTAINS

    subroutine transport_removal_driver(gchem_field)
        implicit none

        real, allocatable, dimension(:) :: average_u_winds          ! Horizontal zonal winds, in m/s
        real, allocatable, dimension(:) :: average_v_winds          ! Horizontal meridional winds, in m/s
        real :: tau_x, tau_y, inverse_tau_total
        real :: i, j, k, n
        real, allocatable, dimension(:,:,:,:) :: gchem_field

        allocate(average_u_winds(nzm), average_v_winds(nzm))
        call calculate_large_scale_horizontal_winds(average_u_winds, average_v_winds)

        do k = 1, nzm
            tau_x = ( nx * dx ) / average_u_winds(k)            ! Transport lifetime in the zonal direction
            tau_y = ( ny * dy ) / average_v_winds(k)            ! Transport lifetime in the meridional direction

            inverse_tau_total = ( 1 / tau_x ) + ( 1 / tau_y )       ! Lifetimes add in parallel; note that inverse of tau_total is just k_total

            gchem_field(:, :, k, :) = ( 1. - ( dtn * inverse_tau_total / ( nx * ny ) ) ) * gchem_field(:, :, k, :)      ! First order loss, but divide by total number of gridboxes in each k
        end do
    end subroutine transport_removal_driver

    subroutine calculate_large_scale_horizontal_winds(average_u_winds, average_v_winds)
        implicit none
        real :: k

        real, allocatable, dimension(:) :: average_u_winds, average_v_winds

        do k = 1, nzm
            average_u_winds(k) = SUM( u(:, :, k) ) / SIZE( u(:, :, k) )
            average_v_winds(k) = SUM( v(:, :, k) ) / SIZE( v(:, :, k) )
        end do

    end subroutine calculate_large_scale_horizontal_winds

end module transport_loss