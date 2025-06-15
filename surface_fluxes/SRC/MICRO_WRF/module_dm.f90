MODULE module_dm
  use grid, only: rank, nsubdomains
  use mpi

  !SAM: Severly truncated version of module_dm from WRFv4.2.1/external/RSL_LITE/

   IMPLICIT NONE

CONTAINS

     SUBROUTINE wrf_dm_gatherv_double ( v, elemsize , km_s, km_e )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e
      REAL*8 v(0:*)
    !SAM#ifndef STUBMPI
    !SAM# ifndef USE_MPI_IN_PLACE
      REAL*8 v_local((km_e-km_s+1)*elemsize)
    !SAM# endif
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs
      INTEGER send_type, myproc, nproc, local_comm, ierr, i
!    INCLUDE 'mpif.h'
      send_type = MPI_DOUBLE_PRECISION
    !SAM    CALL wrf_get_dm_communicator ( local_comm )
    !SAM    CALL wrf_get_nproc( nproc )
    !SAM    CALL wrf_get_myproc( myproc )
    !SAM -- begin
    local_comm = MPI_COMM_WORLD
    nproc = nsubdomains
    myproc = rank
    !SAM -- end
      ALLOCATE( recvcounts(nproc), displs(nproc) )
      i = (km_e-km_s+1)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,local_comm,ierr) ;
      i = (km_s)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,displs,1,MPI_INTEGER,local_comm,ierr) ;
    !SAM# ifdef USE_MPI_IN_PLACE
    !SAM    CALL mpi_allgatherv( MPI_IN_PLACE,                                  &
    !SAM         # else
      DO i = 1,elemsize*(km_e-km_s+1)
        v_local(i) = v(i+elemsize*km_s-1)
      END DO
      CALL mpi_allgatherv( v_local,                                       &
         !SAM         # endif
                           (km_e-km_s+1)*elemsize,                        &
                           send_type,                                     &
                           v,                                             &
                           recvcounts,                                    &
                           displs,                                        &
                           send_type,                                     &
                           local_comm,                                    &
                           ierr )
      DEALLOCATE(recvcounts)
      DEALLOCATE(displs)
    !SAM#endif
      return
     END SUBROUTINE wrf_dm_gatherv_double

     SUBROUTINE wrf_dm_gatherv_single ( v, elemsize , km_s, km_e )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e
      REAL*4 v(0:*)
    !SAM#ifndef STUBMPI
    !SAM# ifndef USE_MPI_IN_PLACE
      REAL*4 v_local((km_e-km_s+1)*elemsize)
    !SAM# endif
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs
      INTEGER send_type, myproc, nproc, local_comm, ierr, i
!    INCLUDE 'mpif.h'
      send_type = MPI_REAL
    !SAM    CALL wrf_get_dm_communicator ( local_comm )
    !SAM    CALL wrf_get_nproc( nproc )
    !SAM    CALL wrf_get_myproc( myproc )
    !SAM -- begin
    local_comm = MPI_COMM_WORLD
    nproc = nsubdomains
    myproc = rank
    !SAM -- end
      ALLOCATE( recvcounts(nproc), displs(nproc) )
      i = (km_e-km_s+1)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,local_comm,ierr) ;
      i = (km_s)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,displs,1,MPI_INTEGER,local_comm,ierr) ;
    !SAM# ifdef USE_MPI_IN_PLACE
    !SAM    CALL mpi_allgatherv( MPI_IN_PLACE,                                  &
    !SAM         # else
      DO i = 1,elemsize*(km_e-km_s+1)
        v_local(i) = v(i+elemsize*km_s-1)
      END DO
      CALL mpi_allgatherv( v_local,                                       &
         !SAM         # endif
                           (km_e-km_s+1)*elemsize,                        &
                           send_type,                                     &
                           v,                                             &
                           recvcounts,                                    &
                           displs,                                        &
                           send_type,                                     &
                           local_comm,                                    &
                           ierr )
      DEALLOCATE(recvcounts)
      DEALLOCATE(displs)
    !SAM#endif
      return
     END SUBROUTINE wrf_dm_gatherv_single

      SUBROUTINE wrf_dm_decomp1d( nt, km_s, km_e )
       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: nt
       INTEGER, INTENT(OUT) :: km_s, km_e
     ! local
       INTEGER nn, nnp,  na, nb
       INTEGER myproc, nproc

    !SAM CALL wrf_get_myproc(myproc)
    !SAM CALL wrf_get_nproc(nproc)
    !SAM -- begin
    myproc = rank 
    nproc = nsubdomains 
    !SAM -- end
       nn = nt / nproc           ! min number done by this task
       nnp = nn
       if ( myproc .lt. mod( nt, nproc ) )   nnp = nnp + 1 ! distribute remainder

       na = min( myproc, mod(nt,nproc) ) ! Number of blocks with remainder that precede this one
       nb = max( 0, myproc - na )        ! number of blocks without a remainder that precede this one
       km_s = na * ( nn+1) + nb * nn     ! starting iteration for this task
       km_e = km_s + nnp - 1             ! ending iteration for this task
      END SUBROUTINE wrf_dm_decomp1d

  !----------------------------------------------------------------------
  
  subroutine task_bcast_fourdim_array_real8(rank_from,array)
    implicit none
    !        include 'mpif.h'

    integer, intent(in) :: rank_from          ! broadcasting task's rank
    real(8), intent(inout) :: array(:,:,:,:) ! array to be broadcast

    integer :: n1,n2,n3,n4        ! dimension lengths
    integer ierr

    real(8), allocatable :: rtmp(:)
    integer :: count, nn, myrank, i,j,k,m
    real(8) :: tmp1, tmp2

    ! get individual dimension lengths
    n1 = SIZE(array,1)
    n2 = SIZE(array,2)
    n3 = SIZE(array,3)
    n4 = SIZE(array,4)

    nn = n1*n2*n3*n4
    allocate(rtmp(nn),STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'Error in allocating array in task_bcast_fourdim_array'
      call task_abort()
    end if

    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    if(myrank.eq.rank_from) then
      count = 1
      do m = 1,n4
        do k = 1,n3
          do j = 1,n2
            do i = 1,n1
              rtmp(count) = array(i,j,k,m)
              count = count + 1
            end do
          end do
        end do
      end do
    end if

    call task_bcast_real8(rank_from,rtmp,nn)

    if(myrank.ne.rank_from) then
      count = 1
      do m = 1,n4
        do k = 1,n3
          do j = 1,n2
            do i = 1,n1
              array(i,j,k,m) = rtmp(count)
              count = count + 1
            end do
          end do
        end do
      end do
    end if

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,992) myrank, tmp1, tmp2
!!$        992 format('Consistency check in task_bcast_fourdim: rank/sum/sum2 = ',I4,2E16.8)


  end subroutine task_bcast_fourdim_array_real8

  !----------------------------------------------------------------------

  subroutine task_bcast_fivedim_array_real4(rank_from,array)
    implicit none
    !        include 'mpif.h'

    integer, intent(in) :: rank_from          ! broadcasting task's rank
    real(4), intent(inout) :: array(:,:,:,:,:) ! array to be broadcast

    integer :: n1,n2,n3,n4,n5        ! dimension lengths
    integer ierr

    real(4), allocatable :: rtmp(:)
    integer :: count, nn, myrank, i1,i2,i3,i4,i5
    real :: tmp1, tmp2

    n1 = SIZE(array,1)
    n2 = SIZE(array,2)
    n3 = SIZE(array,3)
    n4 = SIZE(array,4)
    n5 = SIZE(array,5)

    nn = n1*n2*n3*n4*n5
    allocate(rtmp(nn),STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'Error in allocating array in task_bcast_fivedim_array_real4'
      call task_abort()
    end if

    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    if(myrank.eq.rank_from) then
      count = 1
      do i5 = 1,n5
        do i4 = 1,n4
          do i3 = 1,n3
            do i2 = 1,n2
              do i1 = 1,n1
                rtmp(count) = array(i1,i2,i3,i4,i5)
                count = count + 1
              end do
            end do
          end do
        end do
      end do
    end if

    call task_bcast_float4(rank_from,rtmp,nn)

    if(myrank.ne.rank_from) then
      count = 1
      do i5 = 1,n5
        do i4 = 1,n4
          do i3 = 1,n3
            do i2 = 1,n2
              do i1 = 1,n1
                array(i1,i2,i3,i4,i5) = rtmp(count)
                count = count + 1
              end do
            end do
          end do
        end do
      end do
    end if

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,993) myrank, tmp1, tmp2
!!$        993 format('Consistency check in task_bcast_fivedim: rank/sum/sum2 = ',I4,2E16.8)

  end subroutine task_bcast_fivedim_array_real4

END MODULE module_dm
