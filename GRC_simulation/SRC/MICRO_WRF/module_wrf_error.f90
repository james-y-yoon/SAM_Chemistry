MODULE module_wrf_error
  INTEGER           :: wrf_debug_level = 0
  CHARACTER*256     :: wrf_err_message

CONTAINS
  subroutine wrf_debug(errnum, error_message)
    implicit none
    integer :: errnum
    CHARACTER(LEN=*) error_message

!bloss    write(*,*) 'Error number = ', errnum
    write(*,*) error_message
  end subroutine wrf_debug

  subroutine wrf_message(message_string)
    implicit none
    CHARACTER(LEN=*) message_string

    write(*,*) message_string
  end subroutine wrf_message

  subroutine wrf_error_fatal(error_message)
    implicit none
    CHARACTER(LEN=*) error_message

    write(*,*) error_message
    call task_abort()
  end subroutine wrf_error_fatal
end MODULE module_wrf_error
