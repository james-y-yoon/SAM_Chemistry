  ! Wrapper call to avoid circular dependency between hbuffer and chemistry  modules

subroutine hbuf_chem_init(namelist,deflist,unitlist,status,average_type,count,trcount)
   use chemistry
   implicit none
   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount

   call chem_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

end
