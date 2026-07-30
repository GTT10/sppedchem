program driver_cklink_read
   use chemistry_setup, only: mechdir
   use speedchem, only: ns, nr
   implicit none
   character(len=256) :: env
   integer :: env_len, env_stat

   call get_environment_variable('SC_MECHDIR',env,env_len,env_stat)
   if (env_stat /= 0 .or. env_len == 0) error stop 1
   mechdir = trim(env)
   if (mechdir(len_trim(mechdir):len_trim(mechdir)) /= '/')           &
      mechdir = trim(mechdir)//'/'
   call SCcklink
   write(*,'(a,i0,a,i0)') 'RESULT: PASS - cklink read ns=',ns,' nr=',nr
end program driver_cklink_read
