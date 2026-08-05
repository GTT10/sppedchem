program driver_plog_real_debug
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use working_precision, only: dp
   use chemistry_setup, only: mechdir, use_speedchem
   use speedchem, only: ns, nr, neq, specie, uMW, isumnudiff,        &
                        molefr_to_massfr
   use kinetics_mod, only: reaction_rates, mass_action_productories
   use reacpar, only: nEQREV, iEQREV, n_plog_reactions,             &
                      plog_reaction, plog_kinf_eval, uequilC
   use troepar, only: ntbFALL
   use SCthermodata, only: temperature_array
   use SCmixturethermo, only: SCrho, SCP, rhoY, molar_volumes,       &
                              pressurerhoT
   use speedchem_conV, only: SC_conV
   use sparse_chemistry, only: stoich_r_sp, stoich_p_sp
   use sparse_algebra, only: sparse_value
   implicit none

   character(len=256) :: env
   integer :: env_len, env_stat, i, ieq, ir, ih2, io2, iar, ih, iho2
   real(dp), allocatable :: x(:), ymass(:), state(:), rhs(:), C(:)
   real(dp), allocatable :: kinf(:), k0(:), prod_f(:), prod_b(:), uKc(:)
   real(dp) :: Ta(6), kf, kb, qf, qb, q

   call get_environment_variable('SC_MECHDIR', env, env_len, env_stat)
   if (env_stat /= 0 .or. env_len == 0) error stop 1
   mechdir = trim(env)
   if (mechdir(len_trim(mechdir):len_trim(mechdir)) /= '/')           &
      mechdir = trim(mechdir)//'/'
   use_speedchem = .true.
   call chemistry_input

   if (n_plog_reactions /= 1) error stop 1
   ir = plog_reaction(1)
   ih2=species_index('H2'); io2=species_index('O2'); iar=species_index('AR')
   ih=species_index('H'); iho2=species_index('HO2')
   if (min(ih2,io2,iar,ih,iho2) == 0) error stop 1

   allocate(x(ns),ymass(ns),state(neq),rhs(neq),C(ns),kinf(nr),       &
            k0(ntbFALL),prod_f(nr),prod_b(nr),uKc(nEQREV))
   x=0.0_dp
   x(ih2)=2.0_dp; x(io2)=1.0_dp; x(iar)=7.0_dp
   x(ih)=1.0e-5_dp; x(iho2)=1.0e-8_dp
   x=x/sum(x)
   call molefr_to_massfr(x,ymass)
   state(1)=1000.0_dp
   state(2:neq)=ymass
   SCP=10.0_dp*101325.0_dp
   call rhoY(ymass,state(1))
   call molar_volumes
   Ta=temperature_array(state(1))
   C=1.0e-3_dp*SCrho*uMW*ymass

   call reaction_rates(Ta,k0,kinf)
   call plog_kinf_eval(Ta,pressurerhoT(state(1),ymass),kinf)
   call mass_action_productories(C,prod_f,prod_b)
   uKc=uequilC(Ta)
   ieq=0
   do i=1,nEQREV
      if (iEQREV(i) == ir) ieq=i
   enddo
   if (ieq == 0) error stop 1

   kf=kinf(ir)
   kb=kf*uKc(ieq)
   qf=prod_f(ir)*kf
   qb=prod_b(ir)*kb
   q=qf-qb

   write(*,'(a,i0,a,a)') 'DEBUG reaction=',ir,' equation-index species=', &
      trim(specie(ih))//'+'//trim(specie(io2))//'+'//trim(specie(iar))// &
      '<=>'//trim(specie(iho2))//'+'//trim(specie(iar))
   write(*,'(a,es24.16)') 'DEBUG T=',state(1)
   write(*,'(a,es24.16)') 'DEBUG P=',pressurerhoT(state(1),ymass)
   write(*,'(a,es24.16)') 'DEBUG rho=',SCrho
   write(*,'(a,i0)') 'DEBUG delta_nu=',isumnudiff(ir)
   write(*,'(a,5(1x,es24.16))') 'DEBUG C H O2 AR HO2=',               &
      C(ih),C(io2),C(iar),C(iho2),sum(C)
   write(*,'(a,8(1x,es12.4))') 'DEBUG stoich rH rO2 rAR rHO2 pH pO2 pAR pHO2=', &
      sparse_value(stoich_r_sp,ir,ih),sparse_value(stoich_r_sp,ir,io2), &
      sparse_value(stoich_r_sp,ir,iar),sparse_value(stoich_r_sp,ir,iho2), &
      sparse_value(stoich_p_sp,ir,ih),sparse_value(stoich_p_sp,ir,io2), &
      sparse_value(stoich_p_sp,ir,iar),sparse_value(stoich_p_sp,ir,iho2)
   write(*,'(a,es24.16,1x,l1)') 'DEBUG kinf=',kf,ieee_is_finite(kf)
   write(*,'(a,es24.16,1x,l1)') 'DEBUG uKc=',uKc(ieq),ieee_is_finite(uKc(ieq))
   write(*,'(a,es24.16,1x,l1)') 'DEBUG kb=',kb,ieee_is_finite(kb)
   write(*,'(a,es24.16,1x,l1)') 'DEBUG prod_f=',prod_f(ir),ieee_is_finite(prod_f(ir))
   write(*,'(a,es24.16,1x,l1)') 'DEBUG prod_b=',prod_b(ir),ieee_is_finite(prod_b(ir))
   write(*,'(a,es24.16,1x,l1)') 'DEBUG qf=',qf,ieee_is_finite(qf)
   write(*,'(a,es24.16,1x,l1)') 'DEBUG qb=',qb,ieee_is_finite(qb)
   write(*,'(a,es24.16,1x,l1)') 'DEBUG q=',q,ieee_is_finite(q)

   call SC_conV(neq,0.0_dp,state,rhs)
   write(*,'(a,4(1x,es24.16),4(1x,l1))') 'DEBUG RHS T H O2 HO2=',     &
      rhs(1),rhs(ih+1),rhs(io2+1),rhs(iho2+1),                       &
      ieee_is_finite(rhs(1)),ieee_is_finite(rhs(ih+1)),               &
      ieee_is_finite(rhs(io2+1)),ieee_is_finite(rhs(iho2+1))

contains
   integer function species_index(name) result(idx)
      character(len=*), intent(in) :: name
      integer :: k
      idx=0
      do k=1,ns
         if (trim(specie(k)) == name) then
            idx=k
            return
         endif
      enddo
   end function species_index
end program driver_plog_real_debug
