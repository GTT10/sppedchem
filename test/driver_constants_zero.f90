program constants_zero
 use working_precision, only: dp
 use chemistry_setup, only: mechdir,use_speedchem
 use speedchem, only: ns,nr,specie,Einf,E0,Erev
 use reacpar, only: plog_EoverR
 use sparse_chemistry, only: third_body_sp
 use sparse_algebra, only: sparse_to_dense
 use chemkinII, only: cklen,ckinit,ckrp
 use universal_constants, only: R,Rcal,Rerg,Patm
 implicit none
 real(dp),parameter :: refR=8.31446261815324_dp,ea=10000._dp
 real(dp),allocatable :: dense(:,:),rw(:)
 integer,allocatable :: iw(:)
 character(len=18),allocatable :: cw(:)
 real(dp) :: ru,ruc,pa,expected
 integer :: i,n2,h2,h2o,li,lr,lc,marker,failures,linksize
 character(len=8):: unit
 character(len=:),allocatable :: legacy_blob
 mechdir='./';use_speedchem=.true.
 call get_command_argument(1,unit)
 select case(trim(unit))
 case('KELV');expected=ea
 case('CAL');expected=ea*4.184_dp/refR
 case('KCAL');expected=ea*4184._dp/refR
 case('JOULES');expected=ea/refR
 case('KJOULES');expected=ea*1000._dp/refR
 case default;error stop 3
 end select
 call chemistry_input
 n2=0;h2=0;h2o=0;failures=0
 do i=1,ns
  if(trim(specie(i))=='N2')n2=i
  if(trim(specie(i))=='H2')h2=i
  if(trim(specie(i))=='H2O')h2o=i
 enddo
 if(min(n2,h2,h2o)==0)error stop 3
 allocate(dense(third_body_sp%nr,third_body_sp%nc))
 call sparse_to_dense(third_body_sp,dense)
 call check('explicit zero efficiency',dense(2,n2),0._dp,0._dp)
 call check('default efficiency',dense(2,h2),1._dp,0._dp)
 call check('nondefault efficiency',dense(2,h2o),2._dp,0._dp)
 call check('universal gas constant',R,refR,1e-15_dp)
 call check('ordinary E/R',Einf(1)/Rcal,expected,1e-13_dp)
 call check('explicit REV E/R',Erev(1)/Rcal,expected,1e-13_dp)
 call check('falloff E/R',E0(3)/Rcal,expected,1e-13_dp)
 call check('PLOG E/R node1',plog_EoverR(1),expected,1e-13_dp)
 call check('PLOG E/R node2',plog_EoverR(2),expected,1e-13_dp)
 open(92,file='cklink',form='unformatted',status='old',access='stream')
 read(92)marker
 inquire(unit=92,size=linksize)
 allocate(character(len=linksize-marker-8)::legacy_blob)
 read(92,pos=marker+9)legacy_blob
 close(92)
 open(92,file='legacy_cklink',form='unformatted',status='replace',access='stream')
 write(92)legacy_blob
 close(92)
 open(92,file='legacy_cklink',form='unformatted',status='old')
 open(93,file='ckinit_test.log',status='replace')
 call cklen(92,93,li,lr,lc)
 allocate(iw(li),rw(lr),cw(lc))
 call ckinit(li,lr,lc,92,93,iw,rw,cw)
 call ckrp(iw,rw,ru,ruc,pa)
 close(92);close(93)
 call check('CKRP RU',ru,refR*1e7_dp,1e-15_dp)
 call check('CKRP RUC',ruc,refR/4.184_dp,1e-15_dp)
 call check('CKRP PA',pa,1013250._dp,0._dp)
 if(failures>0)error stop 2
 write(*,'(a)')'RESULT: PASS constants and zero collider'
contains
 subroutine check(name,got,ref,tol)
 character(len=*),intent(in)::name
 real(dp),intent(in)::got,ref,tol
 logical::pass
 pass=abs(got-ref)<=tol*max(abs(ref),1._dp)
 write(*,'(a,1x,l1,2(1x,es25.17e3))')trim(name),pass,got,ref
 if(.not.pass)failures=failures+1
 end subroutine
end program
