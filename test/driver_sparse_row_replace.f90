program test_sparse_row_replace
 use working_precision,only:dp
 use sparse_definitions,only:sparse,allocate
 use sparse_algebra,only:add_line,add_value,sparse_value
 implicit none
 type(sparse)::a
 real(dp)::row(6),expected(6,6)
 integer::i,j,failures,k,nexpected
 failures=0
 call allocate(6,6,0,a);expected=0
 ! Subsequent rows emulate isolated explicit colliders.
 call add_value(a,2,5,1._dp);expected(2,5)=1
 call add_value(a,4,6,1._dp);expected(4,6)=1
 call add_value(a,6,2,1._dp);expected(6,2)=1
 call verify('insertion setup')
 row=[1._dp,0._dp,2._dp,0._dp,3._dp,0._dp]
 call add_line(a,1,row);expected(1,:)=row
 call verify('insert earlier general row')
 ! Replacement larger/equal/smaller while explicit rows follow.
 row=[1._dp,2._dp,3._dp,4._dp,0._dp,5._dp]
 call add_line(a,1,row);expected(1,:)=row
 call verify('replace larger')
 row=[0._dp,2._dp,3._dp,4._dp,5._dp,6._dp]
 call add_line(a,1,row);expected(1,:)=row
 call verify('replace equal')
 row=[0._dp,0._dp,0._dp,0._dp,7._dp,0._dp]
 call add_line(a,1,row);expected(1,:)=row
 call verify('replace smaller')
 ! Existing empty row insertion must not displace explicit colliders.
 row=[0._dp,8._dp,0._dp,0._dp,0._dp,9._dp]
 call add_line(a,3,row);expected(3,:)=row
 call verify('insert existing empty row')
 if(failures>0)then
  print *,'SPARSE_ROW_REPLACE_FAIL',failures
  stop 1
 endif
 print *,'SPARSE_ROW_REPLACE_PASS'
contains
 subroutine verify(label)
 character(len=*),intent(in)::label
 integer::bad
 bad=0;nexpected=count(expected/=0)
 if(a%n/=nexpected)bad=bad+1
 if(a%IA(1)/=1.or.a%IA(a%nr+1)/=a%n+1)bad=bad+1
 if(any(a%IA(2:)<a%IA(:a%nr)))bad=bad+1
 do i=1,6
 do j=1,6
  if(sparse_value(a,i,j)/=expected(i,j))bad=bad+1
 enddo
 enddo
 if(bad>0)then
  failures=failures+bad
  print *,'FAIL ',label,' errors=',bad,' n=',a%n,' expected=',nexpected
 else
  print *,'PASS ',label
 endif
 end subroutine
end program
