subroutine rowmap(n,f,ifcn,t,u,tend,hs,rtol,atol,itol,&
&jacv,ijacv,fdt,ifdt,solout,iout,work,&
&lwork,iwork,liwork,rpar,ipar,idid)


!-----!-----------------------------------------------------------------
!
!   ROWMAP .. A ROW-code with Krylov techniques for large stiff ODEs.
!
!
!   ROWMAP solves the inital value problem for stiff systems of
! first order ODEs
!
!      du
!      -- = f(t,u),    u=(u(1),..,u(n)), f=(f(1),..f(n)).
!      dt
!
!
!   ROWMAP is based on the ROW-methods of order 4 of the code ROS4 of Hairer
! and Wanner (see [4]) and uses Krylov techniques for the solution of linear
! systems. By a special multiple Arnoldi process the order of the basic
! method is preserved with small Krylov dimensions.  Step size control is
! done by embedding with a method of order 3.
!
! Version of January 10, 1997.
!
!   Authors:
!       H. Podhaisky, R. Weiner,
!       Fachbereich Mathematik und Informatik, Universitaet Halle
!       06099 Halle, Germany
!       email: helmut.podhaisky@mathematik.uni-halle.de
!              ruediger.weiner@mathematik.uni-halle.de,
!
!       B.A. Schmitt,
!       Fachbereich Mathematik, Universitaet Marburg
!       35032 Marburg, Germany
!       email: schmitt@mathematik.uni-marburg.de
!
!
!-----!-----------------------------------------------------------------
! References:
!
! [1]  B.A. Schmitt and R. Weiner:
!         Matrix-free W-methods using a multiple Arnoldi Iteration,
!         APNUM 18(1995), 307-320
!
! [2]  R. Weiner, B.A. Schmitt an H. Podhaisky:
!         ROWMAP - a ROW-code with Krylov techniques for large stiff
!         ODEs. Report 39, FB Mathematik und Informatik,
!         Universitaet Halle, 1996
!
! [3]  R.Weiner and B.A.Schmitt:
!         Consistency of Krylov-W-Methods in initial value problems,
!         Tech. Report 14, FB Mathematik und Informatik,
!         Universitaet Halle, 1995
!
! [4]  E. Hairer and G. Wanner:
!         Solving Ordinary Differential Equations II, Springer-Verlag,
!         Second Edition, 1996
!
!
!  For retrieving the latest version visit the ROWMAP World Wide Web
!  homepage at URL :
!   http://www.mathematik.uni-halle.de/institute/numerik/software/
!
!-----!-----------------------------------------------------------------
!
!   INPUT PARAMETERS
!   ----------------
!
!   n      = Number of equations of the ODE system.
!
!   f      = Name (external) of user-supplied subroutine computing the
!            value f(t,u):
!                subroutine f(n,t,u,udot,rpar,ipar)
!                real*8 t,u(n),udot(n),rpar(*)
!                integer n,ipar(*)
!                udot(1)= ...
!            rpar,ipar (user parameters, see below).
!
!   ifcn   = Switch:
!            ifcn = 0:  f(t,u) independent of t (autonomous)
!            ifcn = 1:  f(t,u) may depend on t (nonautonomous), default
!
!   t      = Initial t-value. Changed on exit.
!
!   u      = Initial values for u. Changed on exit.
!
!   tend   = Final t-value.
!
!   hs     = Initial step size guess. Changed on exit.
!            If hs=0d0 the code puts hs to a default value depending on
!            atol, rtol and u.
!
!   rtol, atol = Relative and absolute error tolerances. They
!            can be both scalars or else both vectors of length n.
!            ROWMAP uses a weighted root-mean-square norm to measure the
!            size of error vectors. It is defined by
!
!                wrms=sqrt(sum[(err(i)/scal(i))**2, i=1,n]/n),
!
!            where
!                 scal(i)=atol+rtol*dabs(u(i))              (itol=0)
!            or
!                 scal(i)=atol(i)+rtol(i)*dabs(u(i))        (itol=1).
!
!   itol   = Switch for rtol,atol:
!            itol=0: both are scalars (default).
!            itol=1: both are vectors.
!
!   jacv   = Name (external) of user-supplied subroutine computing
!            the jacobian-vector product.
!            z := Jac(t,u)*v, Jac=df/du:
!                subroutine jacv(n,t,u,v,z,rpar,ipar)
!                real*8 t,u(n),v(n),z(n),rpar(*)
!                integer n,ipar(*)
!                z(1)=...
!            rpar,ipar (see below).
!            This routine is only called if ijacv = 1.
!            Supply a dummy subroutine in the case ijacv = 0
!
!   ijacv  = Switch for the computation of the jacobian-vector products:
!            ijacv = 0: The products are approximated by finite differences.
!                       Subroutine jacv is never called (default).
!            ijacv = 1: The products are computed by the user-supplied
!                       subroutine jacv.
!
!   fdt    = Name (external) of user-supplied subroutine computing the
!            partial derivate of f(t,u) with respect to t.
!            This routine is only called if idft=1 and ifcn=1.
!
!                subroutine fdt(n,t,u,ft,rpar,ipar)
!                real*8 t,u(n),ft(n),rpar(*)
!                integer ipar(*)
!                ft(1)=...
!
!            rpar,ipar (see below).
!
!   ifdt   = Switch for the computation of df/dt:
!            ifdt = 0: df/dt is computed internally by finite differences,
!                      subroutine fdt is never called (default).
!            ifdt = 1: df/dt is supplied by subroutine fdt.
!
!   solout = Name (external) of user-supplied subroutine exporting the
!            numerical solution during integration. If  iout=1:
!            solout it is called after every successful step. Supply
!            a dummy subroutine if iout = 0. It must have the form:
!
!               subroutine solout(n,told,tnew,uold,unew,fold,fnew,ucon,
!              1                  intr,rpar,ipar)
!               integer n,intr,ipar(*)
!               real*8 told,tnew,uold(n),unew(n),fold(n),fnew(n),ucon(n)
!               real*8 rpar(*)
!               ...
!               end
!
!            "intr" serves to interrupt the integration. If solout sets
!             intr .lt. 0, then ROWMAP returns to the calling program.
!
!            solout may produce continuous output by calling the internal
!            subroutine "rowcon" (see below).
!
!   iout   = Gives information on the subroutine solout:
!               iout = 0: subroutine is never called
!               iout = 1: subroutine is used for output (default)
!
!   work   = Array of working space of length "lwork", changed on exit.
!            Serves as working space for all vectors and matrices.
!            The array is used for optional input. For zero input, the
!            default values are set:
!
!            work(1), work(2) -  Parameters for step size selection.
!               The new step size is chosen subject to the restriction
!                  work(1) <= hsnew/hsold <= work(2)
!                Default values: work(1)=0.25d0, work(2)=2d0
!
!            work(3) - The safety factor in step size prediction.
!                Default value: work(3)=0.8d0
!
!            work(4) - UROUND, the machine precision/rounding unit.
!                Default value: work(4)=1d-16
!
!            work(5) - KTOL, tolerance for the iterative solution of
!                the linear equations. The code keeps the weighted
!                root-mean-square norm of the residual of the first stage
!                below KTOL/HS, where HS is the current step size
!                (see also ATOL/RTOL).
!                Default value: work(5)=1d-1.
!
!   lwork  = Length of array work in calling program.
!            "lwork" must be at least: 10+n*(mx+11)+mx*(mx+4),
!            where mx=iwork(3).
!
!   iwork  = Integer working space of length "liwork". The array is used
!            for optional input. For zero input, default values are set:
!            iwork(1)   -  This is the maximal number of allowed steps.
!                Default value is iwork(1)=10000.
!
!            iwork(2)   -  Switch for the choice of integration method:
!                          (see [4], page 110)
!                           1  Method of Shampine
!                           2  Method GRK4T of Kaps-Rentrop
!                           3  Method of Van Veldhuizen (gamma=1/2)
!                           4  Method of Van Veldhuizen ("D-stable")
!                           5  L-stable Method
!                           6  Method GRK4A of Kaps-Rentrop
!                Default value is iwork(2)=2.
!
!            iwork(3)   -  The maximum Krylov dimension allowed (=mx).
!                Default value is iwork(3)=70.
!
!            iwork(4)   -  The logical output unit number "lun" for
!                messages of ROWMAP. Default value is iwork(4)=6.
!
!
!   liwork = Length of array iwork in calling program.
!            "liwork" must be at least mx+20, where mx is iwork(3).
!
!   rpar, ipar = Real and integer parameters (or parameter arrays)
!            that can be used for communication between your calling
!            program and the subroutines f, fdt, fjac, solout.
!
!   OUTPUT PARAMETERS
!   -----------------
!
!   t     = t-Value where the solution is computed (after successful
!           return is t=tend)
!
!   u     = Solution at t.
!
!   hs    = Predicted step size from the last accepted step.
!
!   idid  = Reports on success upon return:
!         idid = 1    computation sucessful
!         idid = 2    computation interrupted by solout
!         idid = -1   stepsize too small, i.e. hs < 10*uround*dabs(t)
!         idid = -2   more than iwork(1) steps
!         idid = -3   input is not consistent
!         idid = -4   internal Krylov matrix is repeatedly singular
!
!   iwork(5) = Number of computed (accepted and rejected) steps.
!   iwork(6) = Number of rejected steps.
!   iwork(7) = Number of function evaluations.
!   iwork(8) = Number of jacobian-times-vector products.
!   iwork(9) = Minimum lenght of array work.
!   iwork(10)= Minimum lenght of array iwork.
!
! ROWMAP uses the following basic linear algebra modules (BLAS):
!       Level 1: DAXPY,DCOPY,DNORM2,DROTG,DROT,DDOT,DSCAL
!       Level 2: DGEMV,DTRSV
!
!-----!-----------------------------------------------------------------
   implicit none
   integer n,lwork,liwork,iwork(liwork)
   integer itol,ijacv,ifdt,iout,ipar(*) ,ifcn,idid
   integer pmq,pmh,pmak1,pmak2,pmak3,pmak4,pmft,pmuu,pmfu0
   integer pmrhs,pml1,pml2,pml3,pml4,pmfm,pmcon,pmscal,pmmax,mx
   double precision rpar(*)
   double precision t,u(n),tend,hs,rtol(*),atol(*),work(lwork),ktol
   external fdt,f,solout,jacv
!
! Globals.
!
   double precision uround
   integer cifcn,citol,cijacv,cifdt,ciout,method,lun
   common /rowmap0/uround,cifcn,citol,cijacv,cifdt,ciout,method,lun

!     For statistics.
   integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3)
   common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim
!
! Parameters for step size control.
!
   integer nstpx
   double precision fac1,fac2,fac3
   common  /rowmap2 /fac1,fac2,fac3,nstpx
!
! Pointers in array "work".
!

!$ OMP THREADPRIVATE(/rowmap0/,/rowmap1/,/rowmap2 /)

   mx=70
   if (iwork(3).gt.19) mx=iwork(3)
   pmq=10
   pmh=pmq+n*mx
   pmak1=pmh+mx*mx
   pmak2=pmak1+n
   pmak3=pmak2+n
   pmak4=pmak3+n
   pmft=pmak4+n
   pmuu=pmft+n
   pmfu0=pmuu+n
   pmrhs=pmfu0+n
   pml1=pmrhs+n
   pml2=pml1+mx
   pml3=pml2+mx
   pml4=pml3+mx
   pmfm=pml4+mx
   pmcon=pmfm+n
   pmscal=pmcon+n
   pmmax=pmscal+n
   iwork(9)=pmmax
!
! Check lwork und liwork.
!
   idid=0
   if (lwork .lt. pmmax) then
      write (lun,9950) lwork, pmmax
9950  format ('Error in ROWMAP: "lwork"= ',i10,' too small. ',&
      &'Minimum is',i10,'.')
      idid=-3
   end if
   iwork(10)=mx+20
   if (liwork.lt.iwork(10)) then
      write (lun,9951) liwork,iwork(10)
9951  format ('Error in ROWMAP: "liwork"=',i6,' too small. ',&
      &'Minimum is',i6,'.')
      idid=-3
   end if
   if (idid.eq.-3) return
!
! Check optional input.
!
   nstpx=10000
   if (1.lt.iwork(1)) nstpx=iwork(1)
   fac1=work(1)
   fac2=work(2)
   fac3=work(3)
   if (fac1.ge.1.or.fac1.lt.1d-2) fac1=0.25d0
   if (fac2.le.1.or.fac2.gt.1d+2) fac2=2d0
   if (fac3.lt.1d-1.or.fac3.ge.1d0) fac3=0.8d0
! The rounding unit.
   uround=1d-16
   if (work(4).gt.1d-40.and.work(4).lt.1d-5) uround=work(4)
   ktol=1d-1
   if (work(5).gt.0) ktol=work(5)
!     default ROW-method is GRK4T (method=2)
   method=2
   if (1.le.iwork(2).and.6.ge.iwork(2)) method=iwork(2)
   cifcn=ifcn
   citol=itol
   cijacv=ijacv
   cifdt=ifdt
   ciout=iout
!     logical output unit number for messages
   lun=6
   if (iwork(4).gt.0) lun=iwork(4)

!
! Call to core integrator.
!

   call rowmapc(n,f,t,u,tend,hs,rtol,atol,jacv,fdt,solout,&
   &ktol,work(pmq),work(pmh),&
   &work(pmak1),work(pmak2),&
   &work(pmak3),work(pmak4), work(pmft),&
   &work(pmuu),work(pmfu0),&
   &work(pmrhs),work(pml1),&
   &work(pml2),&
   &work(pml3),work(pml4),&
   &work(pmfm),work(pmcon),work(pmscal),&
   &iwork(20),mx,rpar,ipar,idid)
!
! Statistics:
!

   iwork(5) = nsteps
   iwork(6) = nstepsr
   iwork(7) = nfeval
   iwork(8) = njacv
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   R O W C O N
!-----!-----------------------------------------------------------------
! This subroutine is for dense output. It can be called by the
! user-supplied subroutine "solout". A Hermite-interpolation
! is used to compute an approximation (order 3) "ucon" to the solution
! of the ODE system at time s. The value s should lie in the
! interval [told,tnew].
!
! Input
! -----
!  n    = Dimension of the ODE system.
!  s    = Point of evaluation.
!  uold = Numerical solution for the ODE at time told, passed by solout.
!  fold = Value of the right hand side at time told, passed by solout.
!  unew = Numerical solution for the ODE at time tnew, passed by solout.
!  fnew = Value of the right hand side at time tnew, passed by solout.
!
! Output
! ------
!  ucon = Approximation at time s.
!
subroutine rowcon(n,s,told,tnew,uold,unew,fold,fnew,ucon)
   implicit none
   integer n
   DOUBLE PRECISION s,told,tnew,uold(n),unew(n),fold(n),fnew(n)
   DOUBLE PRECISION ucon(n),theta,theta2,thetam1,hs,g
   hs=tnew-told
   theta=(s-told)/hs
   theta2=theta**2d0
   thetam1=theta-1d0
   call dcopy(n,0d0,0,ucon,1)
   g=theta2*(3d0-2d0*theta)
   call daxpy(n,1d0-g,uold,1,ucon,1)
   call daxpy(n,g,unew,1,ucon,1)
   g=hs*theta*thetam1**2d0
   call daxpy(n,g,fold,1,ucon,1)
   g=hs*theta2*thetam1
   call daxpy(n,g,fnew,1,ucon,1)
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   R O W M A P C
!-----!-----------------------------------------------------------------
! Here is the core integrator. ROWMAPC prepares the linear stage equations
! and calls subroutine "STAGE" for solving these. ROWMAPC contains the
! step size control.
!
subroutine rowmapc(n,f,t,u,tend,hs,rtol,atol,jacv,fdt,solout,&
&ktol,q,h,ak1,ak2,ak3,ak4,ft,uu,fu0,rhs,l1,l2,&
&l3,l4,fm,ucon,scal,kl,mx,rpar,ipar,idid)
   implicit none
   integer mx,n,kl(mx),asv,mk,ipar(*),idid,i,j,intr,info,nsing
   double precision t,u(n),tend,hs,rtol(*),atol(*),fm(n),ts,ktol,&
   &ktol1
   double precision q(n,mx),h(mx,mx),ak1(n),ak2(n),ak3(n),ak4(n),&
   &ft(n),uu(n), fu0(n),l1(mx),l2(mx),l3(mx),l4(mx),&
   &rhs(n) , dnrm2,unorm,ucon(n),scal(n),&
   &rpar(*),delt,ehg,err,hnew,told
   external f, solout,fdt,jacv

   double precision a21,a31,a32,c21,c31,c32,c41,c42,c43,b1,b2,b3,b4,&
   &e1,e2,e3,&
   &e4,gamma,c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43

   logical reached

!
! Common blocks.
!

!     Globals.
   DOUBLE PRECISION uround
   integer ifcn,itol,ijacv,ifdt,iout,method,lun
   common /rowmap0/uround,ifcn,itol,ijacv,ifdt,iout,method,lun
!     For statistics.
   integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3)
   common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim

!     For step size control.
   integer nstpx
   DOUBLE PRECISION fac1,fac2,fac3
   common  /rowmap2 /fac1,fac2,fac3,nstpx

!     "unorm" is used in subroutine "ROWDQ".
   common /rownorm/unorm

!$ OMP THREADPRIVATE(/rowmap0/,/rowmap1/,/rowmap2 /,/rownorm/)

!
! Initializations.
!
   call rowcoe(method,a21,a31,a32,c21,c31,c32,c41,c42,c43,&
   &b1,b2,b3,b4,e1,e2,e3,e4,gamma,&
   &c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43)
   intr=0
   nfeval=0
   nsteps=0
   nstepsr=0
   njacv=0
   nfdt=0
   nsing=0
   reached=.false.

   do 101 i=1,4
      nkrydim(i,1)=mx
      nkrydim(i,2)=0
      nkrydim(i,3)=0
101 continue
   call f(n,t,u,fm,rpar,ipar)
   nfeval=nfeval+1
   idid=1
!
! Begin integration, loop for successful steps:
!
1  if (t.ge.tend.or.reached) return
   unorm=dnrm2(n,u,1)/dsqrt(dfloat(n))
   if (itol.ne.1) then
!       atol,rtol are scalars
      do 105 i=1,n
         scal(i)=1d0/(atol(1)+rtol(1)*dabs(u(i)))
105   continue
   else
!       atol,rtol are vectors
      do 106 i=1,n
         scal(i)=1d0/(atol(i)+rtol(i)*dabs(u(i)))
106   continue
   end if
   if (nsteps.eq.0.and.hs.le.0d0) then
!          Set initial step size.
      hs=dnrm2(n,scal,1)/dsqrt(dfloat(n))
      hs=(1d0/hs)**0.25*1d-1
   end if
   if ((t+hs*1.005).ge.tend) then
      hs=tend-t
      reached=.true.
   end if
!
! Compute the derivative f_t in the nonautonomous case.
!
   if (ifcn.ne.0) then
      nfdt=nfdt+1
      if (ifdt.ne.1) then
!         finite differences:
         delt=dsqrt(uround*max(1d-5,dabs(t)))
         call f(n,t+delt,u,ft,rpar,ipar)
         nfeval=nfeval+1
         call daxpy(n,-1d0,fm,1,ft,1)
         call dscal(n,1d0/delt,ft,1)
      else
!         user-supplied subroutine:
         call fdt(n,t,u,ft,rpar,ipar)
      end if
   end if
!
! Label for rejected steps:
!
2  ehg=1.D0/(gamma*hs)
   ktol1=ktol/hs
   ts=t
   asv=1
   mk=0
!     Reset the target indices.
   do 102 j=1,mx
      kl(j)=0
102 continue
   if (nsteps.ge.nstpx) then
      write (lun,9990) t
      idid=-2
      write (lun,9992) nstpx
      return
   end if
   nsteps=nsteps+1
!
! 1. stage.
!

   call dcopy (n,fm,1,rhs,1)

   call dscal (n,ehg,rhs,1)
   if (ifcn.ne.0) call daxpy(n,ehg*hs*d1,ft,1,rhs,1)
   call stage(q,h,rhs,l1,kl,asv,mk,mx,n,fm,fu0,u,t,1,ktol1,&
   &ehg,ak1,scal,f,jacv,fdt,rpar,ipar,info)
   if (info.lt.0) goto 3
!
! 2. stage.
!
   call dcopy(n,u,1,uu,1)
   call daxpy(n,hs*a21,ak1,1,uu,1)
   ts=t+c2*hs
   call f(n,ts,uu,rhs,rpar,ipar)
   nfeval=nfeval+1
   call daxpy(n,c21,ak1,1,rhs,1)
   call dscal (n,ehg,rhs,1)
   if (ifcn.ne.0) call daxpy(n,ehg*hs*d2,ft,1,rhs,1)
   call stage(q,h,rhs,l2,kl,asv,mk,mx,n,fm,fu0,u,t,2,ktol1,&
   &ehg,ak2,scal,f,jacv,fdt,rpar,ipar,info)
   if (info.lt.0) goto 3
   call daxpy(n,-c21,ak1,1,ak2,1)
!
! 3. stage.
!
   call dcopy(n,u,1,uu,1)
   call daxpy(n,hs*a31,ak1,1,uu,1)
   call daxpy(n,hs*a32,ak2,1,uu,1)
   ts=t+c3*hs
   call f(n,ts,uu,rhs,rpar,ipar)
   nfeval=nfeval+1
   call dcopy(n,rhs,1,ak4,1)
   call daxpy(n,c31,ak1,1,rhs,1)
   call daxpy(n,c32,ak2,1,rhs,1)
   call dscal(n,ehg,rhs,1)
   if (ifcn.ne.0) call daxpy(n,ehg*hs*d3,ft,1,rhs,1)
   call stage(q,h,rhs,l3,kl,asv,mk,mx,n,fm,fu0,u,t,3,ktol1,&
   &ehg,ak3,scal,f,jacv,fdt,rpar,ipar,info)
   if (info.lt.0) goto 3
   call daxpy(n,-c31,ak1,1,ak3,1)
   call daxpy(n,-c32,ak2,1,ak3,1)
!
! 4. stage.
!
   call dcopy(n,ak4,1,rhs,1)
   call daxpy(n,c41,ak1,1,rhs,1)
   call daxpy(n,c42,ak2,1,rhs,1)
   call daxpy(n,c43,ak3,1,rhs,1)
   call dscal(n,ehg,rhs,1)
   if (ifcn.ne.0) call daxpy(n,ehg*hs*d4,ft,1,rhs,1)
   call stage(q,h,rhs,l4,kl,asv,mk,mx,n,fm,fu0,u,t,4,ktol1,&
   &ehg,ak4,scal,f,jacv,fdt,rpar,ipar,info)
   if (info.lt.0) goto 3
   nsing=0
   call daxpy(n,-c41,ak1,1,ak4,1)
   call daxpy(n,-c42,ak2,1,ak4,1)
   call daxpy(n,-c43,ak3,1,ak4,1)
!
! New solution: uu = u + sum (h* b.i * ak.i,i=1..4).
!
   call dcopy (n,u,1,uu,1)
   call daxpy (n,hs*b1,ak1,1,uu,1)
   call daxpy (n,hs*b2,ak2,1,uu,1)
   call daxpy (n,hs*b3,ak3,1,uu,1)
   call daxpy (n,hs*b4,ak4,1,uu,1)
!
! Embedded solution: fu0 = sum (hs* e.i * ak.i, i=1..4).
!
   call dcopy (n,0d0,0,fu0,1)
   call daxpy (n,hs*e1,ak1,1,fu0,1)
   call daxpy (n,hs*e2,ak2,1,fu0,1)
   call daxpy (n,hs*e3,ak3,1,fu0,1)
   call daxpy (n,hs*e4,ak4,1,fu0,1)
!
! Error estimate, step size control.
!
   err=0d0
   do 550 i=1,n
      err=err+(fu0(i)*scal(i))**2
550 continue
   err=dsqrt(err/n)
   hnew=hs*dmin1(fac2,dmax1(fac1,(1d0/err)**0.25D0 * fac3))
   if (1d-1*hnew.le.dabs(t)*uround) then
      write (lun,9990) t
      idid=-1
      write (lun,9991)
      return
   end if
   if (err .lt. 1d0 ) then
!
! Step is accepted.
!
      told=t
      t=t+hs
      call dcopy(n,fm,1,fu0,1)
      call f(n,t,uu,fm,rpar,ipar)
      nfeval=nfeval+1
      if (iout.ne.0)&
      &call solout(n,told,t,u,uu,fu0,fm,ucon,intr,rpar,ipar)
      call dcopy (n,uu,1,u,1)
      if (intr.lt.0) then
         idid=2
         write (*,9990) t
         write (*,9993)
         return
      end if
      hs=hnew
      goto 1
   else
!
! Step is rejected.
!
      nstepsr=nstepsr+1
      reached=.false.
      hs=hnew
      goto 2
   end if
!
! Matrix is singular.
!
3  nsing=nsing+1
   write (lun,9995)
   if (nsing.ge.3) then
      write (lun,9990) ts
      write (lun,9994)
      idid=-5
      return
   end if
   nstepsr=nstepsr+1
   hs=hs/2d0
   goto 2

9990 format ('Exit of ROWMAP at t=',d10.3)
9991 format ('Step size too small.')
9992 format ('More than ',i7,' steps needed.')
9993 format ('Computation interrupted by "solout".')
9994 format ('Matrix is repeatedly singular')
9995 format ('Warning: Matrix is singular.')
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   S T A G E
!-----!-----------------------------------------------------------------
! The subroutine "STAGE" solves the linear equations by using the multiple
! Arnoldi process, see [1,2,3]. The first step incorporates the new right
! hand side into the Krylov subspace. Then the Krylov subspace is extended
! until the residual is small enough or the maximal dimensions are reached.
! Finally the solution "ak" is computed as a linear combination of
! orthogonalized Krylov vectors.
!
!
subroutine stage(q,h,rhs,l,kl,asv,mk,mx,n,fm,fu0,u,t,st,ktol1,&
&ehg,ak,scal,f,jacv,fdt,rpar,ipar,info)
   implicit none
   integer mx,n,mk,i,asv,m,ikrdim,info
   integer kl(mx),st,krydim(4,4),ipar(*)
   logical done
   DOUBLE PRECISION q(n,mx),h(mx,mx),rhs(n),l(mx),dnrm2,ktol1,nrmv
   DOUBLE PRECISION fm(n),fu0(n),u(n),t,ehg,dfkt,nrmrhs
   DOUBLE PRECISION ak(n),scal(n),rpar(*)
   external f,jacv,fdt
!
!     Maximal Krylov dimension increments per stage, see [3]
!     (autonomous/nonautonomous). For the methods shamp, grk4t,
!     veldd, velds and lstab (b3.eq.0)
!
   data (krydim(i,1),i=1,4) /-12,3,1,3/
   data (krydim(i,2),i=1,4) /-16,5,1,5/
!
!    ... and for grk4a (b3.ne.0):
!
   data (krydim(i,3),i=1,4) /-15,3,4,3/
   data (krydim(i,4),i=1,4) /-19,5,4,5/

!     Globals.
   DOUBLE PRECISION uround
   integer ifcn,itol,ijacv,ifdt,iout,method,lun
   common /rowmap0/uround,ifcn,itol,ijacv,ifdt,iout,method,lun

!     Statistics.
   integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3)
   common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim

!$ OMP THREADPRIVATE(/rowmap0/,/rowmap1/)

   ikrdim=2*(method/6)+ifcn+1
   info=0
   nrmrhs=dnrm2(n,rhs,1)
   if (mk.eq.0) then
      call dcopy(mx,0d0,0,l,1)
   else
      call orthov(n,mk,q,rhs,l,mx)
   end if
   nrmv=dnrm2(n,rhs,1)
   if (nrmv/(1d-12+nrmrhs) .ge. 1d-8.or.st.lt.4) then
!
! Insert the new (orthogonalized) right hand side "rhs" as
! new column in Q.
!
      if (mk.gt.0) asv=asv+1
!        Shift last vectors by 1.
      do 100 i=asv-1,1,-1
         call dcopy(n,q(1,mk+i),1,q(1,mk+i+1),1)
100   continue
      do 101 i=1,mk
         if(kl(i).ge.mk+1) kl(i)=kl(i)+1
101   continue
      l(mk+1)=nrmv
      call dcopy(n,rhs,1,q(1,mk+1),1)
      call dscal(n,1d0/l(mk+1),q(1,mk+1),1)
   end if
   m=0
   done=.false.
!
! Begin Multiple Arnoldi process.
!
1000 m=m+1
   if (m.gt.mk) then
      kl(m)=m+asv
   end if
   call kryarn(n,t,u,mk,m,kl,q,h,mx,ehg,l,fm,fu0,f,jacv,rpar,ipar)
   if (m.gt.mk-st+asv.and.m.ge.asv.and.dabs(h(m,m)).lt.uround) then
!        matrix is singular
      info=-1
      return
   end if
   if (m.gt.mk-st+asv.and.m.ge.asv) then
      call krdfkt(n,m,kl,q,h,mx,l,fu0,dfkt,st,asv,ehg,scal)
   else
      dfkt=1d100
   end if
   if (dfkt.lt.ktol1.and.st.gt.1) done=.true.
   if (dfkt.lt.ktol1.and.m.ge.4) done=.true.
   if(st.gt.1.and.(m-mk).ge.krydim(st,ikrdim)) done=.true.
   if(st.eq.1.and.m.ge.krydim(1,ikrdim)+mx) done=.true.
   if (.not. done) goto 1000
!
! End of multiple Arnoldi process for current stage.
!
   mk=m
!     Used Krylov dimensions.
   nkrydim(st,1)=min0(nkrydim(st,1),mk)
   nkrydim(st,2)=max0(nkrydim(st,2),mk)
   nkrydim(st,3)=nkrydim(st,3)+mk
!     Backsubstitution.
   call dtrsv('u','n','n',m,h,mx,l,1)
!     Compute ak = Q*l.
   call dgemv('n',n,mk,1d0,q,n,l,1,0d0,ak,1)
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   O R T H O V
!-----!-----------------------------------------------------------------
! This subroutine orthogonalizes v with respect to vectors q(.,i),i=1,m,
! by a modified Gram-Schmitt-Algorithm. Dotproducts are stored in "l",
! i.e. l=Q'v.
!
!
subroutine orthov(n,m,q,v,l,mx)
   implicit none
   integer n,m,mx,i
   DOUBLE PRECISION q(n,*), v(n),l(mx),s,ddot
   do 100 i=1,m
      s=ddot(n,q(1,i),1,v,1)
      l(i)=s
      call daxpy(n,-s,q(1,i),1,v,1)
100 continue
   call dcopy(mx-m,0d0,0,l(m+1),1)
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   R E O R T N
!-----!-----------------------------------------------------------------
! The vector q(.,m) is incorporated in the Krylov space, it is scaled to
! norm one with h(m,j), where kl(j)=m. The candidates for further subspace
! extentions have to be re-orthogonalized to q(.,m).
!
subroutine reortn(n,m,kl,q,h,mx)
   implicit none
   integer n,m,kl(*),mx,j0,j
   DOUBLE PRECISION q(n,*),h(mx,*),s,ddot
   j0=m
   do 100 j=m-1,1,-1
100 if (kl(j).ge.m) j0=j
   do 150 j=j0,m-1
      if (kl(j).eq.m.and.h(m,j).gt.1d-14)&
      &call dscal(n,1d0/h(m,j),q(1,m),1)
150 continue
   do 200 j=j0,m-1
      if (kl(j).gt.m) then
         s=ddot(n,q(1,m),1,q(1,kl(j)),1)
         h(m,j)=s
         call daxpy(n,-s,q(1,m),1,q(1,kl(j)),1)
      end if
200 continue
   return
end
!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   K R Y A R N
!-----!-----------------------------------------------------------------
! This subroutine extends the Krylov subspace by one vector and updates
! the QR-decomposition of matrix "h".
!

subroutine kryarn(n,t,u,mk,m,kl,q,h,mx,ehg,l,fm,fu0,f,&
&jacv,rpar,ipar)
   implicit none
   integer n,mk,m,kl(*),mx, i,j,j0

   DOUBLE PRECISION uround
   integer ifcn,itol,ijacv,ifdt,iout,method,lun
   common /rowmap0/uround,ifcn,itol,ijacv,ifdt,iout,method,lun

   integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3),ipar(*)
   common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim
   DOUBLE PRECISION  u(n),t,q(n,*),h(mx,*),ehg,l(*), c,s,fm(n),fu0(n)
   DOUBLE PRECISION  rpar(*)
   external f,jacv

!$ OMP THREADPRIVATE(/rowmap0/,/rowmap1/)


   if (m.gt.mk) call reortn(n,m,kl,q,h,mx)
   l(m) = -l(m)
   j0 = m
   do 50 j=m-1,1,-1
50 if (kl(j).ge.m) j0=j
!     generate new Givens rotations
   do 100 j=j0,m-1
      if (m .gt. mk) then
         call drotg(h(j,j),h(m,j),c,s)
         call drot(m-j-1,h(j,j+1),mx,h(m,j+1),mx,c,s)
      else
!        reuse old rotations
         call ztocs(h(m,j),c,s)
      endif
!       apply to right hand side "l"
      call drot(1,l(j),1,l(m),1,c,s)
100 continue

!
! Multiple Arnoldi-Step with QR-Decomposition.
!
   if (m.le.mk) return
!      Krylov-Step: q(.,kl(m)) = (I - QQ') A q(.,m)
   if (ijacv.ne.1) then
      call rowdq(n,t,u,q(1,m),fu0,fm,q(1,kl(m)),f,rpar,ipar)
      nfeval=nfeval+1
   else
      call jacv(n,t,u,q(1,m),q(1,kl(m)),rpar,ipar)
   end if
   njacv=njacv+1
   call orthov(n,m,q,q(1,kl(m)),h(1,m),mx)
   h(m,m) = h(m,m)-ehg
!     update the QR-decomposition of matrix "h"
   do 210 i=2,m
      do 200 j=i-1,1,-1
200   if (kl(j).ge.i) j0=j
      do 210 j=j0,i-1
         call ztocs(h(i,j),c,s)
         call drot(1,h(j,m),1,h(i,m),1,c,s)
210 continue
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   K R D F K T
!-----!-----------------------------------------------------------------
! This subroutine computes the weighted root-mean-square norm of the residual
!
!               r = (I-gh A) Q l - rhs.
!
subroutine krdfkt(n,m,kl,q,h,mx,r,fu0,dn,st,asv,egh,scal)
   implicit none
   integer n,m,mx,kl(mx),i,j,li,ls,mi,st,asv,j0
   DOUBLE PRECISION q(n,*),h(mx,*),r(m),fu0(*),dn,s,qn,dnrm2,egh,gh
   DOUBLE PRECISION scal(n),wrmsq

   dn = 0D0
   gh=1d0/egh
   j0 = m
   do 10 j=m,1,-1
10 if (kl(j).gt.m) j0=j
   ls = m-j0+1

   do 110 i=0,ls-1
      mi = m-i
      s = r(mi)
      do 100 j=0,i-1
100   s = s-h(mi,m-j)*fu0(ls-j)
      s = s/h(mi,mi)
      fu0(ls-i) = s
      s = dabs(s)
      li = kl(mi)
      if (li.gt.m) then
         qn = dnrm2(n,q(1,li),1)
         h(m+1,mi) = qn
         wrmsq=0d0
         do 101 j=1,n
            wrmsq=wrmsq+(q(j,li)*scal(j))**2
101      continue
         wrmsq=dsqrt(wrmsq/n)
         dn = dn + s*wrmsq
      endif
110 continue
   dn=dn*gh
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   Z T O C S
!-----!-----------------------------------------------------------------
! Reconstructs old Givens rotations computed by DROTG.
!
subroutine ztocs(z,c,s)
   implicit none
   DOUBLE PRECISION z,c,s
   if (dabs(z).lt.1d0) then
      s=z
      c=dsqrt(1d0-s*s)
   else if (dabs(z).gt.1d0) then
      c=1d0/z
      s=dsqrt(1d0-c*c)
   else
      c=0d0
      s=1d0
   end if
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   R O W D Q
!-----!-----------------------------------------------------------------
! This subroutine approximates the Jacobian-vector product
! by finite differences
!
!     v = (f(u+delta y)-f(u))/delta = A y + O(||u||delta^2)
!

subroutine rowdq(n,t,u,y,fu0,fm,v,f,rpar,ipar)
   implicit none
   integer n,ipar(*)
   DOUBLE PRECISION  t,u(n),y(n),fu0(n),fm(n),v(n),delta
   DOUBLE PRECISION  rpar(*),eddelta,unorm
   external f
   common /rownorm/ unorm

!$ OMP THREADPRIVATE(/rownorm/)

   delta=1.d-7*dmax1(1d-5,unorm)
   eddelta=1d0/delta
   call dcopy(n,u,1,fu0,1)
   call daxpy(n,delta,y,1,fu0,1)
   call f(n,t,fu0,v,rpar,ipar)
   call daxpy(n,-1d0,fm,1,v,1)
   call dscal(n,eddelta,v,1)
   return
end

!-----!-----------------------------------------------------------------
! --- S U B R O U T I N E   R O W C O E
!-----!-----------------------------------------------------------------
! This subroutine loads the coefficients of the chosen method.
!
subroutine rowcoe(method,a21,a31,a32,c21,c31,c32,c41,c42,c43,&
&b1,b2,b3,b4,e1,e2,e3,e4,gamma,&
&c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43)
   implicit none
   integer method
   DOUBLE PRECISION a21,a31,a32,c21,c31,c32,c41,c42,c43,b1,b2,b3,b4,&
   &e1,e2,e3,e4,gamma,c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43
   goto (1,2,3,4,5,6),method
10 continue
   c21=g21/gamma
   c31=g31/gamma
   c32=g32/gamma
   c41=g41/gamma
   c42=g42/gamma
   c43=g43/gamma
   c2=a21
   c3=a31+a32
   d1=gamma
   d2=gamma+g21
   d3=gamma+g31+g32
   d4=gamma+g41+g42+g43
   e1=e1-b1
   e2=e2-b2
   e3=e3-b3
   e4=e4-b4
   return

! -- SHAMP
1  a21= 1.00000000000000D0
   a31=0.48000000000000D0
   a32=0.12000000000000D0
   g21=-2.00000000000000D0
   g31=1.32000000000000D0
   g32=0.60000000000000D0
   g41=-0.05600000000000D0
   g42= -0.22800000000000D0
   g43= -0.10000000000000D0
   b1=0.29629629629630D0
   b2=0.12500000000000D0
   b3=0D0
   b4=0.57870370370370D0
   e1=0D0
   e2= -0.04166666666667D0
   e3= -0.11574074074074D0
   e4= 1.15740740740741D0
   GAMMA= 0.50000000000000D0
   goto 10
! --  GRK4T
2  g21=-0.27062966775244D0
   g31=0.31125448329409D0
   g32=0.00852445628482D0
   g41=0.28281683204353D0
   g42=-0.45795948328073D0
   g43=-0.11120833333333D0
   b1=0.21748737165273D0
   b2=0.48622903799012D0
   b3=0.00000000000000D0
   b4=0.29628359035715D0
   a21=0.46200000000000D0
   a31=-0.08156681683272D0
   a32=0.96177515016606D0
   e1=-0.71708850449933D0
   e2=1.77617912176104D0
   e3=-0.05909061726171D0
   e4=0.00000000000000D0
   GAMMA=0.23100000000000D0
   goto 10
! -- VELDS
3  g21=-2.00000000000000D0
   g31=-1.00000000000000D0
   g32=-0.25000000000000D0
   g41=-0.37500000000000D0
   g42=-0.37500000000000D0
   g43= 0.50000000000000D0
   b1=0.16666666666667D0
   b2=0.16666666666667D0
   b3=0D0
   b4=0.66666666666667D0
   a21=1.00000000000000D0
   a31=0.37500000000000D0
   a32=0.12500000000000D0
   e1=1.16666666666667D0
   e2=0.50000000000000D0
   e3=-0.66666666666667D0
   e4=0D0
   GAMMA=0.50000000000000D0
   goto 10
! -- VELDD
4  g21=-0.27170214984937D0
   g31=0.20011014796684D0
   g32=0.09194078770500D0
   g41=0.35990464608231D0
   g42=-0.52236799086101D0
   g43=-0.10130100942441D0
   b1=0.20961757675658D0
   b2=0.48433148684810D0
   b3=0.0D0
   b4=0.30605093639532D0
   a21=0.45141622964514D0
   a31=-0.15773202438639D0
   a32=1.03332491898823D0
   e1=-0.74638173030838D0
   e2=1.78642253324799D0
   e3=-0.04004080293962D0
   e4=0.0D0
   GAMMA=0.22570811482257D0
   goto 10
! -- LSTAB
5  g21=-2.34201389131923D0
   g31=-0.02735980356646D0
   g32=0.21380314735851D0
   g41=-0.25909062216449D0
   g42=-0.19059462272997D0
   g43=-0.22803686381559D0
   b1=0.32453574762832D0
   b2=0.04908429214667D0
   b3=0.00000000000000D0
   b4=0.62637996022502D0
   a21=1.14564000000000D0
   a31=0.52092209544722D0
   a32=0.13429476836837D0
   e1=0.61994881642181D0
   e2=0.19268272217757D0
   e3=0.18736846140061D0
   e4=0D0
   GAMMA=0.57282000000000D0
   goto 10
! -- GRK4A
6  a21=0.43800000000000D0
   a31=0.79692045793846D0
   a32=0.07307954206154D0
   g21=-0.76767239548409D0
   g31=-0.85167532374233D0
   g32=0.52296728918805D0
   g41=0.28846310954547D0
   g42=0.08802142733812D0
   g43=-0.33738984062673D0
   b1=0.19929327570063D0
   b2=0.48264523567374D0
   b3=0.06806148862563D0
   b4=0.25000000000000D0
   e1=0.34632583375795D0
   e2=0.28569317571228D0
   e3=0.36798099052978D0
   e4= 0.0D0
   GAMMA=0.39500000000000D0
   goto 10
end
