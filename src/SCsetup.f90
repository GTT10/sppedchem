!   ***********************************************************************************************
!   **                                                                                           **
!   **         |\   -  -.   ./                                                                   **
!   **         | \./ \/ | ./ /     _____                     __________  __________  ___         **
!   **       __|         /  /     / ___/____  ___  ___  ____/ / ____/ / / / ____/  |/  /         **
!   **       \    .        /.-/   \__ \/ __ \/ _ \/ _ \/ __  / /   / /_/ / __/ / /|_/ /          **
!   **        \   |\.|\/|    /   ___/ / /_/ /  __/  __/ /_/ / /___/ __  / /___/ /  / /           **
!   **         \__\     /___/   /____/ .___/\___/\___/\__,_/\____/_/ /_/_____/_/  /_/            **
!   **                              /_/                                                          **
!   **                                                                                           **
!   **       SpeedCHEM - A fast, portable Fortran library for Chemical Kinetics problems         **
!   **                          http://www.federicoperini.info/speedchem                         **
!   **                 Copyright (C) 2010-2013 Federico Perini. Version 1.0                      **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!   ** License                                                                                   **
!   ** This file is part of SpeedCHEM.                                                           **
!   **                                                                                           **
!   ** SpeedCHEM is free software: you can redistribute it and/or modify                         **
!   ** it under the terms of the GNU General Public License as                                   **
!   ** published by the Free Software Foundation, either version 3 of the                        **
!   ** License, or any later version.                                                            **
!   **                                                                                           **
!   ** SpeedCHEM is distributed in the hope that it will be useful,                              **
!   ** but WITHOUT ANY WARRANTY; without even the implied warranty of                            **
!   ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             **
!   ** GNU General Public License for more details.                                              **
!   **                                                                                           **
!   ** You should have received a copy of the GNU General Public License                         **
!   ** along with SpeedCHEM. If not, see <http://www.gnu.org/licenses/>.                         **
!   **                                                                                           **
!   ** Please send comments, requests, bug reports to federico.perini@gmail.com                  **
!   **                                                                                           **
!   ***********************************************************************************************

!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **           Mechanism and thermodynamic data input            **
!     **                                                             **
!     **   Input files required: SpeedCHEM.dat                       **
!     **                         SCthermo.dat                        **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: tuesday, 25/05/2010                          **
!     **                                                             **
!     *****************************************************************

subroutine SCsetup

!     Current Working precision
   use working_precision

!     Module containing mechanism data
!k2015      use chemistry_setup, only: mechanism
   use chemistry_setup, only: mechanism, mechdir
   use speedchem
   use SCthermodata
   use find_mod
   use kinetics_mod
   use troepar
   use reacpar
   use sparse_definitions
   use sparse_chemistry


   implicit none
!     Local variables
   logical :: present
   integer :: k,j,i,np,irow,icol
   integer, dimension(:), allocatable :: ire
   real (dp)  ::  tmp,tmpa,tmpt1,tmpt2,tmpt3

!     Sparse indexing variables
   real (dp)       , dimension(:), allocatable :: vals, dtmp
   integer, dimension(:),          allocatable :: col_ind, indx
   integer :: nvals, nrows


!     Variables related to thermodynamic database
   integer :: nrecords
   logical, dimension(:), allocatable :: spthermo
   integer, dimension(:), allocatable :: ithermo, count_reacs
   real (dp)        ta,tb,tc,td,te,tf,tg,ts
   integer is, stato
   integer :: ne1, ne2, ne3, ne4
   integer :: max_nreacs, max_nprods

   type(sparseint) :: tmp_isp
   integer,          dimension(:),   allocatable :: tmpip, tmpjp
   real (dp)       , dimension(:),   allocatable :: tmpv

!     Character strings
   character*13 nomefile
   character*18 dummy
   character*2  dummy1,dummy2,dummy3,dummy4
   character*80 mechlabel


!     *****************************************************************

!     Open file for data output
!k2015      open(unit=820,file='SpeedCHEM.out',status='unknown')
   open(unit=820,file=trim(mechdir)//"SpeedCHEM.out",status='unknown')

!     Check presence of mechanism input file
!k2015      nomefile = 'SpeedCHEM.dat'
   nomefile = trim(mechdir)//"SpeedCHEM.dat"
   inquire(file=nomefile,exist=present)

   if (present) then
!k2015         open(unit = 800,file='SpeedCHEM.dat')
      open(unit = 800,file=trim(mechdir)//"SpeedCHEM.dat")

!           Reading mechanism label
      read(800,"(A80)")mechlabel
      mechanism = trim(adjustl(mechlabel))

!           Reading number of species and reactions
      read(800,*)dummy,nr
      read(800,*)dummy,ns

!           Allocate species and reaction indices
      allocate(species(ns), reactions(nr))
      species   = [(j,j=1,ns)]
      reactions = [(j,j=1,nr)]

!           Number of problem equations
      neq = ns + 1
      uneq = 1.e0_dp / neq

!           Reading element number
      read(800,*)
      read(800,*)dummy,nel

!           Allocating variables dimensions
      call SCallocate
!            third_body(:,:) = 0.e0_dp
      call allocate(nr,ns,0,third_body_sp)


!           Reading elements names
      do j=1,nel
         read(800,*)elementi(j)
      end do

!           Reading species names and molecular weights
      read(800,*)
      do j=1,ns
         read(800,"(A18,E17.10)")specie(j),SCMW(j)
!              write(*,"(A18,E17.10)")specie(j),SCMW(j)
      end do
      uMW = 1.e0_dp/SCMW
!           Reading reaction data
      read(800,*)
      read(800,*)
!           Reactants stoichiometric coefficients
      read(800,*)np
!            stoich_r(:,:) = 0.0d0
      do j=1,np
         read(800,*)irow,icol,tmp
         call add_value(stoich_r_sp,irow,icol,tmp)
!              stoich_r(irow,icol) = tmp
      end do

!           Fix sparse matrix representation for unactive species
      stoich_r_sp%nc = ns

!           Maximum number of reactants per reaction
      allocate(count_reacs(nr))
      call sparse_internal_count(stoich_r_sp,count_reacs,dim=2)
      max_nreacs = maxval(count_reacs)
!            max_nreacs = maxval(count(stoich_r/=0.e0_dp, dim=2))


!           Products stoichiometric coefficients
      read(800,*)
      read(800,*)np
!            stoich_p(:,:) = 0.0d0
      do j=1,np
         read(800,*)irow,icol,tmp
!              stoich_p(irow,icol) = tmp
         call add_value(stoich_p_sp,irow,icol,tmp)
      end do

!           Fix sparse matrix representation for unactive species
      stoich_p_sp%nc = ns

!           Assign nudiff
      nudiff_sparse = stoich_p_sp - stoich_r_sp

!           Maximum number of products per reaction
      call sparse_internal_count(stoich_p_sp,count_reacs,dim=2)
      max_nprods = maxval(count_reacs)
      deallocate(count_reacs)

!           Reading reactions collision frequencies
      read(800,*)
      do j=1,nr
         read(800,*)A0(j),Ainf(j),Arev(j)
      end do

!           Reading reactions activation energies
      read(800,*)
      do j=1,nr
         read(800,*)E0(j),Einf(j),Erev(j)
      end do

!           Reading reactions temperature exponents
      read(800,*)
      do j=1,nr
         read(800,*)b0(j),binf(j),brev(j)
      end do

!           Reading HIGH pressure limits explicitly
      read(800,*)
      read(800,*)dummy,np
      allocate(ire(nr))
      HIGH(1:nr) = .FALSE.
      if (np.gt.0) then
         read(800,"(8(1X,I4))")(ire(j),j=1,np)
         do j=1,np
!              third_body(ire(j),:) = 1.e0_dp
            call add_line(third_body_sp,ire(j),[(1.e0_dp,i=1,ns)])
            HIGH(ire(j)) = .TRUE.
         end do
      endif
      deallocate(ire)

!           Reading LOW pressure limits explicitly
      read(800,*)
      read(800,*)dummy,np
      allocate(ire(nr))
      LOW(1:nr) = .FALSE.
      if (np.gt.0) then
         read(800,"(8(1X,I4))")(ire(j),j=1,np)
         do j=1,np
!              third_body(ire(j),:) = 1.e0_dp
            call add_line(third_body_sp,ire(j),[(1.e0_dp,i=1,ns)])
            LOW(ire(j)) = .TRUE.
         end do
      endif
      deallocate(ire)

!           Reading TROE parameters
      read(800,*)dummy
      read(800,*)np
      TROE(1:nr) = .FALSE.
      if(np.gt.0) then
         do j=1,np
            read(800,*)irow,tmpa,tmpt1,tmpt2,tmpt3
            TROE(irow) = .TRUE.
            aTROE(irow) = tmpa
            T1TROE(irow) = tmpt1
            T2TROE(irow) = tmpt2
            T3TROE(irow) = tmpt3
         end do
      endif

!           Check which reactions are reversible
      read(800,*)
      read(800,*)dummy,np

      allocate(ire(nr))
      reversibile(1:nr)=.FALSE.
      if(np.gt.0) then
         read(800,"(8(1X,I4))")(ire(j),j=1,np)
         do j=1,np
            reversibile(ire(j)) = .TRUE.
         end do
      endif
      deallocate(ire)

      if (.not.allocated(inotrev)) then
         call find_indices(.not.reversibile)
         nnotrev = size(indices)
         allocate(inotrev(nnotrev))
         inotrev = indices
      endif

!           Check three-body reactions
      read(800,*)
      read(800,*)dummy,np
      THREE(1:nr) = .FALSE.
!            third_body(:,:) = 0.e0_dp
      allocate(ire(nr))
      if (np.gt.0) then
         read(800,"(8(1X,I4))")(ire(j),j=1,np)
         do j=1,np
            THREE(ire(j)) = .TRUE.
!              third_body(ire(j),:) = 1.e0_dp
            call add_line(third_body_sp,ire(j),[(1.e0_dp,i=1,ns)])
         end do
      endif
      deallocate(ire)

!           Allocate THREE indices
      allocate(count_reacs(nr))
      call sparse_internal_count(third_body_sp,count_reacs,dim=2)
      call find_indices(count_reacs>0)

      allocate(iTHREE(size(indices)))
      iTHREE(:) = indices(:)
      nTHREE    = size(indices)

!           Read third-body efficiencies
      read(800,*)
      read(800,*)dummy,np
!            third_body(:,:) = 0.e0_dp
      if (np.gt.0) then
         do j=1,np
            read(800,*)irow,icol,tmp
!              third_body(irow,icol) = tmp
            call add_value(third_body_sp,irow,icol,tmp)
         end do
      endif

!           Manipulate third_body efficiency matrices for speed
!            third_bodyT = transpose(third_body)

!           If coefficients are diminished by one, the reactor config
!           can have a sparse Jacobian since total concentration is
!           constant
!            third_body_beta = 0.e0_dp
      call allocate(nTHREE, ns, 0, tb_beta_sp   )
      do j=1,size(iTHREE)
         i = iTHREE(j)
!               where(third_body(i,:)/=1.e0_dp)
!     &           third_body_beta(:,i) = 1.e0_dp - third_body(i,:)

         do k = 1, ns
            call        add_value   (tb_beta_sp,    j, k,&
            &1.e0_dp - sparse_value(third_body_sp, i, k) )

         end do

      end do

!            third_body_betaT = transpose(third_body_beta)
!           Fix number of columns
      tb_beta_sp%nc = ns

!           Initi alising nudiff (vectorized operation)
!           and completing variables' definition
!            nudiff = stoich_p - stoich_r

      REV = ((Arev.ne.0.0d0).or.(brev.ne.0.0d0).or.&
      &(Erev.ne.0.0d0))

!           ** REVERSE REACTION RATES - related indices ****************

!           Index for simple, third-body reactions
!           (Effective molecularity but not Troe of Lindemann forms)
      call find_indices( (.not.(HIGH.or.LOW)) .and.&
      &(count_reacs > 0) )
      allocate(iTB(size(indices)))
      iTB = indices
      nTB = size(indices)

!            if (.not.allocated(isumnudiff)) then
!               allocate(isumnudiff(nr))
!               isumnudiff = int(sum(nudiff,2))
!            endif
      tmp_isp = nudiff_sparse
      if (.not.allocated(isumnudiff)) then
         allocate(isumnudiff(nr))
         call sparseint_internal_sum(tmp_isp,isumnudiff,2)
!         isumnudiff = int(sum(nudiff, dim = 2))
      endif


      call spdeallocate(tmp_isp)
      deallocate(count_reacs)


      close(800)

   else
      write(*,*)'Missing SpeedCHEM mech input file, SpeedCHEM.dat. '
      write(*,*)'Provide a SpeedCHEM or cklink mechanism and retry.'
      stop
   endif



!     *****************************************************************
!     Importing SpeedCHEM thermodynamic database

!k2015      nomefile = 'SCthermo.dat'
   nomefile = trim(mechdir)//"SCthermo.dat"
   inquire(file=nomefile,exist=present)

   if (present) then
      open(unit=810,file=nomefile,status='old')

      read(810,*)
      read(810,*)
      read(810,*)dummy,nrecords

!       reading records in the thermodynamic database and associating
!       them to the species
      allocate(spthermo(ns),ithermo(ns))
      do j=1,nrecords

         read(810,"(A18)")dummy
         read(810,*)

         is = 0
         do i=1,ns
!             write(*,"(A18)")specie(i)
            if (specie(i).eq.dummy) is = i
         end do


         read(810,*)ta,tb,tc,td
         read(810,*)te,tf,tg

         if (is.gt.0) then
            spthermo(is) = .TRUE.
            ithermo(is) = 1
            aH(is) = ta
            bH(is) = tb
            cH(is) = tc
            dH(is) = td
            eH(is) = te
            fH(is) = tf
            gH(is) = tg
         endif

         read(810,*)dummy,ts

         read(810,*)ta,tb,tc,td
         read(810,*)te,tf,tg

         if(is.gt.0) then
            tsw(is) = ts
            aL(is) = ta
            bL(is) = tb
            cL(is) = tc
            dL(is) = td
            eL(is) = te
            fL(is) = tf
            gL(is) = tg
         endif

      end do


      if (count(spthermo).lt.ns) then
         call find_indices(.not.spthermo)
         write(*,*)'Missing ',ns-count(spthermo),&
         &' species in thermo database: '
         write(*,"(5(1X,A18))")specie(indices)
         call exit(0)
      endif

      deallocate(spthermo)

!       Reading Element data**********************************

      read(810,"(A1)",IOSTAT=stato)dummy

      if (stato.eq.-1) then
         write(*,*)'Element data not present'
         call exit(0)
      else

!       Allocating arrays
         allocate(nels(ns),el1(ns),el2(ns),el3(ns),el4(ns))
         allocate(nel1(ns),nel2(ns),nel3(ns),nel4(ns))

         nel1(:) = 0
         nel2(:) = 0
         nel3(:) = 0
         nel4(:) = 0

         el1(:) = '  '
         el2(:) = '  '
         el3(:) = '  '
         el4(:) = '  '

!       Storing element data
         do j=1,ns
            read(810,940)nels(j),dummy1,ne1,dummy2,ne2,&
            &dummy3,ne3,dummy4,ne4
!          write(*,940)nels(j),dummy1,ne1,dummy2,ne2,
!     &                         dummy3,ne3,dummy4,ne4


            if(nels(j).lt.4) then
               dummy4 = '  '
               ne4 = 0
               if(nels(j).lt.3) then
                  dummy3 = '  '
                  ne3 = 0
                  if(nels(j).lt.2) then
                     dummy2 = '  '
                     ne2 = 0
                  endif
               endif

            endif

            el1(j) = dummy1
            nel1(j) = ne1

            el2(j) = dummy2
            nel2(j) = ne2

            el3(j) = dummy3
            nel3(j) = ne3

            el4(j) = dummy4
            nel4(j) = ne4


         end do


      endif


      close(810)

   else
      write(*,*)'Missing thermodynamic database SCthermo.dat. '
      call exit(0)

   endif

!     Allocating variables for TROE pressure-dependent reactions
   if (.not.allocated(todotroe)) then

      a1 = 1.e0_dp/6.64385618977472436e0_dp
      a2 = 1.e0_dp/6.64385618977472525e0_dp

      call find_indices(TROE)

      allocate(todotroe(size(indices)))
      allocate(aT2     (size(indices)))
      allocate(uT1T2   (size(indices)))
      allocate(T2T2    (size(indices)))
      allocate(uT3T2   (size(indices)))
!      allocate(expuT3T2(size(indices)))
!      allocate(expuT1T2(size(indices)))
      todotroe = indices
      aT2      = aTROE(indices)
      uT1T2    = 1.e0_dp/T1TROE(indices)
      T2T2     = T2TROE(indices)
      uT3T2    = 1.e0_dp/T3TROE(indices)
!      expuT1T2 = exp(uT1T2)
!      expuT3T2 = exp(uT3T2)

      call find_indices(T2T2.eq.0.e0_dp)
      allocate(zeroT2(size(indices)))
      zeroT2 = indices

   endif

!     ******************************************************************
!     Allocating types of reactions in module reacpar

   call find_indices(.not.LOW.and..not.HIGH.and..not.TROE)
   allocate(Arrhreac(size(indices)))
   Arrhreac = indices

   call find_indices((LOW.or.HIGH).and.TROE)
   allocate(Troereac(size(indices)))
   Troereac = indices

   call find_indices((LOW.or.HIGH).and.(.not.TROE))
   allocate(Lindreac(size(indices)))
   Lindreac = indices

   call find_indices(REV)
   allocate(Revreac(size(indices)))
   Revreac = indices



!     PRODUCTS ********************************************************
!     (data is only needed for reversible reactions)


!     jp, ip = row and column indexes (in stoich_p matrix) referring
!     to non-zero stoichiometric coefficients of species active in
!     reversible reactions
   allocate(ip(stoich_p_sp%n),jp(stoich_p_sp%n),&
   &v_stoich_p(stoich_p_sp%n))
   call extract_rowcol_indices_columnwise(stoich_p_sp,ip,jp,&
   &v_stoich_p)

!     NB: not all the reactions are reversible! save time by
!         reducing the mass action productories arrays to the
!         reversible reactions only (index is ip)
   if (count(reversibile(ip)) < size(ip)) then
      allocate(tmpip(count(reversibile(ip))))
      allocate(tmpjp(count(reversibile(ip))))
      allocate(tmpv (count(reversibile(ip))))

      tmpip = pack(ip,         reversibile(ip))
      tmpjp = pack(jp,         reversibile(ip))
      tmpv  = pack(v_stoich_p, reversibile(ip))

      deallocate(ip,jp,v_stoich_p)
      allocate(ip(size(tmpip)), jp(size(tmpjp)),&
      &v_stoich_p(size(tmpip)))

      ip         = tmpip
      jp         = tmpjp
      v_stoich_p = tmpv

      deallocate(tmpip, tmpjp, tmpv)

   endif

!      call find_indices2D(stoich_p > 0.0e0_dp .and.
!     &                    spread(reversibile, 2, ns))
!      allocate(ip(size(i2D1)),jp(size(i2D2)))
!
!      ip = i2D1
!      jp = i2D2
!
!      if (size(ip) > 0) then
!         allocate(v_stoich_p(size(ip)))
!         v_stoich_p = pack(stoich_p,stoich_p>0.and.
!     &                    spread(reversibile, 2, ns))

!         call find_indices(v_stoich_p == 1.e0_dp)
!         allocate(v_stoich_p1(size(indices)))
!         v_stoich_p1 = indices
!
!         call find_indices(v_stoich_p == 2.e0_dp)
!         allocate(v_stoich_p2(size(indices)))
!         v_stoich_p2 = indices
!         num_vp2 = size(indiceS)
!
!         call find_indices(v_stoich_p.ne.1.e0_dp.and.v_stoich_p.ne.2.e0_dp)
!         allocate(v_stoich_po(size(indices)))
!         v_stoich_po = indices
!         num_vpo = size(indices)

!      endif

   allocate(indice_p(nr,max_nprods))
   indice_p(:,:) = 0

   if (size(ip)>0) then
      do i=1,nr
         call find_indices(ip.eq.i)
         if (size(indices)>0)indice_p(i,1:size(indices)) = indices(:)
      end do
   endif

   call find_indices2D(indice_p > 0)
   allocate(i1p(nr,max_nprods))
   if (count(indice_p>0)>0)allocate(i2p(size(i2D1)))

   i1p = 0
   where (indice_p>0) i1p = 1

!      i1p(:,:) = 0
!      do i=1,size(indice_p,1)
!         do j=1,size(indice_p,2)
!          if (indice_p(i,j).gt.0) i1p(i,j) = 1
!         end do
!      end do


   if (allocated(i2p)) then

      i2p = pack(indice_p,indice_p>0)

!     Reorder v_stoich_p
      allocate(iv_stoich_p(size(v_stoich_p)),ijp(size(v_stoich_p)))

      do j = 1, size(v_stoich_p)
         iv_stoich_p(j) = int(v_stoich_p(i2p(j)))
         ijp(j) = jp(i2p(j))
      end do

   endif

   if (.not.allocated(i2D1p)) then
      call find_indices2D(i1p.eq.1)
      allocate(i2D1p(size(i2D1)))
      i2D1p = i2D1
      allocate(i2D2p(size(i2D2)))
      i2D2p = i2D2
   endif

!     REACTANTS *******************************************************

!      call find_indices2D(stoich_r > 0.0e0_dp)
!      allocate(ir(size(i2D1)),jr(size(i2D2)))
!      ir = i2D1
!      jr = i2D2
!
!      allocate(v_stoich_r(size(ir)))
!      v_stoich_r = pack(stoich_r,stoich_r>0)

   allocate(ir(stoich_r_sp%n),jr(stoich_r_sp%n))
   allocate(v_stoich_r(size(ir)))
   call extract_rowcol_indices_columnwise(stoich_r_sp,ir,jr,&
   &v_stoich_r)


!      call find_indices(v_stoich_r == 1.e0_dp)
!      allocate(v_stoich_r1(size(indices)))
!      v_stoich_r1 = indices
!
!      call find_indices(v_stoich_r == 2.e0_dp)
!      allocate(v_stoich_r2(size(indices)))
!      v_stoich_r2 = indices
!      num_vr2 = size(indices)
!
!      call find_indices(v_stoich_r.ne.1.e0_dp.and.v_stoich_r.ne.2.e0_dp)
!      allocate(v_stoich_ro(size(indices)))
!      v_stoich_ro = indices
!      num_vro = size(indices)


   allocate(indice_r(nr,max_nreacs))
   indice_r(:,:) = 0

   do i=1,stoich_r_sp%nr
      call find_indices(ir.eq.i)
      indice_r(i,1:size(indices)) = indices(:)
   end do

   call find_indices2D(indice_r > 0)
   allocate(i1r(nr,max_nreacs))
   allocate(i2r(size(i2D1)))

   i1r(:,:) = 0
   do i=1,size(indice_r,1)
      do j=1,size(indice_r,2)
         if (indice_r(i,j).gt.0) i1r(i,j) = 1
      end do
   end do

!      do i=1,size(i1r,1)
!      write(*,"(5(I5))")i,indice_r(i,1),indice_r(i,2),i1r(i,1),i1r(i,2)
!      end do
!      pause

   i2r = pack(indice_r,indice_r>0)

!     Reorder v_stoich_r
   allocate(iv_stoich_r(size(v_stoich_r)),ijr(size(v_stoich_r)))

   do j = 1, size(v_stoich_r)
      iv_stoich_r(j) = int(v_stoich_r(i2r(j)))
      ijr(j) = jr(i2r(j))
   end do



   if (.not.allocated(i2D1r)) then
      call find_indices2D(i1r.eq.1)
      allocate(i2D1r(size(i2D1)))
      i2D1r = i2D1
      allocate(i2D2r(size(i2D2)))
      i2D2r = i2D2
   endif



!c     Sort values of ijp, i2D1p, iv_stoich_p for faster computation
!      allocate(dtmp(size(i2D1p)), indx(size(i2D1p)) )
!      dtmp  = dble(i2D1p)
!      indx = [(j,j=1,size(dtmp))]
!      call qsort(dtmp, size(dtmp), indx)
!
!      i2D1p = int(dtmp)
!
!      dtmp = dble(ijp)
!      ijp  = int(dtmp(indx))
!      dtmp = dble(iv_stoich_p)
!      iv_stoich_p = int(dtmp(indx))
!
!      deallocate(dtmp, indx)

!c     Sort values of ijp, i2D1p, iv_stoich_p for faster computation
!      allocate(dtmp(size(i2D1r)), indx(size(i2D1r)) )
!      dtmp  = dble(i2D1r)
!      indx = [(j,j=1,size(dtmp))]
!      call qsort(dtmp, size(dtmp), indx)
!
!      i2D1r = int(dtmp)
!
!      dtmp = dble(ijr)
!      ijr  = int(dtmp(indx))
!      dtmp = dble(iv_stoich_r)
!      iv_stoich_r = int(dtmp(indx))
!
!      deallocate(dtmp, indx)





!     Preparing output file containing mechanism data

   write(820,"(A80)")'***********************************************&
   &********************************'
   write(820,"(A80)")'**                    SpeedCHEM mechanism outpu&
   &t                             **'
   write(820,"(A80)")'***********************************************&
   &********************************'
   write(820,*)
   write(820,*)mechlabel
   write(820,*)
   write(820,900)nr
   write(820,910)ns
   write(820,920)nel
   write(820,*)
   write(820,925)max_nreacs
   write(820,926)max_nprods
   write(820,*)
   write(820,*)'Elements: '
   write(820,"(1X,50(A2,2X))")(elementi(j),j=1,nel)
   write(820,*)
   write(820,*)'Species and molecular weights [g/mol]: '
   do j=1,ns
      write(820,930)specie(j),SCMW(j)
   end do
   write(820,*)



!     ** Writing TSS-solver-related information ***********************

!      write(820,*)'TSS solver-related information*************************
!     &************************************************'

!         write(820,951)'isolver',isolver
!         write(820,951)'ifast',ifast
!         write(820,950)'dt_def',dt_def
!         write(820,950)'threshold',qfqb_threshold
!         write(820,950)'TSSfac',MTSfac
!         write(820,950)'TSSTOL',TSSTOL
!         write(820,950)'RTOL',TOLR
!         write(820,950)'ATOL',TOLA
!         write(820,950)
!         write(820,950)'thugeTS',thugeTS
!         write(820,951)'TSSinj1',TSSinj1
!         write(820,951)'TSSinj2',TSSinj2

!      close(820)
!     *****************************************************************

!      isolver = 0

!      write(*,960)nr,ns
!      write(*,970)

900 format(' Number of reactions:  ',I4)
910 format(' Number of species  :  ',I4)
920 format(' Number of elements :  ',I4)
925 format(' Max number of reactants:  ',I4)
926 format(' Max number of products :  ',I4)
930 format(1X,A18,2X,E16.8)
940 format(22X,I1,1X,4(A2,1X,I2,1X))
950 format(A12,E13.5)
951 format(A12,I13)
960 format(' chemistry link completed: ',I5,' reactions, ',I3,&
   &' species')
970 format(' ------------------------------------------------------')

end subroutine SCsetup




!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **          Build numerical dense Jacobian matrix              **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Last update: sunday, 23/10/2011                           **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

subroutine conV_dense_jac(neq,tempo,yin,jac)

   use working_precision
   use speedchem_conV, only: SC_conV

   implicit none

   integer, intent(in) :: neq
   real (dp)       ,                 intent(inout)   :: tempo
   real (dp)       , dimension(neq), intent(inout)   :: yin
   real (dp)       , dimension(neq,neq), intent(out) :: jac

   real (dp)        :: ysafe, delta
   real (dp)       , dimension(neq) :: y, dydt0, dydt
   integer :: i, j

!     Evaluate derivative at point
   call SC_conV(neq,tempo,yin,dydt0)

   y = yin

!     Evaluate finite differences by small perturbations
   do i=1, neq

      ysafe = y(i)

      if (i==1) then
         delta = 2.e0_dp
      else
         delta = max(1.d-8, ysafe * 1.d-3)
      endif

      y(i)  = y(i) + delta

      call SC_conV(neq,tempo,y,dydt)

      jac(:,i) = (dydt - dydt0)/delta

      y(i) = ysafe

   end do

end subroutine conV_dense_jac


!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **            Evaluate jacobian sparsity pattern               **
!     **                                                             **
!     **                                                             **
!     **   Author:      (C) Federico Perini                          **
!     **   Date created: sunday, 23/10/2011                          **
!     **   Last update : friday, 17/08/2012                          **
!     **                                                             **
!     **   Usage, copying, distribution of the proprietary routines  **
!     **   or any modified versions of them are not allowed without  **
!     **   approval of the author.                                   **
!     **                                                             **
!     *****************************************************************

subroutine jac_sparsity(neq)

   use working_precision
   use speedchem, only: sparse_jac, njac
   use chemistry_setup,  only: print_out_jacobian
   use sparse_chemistry, only: jac_sparse
   use sparse_algebra,   only: dense_to_sparse,&
   &print_sparsity_to_file

   implicit none

   integer, intent(in) :: neq

   integer   :: i, j, n
   real (dp) :: percent

   logical :: verbose    = .true.
   character(len=*), parameter :: fmt_jac = "(10000(1X,I1))",&
   &fmt_pattern = "(' chemistry jacobian sparsity: ',F4.1,'%')"


!     Number of nonzero elements
   njac = jac_sparse%n

   if (print_out_jacobian) then

      call print_sparsity_to_file(jac_sparse,'dat.jacobian')

   endif ! print_file

!     Print to screen
   percent = 100.e0_dp * (1.e0_dp - real(njac,dp)/real(neq**2,dp) )
   if (verbose) write(*,fmt_pattern) percent

!     Done
   sparse_jac = .true.

end subroutine jac_sparsity



!     *****************************************************************
!     **                                                             **
!     **                     SpeedCHEM FORTRAN                       **
!     **                                                             **
!     **              DUMMY SUBROUTINE for ODE solvers               **
!     **                                                             **
!     **                                                             **
!     **   Author:      Federico Perini                              **
!     **   Last update: wedesday, 26/05/2010                         **
!     **                                                             **
!     *****************************************************************

subroutine dummy
end subroutine dummy

































!
!c     *****************************************************************
!c     **                                                             **
!c     **                     SpeedCHEM FORTRAN                       **
!c     **                                                             **
!c     **    Create ek table from thermodynamic properties from       **
!c     **        chemistry link instead than from datahk file         **
!c     **                                                             **
!c     **                                                             **
!c     **   Author:      Federico Perini                              **
!c     **   Last update: monday,  10/10/2011                          **
!c     **                                                             **
!c     *****************************************************************
!
!
!      subroutine create_ek_table
!
!      use working_precision
!      use parametermod,   only: nspl, nsp
!      use propertiesdata, only: ek
!      use particledata,   only: pcriti,tcriti,acentric
!      use speedchem,      only: ns,specie
!      use SCthermodata,   only: aL, bL, cL, dL, eL, fL,
!     &                          aH, bH, cH, dH, eH, fH,   tsw
!      use universal_constants, only: u2, u3, u4, u5, R, kcal_to_joule
!
!      implicit none
!      logical, parameter :: verbose = .false.
!      integer          :: i, j, nsteps = 51
!      real (dp)        :: T, T2, T3, T4, unosuT, lnT
!      real (dp)       , dimension(ns) :: H,H0
!
!
!
!c     Gathering standard enthalpy at 298.15K, H0
!      T      = 298.15e0_dp
!      T2     = T  * T
!      T3     = T2 * T
!      T4     = T3 * T
!      unosuT = 1/T
!      lnT    = log(T)
!
!c     H0 is standard enthalpy in (J/mol)
!      H0 = aL + u2*bL*T + cL*T2*u3 + dL*T3*u4 + u5*eL*T4 + fL*unosuT
!      H0 = H0 * R * T
!
!
!      do i=1,nsteps
!
!         T = dble(i-1) * 100.e0_dp + 1.e-14_dp
!         T2     = T  * T
!         T3     = T2 * T
!         T4     = T3 * T
!         unosuT = 1/T
!         lnT    = log(T)
!
!         where (tsw.gt.T)
!           H = aL + u2*bL*T + cL*T2*u3 + dL*T3*u4 + u5*eL*T4 +
!     &           fL*unosuT
!         else where
!           H = aH + u2*bH*T + cH*T2*u3 + dH*T3*u4 + u5*eH*T4 +
!     &            fH*unosuT
!         end where
!         H = H * R * T
!
!c        Enthalpy tables have already been filled for fuel species in
!c        rinput. Thus, only indices > nspl are computed here.
!
!c        NB: why KIVA computes calorie-to-joule factor as 4.184?!
!         ek(i,nspl+1:nsp) = ( H(nspl+1:nsp) - H0(nspl+1:nsp) )
!     &                    * kcal_to_joule
!
!      end do
!
!c     Optionally, write ek table to external datahk file
!      if (verbose) then
!         open(unit=38,file='SCdatahk.out')
!         do i=1,ns
!            write(38,*  ) adjustl(specie(i))
!            write(38,920) (ek(j,i),j=1,nsteps)
!            write(38,930) tcriti(i),pcriti(i),acentric(i)
!         end do
!         close(38)
!
!         write(*,910)
!
!      endif
!
!
! 910  format(' ek table built.')
! 920  format(1000(e14.6,','))
! 930  format(1x,1000(e14.6,' '))
!
!      end subroutine create_ek_table
