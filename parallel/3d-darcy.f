c ulimit -s unlimited

c ./3ddarcy

c your code still has a segmentation fault, then try these ifort options:

c ifort -g -traceback -check all -heap-arrays

c The ulimit is a command to set your bash shell stack limit to unlimited.

c I just have one last test for you.  Use these ifort compiler options:

c ifort -g -traceback -check all -fp-stack-check

c gfortran ../../3d-darcy.f 3d-darcy-elmt.f -o 3ddarcy -static -march=i386 

c------------------------------------------------------------------------

c This code uses extra arrays to solve the problem and it gives maximum speedup
c Variable information
c isize = maximum size for nodes, supported by this program
c nproblem = problem number to be solved by this program
c numnpmm = total number of nodes in particular problem, read from geom file
c nsdmm = total number of spatial dimensions, is not currently used in the code.
c numgp = Number of generation point (0 - No generation, > 0 - generate nodal points)
c numelmm = total number of elements in particular problem, read from geom file
c nenmm = total number nodes in one element, is not currently used in the code.
c itrmm = number of iteration to converge solution, need to confirm.
c epsmm = 1 why? not cleared, epsilon (user defined penality parameter) 
c alpha = tolerance parameter
c x = stores the nodal values
c xn = contains the nodal values of previous mesh, it may be removed
c ienmm = stores the information of each element that contains nodes
c xinit = initial value and it may be removed
c ibcx = contains the information of boundary condition in x-axis
c ibcy = contains the information of boundary condition in y-axis
c lstep = load step number or current iteration number
c a = stiffness matrix
c c = damping matrix (SA)

c------------------------------------------------------------------------

      program DRIVER
c
c------------------------------------------------------------------------
c
c Arif Masud, Summer 1991.
c Modified for 3d (SA)
c------------------------------------------------------------------------
c
        implicit none

        integer, parameter :: isize=1200000

	real*8 xinit(3,isize),x(3,isize),delt
        real*8 sDarcyx(4,4,isize),sDarcyy(4,4,isize),
     &         sDarcyz(4,4,isize),sDarcyn(4,4,isize)
        real*8 pDarcyx(isize),pDarcyy(isize),
     &         pDarcyz(isize),pDarcyn(isize)
        real*8 uxn(4,isize),uyn(4,isize),uzn(4,isize),
     &         unn(4,isize)
        real*8 starttime,endtime,omp_get_wtime
        real*8 time1,time2,time3,time
	integer nsdmm,nenmm,numnpmm,numelmm
	integer ienmm(4,isize),iBCmm(isize),iBCxmm(isize),
     &          iBCymm(isize),iBCzmm(isize),omp_get_max_threads
        integer i,j,maxtstep,tstep,reorder,trans,maxtr

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

c
c.... initialization for transient analysis (SA)
c
        reorder = 0
	trans = 0
        time1 = 0.d0
        time2 = 0.d0
        time3 = 0.d0

	if (reorder.eq.1) then
	  print*,'mesh reordering used'
	else
	  print*,'mesh reordering not used'
	end if

	if (trans.eq.1) then
          maxtstep = 5
          delt = 1.d0
	  print*,'running transient algorithm'
	else
          maxtstep = 1
	  print*,'running steady-state algorithm'
	end if

        maxtr = omp_get_max_threads()

	if (trans.eq.1) then
!$omp parallel do num_threads(maxtr) private (i, j)
!$omp& shared (numelmm, nenmm, uxn, uyn, uzn, unn)
          do i = 1, numelmm
            do j = 1, nenmm
              uxn(j,i) = 0.d0
              uyn(j,i) = 0.d0
              uzn(j,i) = 0.d0
              unn(j,i) = 0.d0
            end do
          end do
	end if
c
c.... reading mesh data
c
        call inputMMove(xinit, ienmm, iBCmm, iBCxmm, iBCymm, iBCzmm)

      do i = 1, 10 !for average time
        print*, 'iteration = ', i
c
c.... reordering mesh
c
        if (reorder.eq.1) then
          starttime = omp_get_wtime()
          call reorderMesh(xinit,ienmm,iBCmm,iBCxmm,iBCymm,iBCzmm,x)
          endtime = omp_get_wtime()
          time1 = time1 + (endtime - starttime)
cs          print *,'Mesh Reordering Time: ',time1
        end if
c
c.... computing stiffness matrix
c
          starttime = omp_get_wtime()
          call stiffMDarcy(xinit,ienmm,sDarcyx,sDarcyy,sDarcyz,
     &                     sDarcyn,pDarcyx,pDarcyy,pDarcyz,pDarcyn,
     &                     uxn,uyn,uzn,unn,delt)
          endtime = omp_get_wtime()
          time2 = time2 + (endtime - starttime)
cs          print *,'Stiffness Matrix Time: ',time2

        do tstep = 1, maxtstep
c
c.... solving matrix equation
c
          starttime = omp_get_wtime()
          call stableMDarcy(xinit,ienmm,sDarcyx,sDarcyy,sDarcyz,sDarcyn,
     &                      pDarcyx,pDarcyy,pDarcyz,pDarcyn,iBCmm,
     &                      iBCxmm,iBCymm,iBCzmm,uxn,uyn,uzn,unn,tstep)
          endtime = omp_get_wtime()
          time = endtime - starttime
cs          print *,'Equation Solving Time: ',time3
cs          print *,'Total Time: ',time1+time2+time3
        end do
        time3 = time3 + time

      end do !for average time

      print *,'Mesh Reordering Time: ',time1/10.d0
      print *,'Stiffness Matrix Time: ',time2/10.d0
      print *,'Equation Solving Time: ',time3/10.d0
      print *,'Total Time: ',(time1+time2+time3)/10.d0

	end program DRIVER
c
c-----------------------------end DRIVER--------------------------------
c
c
c-----------------------------inputMMove-------------------------------
c
        subroutine inputMMove(x, ien, iBC, iBCx, iBCy, iBCz)

c       This routine reads coordinates and boundary conditions.
c       SA, 12,2014.
c       Modified for 3D by SA on 6, 2015.
c---------------------------------------------------------------------
c
        implicit none

        real*8 x(3,1)
        integer ien(4,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
        integer nsdmm, nenmm, numelmm, numnpmm
        integer n, i, maxtr, omp_get_max_threads

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

	maxtr = omp_get_max_threads()

        open(unit=1,file='geom_Darcy.dat')
        read(1,*),numnpmm,nsdmm

        do i = 1,numnpmm
          read(1,*),n,x(1,i),x(2,i),x(3,i)
        end do

        read(1,*),numelmm,nenmm

        do i = 1,numelmm
          read(1,*),n,ien(1,i),ien(2,i),ien(3,i),ien(4,i)
        end do
        close(1)

!$omp parallel do num_threads(maxtr) private(i)
!$omp& shared (numnpmm, iBC, iBCx, iBCy, iBCz)
        do i = 1,numnpmm
          iBC(i) = 0
          iBCx(i) = 0
          iBCy(i) = 0
          iBCz(i) = 0
        end do

        open(unit=1,file='BC_Darcy.dat')
        do
          read(1,*),n,iBC(i),iBCx(i),iBCy(i),iBCz(i)
          if (n.eq.0) exit
        end do
        close(1)

        return
        end
c
c---------------------------end inputMMove-----------------------------
c
c
c----------------------------stiffMDarcy-------------------------------
c
c       Input: x and ien
c       Output: a (contains the values of stiffness matrix)

	subroutine stiffMDarcy (x, ien, ax, ay, az, an, px, py, pz, pn,
     &                          uxn, uyn, uzn, unn, delt)
c
c------------------------------------------------------------------------
c
c  This routine calculates the RHS and LHS matrices.
c
c------------------------------------------------------------------------
c
	implicit none

	real*8 x(3,1), ax(4,4,1), ay(4,4,1), az(4,4,1), an(4,4,1), 
     &         px(1), py(1), pz(1), pn(1)
	real*8 xl(3,4), s(16,16), p(16)
        real*8 uxn(4,1), uyn(4,1), uzn(4,1), unn(4,1)
        real*8 c(16,16), delt, st !size 16 is more efficient???
	integer nsdmm, nenmm, numnpmm, numelmm
	integer i, j, n, ien(4,1)
	integer ndf, nst, trans, maxtr, omp_get_max_threads

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm
c
c.... initialization
c
        trans = 0
	ndf = 4
	nst = ndf * nenmm
        st = 0.7d0

	maxtr = omp_get_max_threads()

!$omp parallel do num_threads(maxtr) private (n)
!$omp& shared (numnpmm, px, py, pz, pn)

	do n = 1, numnpmm
	  px (n) = 0.d0
	  py (n) = 0.d0
	  pz (n) = 0.d0
	  pn (n) = 0.d0
	end do
c
c.... loop through elements
c
!$omp parallel num_threads(maxtr) private (n, i, j, xl, p, s)
!$omp& shared (numelmm, nenmm, nsdmm, nst, px, py, pz, pn)
!$omp& shared (x, ien, ax, ay, az, an, trans)

!$omp do
	do n = 1, numelmm

	  do j = 1, nenmm
            do i = 1, nsdmm
	      xl(i,j) = x(i,ien(j,n))
            end do
	  end do

	  do j = 1, nst
	    p(j) = 0.d0
	    do i = 1, nst
	      s(i,j) = 0.d0
  	    end do
	  end do 
c
c....  intializing damping matrix
c
          do i = 1, nenmm
            do j = 1, nenmm
              c(i,j) = 0.d0
            end do
          end do
c
c....  computing L.H.S. matrix and R.H.S vector
c
          call elmt30(xl, s, p, c)

cs          do i = 1, nenmm
cs            do j = 1, nenmm
cs              c(i,j) = (st*c(i,j))/delt
cs            end do
cs          end do

	  do j = 1, nenmm
            do i = 1, nenmm
              if (trans.eq.1) then
c
c.... residue r0 = Ax0 - b
c
	        px(ien(j,n)) = px(ien(j,n)) + p(4*j-3) -
     &                          st*c(i,j)*uxn(j,n)/delt
	        py(ien(j,n)) = py(ien(j,n)) + p(4*j-2) -
     &                          st*c(i,j)*uyn(j,n)/delt
	        pz(ien(j,n)) = pz(ien(j,n)) + p(4*j-1) -
     &                          st*c(i,j)*uzn(j,n)/delt
	        pn(ien(j,n)) = pn(ien(j,n)) + p(4*j) -
     &                          st*c(i,j)*unn(j,n)/delt

c
c.... stiffness matrix
c
                ax(i,j,n) = s(4*i-3,4*j-3) + st*c(i,j)/delt
                ay(i,j,n) = s(4*i-2,4*j-2) + st*c(i,j)/delt
                az(i,j,n) = s(4*i-1,4*j-1) + st*c(i,j)/delt
                an(i,j,n) = s(4*i,4*j) + st*c(i,j)/delt
              else
c
c.... residue r0 = Ax0 - b
c
                px(ien(j,n))=px(ien(j,n))+p(4*j-3)
                py(ien(j,n))=py(ien(j,n))+p(4*j-2)
                pz(ien(j,n))=pz(ien(j,n))+p(4*j-1)
                pn(ien(j,n))=pn(ien(j,n))+p(4*j)
c
c.... stiffness matrix
c
                ax(i,j,n) = s(4*i-3,4*j-3)
                ay(i,j,n) = s(4*i-2,4*j-2)
                az(i,j,n) = s(4*i-1,4*j-1)
                an(i,j,n) = s(4*i,4*j)
              end if

	    end do !loop i
	  end do !loop j
        end do !loop n

!$omp end do
!$omp barrier
!$omp end parallel

        return
        end

c----------------------end stiffMDarcy----------------------------


c------------------------stableMDarcy-----------------------------

c  Subroutine: A subroutine is referenced by a CALL statement.

 	subroutine stableMDarcy(x,ien,ax,ay,az,an,px,py,pz,pn,
     &                          iBC,iBCx,iBCy,iBCz,uxn,uyn,uzn,unn,t)
c
c-----------------------------------------------------------------
c
	implicit none

        integer, parameter :: isize=1200000

	integer ien(4,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
	integer ndf, nst, i, j, n, t, trans
	integer nsdmm, nenmm, numnpmm, numelmm
	real*8 x(3,1), ax(4,4,1), ay(4,4,1), az(4,4,1), an(4,4,1) 
        real*8 px(1), py(1), pz(1), pn(1)
        real*8 ux(isize), uy(isize), uz(isize), un(isize)
	real*8 ul(4,4), s(16,16), xl(3,4)
        real*8 unn(4,1),uxn(4,1),uyn(4,1),uzn(4,1)
        character (len=10) :: fname

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

	ndf = 4
	nst = ndf * nenmm
	trans = 0

        call solveMMoveCG(ien, px, ax, ux, iBCx)
        call solveMMoveCG(ien, py, ay, uy, iBCy)
        call solveMMoveCG(ien, pz, az, uz, iBCz)
        call solveMMoveCG(ien, pn, an, un, iBC)
c
c.... D E B U G G I N G
c        
cs        print *,'ux = ',ux(4)
cs        print *,'uy = ',uy(4)
cs        print *,'uz = ',uz(4)
cs        print *,'un = ',un(4)
c
c.... writing results in files
c
c        open(unit=1,file='ux1.txt');
c        do n = 1, numnpmm
c          write(1,*),x(1,n),x(2,n),x(3,n),ux(n)
c        end do
c        close(1)
c
c        open(unit=1,file='uy1.txt');
c        do n = 1, numnpmm
c          write(1,*),x(1,n),x(2,n),x(3,n),uy(n)
c        end do
c        close(1)
c
c        open(unit=1,file='uz1.txt');
c        do n = 1, numnpmm
c          write(1,*),x(1,n),x(2,n),x(3,n),uz(n)
c        end do
c        close(1)
c
c        open(unit=1,file='un1.txt');
c        do n = 1, numnpmm
c          write(1,*),x(1,n),x(2,n),x(3,n),un(n)
c        end do
c        close(1)
c
c....  saving ux and un for next time step (SA)
c
	if (trans.eq.1) then
!$omp parallel do private (n, j)
!$omp& shared (numelmm, nenmm, uxn, uyn, uzn, unn)
!$omp& shared (ux, uy, uz, un)

          do n = 1, numelmm
            do j = 1, nenmm
              uxn(j,n) = ux(ien(j,n))
              uyn(j,n) = uy(ien(j,n))
              uzn(j,n) = uz(ien(j,n))
              unn(j,n) = un(ien(j,n))
            end do
          end do
	end if
c
c....   need to copy ul, xl, s
c
c        do n = 1, numelmm
c
c	  do j = 1, nenmm
c	    ul(1,j) = ux(ien(j,n))
c	    ul(2,j) = uy(ien(j,n))
c	    ul(3,j) = uz(ien(j,n))
c	    ul(4,j) = un(ien(j,n))
c
c	    xl(1,j) = x(1,ien(j,n))
c	    xl(2,j) = x(2,ien(j,n))
c	    xl(3,j) = x(3,ien(j,n))
c
c	    do i = 1, nenmm
c              s(4*i-3,4*j-3) = ax(i,j,n) 
c              s(4*i-2,4*j-2) = ay(i,j,n) 
c              s(4*i-1,4*j-1) = az(i,j,n) 
c              s(4*i,4*j) = an(i,j,n)
c	    end do
c
c	  end do	
c          call elmtcheck(ul,xl,s)
c
c	end do

	return
	end

c----------------------end stableMDarcy--------------------------


c------------------------solveMMoveCG----------------------------

c Input: ien, a, iBC, r (also modified in this function)
c output: x

	subroutine solveMMoveCG (ien, r, a, x, iBC)

c----------------------------------------------------------------
c
c  This subroutine uses Conjugate Gradient method
c  to solve the system of algebraic quations.
c
c---------------------------------------------------------------
c
	implicit none

        integer, parameter :: isize=1200000
        integer, parameter :: maxp=32

	integer ien(4,1), iBC(1)
	integer nsdmm, nenmm, numnpmm, numelmm, itrmm
        integer tn, tid, omp_get_thread_num, omp_get_max_threads
	integer iter, n, i, j
        integer sharedn(maxp,isize), counter(isize), countn
	real*8 pap, rbr, alpha, beta, rbr0, tmp
	real*8 r(1), a(4,4,1), x(1) !, bin(isize,maxp)
	real*8 p(isize), q(isize), B(isize), tol1

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

	data tol1 /1.d-2/
cs        itrmm = 200
        itrmm = 1

c
c.... initialize the variables
c
cs      tn = omp_get_max_threads()

!$omp parallel do private (n, i)
!$omp&  shared (numnpmm, x, p)
cs!$omp&  shared (numnpmm, B, x, p, bin, tn)

      do n = 1, numnpmm
cs        B(n) = 0.d0
	x(n) = 0.d0
	p(n) = 0.d0
cs        do i = 1, tn
cs          bin(n,i) = 0.d0
cs        end do
      end do
c
c.... computing B matrix for Dirichlet BCs
c
cs!$omp parallel private (n, i, tid)
cs!$omp& shared (numelmm, numnpmm, nenmm, a)
cs!$omp& shared (numelmm, numnpmm, nenmm, B, a, bin, tn)

      tid = omp_get_thread_num() + 1

cs!$omp do
cs      do n = 1, numelmm
cs        do i = 1, nenmm
cs          bin(ien(i,n),tid)=bin(ien(i,n),tid)+a(i,i,n)
cs        enddo
cs      end do
cs!$omp end do

cs!$omp do
cs      do n = 1, numnpmm
cs        do i = 1, tn
cs          B(n) = B(n) + bin(n,i)
cs        enddo
cs      enddo
cs!$omp end do
cs!$omp end parallel

c
c.... computing rbr
c
      rbr = 0.

!$omp parallel do private (n)
!$omp& shared (numnpmm, r, p)
!$omp& reduction (+:rbr)
      do n = 1, numnpmm
        p(n) = -r(n) !initializing p
        rbr  = rbr + r(n) * r(n)
      end do

      rbr0 = tol1 * rbr
c
c.... iterate for convergence
c
      do 5000 iter = 1, itrmm
c
c.... perform the pAp product
c
!$omp parallel private (n, i, j, tid)
!$omp& shared (numnpmm, numelmm, nenmm, ien, p, a, sharedn)

!$omp do
        do n = 1, numnpmm
          q(n) = 0.d0
        end do
!$omp end do

!$omp do
        do n = 1, numelmm
          do i = 1, nenmm
            do j = 1, nenmm
              q(ien(i,n)) = q(ien(i,n)) + a(i,j,n)*p(ien(j,n))
            end do
          end do
        end do
!$omp end do

c*********************************************
c.... counting shared nodes by SA
c
c        tid = omp_get_thread_num() + 1
c
c        do i = 1, tn
c          do j = 1, numnpmm
c            sharedn(i,j) = 0
c          end do
c        end do
c
c!$omp do 
c        do n = 1, numelmm
cc          print *, tid, n
c          do i = 1, nenmm
c            sharedn(tid,ien(i,n)) = 1
cc            print *,tid,ien(i,n),sharedn(tid,ien(i,n))
c          end do
c        end do
c!$omp end do
c
cc        if (tid.eq.1) then
cc          do i = 1, tn
cc            do j = 1, numnpmm
cc              print *,i,j,sharedn(i,j)
cc            end do
cc          end do
cc        end if

!$omp end parallel

c        do i = 1, numnpmm
c          counter(i) = 0
c        end do
c
c        do i = 1, tn
c          do j = 1, numnpmm
cc            print *,sharedn(i,j)
c            if (sharedn(i,j).gt.0) then 
c              counter(j) = counter(j) + 1
c            end if
c          end do
c        end do
c
c        countn = 0
c        do i = 1, numnpmm
c          if (counter(i).gt.1) then
c            countn = countn + 1
c          end if
c        end do
c        
c        print *, 'Shared Nodes = ', countn
cc
cc**********************************************
cc
c
c.... calculate alpha
c
        pap = 0.d0

!$omp parallel do  private (n)
!$omp& shared (numnpmm, q, p)
!$omp& reduction (+:pap)

	do n = 1, numnpmm
	  pap = pap + p(n) * q(n)
        end do

	alpha = rbr / pap 
c
c.... compute solution
c
        tmp = rbr

	rbr = 0.d0

!$omp parallel do private(n)
!$omp& shared(numnpmm, x, alpha, p, r, q)
!$omp& reduction (+:rbr)
cs!$omp& shared(numnpmm, x, alpha, p, r, q, B)

	do n = 1, numnpmm
          if (iBC(n) .gt. 0) then 
cs	    x(n) = r(n)/B(n) !seems dirichlet BCs only (SA)
            x(n) = iBC(n)
          else
	    x(n) = x(n) + alpha * p(n)
          end if
	  r(n) = r(n) + alpha * q(n)
	  rbr = rbr + r(n) * r(n)
        end do

        beta = rbr / tmp
c
c.... check for convergence
c
	if (rbr .le. rbr0) goto 6000
c
c.... calculate a new search direction
c
!$omp parallel do private (n)
!$omp& shared (numnpmm, p, r, beta)

        do n = 1, numnpmm
          p(n) = -r(n) + beta * p(n)
        end do

c   here p is the q of the algorithm
c.... end of iteration
c
5000	continue
c
c.... if converged
c
6000	continue

	return
	end
c
c--------------------------end solveMMoveCG--------------------------
c
c
c-----------------------------reorderMesh----------------------------
c
        subroutine reorderMesh(x, ien, iBC, iBCx, iBCy, iBCz, xtemp)
c
c--------------------------------------------------------------------
c
c  This subroutine reorder the mesh using tetra tree
c  and it adjust all the nodes and their boundary conditions
c
c
c Masroor Hussain, 2008
c--------------------------------------------------------------------
c
        implicit none
        integer, parameter:: isize=1200000
        integer, parameter:: maxp=32

        real*8 x(3,1), xtemp(3,1), tmin(3), tmax(3), tx(3,isize)
c        real*8 sp(3,128) ! splitter
        real*8 sp(3,maxp*maxp) ! splitter
        integer ien(4,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
        integer eleid(isize), h, l, r, i, j
c        integer neleid(isize), sid(64), part(64), cid(8,8)
c        integer gtsp(8,8), gtpart(8,8), gsize(8)
        integer neleid(isize), sid(maxp*maxp),
     &          part(maxp*maxp), cid(maxp,maxp)
        integer gtsp(maxp,maxp), gtpart(maxp,maxp), gsize(maxp)
        integer chunksize, osize, n, lr, rr, k
        integer nsdmm, nenmm, numelmm, numnpmm
        integer tn, tid, omp_get_thread_num, omp_get_max_threads

        common / infmMM / nsdmm, nenmm, numelmm, numnpmm

        call makePartitioningArrays(x, ien, xtemp, eleid)
        tn = omp_get_max_threads()
cs        chunksize = 64
        chunksize = 4096 !for one million elements

        if (tn .gt. 1) then

!$omp parallel do private (i) shared (sp)
c          do i = 1, 128
          do i = 1, maxp*maxp
            sp(1,i) = 0.d0
            sp(2,i) = 0.d0
            sp(3,i) = 0.d0
          end do

!$omp parallel do private (i, j) 
!$omp& shared (cid)
          do i = 1, maxp
            gsize(i) = 0
            do j = 1, maxp
              cid(i,j) = 0
            end do
          end do
          osize = numelmm / tn

!$omp parallel private (tmin, tmax, tid, h, l, r, i, j)
!$omp& shared (osize,numelmm,nsdmm,xtemp,eleid,sp,tn,chunksize,sid)
          tid = omp_get_thread_num ()
          l = tid * osize + 1
          r = l + osize - 1
          if(tid .eq. (tn - 1)) r = numelmm
          h = 0
          call calcMinMax(xtemp, l, r, tmin, tmax)

          call recursiveCreateTree (xtemp, eleid,
     &             l, r, tmin, tmax, h, chunksize)

          tid = omp_get_thread_num ()
          l = tid * osize + 1
          r = l + osize - 1
          do i = 1, 8
            j = l+i*osize/8 - 1
            if (i .eq. 8 .and. (tid .eq. (tn - 1))) j = numelmm
            sp(1,tid*8+i) = xtemp(1,j)
            sp(2,tid*8+i) = xtemp(2,j)
            sp(3,tid*8+i) = xtemp(3,j)
            sid(tid*8+i) = j
          end do
!$omp end parallel

          l = 1
          r = 8*tn
          h = 0
cs          chunksize = 64
          chunksize = 4096

          call calcMinMax(sp, l, r, tmin, tmax)

          call recursiveCreateTree(sp, sid,
     &             l, r, tmin, tmax, h, chunksize)

!$omp parallel private (i, j, tid)
!$omp& shared (part, sid, tn, osize, cid)
!$omp do
          do i = 1, 8*tn
            j = (sid(i)-1)/osize
            if (j .ge. tn) j = tn - 1
            part(i) = j + 1
          end do
!$omp end do

          tid = omp_get_thread_num()
          i = tid*8 + 1
          j = 8

          call countingSort(part(i), sid(i), j)

          i = tid + 1

          do j = 1, 8
            if (sid(j+tid*8) .gt. cid(i, part(j+tid*8))) then
              cid(i, part(j+tid*8)) = sid(j+tid*8)
            end if
          end do
!$omp end parallel
c
c.... copy data
c
!$omp parallel private (i, j)
!$omp& shared (gtsp, cid, gtpart, tn)
!$omp do
          do i = 1, tn
            do j = 1, tn
              gtsp(j,i) = cid(j, i)
              gtpart(j,i) = j
            end do
          end do
!$omp end do

!$omp do
          do i = 1, tn
            call insertionSort(gtsp(1,i),gtpart(1,i), tn)
          end do
!$omp end do
!$omp end parallel
 
          k = 1
          n = 1
          do i = 1, tn
            l = (i-1) * osize
            do j = 1, tn
              if (j .eq. 1) then
                k = gtsp(j,i) - l
              else
                if (gtsp(j-1,i) .eq. 0) then
                  k = gtsp(j,i) - l
                else
                  k = gtsp(j,i) - gtsp(j-1,i)
                end if
              end if
              if (k .gt. 0) then
                gsize(gtpart(j,i)) = gsize(gtpart(j,i)) + k
              end if
            end do
          end do

          do i = 2, tn
            gsize(i) = gsize(i) + gsize(i-1)
          end do

!$omp parallel private (i, j, k, l, r, h, lr, rr, tid, tmin, tmax)
!$omp& shared (gsize, xtemp, tx, eleid, neleid, gtsp, tn, chunksize)
!$omp& shared (gtpart)
          tid = omp_get_thread_num ()
          if (tid .eq. 0) then
            l = 1
          else
            l = gsize(tid) + 1
          endif
          r = gsize(tid+1)
          do i = 1, tn
            lr = (i-1) * osize + 1
            do j = 1, tn
              if (gtsp(j,i) .eq. 0) then
                rr = 0
              else
                rr = gtsp(j,i)
                if (j .gt. 1 .and. gtsp(j-1,i) .gt. 0) then
                  lr = gtsp(j-1,i) + 1
                end if
              end if
              if ((gtpart(j,i) .eq. (tid+1)) .and. (rr .gt. 0))then
                do k = lr, rr
                  tx(1,l) = xtemp(1,k)
                  tx(2,l) = xtemp(2,k)
                  tx(3,l) = xtemp(3,k)
                  neleid(l) = eleid(k)
                  l = l + 1
                end do
              end if
            end do
          end do
!$omp end parallel

          if (tid .eq. 0) then
            l = 1
          else
            l = gsize(tid) + 1
          endif
          r = gsize(tid+1)
          h = 0
c          chunksize = 1
cs          chunksize = 64
          chunksize = 4096

          call calcMinMax(tx, l, r, tmin, tmax)
          call recursiveCreateTree (tx, neleid,
     &             l, r, tmin, tmax, h, chunksize)


          call rearrangeData(ien, x, iBC, iBCx, iBCy,
     &                          iBCz, neleid, tx)

        else ! if (tn .eq. 1) then

          h = 0
          l = 1
          r = numelmm
cs          chunksize = 64
          chunksize = 4096

          call calcMinMax(xtemp, l, r, tmin, tmax)

c        print *,'recursive 4'
          call recursiveCreateTree(xtemp, eleid,
     &             l, r, tmin, tmax, h, chunksize)
          call rearrangeData(ien, x, iBC, iBCx, iBCy,
     &                          iBCz, eleid, xtemp)

        end if

        return
        end
c
c---------------------end of reorderMesh subroutine ------------------
c
c
c----------------------------insertionSort-----------------------------
c
        subroutine insertionSort(idata, sdata, dsize)

        integer idata(1), sdata(1), dsize, i, j, ikey, skey

        do j = 2, dsize
          ikey = idata(j)
          skey = sdata(j)
          i = j - 1
          do while ((i .gt. 0) .and. (idata(i) .gt. ikey))
            idata(i+1) = idata(i)
            sdata(i+1) = sdata(i)
            i = i - 1
          end do
          idata(i+1) = ikey
          sdata(i+1) = skey
        end do

        end
c
c---------------------------end insertionSort-------------------------
c
c
c--------------------------------calcMinMax--------------------------
c
        subroutine calcMinMax(data, l, r, tmin, tmax)

        implicit none

        real*8 data(2,1), tmin(1), tmax(1)
        integer l, r, i, j
        integer nsdmm, nenmm, numelmm, numnpmm

        common / infmMM / nsdmm, nenmm, numelmm, numnpmm
c
c.... get the maximum and minimum value
c
        do i = l, r
          if (i .eq. l) then
            do j = 1, nsdmm
              tmin(j) = data(j,i)
              tmax(j) = data(j,i)
            end do
          else
            do j = 1, nsdmm
              if (tmin(j) .gt. data(j,i)) then
                tmin(j) = data(j,i)
              end if
              if (tmax(j) .lt. data(j,i)) then
                tmax(j) = data(j,i)
              end if
            end do
          endif
        end do

        end
c
c------------------------------end calcMinMax----------------------
c
c
c---------------------Recursively create tetra tree--------------------
c
c  Input:       eleval, eleid
c  Output:      eleval, eleid, tree
c  this function is only written for 2D have to generalized to 3D later
c  first divide the problem into eight spaces for 3D
c  Abberivation
c  F -> Front, B -> Back, L -> Left, R -> Right, U -> Up and D -> Down
c  Eight spaces are FLU, FRU, BLU, BRU, FLD, FRD, BLD and BRD for 3D

        RECURSIVE subroutine recursiveCreateTree(eleval, eleid,
     &                         l, r, tmin, tmax, h, chunksize)

        implicit none

        real*8 tmin(1), tmax(1), cmin(3), cmax(3)
        real*8 eleval(3,1)
        integer indices(9), chunksize
        integer l, r, size, h !height
        integer i, j, eleid(1)

        size = r - l + 1
          
        if (size .gt. chunksize) then
          call rearrange(eleval, eleid, l, r, tmin, tmax, indices)
          l = indices(1)
          r = indices(2) - 1
          do i = 1, 3
            cmin(i) = tmin(i)
            cmax(i) = (tmin(i) + tmax(i))/2.0
          end do

!         print *, 'calling 1'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices(2)
          r = indices(3) - 1
          do i = 1, 3
            cmin(i) = tmin(i)
            cmax(i) = (tmin(i) + tmax(i))/2.0
          end do
          cmin(3) = cmax(3)
          cmax(3) = tmax(3)

!         print *, 'calling 2'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices(3)
          r = indices(4) - 1
          do i = 1, 3
            cmin(i) = tmin(i)
            cmax(i) = (tmin(i) + tmax(i))/2.0
          end do
          cmin(2) = cmax(2)
          cmax(2) = tmax(2)

!         print *, 'calling 3'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices(4)
          r = indices(5) - 1
          do i = 1, 3
                  cmin(i) = tmin(i)
                  cmax(i) = (tmin(i) + tmax(i))/2.0
          end do
          cmin(2) = cmax(2)
          cmax(2) = tmax(2)
          cmin(3) = cmax(3)
          cmax(3) = tmax(3)

!         print *, 'calling 4'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices(5)
          r = indices(6) - 1
          do i = 1, 3
            cmin(i) = tmin(i)
            cmax(i) = (tmin(i) + tmax(i))/2.0
          end do
          cmin(1) = cmax(1)
          cmax(1) = tmax(1)

!         print *, 'calling 5'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices(6)
          r = indices(7) - 1
          do i = 1, 3
            cmin(i) = tmin(i)
            cmax(i) = (tmin(i) + tmax(i))/2.0
          end do
          cmin(1) = cmax(1)
          cmax(1) = tmax(1)
          cmin(3) = cmax(3)
          cmax(3) = tmax(3)

!         print *, 'calling 6'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices (7)
          r = indices (8) - 1
          do i = 1, 3
            cmin(i) = tmin(i)
            cmax(i) = (tmin(i) + tmax(i))/2.0
          end do
          cmin(1) = cmax(1)
          cmax(1) = tmax(1)
          cmin(2) = cmax(2)
          cmax(2) = tmax(2)

!         print *, 'calling 7'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)
          l = indices(8)
          r = indices(9)

          do i = 1, 3
            cmin(i) = (tmin(i) + tmax(i))/2.0
            cmax(i) = tmax(i)
          end do

!         print *, 'calling 8'
          call recursiveCreateTree(eleval, eleid, l, r, cmin,
     &                             cmax, h, chunksize)

        end if   ! chunk size

        return
      end
c
c-------------------------end of recursive subroutine-------------------
c
c
c--------------------make array for temporary processing---------------
c computing element centroid
c done for 3D
        subroutine makePartitioningArrays(x, ien, xtemp, eleid)

        implicit none

        real*8 x(3,1), xtemp(3,1)
        integer eleid(1), ien(4,1), i, j, k, nn
        integer nsdmm, nenmm, numelmm, numnpmm

        common / infmMM / nsdmm, nenmm, numelmm, numnpmm

!$omp parallel do private (i, j) shared (numelmm, nsdmm, eleid, xtemp)
        do i = 1, numelmm
          eleid(i) = i
          do j = 1, nsdmm
            xtemp(j,i) = 0
          end do
        end do

!$omp parallel private (i, j, k, nn)
!$omp& shared (ien, x, xtemp, numelmm, nenmm, nsdmm)
!$omp do
        do i = 1, numelmm
          do j = 1, nenmm
            nn = ien(j,i)
            do k = 1, nsdmm
              xtemp(k,i) = xtemp(k,i) + x(k,nn)
            end do
          end do
        end do
!$omp end do

!$omp do
        do i = 1, numelmm
          do j = 1, nsdmm
            xtemp(j,i) = xtemp(j,i) / nenmm
          end do
        end do
!$omp end do
!$omp end parallel

        return
      end
c
c------------------end of makePartitioningArrays-----------------------
c
c
c---------------------rearrangeData function---------------------------
c
        subroutine rearrangeData(ien, x, iBC, iBCx, iBCy, iBCz,
     &                                  id, xtemp)

        implicit none
        integer, parameter :: isize=1200000

        real*8 x(3,1), xtemp(3,1)
        integer temp(8,isize), i, j, k, ni, id(1)
        integer ien(4,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
        integer nsdmm, nenmm, numelmm, numnpmm

        common / infmMM / nsdmm, nenmm, numelmm, numnpmm
c
c.... copy the elements
c
c!$omp parallel private (i, j)
c!$omp& shared (numelmm, numnpmm, nenmm, temp, ien, id)
c!$omp do
        do i = 1, numelmm
          do j = 1, nenmm
            temp(j,i) = ien(j,id(i))
          enddo
        enddo
c!$omp end do

c!$omp do
        do i = 1, numelmm
          do j = 1, nenmm
            ien(j,i) = temp(j,i)
          enddo
        enddo
c!$omp end do

c!$omp do
        do i = 1, numnpmm
          id(i) = 0
        end do
c!$omp end do
c
c.... re-numbering nodes
c
        ni = 1
c!$omp do
        do i = 1, numelmm
          do j = 1, nenmm
            if (id(ien(j,i)) .eq. 0) then
              id(ien(j,i)) = ni
              ien(j,i) = ni
              ni = ni + 1
            else
              k = ien(j,i)
              ien(j,i) = id(k)
            end if
          end do
        end do
c!$omp end do
c!$omp end parallel
c
c.... copy nodes
c
c!$omp parallel private (i, j)
c!$omp& shared (iBC, iBCx, iBCy, iBCz, x, xtemp)
c!$omp& shared (temp, id, numnpmm, nsdmm)
c!$omp do
        do i = 1, numnpmm
          do j = 1, nsdmm
            xtemp(j,id(i)) = x(j,i)
          enddo
        enddo
c!$omp end do

c!$omp do
        do i = 1, numnpmm
          do j = 1, nsdmm
            x(j,i) = xtemp(j,i)
          enddo
        enddo
c!$omp end do

c!$omp do
        do i = 1, numnpmm
          temp(1,id(i)) = iBC(i)
          temp(2,id(i)) = iBCx(i)
          temp(3,id(i)) = iBCy(i)
          temp(4,id(i)) = iBCz(i)
        end do
c!$omp end do

c!$omp do
        do i = 1, numnpmm
          iBC(i) = temp(1,i)
          iBCx(i) = temp(2,i)
          iBCy(i) = temp(3,i)
          iBCz(i) = temp(4,i)
        end do
c!$omp end do
c!$omp end parallel

        return
        end
c
c---------------------end of rearrangeData ---------------------------
c
c
c----------------------rearrange subroutine---------------------------
c       done for 3D
        
        subroutine rearrange(eleval, eleid, l, r, tmin, tmax, indices)
        
        implicit none
        
        real*8 eleval(3,1), tmin(1), tmax(1)
        integer l, r, indices(9), mid, left, right
        integer eleid(1), c_axis

        indices(1) = l
        indices(9) = r
        left = l
        right = r
c
c.... partition for +ve x-axis
        c_axis = 0
        call partition(eleval, eleid, left, right, tmin(1),
     &                  tmax(1), c_axis, mid)
        indices(5) = mid
        left = indices(1)
        right = mid - 1
c
c.... partition for -ve x-axis
c
        c_axis = 1
        call partition(eleval, eleid, left, right, tmin(2),
     &                  tmax(2), c_axis, mid)
        indices (3) = mid

        left = indices(5)
        right = r
        call partition(eleval, eleid, left, right, tmin(2),
     &                  tmax(2), c_axis, mid)
        indices(7) = mid
c
c.... partition for z-axis
c
        c_axis = 2
        left = indices(1)
        right = indices(3) - 1
        call partition(eleval, eleid, left, right, tmin(3),
     &                  tmax(3), c_axis, mid)

        indices(2) = mid

        left = indices(3)
        right = indices(5) - 1
        call partition(eleval, eleid, left, right, tmin(3),
     &                  tmax(3), c_axis, mid)

        indices(4) = mid

        left = indices(5)
        right = indices(7) - 1
        call partition(eleval, eleid, left, right, tmin(3),
     &                  tmax(3), c_axis, mid)
        indices(6) = mid

        left = indices(7)
        right = indices(9)
        call partition(eleval, eleid, left, right, tmin(3),
     &                  tmax(3), c_axis, mid)
        indices(8) = mid

      end
c
c-------------------------end of rearrange subroutine----------------------
c
c
c-------------------------partition algorithm---------------------------
c       done for 3D
        subroutine partition(eleval, eleid, l, r, tmin, tmax, b, mid)
 
        implicit none

        real*8 eleval(3,1), tmin, tmax, pivotvalue, v1
        integer l, r, b, i, j
        integer eleid(1), mid, temp

        i = l - 1
        pivotvalue = (tmin + tmax) / 2.0

        do j = l, r
          if (eleval(b+1,j) .lt. pivotvalue) then
            if (i .ne. j) then
              i = i + 1
              v1 = eleval(1,j)
              eleval(1,j) = eleval(1,i)
              eleval(1,i) = v1

              v1 = eleval(2,j)
              eleval(2,j) = eleval(2,i)
              eleval(2,i) = v1

              v1 = eleval(3,j)
              eleval(3,j) = eleval(3,i)
              eleval(3,i) = v1

              temp = eleid(i)
              eleid(i)= eleid(j)
              eleid(j)= temp
            endif
          endif
        end do

        mid = i + 1

        return
        end
c
c--------------------------end of partition---------------------
c
c
c---------------------------metis partition---------------------
c
c
c        subroutine testmetis()
c
c        implicit none
c
c        integer, parameter :: nels=2, nnds=6, npel=4
c        integer eptr(nels+1), nodes(nels*npel), epart(nels),
c     &          npart(nnds), n
c        integer, pointer :: vwgt=>null(), vsize=>null(), mopts=>null()
c        real*8, pointer :: tpwgts=>null()
c
c        eptr=(/0,4,8/)
c        nodes=(/0,1,2,3,1,4,5,2/) ! Element 1 has nodes 0 1 2 3
c                            ! Element 2 has nodes 1 4 5 2
c
c        call METIS_PartMeshNodal(nels,nnds,eptr,nodes,vwgt,vsize,2,
c     &                                  tpwgts,mopts,n,epart,npart) 
c        
c        print*, npart; print*, epart
c        
c        return
c        end 
c
c-------------------------end metis partition--------------------
c
