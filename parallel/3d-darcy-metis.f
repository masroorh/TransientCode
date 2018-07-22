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

        integer, parameter :: isize=30000

	real*8 xinit(3,isize), sDarcyx(4,4,isize), delt,
     &         sDarcyy(4,4,isize), sDarcyz(4,4,isize),
     &         sDarcyn(4,4,isize), pDarcyx(isize), 
     &         pDarcyy(isize), pDarcyz(isize), pDarcyn(isize)
        real*8 uxn(4,isize),uyn(4,isize),uzn(4,isize),
     &         unn(4,isize)
        real*8 starttime,endtime,omp_get_wtime
        real*8 time,time2,time3
	integer ienmm(4,isize), iBCmm(isize), iBCxmm(isize),
     &          iBCymm(isize), iBCzmm(isize)
	integer nsdmm, nenmm, numnpmm, numelmm
        integer i,j,maxtstep,tstep,reorder

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

c
c.... initialization for transient analysis (SA)
c
        delt = 1.d-3
        maxtstep = 5
        time2 = 0.d0
        time3 = 0.d0

!$omp parallel do private (i, j)
!$omp& shared (numelmm, nenmm, uxn, uyn, uzn, unn)
        do i = 1, numelmm
          do j = 1, nenmm
            uxn(j,i) = 0.d0
            uyn(j,i) = 0.d0
            uzn(j,i) = 0.d0
            unn(j,i) = 0.d0
          end do
        end do
c
c.... reading mesh data
c
        call inputMMove(xinit, ienmm, iBCmm, iBCxmm, iBCymm, iBCzmm)

      do i = 1, 10 !for average time
        print*, 'iteration = ', i
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
c
c.... time steps through transient phase
c
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
cs          print *,'Equation Solving Time: ',time2
cs          print *,'Total Time: ',time1+time2
        end do
	time3 = time3 + time
      end do !for average time

      print *,'Stiffness Matrix Time: ',time2/10.d0
      print *,'Equation Solving Time: ',time3/10.d0
      print *,'Total Time: ',(time2+time3)/10.d0

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

        integer, parameter :: isize=30000

        real*8 x(3,1), x_temp(3,isize)
        integer ien(4,1), ien_temp(4,isize)
        integer iBC(1), iBCx(1), iBCy(1), iBCz(1)
        integer nsdmm, nenmm, numelmm, numnpmm
        integer n, i, j, k

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

c        open(unit=1,file='geom_Darcy.dat')
c        read(1,*),numnpmm,nsdmm
c        do i = 1,numnpmm
c          read(1,*),n,x(1,i),x(2,i),x(3,i)
c        end do
c
c        read(1,*),numelmm,nenmm
c        do i = 1, numelmm
c          read(1,*),n,ien(1,i),ien(2,i),ien(3,i),
c     &                  ien(4,i)
c        end do
c        close(1)

        open(unit=1,file='geom_Darcy.dat')
        read(1,*),numnpmm,nsdmm
        do i = 1,numnpmm
          read(1,*),n,x_temp(1,i),x_temp(2,i),x_temp(3,i)
        end do

        read(1,*),numelmm,nenmm
        do i = 1, numelmm
          read(1,*),n,ien_temp(1,i),ien_temp(2,i),ien_temp(3,i),
     &                  ien_temp(4,i)
        end do
        close(1)

        open(unit=1,file='metis_geom_Darcy.dat')
        do i = 1, numnpmm
          read(1,*),n
          do j = 1, nsdmm
            x(j,i) = x_temp(j,n+1)
          end do
          
          do k = 1, numelmm
            do j = 1, nenmm
              if (ien_temp(j,k).eq.i) then
                ien(j,k) = n + 1 !starts with 1
cs		print*,ien(j,k)
              end if
            end do
          end do
        end do
        close(1)

c        open(unit=1,file='x.txt')
c        do n = 1, numnpmm
c          write(1,*),x(1,n),x(2,n),x(3,n)
c        end do
c        close(1)
c
c        open(unit=1,file='ien.txt')
c        do n = 1, numelmm
c          write(1,*),n,ien(1,n),ien(2,n),ien(3,n),
c     &                  ien(4,n)
c        end do
c        close(1)

!$omp parallel do private(i)
!$omp& shared (numnpmm, iBC, iBCx, iBCy, iBCz)
        do i = 1, numnpmm
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
        real*8 c(4,4), delt, st
	integer nsdmm, nenmm, numnpmm, numelmm
	integer i, j, n, ien(4,1)
	integer ndf, nst, trans

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm
c
c.... initialization
c
        trans = 1
	ndf = 4
	nst = ndf * nenmm
        st = 0.7d0

!$omp parallel do private (n)
!$omp& shared (numnpmm, px, py, pz, pn)

	do n = 1, numnpmm
	   px(n) = 0.d0
	   py(n) = 0.d0
	   pz(n) = 0.d0
	   pn(n) = 0.d0
	end do
c
c.... loop through elements
c
!$omp parallel private (n, i, j, xl, p, s)
!$omp& shared (numelmm, nenmm, nsdmm, nst, px, py, pz, pn)
!$omp& shared (ien, x, ax, ay, az, an)

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

        integer, parameter :: isize=30000

	integer ien(4,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
	integer ndf, nst, i, j, n, t
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
c	      s(4*i,4*j) = an(i,j,n)
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

        integer, parameter :: isize=30000
        integer, parameter :: maxp=32

	integer ien(4,1), iBC(1)
	integer nsdmm, nenmm, numnpmm, numelmm, itrmm
        integer tn, tid, omp_get_thread_num, omp_get_max_threads
	integer iter, n, i, j
	real*8 pap, rbr, alpha, beta, rbr0, tmp
	real*8 r(1), a(4,4,1), x(1), bin(isize,maxp)
	real*8 p(isize), q(isize), B(isize), tol1

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

	data tol1 / 1.d0 /
cs        itrmm = 200
        itrmm = 10

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

cs!$omp parallel private (n, i, tid)
cs!$omp& shared (numelmm, numnpmm, nenmm, B, a, bin, tn)
cs
cs      tid = omp_get_thread_num() + 1
cs
cs!$omp do
cs      do n = 1, numelmm
cs        do i = 1, nenmm
cs          bin(ien(i,n),tid)=bin(ien(i,n),tid)+a(i,i,n)
cs        enddo
cs      end do
cs!$omp end do
cs
cs!$omp do
cs      do n = 1, numnpmm
cs        do i = 1, tn
cs          B(n) = B(n) + bin(n,i)
cs        enddo
cs      enddo
cs!$omp end do
cs
cs!$omp end parallel

c
c.... computing rbr
c
      rbr = 0.

!$omp parallel do private (n) reduction (+:rbr)
!$omp&  shared (numnpmm, r, p)

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
!$omp parallel private (n, i, j)
!$omp& shared (numnpmm, nenmm, ien, p, a, q)
cs!$omp parallel private (n, i, j, tid)
cs!$omp& shared (numnpmm, nenmm, ien, p, a, tn, bin, q)

cs        tid = omp_get_thread_num() + 1

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
!$omp end parallel

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
	    x(n) = iBC(n) !seems dirichlet BCs only (SA)
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
c----------------------end solveMMoveCG-----------------------
c
