c ulimit -s unlimited

c ./3ddarcy

c your code still has a segmentation fault, then try these ifort options:

c ifort -g -traceback -check all -heap-arrays

c The ulimit is a command to set your bash shell stack limit to unlimited.

c I just have one last test for you.  Use these ifort compiler options:

c ifort -g -traceback -check all -fp-stack-check

c gfortran ../../3d_darcy_trans.f -o 3ddarcytrans -static -march=i386 

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

        integer, parameter :: isize=200000

	real*8 xinit(3,isize),epsmm,
     &         sDarcyx(8,8,isize),sDarcyy(8,8,isize),sDarcyz(8,8,isize),
     &         sDarcyn(8,8,isize),pDarcyx(isize), 
     &         pDarcyy(isize),pDarcyz(isize),pDarcyn(isize)
	integer nsdmm,nenmm,numnpmm,numelmm
	integer ienmm(8,isize),iBCmm(isize),iBCxmm(isize),
     &          iBCymm(isize),iBCzmm(isize)

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

c
c.... reading mesh data
c
        call inputMMove(xinit, ienmm, iBCmm, iBCxmm, iBCymm, iBCzmm)
c
c.... computing stiffness matrix
c
        call stiffMDarcy(xinit,ienmm,sDarcyx,sDarcyy,sDarcyz,sDarcyn,
     &                         pDarcyx,pDarcyy,pDarcyz,pDarcyn)
c
c.... solving matrix equation
c
	call stableMDarcy(xinit,ienmm,sDarcyx,sDarcyy,sDarcyz,sDarcyn,
     &                          pDarcyx,pDarcyy,pDarcyz,pDarcyn,iBCmm,
     &                          iBCxmm,iBCymm,iBCzmm)

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
        integer ien(8,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
        integer nsdmm, nenmm, numelmm, numnpmm
        integer n, i

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

        open(unit=1,file='geom_Darcy.dat')
        read(1,*),numnpmm,nsdmm
        do i = 1,numnpmm
          read(1,*),n,x(1,i),x(2,i),x(3,i)
        end do

        read(1,*),numelmm,nenmm
        do i = 1,numelmm
          read(1,*),n,ien(1,i),ien(2,i),ien(3,i),ien(4,i),
     &                  ien(5,i),ien(6,i),ien(7,i),ien(8,i)
        end do
        close(1)

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

	subroutine stiffMDarcy (x, ien, ax, ay, az, an, px, py, pz, pn)
c
c------------------------------------------------------------------------
c
c  This routine calculates the RHS and LHS matrices.
c
c------------------------------------------------------------------------
c
	implicit none

	real*8 x(3,1), ax(8,8,1), ay(8,8,1), az(8,8,1), an(8,8,1), 
     &         px(1), py(1), pz(1), pn(1)
	real*8 xl(3,8), s(32,32), p(32)
	integer nsdmm, nenmm, numnpmm, numelmm
	integer i, j, n, ien(8,1)
	integer ndf, nst

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm
c
c.... initialization
c
	ndf = 4
	nst = ndf * nenmm

	do n = 1, numnpmm
	   px (n) = 0.d0
	   py (n) = 0.d0
	   pz (n) = 0.d0
	   pn (n) = 0.d0
	end do
c
c.... loop through elements
c
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
c....  computing L.H.S. matrix and R.H.S vector
c
          call elmt30(xl,s,p)

	  do j = 1, nenmm
c
c.... residue r0 = Ax0 - b
c
	    px(ien(j,n)) = px(ien(j,n)) + p(4*j-3)
	    py(ien(j,n)) = py(ien(j,n)) + p(4*j-2)
	    pz(ien(j,n)) = pz(ien(j,n)) + p(4*j-1)
	    pn(ien(j,n)) = pn(ien(j,n)) + p(4*j)
c
c.... stiffness matrix
c
	    do i = 1, nenmm
              ax(i,j,n) = s(4*i-3,4*j-3)
              ay(i,j,n) = s(4*i-2,4*j-2)
              az(i,j,n) = s(4*i-1,4*j-1)
              an(i,j,n) = s(4*i,4*j)

	    end do !loop i
	  end do !loop j
        end do !loop n

        return
        end

c----------------------end stiffMDarcy----------------------------


c------------------------stableMDarcy-----------------------------

c  Subroutine: A subroutine is referenced by a CALL statement.

 	subroutine stableMDarcy(x, ien, ax, ay, az, an, px, py, pz,
     &                                  pn, iBC, iBCx, iBCy, iBCz)
c
c-----------------------------------------------------------------
c
	implicit none

        integer, parameter :: isize=200000

	integer ien(8,1), iBC(1), iBCx(1), iBCy(1), iBCz(1)
	integer ndf, nst
	integer nsdmm, nenmm, numnpmm, numelmm
        integer i, j, n
	real*8 x(3,1), ax(8,8,1), ay(8,8,1), az(8,8,1), an(8,8,1), 
     &         px(1), py(1), pz(1), pn(1),
     &         ux(isize), uy(isize), uz(isize), un(isize)
	real*8 ul(8,8), s(32,32), xl(3,8)

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
        print *,'ux = ',ux(4)
        print *,'uy = ',uy(4)
        print *,'uz = ',uz(4)
        print *,'un = ',un(4)
c
c.... writing results in files
c
        open(unit=1,file='ux1.txt');
        do n = 1, numnpmm
          write(1,*),x(1,n),x(2,n),x(3,n),ux(n)
        end do
        close(1)

        open(unit=1,file='uy1.txt');
        do n = 1, numnpmm
          write(1,*),x(1,n),x(2,n),x(3,n),uy(n)
        end do
        close(1)

        open(unit=1,file='uz1.txt');
        do n = 1, numnpmm
          write(1,*),x(1,n),x(2,n),x(3,n),uz(n)
        end do
        close(1)

        open(unit=1,file='un1.txt');
        do n = 1, numnpmm
          write(1,*),x(1,n),x(2,n),x(3,n),un(n)
        end do
        close(1)
c
c....   need to copy ul, xl, s
c
        do n = 1, numelmm

	  do j = 1, nenmm
	    ul(1,j) = ux(ien(j,n))
	    ul(2,j) = uy(ien(j,n))
	    ul(3,j) = uz(ien(j,n))
	    ul(4,j) = un(ien(j,n))

	    xl(1,j) = x(1,ien(j,n))
	    xl(2,j) = x(2,ien(j,n))
	    xl(3,j) = x(3,ien(j,n))

	    do i = 1, nenmm
              s(4*i-3,4*j-3) = ax(i,j,n) 
              s(4*i-2,4*j-2) = ay(i,j,n) 
              s(4*i-1,4*j-1) = az(i,j,n) 
	      s(4*i,4*j) = an(i,j,n)
	    end do

	  end do	
c          call elmtcheck(ul,xl,s)

	end do

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

        integer, parameter :: isize=200000

	integer ien(8,1), iBC(1)
	integer nsdmm, nenmm, numnpmm, numelmm, itrmm
	integer iter, n, i, j
	real*8 pap, rbr, alpha, beta, rbr0, tmp
	real*8 r(1), a(8,8,1), x(1)
	real*8 p(isize), q(isize), 
     &         B(isize), tol1

        common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

	data tol1 / 1.d0 /
        itrmm = 200
c        itrmm = 1

c
c.... initialize the variables
c
      do n = 1, numnpmm
        B(n) = 0.d0
	x(n) = 0.d0
      end do
      rbr = 0.
c
c.... computing B matrix for Dirichlet BCs
c
      do n = 1, numelmm
        do i = 1, nenmm
          B(ien(i,n)) =  B(ien(i,n)) + a(i,i,n)
        end do
      end do

c
c.... computing rbr
c
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
        do n = 1, numnpmm
          q(n) = 0.d0
        end do

        do n = 1, numelmm        
          do i = 1, nenmm
            do j = 1, nenmm
              q(ien(i,n)) = q(ien(i,n)) + a(i,j,n)*p(ien(j,n))
            end do
          end do
        end do

c        do n = 1, numnpmm
c          print *, n, q(n)
c        end do
c
c.... calculate alpha
c
        pap = 0.d0
	do n = 1, numnpmm
	  pap = pap + p(n) * q(n)
        end do

	alpha = rbr / pap 

c        print *, alpha,rbr,pap
c
c.... compute solution
c
        tmp = rbr

	rbr = 0.
	do n = 1, numnpmm
          if (iBC(n) .gt. 0) then 
	    x(n) = r(n)/B(n) !seems dirichlet BCs only (SA)
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

c----------------------end solveMMoveCG-----------------------


