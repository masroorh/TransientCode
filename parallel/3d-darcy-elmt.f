c
      subroutine elmt30(xl, s, p, c)
c
c.... 2 Node Linear Element
c 
c....    DARCYS LAW + CONSERVATION OF MASS
c
c....With stabilization scheme.
c
c....Copyright(c): A Masud (4/98)  
c....Copyright(c): A Masud (3/01)
c---------------------------------------------------------

        implicit none
c
c.... Subroutine Argument Definitions and Declarations:
c         d   = material set vector (material properties and constants)
c         xl  = nodal coordinate array
c         ix  = global node number array
c         s   = System Stiffness Matrix
c         c   = Damping matrix
c         p   = System Right-Hand-Side Matrix
c         ndf = 4 = no of degrees of freedom per node: v, p
c         ndm = 1 (x,y cartesian coordinates at the nodes)
c         nst = number of element equations (ndf*number of nodes in element)
c         isw = task control number from calling program
c
c....Subroutine Variable Definitions and Declarations (Non-shared):
c         lint  = number of quadrature pts (numerical int)
c         w     = quadrature point weight 
c         det   = jacobian determinant
c         shl   = ref element shape fnctn and derivatives 
c         shg   = global shape functions and derivatives 
c         xlq   = coord array of quadr pts

      real*8    d(8),xl(3,1),s(16,1),p(1)
      real*8    shp(4,4),shg3l(4,4),shg2l(4,4),shg1l(4,4)
      real*8    shg(4,4),w(4) !wts for gauss quad (SA)
      real*8    ue(4),duex(4),duey(4),duez(4)
      real*8    un,upnx,upny,upnz,eps1,eps2
      real*8    c1,c2,c3,c4,c5,c6,det
      real*8    h2,vis,del1,alpha,gf1,gf2,djx,djy,djz,djn
      real*8    dix,diy,diz,din
      real*8    kap,gc,rho,pi,phi
      real*8    xint,cx,cy,cz,xc,rc,radius,diag
      real*8    xs11,xs12,xs13,xs21,xs22,xs23,
     &          xs31,xs32,xs33,c(4,1)
      integer   ndf,nst,i,j,l,lint
      integer   nsdmm,nenmm,numnpmm,numelmm

      common /infmMM/ nsdmm, nenmm, numelmm, numnpmm

      data pi /3.14159265359d0/

      data shp /0.25d0, 0.25d0, 0.25d0, 0.25d0,
     &          0.25d0, 0.25d0, 0.25d0, 0.25d0,
     &          0.25d0, 0.25d0, 0.25d0, 0.25d0,
     &          0.25d0, 0.25d0, 0.25d0, 0.25d0/

      data shg1l /1.d0, 0.d0, 0.d0, -1.d0,
     &            1.d0, 0.d0, 0.d0, -1.d0,
     &            1.d0, 0.d0, 0.d0, -1.d0,
     &            1.d0, 0.d0, 0.d0, -1.d0/

      data shg2l /0.d0, 1.d0, 0.d0, -1.d0,
     &            0.d0, 1.d0, 0.d0, -1.d0,
     &            0.d0, 1.d0, 0.d0, -1.d0,
     &            0.d0, 1.d0, 0.d0, -1.d0/

      data shg3l /0.d0, 0.d0, 1.d0, -1.d0,
     &            0.d0, 0.d0, 1.d0, -1.d0,
     &            0.d0, 0.d0, 1.d0, -1.d0,
     &            0.d0, 0.d0, 1.d0, -1.d0/

      data w /1.0d0,1.0d0,1.0d0,1.0d0/

      ndf = 4
      nst = ndf * nenmm
       
c--------------------------------------------------------------------
c.... MATERIAL PROPERTY SET
c....     Map of {d}:
c             1 --- Dynamic Viscosity 
c             2 --- Permeability
c             3 --- Density
c             4 --- Epsilon 1
c             5 --- Epsilon 2
c             6 --- alpha
c             7 -- 1-gravity
c             8 -- 2-gravity
c--------------------------------------------------------------------
           d(1) = 1.d0 !dynamic viscosity of water (poise @ 20deg)
           d(2) = 2.d0 !permeability of silty gravels (cm/s)
cs           d(1) = 0.001002d0 !dynamic viscosity of water (poise @ 20deg)
cs           d(2) = 5.0d-6 !permeability of silty gravels (cm/s)
           d(3) = 0.9982071 !density of water (g/cm^3 @ 20deg)
           d(4) = 1.d0 !epsilon 1
           d(5) = 1.d0 !epsilon 2
           d(6) = 1.d-2 !alpha
           d(7) = 40.d0 !gravity 1 (cm/sec^2)
           d(8) = 980.d0 !gravity 2 (cm/sec^2)
c--------------------------------------------------------------------
c     
c.... Specify Integration Pt Loop Control for Full Integration
c
      lint = 4
c
c.... Compute Element Geometry Factor h
c
      h2 = 0.d0
      do i = 1,2
        h2 = h2+(xl(1,i)-xl(1,i+3))**2+(xl(2,i)-xl(2,i+3))**2
     &                          *(xl(3,i)-xl(3,i+3))**2
      end do
      h2 = dsqrt(h2)/2.d0 !average of two sides
      h2 = h2*h2

c.... Set Up Material Properties 
      vis = d(1)
      kap = d(2)
      rho = d(3)
      eps1 = d(4)
      eps2 = d(5)
      alpha = d(6)
      gf1 = 0.d0
      gf2 = 0.d0
      gc = 1.0d0 

      c2 = vis/kap !used in stiff matrix
      c3 = kap/vis !used in stiff matrix
      c4 = h2*c2 !used in stiff matrix
      c5 = 1.d0/vis !used in source term
      c6 = rho/gc !used in source term
c
c.... numerical integration
c
      do l = 1, lint

	do i = 1, nenmm
	  shg(1,i) = shg1l(i,l)
	  shg(2,i) = shg2l(i,l)
	  shg(3,i) = shg3l(i,l)
	  shg(4,i) = shp(i,l)
	end do
c
c.... get the jacobian
c
        xs11 = shg1l(1,l)*xl(1,1) + shg1l(2,l)*xl(1,2)
     &           + shg1l(3,l)*xl(1,3) + shg1l(4,l)*xl(1,4)
        xs12 = shg2l(1,l)*xl(1,1) + shg2l(2,l)*xl(1,2)
     &           + shg2l(3,l)*xl(1,3) + shg2l(4,l)*xl(1,4)
        xs13 = shg3l(1,l)*xl(1,1) + shg3l(2,l)*xl(1,2)
     &           + shg3l(3,l)*xl(1,3) + shg3l(4,l)*xl(1,4)
        xs21 = shg1l(1,l)*xl(2,1) + shg1l(2,l)*xl(2,2)
     &           + shg1l(3,l)*xl(2,3) + shg1l(4,l)*xl(2,4)
        xs22 = shg2l(1,l)*xl(2,1) + shg2l(2,l)*xl(2,2)
     &           + shg2l(3,l)*xl(2,3) + shg2l(4,l)*xl(2,4)
        xs23 = shg3l(1,l)*xl(2,1) + shg3l(2,l)*xl(2,2)
     &           + shg3l(3,l)*xl(2,3) + shg3l(4,l)*xl(2,4)
        xs31 = shg1l(1,l)*xl(3,1) + shg1l(2,l)*xl(3,2)
     &           + shg1l(3,l)*xl(3,3) + shg1l(4,l)*xl(3,4)
        xs32 = shg2l(1,l)*xl(3,1) + shg2l(2,l)*xl(3,2)
     &           + shg2l(3,l)*xl(3,3) + shg2l(4,l)*xl(3,4)
        xs33 = shg3l(1,l)*xl(3,1) + shg3l(2,l)*xl(3,2)
     &           + shg3l(3,l)*xl(3,3) + shg3l(4,l)*xl(3,4)

        det = xs11 * (xs22 * xs33 - xs23 * xs32) -
     &        xs12 * (xs21 * xs33 - xs23 * xs31) +
     &        xs13 * (xs21 * xs32 - xs22 * xs31)
cs	print*, det
cs        if (det.eq.0.d0) then
c          print*,'shahab'
cs          det = 0.001 !done for metis
cs        end if

        c1 = det*w(l)
c
c.... Transform to x, y, z derivatives (feap:shp3d) (SA)
c
cs        do i = 1, nenmm
cs          shg(1,i)=(xs11*shg1l(i,l)+xs21*shg2l(i,l)+xs31*shg3l(i,l))/det
cs          shg(2,i)=(xs12*shg1l(i,l)+xs22*shg2l(i,l)+xs32*shg3l(i,l))/det
cs          shg(3,i)=(xs13*shg1l(i,l)+xs23*shg2l(i,l)+xs33*shg3l(i,l))/det
cs        end do
c       
c.... Computing right hand side contribution of gravity & body forces
c
        cx = 0.d0
        cy = 0.d0
        cz = 0.d0
c
c.... Computing the location of the integration point to be used in applying phi
c
        do j = 1, nenmm
          cx = cx + shg(4,j)*xl(1,j)
          cy = cy + shg(4,j)*xl(2,j)
          cz = cz + shg(4,j)*xl(3,j)
        end do
c
c.... Added for tetrahedron
c
cs        cx = cx/nenmm
cs        cy = cy/nenmm
cs        cz = cz/nenmm
c
c.... form constants needed for stiffness mtx terms and right hand side
c
        do j = 1, nenmm
          djx = shg(1,j)*c1
          djy = shg(2,j)*c1
          djz = shg(3,j)*c1
          djn = shg(4,j)*c1
c     
c.... source terms      
c
	  call uexact_sin_2pi(cx,cy,cz,ue,duex,duey,duez,d)
cs	  p(ndf*j) = p(ndf*j) - djn*ue(4)
          phi = 8.d0*pi*pi*ue(4)
          p(ndf*j) = p(ndf*j)-c3*djn*phi
c
c.... added the following for vx, vy and vz (SA)
c
cs          p(ndf*j-3) = p(ndf*j-3) - djx*ue(4)
cs          p(ndf*j-2) = p(ndf*j-2) - djy*ue(4)
cs          p(ndf*j-1) = p(ndf*j-1) - djz*ue(4)
          p(ndf*j-3)=p(ndf*j-3)-0.5d0*c2*h2*alpha*djx*phi
          p(ndf*j-2)=p(ndf*j-2)-0.5d0*c2*h2*alpha*djy*phi
          p(ndf*j-1)=p(ndf*j-1)-0.5d0*c2*h2*alpha*djz*phi
c
c.... forming element stiffness matrix
c
          do i = 1, nenmm
            dix = shg(1,i)
            diy = shg(2,i)
            diz = shg(3,i)
            din = shg(4,i)

cs            s(3*i-2,3*j-2) = s(3*i-2,3*j-2) + 0.5d0*c2*dix*djx*eps1
cs     &                            + 0.5d0*c3*dix*djx - din*djn
cs            s(3*i-2,3*j-1) = s(3*i-2,3*j-1) + 0.d0
cs            s(3*i-2,3*j) = s(3*i-2,3*j) + 0.d0
cs
cs            s(3*i-1,3*j-2) = s(3*i-1,3*j-2) + 0.d0
cs            s(3*i-1,3*j-1) = s(3*i-1,3*j-1) + 0.5d0*c2*diy*djy*eps1
cs     &                            + 0.5d0*c3*diy*djy - din*djn
cs            s(3*i-1,3*j) = s(3*i-1,3*j) + 0.d0
cs
cs            s(3*i,3*j-2) = s(3*i,3*j-2) + 0.d0
cs            s(3*i,3*j-1) = s(3*i,3*j-1) + 0.d0
cs            s(3*i,3*j) = s(3*i,3*j) + eps2*din*djn/eps1

            s(4*i-3,4*j-3) = s(4*i-3,4*j-3) + 0.5d0*c2*din*djn
     &                             + 0.5d0*alpha*c2*h2*dix*djx
            s(4*i-3,4*j-2) = s(4*i-3,4*j-2) + 0.d0
            s(4*i-3,4*j-1) = s(4*i-3,4*j-1) + 0.d0
            s(4*i-3,4*j) = s(4*i-3,4*j) + 0.d0

            s(4*i-2,4*j-3) = s(4*i-2,4*j-3) + 0.d0
            s(4*i-2,4*j-2) = s(4*i-2,4*j-2) + 0.5d0*c2*din*djn
     &                             + 0.5d0*alpha*c2*h2*diy*djy
            s(4*i-2,4*j-1) = s(4*i-2,4*j-1) + 0.d0
            s(4*i-2,4*j) = s(4*i-2,4*j) + 0.d0

            s(4*i-1,4*j-3) = s(4*i-1,4*j-3) + 0.d0
            s(4*i-1,4*j-2) = s(4*i-1,4*j-2) + 0.d0
            s(4*i-1,4*j-1) = s(4*i-1,4*j-1) + 0.5d0*c2*din*djn
     &                             + 0.5d0*alpha*c2*h2*diz*djz
            s(4*i-1,4*j) = s(4*i-1,4*j) + 0.d0

            s(4*i,4*j-3) = s(4*i,4*j-3) + 0.d0
            s(4*i,4*j-2) = s(4*i,4*j-2) + 0.d0
            s(4*i,4*j-1) = s(4*i,4*j-1) + 0.d0
            s(4*i,4*j) = s(4*i,4*j) + 0.5d0*c3*(dix*djx+diy*djy)
     &                                  + eps1*0.5d0*c2*din*djn
c
c.... forming damping matrix (SA)
c
            c(i,j) = c(i,j) + din*djn

         end do !i loop
       end do !j loop
      end do !l loop

      return

      end !elmt30

c**********************************************************************
	subroutine uexact_cos_pi(cx,cy,ue,duex,duey,d)
c....	Written by Arif Masud.   (Fall 2001)
c....   This solution is valid over unit (square) domain
c....	It is a modification of the original problem designed by T. Hughes.
c**************************************************************************
	implicit none
        real*8 d(*)  
	real*8 ue(4),duex(4),duey(4),duez(4)
	real*8 p,vx,vy,vz,vxx,vyy,vzz,vxy,vyx
        real*8 vxz,vzx,vyz,vzy,cx,cy,cz
	real*8 pi,c3
 	Data pi /3.14159265359d0/
c	
c....	----------------> Calculate the Exact Solution <------------------
c
c     NOTE: lenght: l = 1.
      c3 = d(2)/d(1)

      p=Dcos(pi*cx)*Dcos(pi*cy)*Dcos(pi*cz)
      vx=-c3*pi*Dsin(pi*cx)*Dcos(pi*cy)*Dcos(pi*cz)
      vy=-c3*pi*Dcos(pi*cx)*Dsin(pi*cy)*Dcos(pi*cz)
      vz=-c3*pi*Dcos(pi*cx)*Dcos(pi*cy)*Dsin(pi*cz)

      vxx=-c3*pi*pi*p
      vyy=vxx
      vzz=vxx

      vxy=c3*pi*pi*Dsin(pi*cx)*Dsin(pi*cy)*Dcos(pi*cz)
      vyx=vxy

      vxz=c3*pi*pi*Dsin(pi*cx)*Dcos(pi*cy)*Dsin(pi*cz)
      vzx=vxz

      vyz=c3*pi*pi*Dcos(pi*cx)*Dsin(pi*cy)*Dsin(pi*cz)
      vzy=vyz

      ue(1) = vx
      ue(2) = vy
      ue(3) = vz
      ue(4) = p

      duex(1) = vxx
      duex(2) = vxy
      duex(3) = vxz
      duex(4) = -pi*Dsin(pi*cx)*Dcos(pi*cy)*Dcos(pi*cz)

      duey(1) = vyx
      duey(2) = vyy
      duey(3) = vyz
      duey(4) = -pi*Dcos(pi*cx)*Dsin(pi*cy)*Dcos(pi*cz)

      duez(1) = vzx
      duez(2) = vzy
      duez(3) = vzz
      duez(4) = -pi*Dcos(pi*cx)*Dcos(pi*cy)*Dsin(pi*cz)

      return
      end

c**********************************************************************
	subroutine uexact_cos_2pi(cx,cy,ue,duex,duey,d)
c....	Written by Arif Masud.   (Fall 2001)
c....   This solution is valid over unit (square) domain
c....	It is a modification of the original problem designed by T. Hughes.
c**************************************************************************
	implicit none
        real*8 d(1)  
	real*8 ue(4),duex(4),duey(4),duez(4)
	real*8 p,vx,vy,vz,vxx,vyy,vzz,vxy,vyx
        real*8 vxz,vzx,vyz,vzy,cx,cy,cz,pi,c3
 	Data pi/3.14159265359d0/
c	
c....	----------------> Calculate the Exact Solution <------------------
c
c     NOTE: lenght: l = 1.
      c3 = d(2)/d(1)

      p=Dcos(pi*cx)*Dcos(pi*cy)*Dcos(pi*cz)
      vx=-c3*pi*Dsin(pi*cx)*Dcos(pi*cy)*Dcos(pi*cz)
      vy=-c3*pi*Dcos(pi*cx)*Dsin(pi*cy)*Dcos(pi*cz)
      vz=-c3*pi*Dcos(pi*cx)*Dcos(pi*cy)*Dsin(pi*cz)

      vxx=c3*pi*pi*p
      vyy=vxx
      vzz=vxx

      vxy=c3*pi*pi*Dsin(pi*cx)*Dsin(pi*cy)*Dcos(pi*cz)
      vyx=vxy

      vxz=c3*pi*pi*Dsin(pi*cx)*Dcos(pi*cy)*Dsin(pi*cz)
      vzx=vxz

      vyz=c3*pi*pi*Dcos(pi*cx)*Dsin(pi*cy)*Dsin(pi*cz)
      vzy=vyz

      ue(1) = vx
      ue(2) = vy
      ue(3) = vz
      ue(4) = p

      duex(1) = vxx
      duex(2) = vxy
      duex(3) = vxz
      duex(4) =-pi*Dsin(pi*cx)*Dcos(pi*cy)*Dcos(pi*cz)

      duey(1) = vyx
      duey(2) = vyy
      duey(3) = vyz
      duey(4) =-pi*Dcos(pi*cx)*Dsin(pi*cy)*Dcos(pi*cz)

      duez(1) = vzx
      duez(2) = vzy
      duez(3) = vzz
      duez(4) =-pi*Dcos(pi*cx)*Dcos(pi*cy)*Dsin(pi*cz)

      return
      end

c**********************************************************************
	subroutine uexact_sin_2pi(cx,cy,cz,ue,duex,duey,duez,d)
c
c....	Written by Arif Masud.   (Fall 2000)
c....   This solution is valid over unit (square) domain
c
c**************************************************************************
	implicit none
        real*8 d(*)  
	real*8 ue(4),duex(4),duey(4),duez(4)
	real*8 p,vx,vy,vz,vxx,vxy,vxz,vyx,vyy,vyz,cx,cy,cz
        real*8 vzx,vzy,vzz
	real*8 pi,c3,c
 	Data pi/3.14159265359d0/
c
c....	----------------> Calculate the Exact Solution <------------------
c
c     NOTE: lenght: l = 1.
cs      c3 = -d(2)/d(1)
      c3 = d(2)/d(1)
      c = 2.d0*pi
      p = Dsin(c*cx)*Dsin(c*cy)*Dsin(c*cz)
      vx = -c*c3*Dcos(c*cx)*Dsin(c*cy)*Dsin(c*cz)
      vy = -c*c3*Dsin(c*cx)*Dcos(c*cy)*Dsin(c*cz)
      vz = -c*c3*Dsin(c*cx)*Dsin(c*cy)*Dcos(c*cz)

      vxx = -4.d0*c3*pi*pi*p
      vyy = vxx
      vzz = vxx

      vxy = -c*c3*pi*Dcos(c*cx)*Dcos(c*cy)*Dsin(c*cz)
      vyx = vxy

      vxz = -c*c3*pi*Dcos(c*cx)*Dsin(c*cy)*Dcos(c*cz)
      vzx = vxz

      vyz = -c*c3*pi*Dsin(c*cx)*Dcos(c*cy)*Dcos(c*cz)
      vzy = vyz

      ue(1) = vx
      ue(2) = vy
      ue(3) = vz
      ue(4) = p

      duex(1) = vxx
      duex(2) = vxy
      duex(3) = vxz
      duex(4) = -c*Dcos(c*cx)*Dsin(c*cy)*Dsin(c*cz)

      duey(1) = vyx
      duey(2) = vyy
      duey(3) = vyz
      duey(4) = -c*Dsin(c*cx)*Dcos(c*cy)*Dsin(c*cz)

      duez(1) = vzx
      duez(2) = vzy
      duez(3) = vzz
      duez(4) = -c*Dsin(c*cx)*Dsin(c*cy)*Dcos(c*cz)

      return
      end
c
c**************************************************************************
c
	subroutine uexact_rough(XI,YI,ue,duex,duey)
c
c....	Written by Kaiming Xia.   (Fall 2000)
c....   This solution is valid over unit (square) domain
c
c**************************************************************************
        implicit none
	real*8 ue(4),duex(4),duey(4)
	real*8 K2,PI,XI,YI,XP,YP,CXI,CXP,CYI,CYP,FI1,FI2,F,SXI,SXP,
     &       SYI,SYP,VXI,VXP,VYI,VYP,VX,VY
	real*8 VXX,VXXI,VXXP,VXY,VXYI,VXYP,VYY,VYYI,VYYP,
     &       VYX,VYXI,VYXP

      K2=0.5d0
      PI=22.d0/7.d0
c
c....	-----------------> Set up material properties <------------------
c

c	
c....	----------------> Calculate the Exact Solution <------------------
c
c     NOTE: lenght: l = 1.
	XP=-1.0+XI
	YP=-1.0+YI

     	CXI=1.d0-XI**2/2.d0+(1.d0+4.d0*K2)*XI**4/24.d0-
     1(1.d0+44.d0*K2+16.d0*K2**2)*XI**6/720.d0

	CXP=1.d0-XP**2/2.d0+(1.d0+4.d0*K2)*XP**4/24.d0-
     1(1.d0+44.d0*K2+16.d0*K2**2)*XP**6/720.d0
	 
	CYI=1.d0-YI**2/2.d0+(1.d0+4.d0*K2)*YI**4/24.d0-
     1(1.d0+44.d0*K2+16.d0*K2**2)*YI**6/720.d0 

	CYP=1.d0-YP**2/2.d0+(1.d0+4.d0*K2)*YP**4/24.d0-
     1(1.d0+44.d0*K2+16.d0*K2**2)*YP**6/720.d0

c ... Modification by Loula
      FI1=-dlog((1.d0-(CXI**2)*(CYI**2))/(CXI**2*CYI**2))/(4.d0*PI) 
      FI2= dlog((1.d0-(CXP**2)*(CYP**2))/(CXP**2*CYP**2))/(4.d0*PI)

	SXI=XI-(1.d0+K2)*XI**3/6.d0+(1.d0+14.d0*K2+16.d0*k2**2)*XI**5
     & /120.d0 
	SXP=XP-(1.d0+K2)*XP**3/6.d0+(1.d0+14.d0*K2+16.d0*k2**2)*XP**5
     & /120.d0
	SYI=YI-(1.d0+K2)*YI**3/6.d0+(1.d0+14.d0*K2+16.d0*k2**2)*YI**5
     & /120.d0
	SYP=YP-(1.d0+K2)*YP**3/6.d0+(1.d0+14.d0*K2+16.d0*k2**2)*YP**5
     &  /120.d0
	     
      VXI=1.d0/(2.d0*PI)*(CXI**2*CYI**2)/(1.d0-CXI**2*CYI**2)*(CXI*SXI*
     1  CYI**2*(CXI**2*CYI**2)+CXI*SXI*CYI**2*(1.d0-CXI**2*CYI**2))/
     2 (CXI**2*CYI**2)**2

      VXP=-1.d0/(2.d0*PI)*(CXP**2*CYP**2)/(1.d0-CXP**2*CYP**2)*(CXP*SXP*
     1  CYP**2*(CXP**2*CYP**2)+CXP*SXP*CYP**2*(1.d0-CXP**2*CYP**2))/
     2  (CXP**2*CYP**2)**2 

      VYI=1.d0/(2.d0*PI)*(CXI**2*CYI**2)/(1.d0-CXI**2*CYI**2)*(CYI*SYI*
     1 CXI**2*(CXI**2*CYI**2)+CYI*SYI*CXI**2*(1.d0-CXI**2*CYI**2))/
     2  (CXI**2*CYI**2)**2
	 
      VYP=-1.d0/(2.d0*PI)*(CXP**2*CYP**2)/(1.d0-CXP**2*CYP**2)*(CYP*SYP*
     1    CXP**2*(CXP**2*CYP**2)+CYP*SYP*CXP**2*(1.d0-CXP**2*CYP**2))/
     2   (CXP**2*CYP**2)**2

     	F=FI1+FI2
      VX=VXI+VXP
	VY=VYI+VYP

c THE FOLLOWING NEEDS TO BE FIXED. Arif, Sep 30,2001     	
c      differential of velocity

      VXXI=((-SXI**2+CXI**2)*(1.0+CYI**4)*(1.0-CXI**2*CYI**2)*(CXI**2+
     1 CYI**2)-CXI*SXI*(1.0+CYI**4)*(2.0*CXI*SXI*CYI**2*(CXI**2+CYI**2)
     2 +(1.-CXI**2*CYI**2)*(-2.0*CXI*SXI)))/
     3 (2*PI*((1.-CXI**2*CYI**2)*(CXI**2+CYI**2))**2) 

      VXYI=(CXI*SXI*(-4.0*CYI**3*SYI)*(1.0-CXI**2*CYI**2)*
     1 (CXI**2+CYI**2)-CXI*SXI*(1.0+CYI**4)*(2.0*CXI**2*CYI*SYI*
     2 (CXI**2+CYI**2)+(1.0-CXI**2*CYI**2)*(-2.0*CYI*SYI)))/
     3 (2*PI*((1.-CXI**2*CYI**2)*(CXI**2+CYI**2))**2) 

      VYYI=((-SYI**2+CYI**2)*(1.0+CXI**4)*(1.0-CXI**2*CYI**2)*
     1 (CXI**2+CYI**2)-CYI*SYI*(1.0+CXI**4)*(2.0*CYI*SYI*CXI**2*(CXI**2
     2 +CYI**2)+(1.-CXI**2*CYI**2)*(-2.0*CYI*SYI)))/
     3 (2*PI*((1.-CXI**2*CYI**2)*(CXI**2+CYI**2))**2) 
      
       VYXI=(CYI*SYI*(-4.0*CXI**3*SXI)*(1.0-CXI**2*CYI**2)*
     1 (CXI**2+CYI**2)-CYI*SYI*(1.0+CXI**4)*(2.*CYI**2*CXI*SXI*
     2 (CXI**2+CYI**2)+(1.0-CXI**2*CYI**2)*(-2.0*CXI*SXI)))/
     3 (2*PI*((1.-CXI**2*CYI**2)*(CXI**2+CYI**2))**2)


      VXXP=-((-SXP**2+CXP**2)*(1.0+CYP**4)*(1.0-CXP**2*CYP**2)*(CXP**2+
     1CYP**2)-CXP*SXP*(1.0+CYP**4)*(2.0*CXP*SXP*CYP**2*(CXP**2+CYP**2)+
     2(1.-CXP**2*CYP**2)*(-2.0*CXP*SXP)))/
     3(2*PI*((1.-CXP**2*CYP**2)*(CXP**2+CYP**2))**2)
	
      VXYP=-(CXP*SXP*(-4.0*CYP**3*SYP)*(1.0-CXP**2*CYP**2)*
     1(CXP**2+CYP**2)-CXP*SXP*(1.0+CYP**4)*(2.*CXP**2*CYP*SYP*
     2(CXP**2+CYP**2)+(1.0-CXP**2*CYP**2)*(-2.0*CYP*SYP)))/
     3(2*PI*((1.-CXP**2*CYP**2)*(CXP**2+CYP**2))**2)
      
      VYYP=-((-SYP**2+CYP**2)*(1.0+CXP**4)*(1.0-CXP**2*CYP**2)*
     1(CXP**2+CYP**2)- CYP*SYP*(1.0+CXP**4)*(2.0*CYP*SYP*CXP**2*(CXP**2
     2+CYP**2)+(1.-CXP**2*CYP**2)*(-2.0*CYP*SYP)))/
     3(2*PI*((1.-CXP**2*CYP**2)*(CXP**2+CYP**2))**2)

      VYXP=-(CYP*SYP*(-4.0*CXP**3*SXP)*(1.0-CXP**2*CYP**2)*
     1(CXP**2+CYP**2)-CYP*SYP*(1.0+CXP**4)*(2.*CYP**2*CXP*SXP*
     2(CXP**2+CYP**2)+(1.0-CXP**2*CYP**2)*(-2.0*CXP*SXP)))/
     3(2*PI*((1.-CXP**2*CYP**2)*(CXP**2+CYP**2))**2)

        VXX=VXXI+VXXP
        VXY=VXYI+VXYP
        VYY=VYYI+VYYP
        VYX=VYXI+VYXP

	ue(1)  = VX
	ue(2)  = VY
	ue(3)  = F

	duex(1) = VXX
	duex(2) = VYX
	duex(3) = -VX

	duey(1) = VXY
	duey(2) = VYY
	duey(3) = -VY
	return
	end


c*********************************************************************
      subroutine stcn30(ix,xl,ul,lint,s,dt,st,ndf,ndm,d,iprob)
c
c.... Program to take displacements evaluated at quadrature points 
c     and project it onto the nodes of the current element
c     storing it for later global assembly.
c
      implicit none
      
      integer         numnp,nummat,nen,neq,ipr
      common /cdata/  numnp,nummat,nen,neq,ipr

      integer         n,nel, numel
      common /eldata/ n,nel,numel

      integer   ndf,lint,ndm,j,k,l,ll,iprob
      integer   ix(nel)
      real*8    xl(ndm,*), ul(ndf,*), dt(numnp), st(numnp,*)
      real*8    s(nel,*)
      real*8    c1,xg,ue(3),duex(3),duey(3),w,det,cx,cy
      real*8    shg(3,9)!,shl(3,9), shls(3,9), shgs(3,9)
      real*8    d(*)  

      save

      call pzero(s,nel*nel)
      do l = 1,nel
          ll = iabs(ix(l))
            dt(ll) = dt(ll)+1.d0
            st(ll,1) = st(ll,1)+ul(1,l)
            st(ll,2) = st(ll,2)+ul(2,l)
            st(ll,3) = st(ll,3)+ul(3,l)
c     resultant velocity field
            st(ll,4) = dsqrt(st(ll,1)**2+st(ll,2)**2)
      end do
c
      if (iprob .ge. 2) then
      do l = 1,nel
      if(iprob.eq.2)call uexact_cos_pi (xl(1,l),xl(2,l),ue,duex,duey,d)
      if(iprob.eq.3)call uexact_cos_2pi(xl(1,l),xl(2,l),ue,duex,duey,d)
      if(iprob.eq.4)call uexact_sin_2pi(xl(1,l),xl(2,l),ue,duex,duey,d)
          ll = iabs(ix(l))
            dt(ll) = dt(ll)+1.d0
            st(ll,6) = st(ll,6)+ue(1)
            st(ll,7) = st(ll,7)+ue(2)
            st(ll,8) = st(ll,8)+ue(3)
c     resultant velocity field
            st(ll,9) = dsqrt(st(ll,6)**2+st(ll,7)**2)
      end do
	endif

      if (iprob .eq. 1) then
c.... GET THE EXACT SOLUTION PROJECTED
c.... Loop Over Quadrature Pts
      do l = 1,lint
c.... Compute Local & Global Element Shape Functions
c ... need to calculate by masroor
!         if(nel.eq.3.or.nel.eq.6)then
!           call shlt(l,lint,nel,w,shl,shls)
!           call shgt(xl,nel,shl,shls,nummat,nen,det,shg,shgs) 
!         elseif(nel.eq.4.or.nel.eq.9)then
!           call shlq(l,lint,nel,w,shl,shls)
!           call shgq(xl,nel,shl,shls,nummat,nen,det,shg,shgs) 
!         endif
! 	
        c1 = det*w

c     Computing the location of the integration point to be used in applying phi
        cx = 0.d0
	  cy = 0.d0
          do 150 j=1,nel
            k=ix(j)
            cx = cx + shg(3,j)*xl(1,j)
	      cy = cy + shg(3,j)*xl(2,j)
 150      continue

c.... Compute the exact quantities
          call uexact_rough  (cx,cy,ue,duex,duey)
c.... compute consistent projection operator
c	do i = 1,nel
c	  xg     = shg(3,i)*c1
c	  do j = 1,nel
c	    s(i,j) = s(i,j) + xg*shg(3,j)
c	  end do
c	end do
c.... compute lumped projection and assemble the stress integrals
	do j = 1,nel
	  ll = iabs(ix(j))
	  if(ll.gt.0) then
	    xg     = c1*shg(3,j)
	    dt(ll) = dt(ll) + xg
	    st(ll,6) = st(ll,6) + ue(1)*xg
	    st(ll,7) = st(ll,7) + ue(2)*xg
	    st(ll,8) = st(ll,8) + ue(3)*xg
c     resultant exact velocity field
          st(ll,9) = dsqrt(st(ll,6)**2+st(ll,7)**2)
	  endif
	end do
      end do
	endif
c
      return
      end

c--------------------------------------------------------------------

      subroutine pzero (v, nn)

       implicit  none

       integer   n,nn
       real*8    v(nn)

      save

       do n = 1,nn
         v(n) = 0.0d0
       end do

      end

c*********************************************************************
      subroutine elmtcheck(ul,xl)
c
c.... 2 Node Linear Element
c 
c....    DARCYS LAW + CONSERVATION OF MASS
c
c....With stabilization scheme.
c
c....Copyright(c): A Masud (4/98)  
c....Copyright(c): A Masud (3/01) 
c---------------------------------------------------------

        implicit none

c....Subroutine Argument Definitions and Declarations:
c         d   = material set vector (material properties and constants)
c         ul  = displacement vector (solution)
c         xl  = nodal coordinate array
c         ndf = 4 = no of degrees of freedom per node: v1, v2, p
c         nst = number of element equations (ndf*number of nodes in element)
c
c....Subroutine Variable Definitions and Declarations (Non-shared):
c         lint  = number of quadrature pts (numerical int)
c         w     = quadrature point weight 
c         det   = jacobian determinant
c         shl   = ref element shape fnctn and derivatives 
c         shg   = global shape functions and derivatives 

      real*8 un,upnx,upny,upnz,eps1,eps2
      real*8 shp(4,4),shg3l(4,4),shg2l(4,4),shg1l(4,4)
      real*8 d(8),ul(4,*),xl(3,*)  
      real*8 s(16,16),det
      real*8 c1,c2,c3,c4,c5,c6
      real*8 xs11,xs12,xs13,xs21,xs22,xs23,xs31,xs32,xs33	
      real*8 shg(4,4), w(4)
      real*8 h1,h2,vis,del1,alpha,gf1,gf2,djx,djy,djn
      real*8 dix,diy,diz,din
      real*8 kap,gc,rho,qtild,pi
      real*8 u(4),dux(4),duy(4),duz(4),ue(4),duex(4),duey(4),duez(4)
      real*8 el2(4),eprix(4),epriy(4),epriz(4)
      real*8 el2el(4),eprixel(4),epriyel(4),eprizel(4)
      real*8 rsdo,al2u,al2p,al2div,h1u,h1p,dl10
      real*8 xint,yint,zint,cx,cy,cz,xc,yc,zc,rc,radius,diag
      integer ndf,nst,iprob
      integer i,j,l,lint
      integer il,l1
      integer nsdmm,nenmm,numnpmm,numelmm

      common /infmMM/ nsdmm,nenmm,numelmm,numnpmm
      
        Data pi/3.14159265359d0/

      data shp /0.25, 0.25, 0.25, 0.25,
     &          0.25, 0.25, 0.25, 0.25,
     &          0.25, 0.25, 0.25, 0.25,
     &          0.25, 0.25, 0.25, 0.25/

      data shg1l /-1, 1, 0, 0,
     &            -1, 1, 0, 0,
     &            -1, 1, 0, 0,
     &            -1, 1, 0, 0/

      data shg2l /-1, 0, 1, 0,
     &            -1, 0, 1, 0,
     &            -1, 0, 1, 0,
     &            -1, 0, 1, 0/

      data shg3l /-1, 0, 0, 1,
     &            -1, 0, 0, 1,
     &            -1, 0, 0, 1,
     &            -1, 0, 0, 1/

      data w /1.00d0,1.00d0,1.00d0,1.00d0/

      ndf = 4
      nst = ndf * nenmm
      iprob = 4

c--------------------------------------------------------------------
c.... MATERIAL PROPERTY SET
c....     Map of {d}:
c             1 --- Dynamic Viscosity
c             2 --- Permeability
c             3 --- Density
c             4 --- Epsilon 1
c             5 --- Epsilon 2
c             6 --- alpha
c             7 -- 1-gravity
c             8 -- 2-gravity
c--------------------------------------------------------------------
cs           d(1) = 0.001002 !dynamic viscosity of water (poise @ 20deg)
cs           d(2) = 5.0d-6 !permeability of silty gravels (cm/s)
           d(1) = 1.d0 !dynamic viscosity of water (poise @ 20deg)
           d(2) = 1.d0 !permeability of silty gravels (cm/s)
           d(3) = 0.9982071 !density of water (g/cm^3 @ 20deg)
           d(4) = 5.d0 !epsilon 1
           d(5) = 5.d0 !epsilon 2
           d(6) = 1.d1 !alpha
           d(7) = 40.d0 !gravity 1 (cm/sec^2)
           d(8) = 980.d0 !gravity 2 (cm/sec^2)
c--------------------------------------------------------------------

c
c.... Specify Integration Pt Loop Control for Full Integration
c
      lint = 4
c
c.... Compute Element Geometry Factor h
c
cs      h2 = 0.d0
cs      do i = 1,2
cs        h2 = h2+(xl(1,i)-xl(1,i+2))**2+(xl(2,i)-xl(2,i+2))**2
cs     &                          *(xl(3,i)-xl(3,i+2))**2
cs      end do
c
c.... Set Up Material Properties
c
      vis = d(1)
      kap = d(2)
      rho = d(3)
      eps1 = d(4)
      eps2 = d(5)
      alpha = d(6)
      gf1 = 0.d0
      gf2 = 0.d0
      gc = 1.0d0

      c2 = vis/kap !used in stiff matrix
      c3 = kap/vis !used in stiff matrix
      c4 = h2*c2 !used in stiff matrix
      c5 = 1.d0/vis !used in source term
      c6 = rho/gc !used in source term

c
c	--------------------> Define Variables <----------------------
c
c     Some variables used in this subroutine
c     j: degree of freedom (1,...,ndof)
c     u(j): finite element solution 
c     du(j): derivative of finite element solution 
c     ue(j): exact solution
c     due(j): derivative of exact solution
c     el2(j): error in L2
c     epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c     epriel(j): error in the seminorm of H1 in the element domain
   
c
c.... clear the global arrays 
c
	call pzero (el2, ndf)
	call pzero (eprix, ndf)
	call pzero (epriy, ndf)
c
c.... clear the element arrays
c
	call pzero (el2el, ndf)
	call pzero (eprixel, ndf)
	call pzero (epriyel, ndf)
c
c.... Loop Over Quadrature Pts
c
      do l = 1, lint

        do i = 1, nenmm
          shg(1,i) = shg1l(i,l)
          shg(2,i) = shg2l(i,l)
          shg(3,i) = shg3l(i,l)
          shg(4,i) = shp(i,l)
        end do
c
c.... get the jacobian
c
        xs11 = shg1l(1,l)*xl(1,1) + shg1l(2,l)*xl(1,2)
     &           + shg1l(3,l)*xl(1,3) + shg1l(4,l)*xl(1,4)
        xs12 = shg2l(1,l)*xl(1,1) + shg2l(2,l)*xl(1,2)
     &           + shg2l(3,l)*xl(1,3) + shg2l(4,l)*xl(1,4)
        xs13 = shg3l(1,l)*xl(1,1) + shg3l(2,l)*xl(1,2)
     &           + shg3l(3,l)*xl(1,3) + shg3l(4,l)*xl(1,4)
        xs21 = shg1l(1,l)*xl(2,1) + shg1l(2,l)*xl(2,2)
     &           + shg1l(3,l)*xl(2,3) + shg1l(4,l)*xl(2,4)
        xs22 = shg2l(1,l)*xl(2,1) + shg2l(2,l)*xl(2,2)
     &           + shg2l(3,l)*xl(2,3) + shg2l(4,l)*xl(2,4)
        xs23 = shg3l(1,l)*xl(2,1) + shg3l(2,l)*xl(2,2)
     &           + shg3l(3,l)*xl(2,3) + shg3l(4,l)*xl(2,4)
        xs31 = shg1l(1,l)*xl(3,1) + shg1l(2,l)*xl(3,2)
     &           + shg1l(3,l)*xl(3,3) + shg1l(4,l)*xl(3,4)
        xs32 = shg2l(1,l)*xl(3,1) + shg2l(2,l)*xl(3,2)
     &           + shg2l(3,l)*xl(3,3) + shg2l(4,l)*xl(3,4)
        xs33 = shg3l(1,l)*xl(3,1) + shg3l(2,l)*xl(3,2)
     &           + shg3l(3,l)*xl(3,3) + shg3l(4,l)*xl(3,4)

        det = xs11 * (xs22 * xs33 - xs23 * xs32) -
     &        xs12 * (xs21 * xs33 - xs23 * xs31) +
     &        xs13 * (xs21 * xs32 - xs22 * xs31)

        c1 = det*w(l)
c
c....	clear the finite element solutions
c
       call pzero (u, ndf)
       call pzero (dux, ndf)
       call pzero (duy, ndf)
       call pzero (duz, ndf)
 
       xint = 0.d00
       yint = 0.d00
       zint = 0.d00

       do l1=1,nenmm

         xint = xint + shg(3,l1)*xl(1,l1)
         yint = yint + shg(3,l1)*xl(2,l1)
         zint = zint + shg(3,l1)*xl(3,l1)

         do j=1,ndf
	   dux(j) = dux(j) + shg(1,l1)*ul(j,l1)
	   duy(j) = duy(j) + shg(2,l1)*ul(j,l1)
	   duz(j) = duz(j) + shg(3,l1)*ul(j,l1)
           u(j) = u(j) + shg(3,l1)*ul(j,l1)
         end do
       end do
c
c	---------------------> Exact Solution <-----------------------
c
	if (iprob .eq. 1) then
          call uexact_rough(xint,yint,zint,ue,duex,duey,duez)
        end if
	if (iprob .eq. 2) then
          call uexact_cos_pi(xint,yint,zint,ue,duex,duey,duez,d)
        end if
	if (iprob .eq. 3) then
          call uexact_cos_2pi(xint,yint,zint,ue,duex,duey,duez,d)
        end if
	if (iprob .eq. 4) then
          call uexact_sin_2pi(xint,yint,zint,ue,duex,duey,duez,d)
        end if

c	---------------------> Error Evaluation <---------------------

c....	loop over nodal vector
c
	 do j = 1, ndf
	   un   = c1 * ((u(j)-ue(j))**2)
	   upnx = c1 * ((dux(j)-duex(j))**2)
	   upny = c1 * ((duy(j)-duey(j))**2)
	   upnz = c1 * ((duz(j)-duez(j))**2)
	   el2el(j)   = el2el(j)   + un
	   eprixel(j) = eprixel(j) + upnx
	   epriyel(j) = epriyel(j) + upny
	   eprizel(j) = eprizel(j) + upnz
         end do
      end do     ! End of integration loop

      do j = 1,ndf
        el2(j)   = el2(j)   + el2el(j)
	eprix(j) = eprix(j) + eprixel(j)
	epriy(j) = epriy(j) + epriyel(j)
	epriz(j) = epriz(j) + eprizel(j)
      end do
c
c	--------------------> Error in all domain <---------------------
c
c.... some variables
c

c	al2u    = || e ||       ( L2 norm of u )

c	al2p    = || psi ||     ( L2 norm of p )

c	al2div  = |  e  |       ( L2 norm of div(u) )

c	h1u     = |  e  |       ( H1 seminorm of u )

c	h1p     = |  e  |       ( H1 seminorm of p )

        rsdo  = 1.d-30
        dl10  = dlog(10.d0)
        al2u  = dsqrt( el2(1) + el2(2) + el2(3) )
	al2p  = dsqrt( el2(4) )
	h1u   = dsqrt( eprix(1) + eprix(2) + epriy(1) +  epriy(2) 
     &                                  + epriz(1) + epriz(2) )
	al2div= dsqrt( eprix(1) + epriy(2) + epriz(3) )
	h1p   = dsqrt( eprix(3) + epriy(3) + epriz(3) )
c
c.... calculate the log
c
        if(al2u.gt.rsdo)  al2u  = dlog(al2u) / dl10
        if(al2p.gt.rsdo)  al2p  = dlog(al2p) / dl10
        if(h1u.gt.rsdo)   h1u   = dlog(h1u)  / dl10
        if(al2div.gt.rsdo)al2div= dlog(al2div) / dl10
        if(h1p.gt.rsdo)   h1p   = dlog(h1p)  / dl10
c
c	-------------------------> Write <--------------------------
c
        open(unit=1,file='norm.txt',access='append')
	write(1,*), al2u,al2p,al2div,h1p
        close(1)
	
      return 
c
c     Exit Subroutine
c
      end
c	--------------------------> End <----------------------------
