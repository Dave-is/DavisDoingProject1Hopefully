!---------------------------------------------------
      program Project1
!---------------------------------------------------
      implicit none
      integer LE,UE,j,w,r,de,db,z,iu,il,jpp,jp,jm,jmm,t
      real*8 h,UB,LB,i,rx,dx,fl,fe,k1,k2,k3,k4
      character*80 fname
      real*8:: X(10000,10000)






!....set delta x
      dx = .01


!....set x interval (ib is low end, ie is high end)
      fl = 0.0
      fe = 5.+dx !....the +dx is to account for Fortran indexing arrays from 1
                 !....similar reasoning is used below with the UB variable



!....set timestep
      h = 0.01


!....set time interval (lowerbound/upperbound)
      UB = 30.+h
!....the +h is due to fortran indexing arrays from 1. Thus the entry X(j,t) is
!....actually for X(j,t-1), so we need 1 "additional" timestep to makeup for it.
      LB = 0.
      w = int((UB-LB)/h)
!....initialize x
      ! allocate(X(iu,w))

!.....Fortran doesn't like noninteger do loop conditions so we
!.....must provide it some integers. UpperEnd/LowerEnd
      UE = int(UB/h)
      LE = int((LB+h)/h)
      il = int(fl)
      iu = int((fe-fl)/dx)


!.....set up initial conditions/initialize our function
!....TESTING: Write IC to a file to make sure everything works.
      open(unit = 4, file = "IC.dat")
      do j = 1, iu
        X(j,1) = 1./((cosh(fl/2))**2)
        fl = fl+dx
        write(4,*)X(j,1)
      enddo
      ! open(unit = 11, file = "Solution.dat")
      ! do j = 1, iu
      !   X(j,50) = (-1./2.)*1./((cosh((fl-30)/2))**2)
      !   fl = fl+dx
      !   write(11,*)X(j,1)
      ! enddo
      ! close(11)
      close(4)
!....TESTING
!....iu should equal 5000, w, should equal 1000
! write(6,*) iu,w,il,iu,X(2000,1)


!....Integrate using Adams-Bashforth, a two-step linear multistep method
!....It requires two steps to begin, so we will compute the first unknown
!....time using Euler forward
      do j = 1, iu
        t=2
        !....Initializing these variables will make filling ghost cells more straightforward.
        jpp = j+2
        jp = j+1
        jm = j-1
        jmm = j-2
        !....Sech is "periodic" in the sense that sech(-x) = sech(x)
        !....So, in order to fill ghost cells, we will determine if
        !....any of the j-minus vars are out of bounds and fix if so.
        !....This is a bit strange with Fortran's indexing, but if jm is 0
        !....(which is really -1) then it should be 2 (which is really i=1) (so add 2)
        !....If jmm is -1 (or what should be -2), then the adjustment is to make
        !....it 3 (so add 4)
        if(jm.lt.1)jm = jm+2
        if(jmm.lt.1)jmm = jmm+4
        X(j,2) = X(j,1) + (-(X(jpp,t-1)+2*X(jp,t-1)+2*X(jm,t-1)-X(jmm,t-1))+ &
          (6*X(j,t-1)*((-X(jp,t-1)+(X(jm,t-1))))))*h
          write(6,*)X(j,2)
      enddo



       do t = 3, 150, 1
           write(fname,'(i4.4)'),t
           open(11, file = fname)
        do j = 1, iu
          !....All of this business is explained above
          jpp = j+2
          jp = j+1
          jm = j-1
          jmm = j-2
          if(jm.lt.1)jm = jm+2
          if(jmm.lt.1)jmm = jmm+4
          X(j,t) = X(j,t-1) + h*((3/2)*X(j,t-1)-(1/2)*(X(j,t-2)))
          ! write(6,*)X(j,t-2)
          write(11,*)X(j,t)

        enddo
        close(11)
       enddo



return
end
