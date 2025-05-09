      PROGRAM gravothermal
c**********************************************************************
c  + Version 2.0: Jan 2024
c----------------------------------------------------------------------
c  Program to solve the equations of gravothermal fluids to
c  study the gravothermal collapse.
c----------------------------------------------------------------------
c  Author: Frank van den Bosch                        Yale University  
c**********************************************************************

      INCLUDE 'paramfile.h'

      INTEGER   k,j,n,iN1,iN2,iN3,Niter1,Niter2,Niter3
      REAL*8    drho1,drho2,rho_old1,rho_old2,rho_orig,rmed,sigma_orig
      REAL*8    tcc,t0cmax,rho0min
      REAL      ai1,aver_iter1,aver_iter2,aver_iter3
      CHARACTER outfil1*60,outfil2*60,outfil3*60
      
      INTEGER   lblnk
      EXTERNAL  lblnk
      
c**********************************************************************
 
c---initialize elapsed time

      CALL get_time

c---read parameters
      
      CALL read_param
      
c---set parameters and specify characteristic quantities
      
      CALL set_param
      
      IF (imode.EQ.3) CALL integrate_potential

c---open output file to log progress
 
      outfil1='/time_evolution.dat'
      outfil2='/timestep.log'
      outfil3='/logfile'

      OPEN(98,file=moddir(1:lblnk(moddir))//outfil1(1:lblnk(outfil1)),
     &     status='UNKNOWN')
      OPEN(99,file=moddir(1:lblnk(moddir))//outfil2(1:lblnk(outfil2)),
     &     status='UNKNOWN')
      OPEN(97,file=moddir(1:lblnk(moddir))//outfil3(1:lblnk(outfil3)),
     &     status='UNKNOWN')

c---set up radial grid
       
      CALL setup_grid
      
c---set up radial grid, and compute initial masses from the
c   initial density profile
      
      CALL initialize_grid
      
c---write out initial profiles on a much more extended grid.
c   this is only for testing/plotting purposes (not used in computation)      

      CALL initial_profiles

c---compute various properties from density profile

      CALL get_properties

c---output column captions for interpretation

      WRITE(*,*)'           >>> Gravo-Thermal Collapse <<< '
      WRITE(*,*)' '
      WRITE(*,*)'   istep        time           dt         rho0       
     &vmax         CPU     lg[epsdt]  lg[drmax]  lg[dumax]  <Niter1>  <N
     &iter2>  <Niter3>'
      WRITE(*,*)'-------------------------------------------------------
     &------------------------------------------------------------------
     &----------------'

      WRITE(97,*)' '
      WRITE(97,*)'   istep        time           dt         rho0       
     & vmax         CPU     lg[epsdt]  lg[drmax]  lg[dumax]  <Niter1>  <
     &Niter2>  <Niter3>'
      WRITE(97,*)'------------------------------------------------------
     &------------------------------------------------------------------
     &-----------------'
 
c---integrate forward in time (stop when central density becomes
c   three times the original central density, or max time is reached)
 
      k = 0           ! global integration step counter (never reset)
      n = 0           ! counts steps since last stats write
      j = 0           ! counts profile output snapshots
      t = 0.0d0       ! current simulation time
      dt = 1.0d-6     ! initial time step (will be updated adaptively)
      dumax = eps_du  ! upper limit on relative change in u
      drmax = eps_dr  ! upper limit on relative change in radius
      rmed = r(1)/2.0d0 ! initial estimate of median radius
      ai1=1.0         ! UNCLEAR WHAT THIS IS
      
      WRITE(98,*)k,t,rho(1),maxvel,t*t0_Gyr,rho(1)*rho_s_Msunpc3,
     &           maxvel*v0,Mtotal,rc1,rc2,rc3,rhochar1,rhochar2,
     &           rhochar3,rhochar4,rhochar5,rhochar6,rhochar7,rhochar8,
     &           rmed,minKn,r_BH,M_BH
      WRITE(*,88) k,t,dt,rho(1),maxvel,dt_tot,DLOG10(eps_dt),
     &            DLOG10(drmax),DLOG10(dumax),ai1,ai1,ai1
      WRITE(97,88)k,t,dt,rho(1),maxvel,dt_tot,DLOG10(eps_dt),
     &            DLOG10(drmax),DLOG10(dumax),ai1,ai1,ai1

c---initialize few parameters that quantify specific epochs, and 
c   parameters that are used to gauge evolution in central density

      tcc = -1.0d0
      t0cmax = 0.0d0
      rho_old1 = rho(1)
      rho_old2 = rho(1)
      rho_orig = rho(1)
      rho0min = rho(1)
      apply_impulse_iter = .FALSE.
      applied_impulse = .FALSE.

c---integrate time steps until stopping criterion is reached
      
      Niter1 = 0.0
      Niter2 = 0.0
      Niter3 = 0.0
      DO WHILE (t.LT.tstop)

        CALL set_time_step
        IF (k.EQ.0) dt=1.0E-7

c---ensure that an integration is performed at tfly

        IF (flyby_on .AND. .NOT. applied_impulse .AND.
     &      t .LT. tfly .AND. t + dt .GE. tfly) THEN
          dt = tfly - t
          apply_impulse_iter = .TRUE.
        END IF

c---increment counters

        k = k + 1
        n = n + 1
        t = t + dt

c---check if we are at the flyby time

        ! IF ( apply_impulse_iter ) THEN
        !   CALL apply_impulse
        !   applied_impulse = .TRUE.
        ! END IF

        CALL integrate_time_step(k,iN1,iN2,iN3)

        drho1 = ABS(rho(1)-rho_old1)/rho_old1
        drho2 = ABS(rho(1)-rho_old2)/rho_old2
        rmed = r(1)/2.0d0

        Niter1 = Niter1 + iN1
        Niter2 = Niter2 + iN2
        Niter3 = Niter3 + iN3
        
c---keep track of density and time at maximum core

        IF (rho(1).LT.rho0min) THEN
          rho0min = rho(1)
          t0cmax = t
        END IF

c---stop criterion (we have reached well into core collapse)

        IF (rho(1).GT.1500.0.AND.t.GT.50.0) GOTO 101
        
c---if density becomes NaN, terminate

        IF (isnan(rho(1))) GOTO 102

c---keep track of core collapse time, defined as time when central
c   density becomes 1000x minimum value over its past history

        IF (rho(1).GT.(1000.0d0*rho0min).AND.tcc.LT.0.0d0) THEN
          tcc = t
        END IF
          
c---write certain statistics to file (everytime when log of central 
c    density has changed by >0.02)

        IF (drho1.GT.0.02) THEN 
          n = 0  
          rho_old1 = rho(1)
          CALL get_properties
          WRITE(98,*)k,t,rho(1),maxvel,t*t0_Gyr,rho(1)*rho_s_Msunpc3,
     &           maxvel*v0,Mtotal,rc1,rc2,rc3,rhochar1,rhochar2,
     &           rhochar3,rhochar4,rhochar5,rhochar6,rhochar7,rhochar8,
     &           rmed,minKn,r_BH,M_BH
        END IF  

c---write profiles to file (everytime when log of central 
c    density has changed by >0.1)

        IF (drho2.GT.0.1 .OR. apply_impulse_iter) THEN 
          j = j + 1
          rho_old2 = rho(1)
          CALL write_output(j)
          WRITE(99,*)j,t,t*t0_Gyr
        END IF  

c---update logfile

        IF (MOD(k,100000).EQ.0 .OR. apply_impulse_iter) THEN 
          CALL get_time
          aver_iter1 = FLOAT(Niter1)/100000.0
          aver_iter2 = FLOAT(Niter2)/100000.0
          aver_iter3 = FLOAT(Niter3)/100000.0
          Niter1=0
          Niter2=0
          Niter3=0
          WRITE(*,88) k,t,dt,rho(1),maxvel,dt_tot,
     &      DLOG10(dtmax),DLOG10(drmax),DLOG10(dumax),
     &      aver_iter1,aver_iter2,aver_iter3
          WRITE(97,88) k,t,dt,rho(1),maxvel,dt_tot,
     &      DLOG10(dtmax),DLOG10(drmax),DLOG10(dumax),
     &      aver_iter1,aver_iter2,aver_iter3
        END IF

        IF (apply_impulse_iter) apply_impulse_iter = .FALSE.

      END DO
  101 CONTINUE
  
c---write details of final state of system to files

      CALL get_properties
      
      WRITE(98,*)k,t,rho(1),maxvel,t*t0_Gyr,rho(1)*rho_s_Msunpc3,
     &  maxvel*v0,Mtotal,rc1,rc2,rc3,rhochar1,rhochar2,
     &  rhochar3,rhochar4,rhochar5,rhochar6,rhochar7,rhochar8,
     &  rmed,minKn,r_BH,M_BH
      WRITE(97,88) k,t,dt,rho(1),maxvel,dt_tot,
     &  DLOG10(dtmax),DLOG10(drmax),DLOG10(dumax),
     &  aver_iter1,aver_iter2,aver_iter3

  102 CONTINUE
      CLOSE(98)
      CLOSE(99)
     
      CALL write_time(1)

c---output final stats

      WRITE(97,*)' '
      WRITE(97,*)'-------------------------------------------------------
     &-----------'
      WRITE(97,87)cvir,sigma_m,sigma_m*sigma0,rho_orig,rho0min,
     &           t0cmax,tcc,t0_Gyr,rho_s_Msunpc3,r_BH,M_BH      
      WRITE(97,*)'-------------------------------------------------------
     &-----------'
      CLOSE(97)

c---
           
 87   FORMAT(F7.3,2X,F8.4,2X,F9.5,2X,2(E12.7,2X),6(F13.6,2X))           
 88   FORMAT(I9,2X,F10.4,3(2X,E11.3),2X,4(F10.4,2X),3(F8.2,2X))

      STOP
      END

c***********************************************************************

      SUBROUTINE initialize_grid
c-----------------------------------------------------------------------
c  Initialize physical quantities on the radial grid:
C  - Enclosed mass M(i)
C  - Velocity dispersion v2(i)
C  - Density rho(i)
C  - Pressure P(i)
C  - Internal energy u(i)
C  Also outputs a diagnostic file with hydrostatic equilibrium checks.
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'
      
      INTEGER   i,ierr
      REAL*8    rbuf,rrr,rmed,dP,drho,dr,Yvec,rho1,Kn,SS,SSerr
      CHARACTER outfile*60

      INTEGER   lblnk
      EXTERNAL  lblnk

      REAL*8    sigr_NFW,Menc,toint3
      EXTERNAL  sigr_NFW,Menc,toint3
       
c-----------------------------------------------------------------------
C Begin grid initialization
c-----------------------------------------------------------------------

      WRITE(*,*)' ...Initializing Grid...'
      WRITE(*,*)' '

c--- If using truncated NFW (imode = 2), precompute ftr from integral
      IF (imode.EQ.2) THEN
        CALL dqagi(toint3,rtr,1,1.0D-6,1.0D-4,SS,SSerr,Neval,
     &          ierr,Nlimit,Nlenw,last,iwork,work)
        ftr = SS
      END IF
      
c--- If using flat-topped core profile (imode = 5), define central density
      IF (imode.EQ.5) THEN
        rho(1) = 1.0d0/(r(1) * (1.0d0+r(1))**2)  
      END IF  

c-----------------------------------------------------------------------
C Loop over radial grid and compute physical quantities
c-----------------------------------------------------------------------

      M(0) = 0.0d0
      maxvel = 0.0d0
      minKn = 9.9D+19

      DO i=1,Ngrid
        rrr = r(i)                                 ! Outer radius of shell
        rmed=(r(i)+r(i-1))/2.0d0                   ! Midpoint of shell
        IF (imode.EQ.5) rmed = MAX(rmed,r(1))      ! Avoid rmed < r(1) if flat core

        M(i)   = Menc(rrr)                         ! Enclosed mass at r(i)
        v2(i)  = sigr_NFW(rmed)                    ! Velocity dispersion squared
        rho(i) = 3.0d0 * ((M(i)-M(i-1)) /          ! Shell-averaged density
     &                 (r(i)**3 - r(i-1)**3))
        P(i)   = rho(i) * v2(i)                    ! Pressure (isotropic, ideal gas)
        u(i)   = 1.5d0 * v2(i)                     ! Internal energy (3/2 * v^2)

        Kn = (1.0d0/sigma_m) / DSQRT(P(i))         ! Knudsen number
        maxvel = MAX(DSQRT(v2(i)),maxvel)
        minKn  = MIN(Kn,minKn)
      END DO

c-----------------------------------------------------------------------
C Apply central smoothing if using regular NFW profile (imode = 1)
C Helps reduce artificial gradients in innermost cell
c-----------------------------------------------------------------------

      IF (imode.EQ.1) THEN
        rho1 = 1.0d0/(r(1) * (1.0d0+r(1))**2)      ! Idealized central value
        rho(1) = 2.0d0 * rho1 - rho(2)             ! Extrapolated central density
        P(1) = P(2) - ((r(2)-r(0))/(r(3)-r(1))) *   ! Linear extrapolated pressure
     &         (P(3)-P(2))
        v2(1) = P(1)/rho(1)                        ! Updated dispersion
        u(1) = 1.5d0 * v2(1)                       ! Updated energy
      END IF

c-----------------------------------------------------------------------
C Write diagnostics to initial_grid.dat:
C Contains values needed to assess hydrostatic equilibrium
c-----------------------------------------------------------------------

      outfile='/initial_grid.dat'
      OPEN(13,file=moddir(1:lblnk(moddir))//outfile(1:lblnk(outfile)),
     &     status='UNKNOWN')

      DO i=1,Ngrid-1
        rrr = r(i)
        rmed=(r(i)+r(i-1))/2.0d0
        IF (imode.EQ.5) rmed = MAX(rmed,r(1))

        dP   = P(i+1)  - P(i)
        drho = rho(i+1)+ rho(i)
        dr   = r(i+1)  - r(i-1)
        IF (imode.EQ.5.AND.i.EQ.1) dr = r(i+1) - r(i)

        Yvec = -((4.0d0/M(i)) * ((r(i)**2)/dr) * (dP/drho)) - 1.0d0
        Kn   = (1.0d0/sigma_m) / DSQRT(P(i))

        WRITE(13,99)i,rrr,rmed,M(i),v2(i),rho(i),P(i),u(i),Yvec,Kn
      END DO
      CLOSE(13)

c-----------------------------------------------------------------------
C Write full initial grid state to standard output files
c-----------------------------------------------------------------------

      CALL write_output(0) 
              
 99   FORMAT(I4,2X,9(E14.6,2X))

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE initial_profiles
c-----------------------------------------------------------------------
c  Compute enclosed mass, density, velocity dispersion, etc. on a
c  large grid spanning much wider range than used in computation.
c  ONLY USED FOR PLOTTING PURPOSES
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'
      
      INTEGER   i
      REAL*8    xlgr,rrr,yy1,yy2,yy3,yy4,yy5,yy6,yy7,yy8,yy4a,yy4b
      CHARACTER outfile*60

      INTEGER   lblnk
      EXTERNAL  lblnk

      REAL*8    sigr_NFW,Menc,Mbaryon,denss
      EXTERNAL  sigr_NFW,Menc,Mbaryon,denss
       
c---

      WRITE(*,*)' ...Computing Initial Profiles...'
      WRITE(*,*)' '

c---compute enclosed mass, velocity dispersion, density, pressure
c   and internal energy (in dimensionless units)
c   NOTE:  r(i) is already dimensionless  (r/rs)

      outfile='/initial_profiles.dat'
      OPEN(13,file=moddir(1:lblnk(moddir))//outfile(1:lblnk(outfile)),
     &     status='UNKNOWN')

      DO i=1,500
        xlgr = -4.9 + FLOAT(i-1)/499.0 * 7.0
        rrr = 10.0**xlgr
        yy1 = Menc(rrr)
        yy2 = Mbaryon(rrr)
        yy3 = denss(rrr)
        yy4 = sigr_NFW(rrr)
        yy4a= sigr_NFW(0.99*rrr)
        yy4b= sigr_NFW(1.01*rrr)
        yy5 = (yy4b-yy4a) / (0.02*rrr)
        yy6  = yy3 * yy4
        yy7  = 1.5d0  * yy4
        yy3  = DLOG10(yy3)
        yy8 = (1.0d0/sigma_m) / DSQRT(yy6)
        WRITE(13,91)i,xlgr,rrr,yy1,yy2,yy3,yy4,yy5,yy6,yy7,yy8
      END DO
      CLOSE(13)
            
 91   FORMAT(I4,2X,10(E14.6,2X))

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE setup_grid
c-----------------------------------------------------------------------
c  Set up a radial grid
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'
      
      INTEGER  i,j,imax
      REAL*8   xlgr,testrho,rtest(Ngrid_max)
      
      EXTERNAL rhotilde4
      REAL*8   rhotilde4
      
c---

      
c---in case of exponential truncation, adjust the maximum radius of the
c   grid (if needed)

      IF (imode.EQ.4) THEN
        imax = 0
        DO i=1,Ngrid
          xlgr = xlgrmin + DBLE(i-1)/DBLE(Ngrid-1) * (xlgrmax-xlgrmin)
          rtest(i) = 10.0d0**xlgr
          testrho = DLOG10(rhotilde4(rtest(i)))
          IF (testrho.LT.-9.0.AND.imax.EQ.0) imax = i-1
        END DO
        IF (imax.GT.0) xlgrmax = DLOG10(rtest(imax))
      END IF
      
c---make radial grid to be used

      r(0) = 0.0d0
      DO i=1,Ngrid
        xlgr = xlgrmin + DBLE(i-1)/DBLE(Ngrid-1) * (xlgrmax-xlgrmin)
        r(i) = 10.0d0**xlgr
      END DO
      
c---write to screen important quantities to screen & logfile
      
      CALL OUTPUT_PARAMS
      
      RETURN
      END

c***********************************************************************

      SUBROUTINE apply_impulse
c-----------------------------------------------------------------------
c  Apply an impulsive velocity kick to each shell due to a flyby
c  Formula: Δv² = (2/3) * 4 * G² * M_p² * r² / (v_p² * b⁴)
C  Adds to the velocity dispersion squared array v2(i)
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER i
      REAL*8  rmed_i,vkick2

      DO i = 1, Ngrid
         rmed_i = 0.5d0 * (r(i) + r(i-1))
         vkick2 = (8.d0 * Mfly**2 * rmed_i**2) / 
     &            (3.d0 * vfly**2 * bfly**4)
         v2(i) = v2(i) + vkick2
         u(i) = 1.5d0 * v2(i)
         P(i) = rho(i) * v2(i)
      END DO

      RETURN
      END

c***********************************************************************

      SUBROUTINE integrate_time_step(istep,iter1,iter2,iter3)
c-----------------------------------------------------------------------
c  Perform a single time step t --> t + dt
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i,istep,iter1,iter2,iter3
      REAL*8   a0(Ngrid),a1(Ngrid),a2(Ngrid),a3(Ngrid)
      REAL*8   a4(Ngrid),a5(Ngrid)

c---

c---store current state in case we need to revert back to it
      DO i=1,Ngrid
        a0(i) = M(i)
        a1(i) = v2(i)
        a2(i) = r(i)
        a3(i) = P(i)
        a4(i) = u(i)
        a5(i) = rho(i)
      END DO
              
c---compute the luminosities for all grid points

      CALL compute_luminosities

c----conduct heat across neighboring cells. 

      iter1 = 1
      iter2 = 1
      iter3 = 1
 33   CALL conduct_heat
 
c---if maximum conduction too large, revert back and use smaller
c   timestep. Iterate if needed

      IF (dumax.GT.eps_du) THEN
        t = t - dt
        dt = dt * 0.95d0 * (eps_du/dumax)
        t = t + dt
        DO i=1,Ngrid
          M(i) = a0(i)
          v2(i)= a1(i)
          r(i) = a2(i)
          P(i) = a3(i)
          u(i) = a4(i)
          rho(i)=a5(i)
        END DO
        iter1 = iter1 + 1
        IF (iter1.EQ.10) CALL Terminate('No convergence achieved iter1')
        GOTO 33
      END IF

c---evaporate mass due to collisions with host particles
 
      IF (Gamma_evap.NE.0.0d0) CALL evaporate

C---Heating due to impulsive encounter
      IF ( apply_impulse_iter ) THEN
        CALL apply_impulse
        applied_impulse = .TRUE.
      END IF

c---revirialize

 34   CONTINUE
      CALL revirialize
 
      IF (Failed) THEN
        t = t - dt
        dt = dt / 2.0d0
        t = t + dt
        DO i=1,Ngrid
          M(i) = a0(i)
          v2(i)= a1(i)
          r(i) = a2(i)
          P(i) = a3(i)
          u(i) = a4(i)
          rho(i)=a5(i)
        END DO
        iter2 = iter2 + 1
        IF (iter2.EQ.10) CALL Terminate('No convergence achieved iter3')
        WRITE(*,*)' TimeStep Iteration:',iter2,dt
        GOTO 33
      END IF

c---check to see if relaxation step was successful 
c   NOTE: in first time step we accept larger dr in order to allow 
c         initial numerical virialization

      IF (drmax.GT.eps_dr.AND.istep.NE.1) THEN
        iter3 = iter3 + 1
        IF (iter3.LE.20000) THEN
          GOTO 34
        ELSE
          WRITE(*,*)'WARNING: Max revirialization iterations reached'
cc          CALL Terminate('Maximum revirialization iterations reached')
        END IF
      END IF

c---update the actual eps_dt values used
      
      dtmax = dt/trelmin
      
      RETURN
      END
      
c***********************************************************************

      SUBROUTINE compute_luminosities
c-----------------------------------------------------------------------
c  explain
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   dTdr,fac1,fac2,Pmed,vmed
      REAL*8   y1,y2,x1,x2
      
c---

c---New method of computing the temperature gradients [CORRECT]

      DO i=1,Ngrid-1
        IF (imode.EQ.5.AND.i.EQ.1) THEN
          dTdr = (v2(2) - v2(1)) / (r(2) - r(1))
          vmed = DSQRT(v2(1))
          Pmed = P(1)
        ELSE 
          dTdr = (v2(i+1) - v2(i)) / (r(i+1) - r(i-1))
          vmed = DSQRT((v2(i+1) + v2(i)) / 2.0d0)  
          Pmed = (P(i+1) + P(i)) / 2.0d0  
        END IF
        fac1 = -3.0d0 * vmed * r(i)**2
        fac2 = (apar/bpar) * (sigma_m)**2 + ((1.0d0/cpar) / Pmed)
        L(i) = (fac1/fac2) * dTdr 
      END DO

c---boundary conditions

      L(0) = 0.0d0

c---use linear extrapolation to compute the flux out of the grid
      
      y2 = L(Ngrid-2)
      y1 = L(Ngrid-1)
      x2 = r(Ngrid-2)
      x1 = r(Ngrid-1)
      
c---two options for grid-boundary. I have experimented with both
c   and find no significant differences
      
cc      L(Ngrid) = y2 + ((r(Ngrid)-x2)/(x1-x2)) * (y1-y2)
      L(Ngrid) = 0.0d0
            
      RETURN
      END

c***********************************************************************

      SUBROUTINE conduct_heat
c-----------------------------------------------------------------------
c  Conduct heat and adjust internal energies accordingly. During
c  this step we ignore PdV work and we assume that densities of the
c  cells remain fixed. As a consequence of the change in u, the 
c  pressure changes are computed as well, which are used in the next
c  step to find the new hydrostatic equilibrium 
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   dudt,du
            
c---

c---calculate the temporal change in internal energy for each cell,
c   and update internal energy and corresponding pressure

      dumax = 0.0d0 
      DO i=1,Ngrid
        dudt = -(L(i) - L(i-1)) / (M(i) - M(i-1))
        du = dudt * dt
        u(i) = u(i) + du
        IF (DABS(du/u(i)).GT.dumax) THEN
          dumax=DABS(du/u(i))
          idumax = i
        END IF  
        P(i) = (2.0d0/3.0d0) * rho(i) * u(i)
      END DO

      RETURN
      END

c***********************************************************************

      SUBROUTINE evaporate
c-----------------------------------------------------------------------
c  Evaporate mass using specified evaporation rate
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   Mcum,Mcell
      REAL*8   Mnew(0:Ngrid)
      
c---

c---first lower masses in each bin and compute cumulative mass and
c   update mass grid. Also adjust densities, such that these are
c   in agreement. 

      Mcum = 0.0d0
      Mnew(0) = 0.0D0

      DO i=1,Ngrid
        Mcell = M(i) - M(i-1)
        Mcell = Mcell - Gamma_evap * Mcell * dt
        rho(i)= rho(i)- Gamma_evap * rho(i)* dt
        Mcum  = Mcum + Mcell
        Mnew(i)  = Mcum
      END DO

      DO i=1,Ngrid
         M(i) = Mnew(i)
      END DO
      
      Mtotal = Mcum
            
      RETURN
      END

c***********************************************************************

      SUBROUTINE revirialize
c-----------------------------------------------------------------------
c  Solve tridiagonal set of equations to re-establish hydrodynamic
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

c---

      IF (iRevir.EQ.1) CALL revirialize1
      IF (iRevir.EQ.2) CALL revirialize2
      IF (iRevir.EQ.3) CALL revirialize3

      END

c***********************************************************************

      SUBROUTINE revirialize1
c-----------------------------------------------------------------------
c  Solve tridiagonal set of equations to re-establish hydrodynamic
c  Iteratively re-establish hydrostatic equilibrium until convergence
c  criterion is met. Each step we solve a tridiagonal set of equations.
c
c  My old method
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER i
      REAL*8  dP,dr,drho,dd,r3a,r3b,r3c,r3d,q1,q2,c1,c2,dV_over_V,Kn,MM
      REAL*8  avec(Ngrid-1),bvec(Ngrid-1),cvec(Ngrid-1)
      REAL*8  Xvec(Ngrid-1),Yvec(Ngrid-1)
          
      REAL*8  Mbaryon
            
c---
     
      Failed =.FALSE.

c---determine elements of tridiagonal matrix

      DO i=1,Ngrid-1
        dP   = P(i+1)  - P(i)
        drho = rho(i+1)+ rho(i)
        dr   = r(i+1)  - r(i-1)
         
        drho = MAX(1.0D-20,drho)
        IF (dP.EQ.0.0d0)   CALL Terminate('dP zero') 
        IF (drho.EQ.0.0d0) CALL Terminate('drho zero') 

        r3a = r(i+1)**3 / (r(i+1)**3 - r(i)**3)
        r3c = r(i)**3   / (r(i)**3 - r(i-1)**3)
        r3b = r3a - 1.0d0
        r3d = r3c - 1.0d0
        
        q1 = r(i+1) / dr
        q2 = q1 - 1.0d0

        MM = M(i) + Mbaryon(r(i)) 
        dd = -((4.0d0/MM) * (r(i)**2/dr) * (dP/drho))
        
        c1 = 5.0d0 + 5.0d0 * dd * (P(i)/dP) - 3.0d0 * (rho(i+1)/drho)
        c2 = c1 - 2.0d0   

        Yvec(i) = dd - 1.0d0
        
        avec(i) = r3d * c2 - q2
        bvec(i) =-(r3b * c1 + r3c * c2 + 2.0d0)
        cvec(i) = r3a * c1 + q1
      END DO
        
c---solve tridiagonal equation

      CALL tridag(avec,bvec,cvec,Yvec,Xvec,Ngrid-1)

c---update radii, and corresponding properties

      drmax = 0.0d0
      maxvel = 0.0d0
      minKn = 9.9d+19
      DO i=1,Ngrid
        r3c = r(i)**3   / (r(i)**3 - r(i-1)**3)

        IF (i.NE.Ngrid) THEN
          IF (DABS(Xvec(i)).GT.drmax) THEN
            drmax = DABS(Xvec(i))
            idrmax= i
          END IF  
          r(i) = r(i) * (1.0d0 + Xvec(i))
        END IF
          
        IF (i.EQ.1) THEN
          dV_over_V = 3.0 * r3c * Xvec(i)
        ELSE
          IF (i.EQ.Ngrid) THEN
            dV_over_V = -3.0 * (r3c-1.0d0) * Xvec(i-1)
          ELSE
            dV_over_V = 3.0 * (r3c*Xvec(i) - (r3c-1.0d0)*Xvec(i-1))
          END IF  
        END IF        

        P(i) = P(i) * (1.0d0 - (5.0d0/3.0d0) * dV_over_V)
        rho(i) = rho(i) * (1.0d0 - dV_over_V)
        v2(i) = P(i)/rho(i)

        IF (v2(i).LT.0.0d0) THEN
          WRITE(*,*)' Failed ',i,t,dt,P(i),rho(i),dumax
          Failed=.TRUE.
          RETURN
        END IF
          
        u(i) = 1.5d0  * v2(i)
        maxvel = MAX(DSQRT(v2(i)),maxvel)

        Kn = (1.0d0/sigma_m) / DSQRT(P(i))     
        minKn = MIN(Kn,minKn)
      END DO
             
      RETURN
      END
      
c***********************************************************************

      SUBROUTINE revirialize2
c-----------------------------------------------------------------------
c  Solve tridiagonal set of equations to re-establish hydrodynamic
c  Iteratively re-establish hydrostatic equilibrium until convergence
c  criterion is met. Each step we solve a tridiagonal set of equations.
c
c  My new method
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER i
      REAL*8  dP,dr,drho,dd,r3a,r3b,r3c,r3d,q1,q2,c1,c2,dV_over_V,Kn,MM
      REAL*8  avec(Ngrid-1),bvec(Ngrid-1),cvec(Ngrid-1)
      REAL*8  Xvec(Ngrid-1),Yvec(Ngrid-1)
          
      REAL*8  Mbaryon
            
c---
     
      Failed =.FALSE.

c---determine elements of tridiagonal matrix

      DO i=1,Ngrid-1
        dP   = P(i+1)  - P(i)
        drho = rho(i+1)+ rho(i)
        dr   = r(i+1)  - r(i-1)
         
        drho = MAX(1.0D-20,drho)
        IF (dP.EQ.0.0d0)   CALL Terminate('dP zero') 
        IF (drho.EQ.0.0d0) CALL Terminate('drho zero') 

        r3a = r(i+1)**3 / (r(i+1)**3 - r(i)**3)
        r3c = r(i)**3   / (r(i)**3 - r(i-1)**3)
        r3b = r3a - 1.0d0
        r3d = r3c - 1.0d0
        
        q1 = r(i+1) / dr
        q2 = q1 - 1.0d0

        MM = M(i) + Mbaryon(r(i)) 
        dd = -((4.0d0/MM) * (r(i)**2/dr) * (dP/drho))
        
        c1 = 5.0d0 * dd * (P(i+1)/dP) - 3.0d0 * (rho(i+1)/drho)
        c2 = 5.0d0 * dd * (P(i)/dP)   + 3.0d0 * (rho(i)/drho)

        Yvec(i) = dd - 1.0d0
        
        avec(i) = r3d * c2 - q2
        bvec(i) =-2.0d0 - r3b * c1 - r3c * c2
        cvec(i) = r3a * c1 + q1
      END DO
        
c---solve tridiagonal equation

      CALL tridag(avec,bvec,cvec,Yvec,Xvec,Ngrid-1)

c---update radii, and corresponding properties

      drmax = 0.0d0
      maxvel = 0.0d0
      minKn = 9.9d+19
      DO i=1,Ngrid
        r3c = r(i)**3   / (r(i)**3 - r(i-1)**3)

        IF (i.NE.Ngrid) THEN
          IF (DABS(Xvec(i)).GT.drmax) THEN
            drmax = DABS(Xvec(i))
            idrmax= i
          END IF  
          r(i) = r(i) * (1.0d0 + Xvec(i))
        END IF
          
        IF (i.EQ.1) THEN
          dV_over_V = 3.0 * r3c * Xvec(i)
        ELSE
          IF (i.EQ.Ngrid) THEN
            dV_over_V = -3.0 * (r3c-1.0d0) * Xvec(i-1)
          ELSE
            dV_over_V = 3.0 * (r3c*Xvec(i) - (r3c-1.0d0)*Xvec(i-1))
          END IF  
        END IF        

        P(i) = P(i) * (1.0d0 - (5.0d0/3.0d0) * dV_over_V)
        rho(i) = rho(i) * (1.0d0 - dV_over_V)
        v2(i) = P(i)/rho(i)

        IF (v2(i).LT.0.0d0) THEN
          WRITE(*,*)' Failed ',i,t,dt,P(i),rho(i),dumax
          Failed=.TRUE.
          RETURN
        END IF
          
        u(i) = 1.5d0  * v2(i)
        maxvel = MAX(DSQRT(v2(i)),maxvel)

        Kn = (1.0d0/sigma_m) / DSQRT(P(i))     
        minKn = MIN(Kn,minKn)
      END DO
             
      RETURN
      END
      
c***********************************************************************

      SUBROUTINE revirialize3
c-----------------------------------------------------------------------
c  Solve tridiagonal set of equations to re-establish hydrodynamic
c  Iteratively re-establish hydrostatic equilibrium until convergence
c  criterion is met. Each step we solve a tridiagonal set of equations.
c
c  Method of Zhong, Yang & Yu (eq. A5).
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER i
      REAL*8  dP,dr,drho,dd,r3a,r3b,r3c,r3d,q1,q2,c1,c2,dV_over_V,Kn,MM
      REAL*8  avec(Ngrid-1),bvec(Ngrid-1),cvec(Ngrid-1)
      REAL*8  Xvec(Ngrid-1),Yvec(Ngrid-1)
          
      REAL*8  Mbaryon
            
c---
     
      Failed =.FALSE.

c---determine elements of tridiagonal matrix

      DO i=1,Ngrid-1
        dP   = P(i+1)  - P(i)
        drho = rho(i+1)+ rho(i)
        dr   = r(i+1)  - r(i-1)
         
        drho = MAX(1.0D-20,drho)
        IF (dP.EQ.0.0d0)   CALL Terminate('dP zero') 
        IF (drho.EQ.0.0d0) CALL Terminate('drho zero') 

        r3a = r(i+1)**3 / (r(i+1)**3 - r(i)**3)
        r3c = r(i)**3   / (r(i)**3 - r(i-1)**3)
        r3b = r3a - 1.0d0
        r3d = r3c - 1.0d0
        
        q1 = r(i+1) / dr
        q2 = q1 - 1.0d0

        MM = M(i) + Mbaryon(r(i)) 
        dd = -((4.0d0/MM) * (r(i)**2/dr) * (dP/drho))
        
        c1 = 5.0d0 * dd * (P(i+1)/dP) - 3.0d0 * (rho(i+1)/drho)
        c2 = 5.0d0 * dd * (P(i)/dP)   + 3.0d0 * (rho(i)/drho)

        Yvec(i) = dd - 1.0d0
        
        avec(i) = r3d * c2 - q2
        bvec(i) =-2.0d0 * dd - r3b * c1 - r3c * c2
        cvec(i) = r3a * c1 + q1
      END DO
        
c---solve tridiagonal equation

      CALL tridag(avec,bvec,cvec,Yvec,Xvec,Ngrid-1)

c---update radii, and corresponding properties

      drmax = 0.0d0
      maxvel = 0.0d0
      minKn = 9.9d+19
      DO i=1,Ngrid
        r3c = r(i)**3   / (r(i)**3 - r(i-1)**3)

        IF (i.NE.Ngrid) THEN
          IF (DABS(Xvec(i)).GT.drmax) THEN
            drmax = DABS(Xvec(i))
            idrmax= i
          END IF  
          r(i) = r(i) * (1.0d0 + Xvec(i))
        END IF
          
        IF (i.EQ.1) THEN
          dV_over_V = 3.0 * r3c * Xvec(i)
        ELSE
          IF (i.EQ.Ngrid) THEN
            dV_over_V = -3.0 * (r3c-1.0d0) * Xvec(i-1)
          ELSE
            dV_over_V = 3.0 * (r3c*Xvec(i) - (r3c-1.0d0)*Xvec(i-1))
          END IF  
        END IF        

        P(i) = P(i) * (1.0d0 - (5.0d0/3.0d0) * dV_over_V)
        rho(i) = rho(i) * (1.0d0 - dV_over_V)
        v2(i) = P(i)/rho(i)

        IF (v2(i).LT.0.0d0) THEN
          WRITE(*,*)' Failed ',i,t,dt,P(i),rho(i),dumax
          Failed=.TRUE.
          RETURN
        END IF
          
        u(i) = 1.5d0  * v2(i)
        maxvel = MAX(DSQRT(v2(i)),maxvel)

        Kn = (1.0d0/sigma_m) / DSQRT(P(i))     
        minKn = MIN(Kn,minKn)
      END DO
             
      RETURN
      END
      
c***********************************************************************

      SUBROUTINE get_properties
c-----------------------------------------------------------------------
c  Compute various instantaneous properties from density profile
c-----------------------------------------------------------------------

c---compute differently defined core radii
c   1: radius where density is half the `central' density
c   2: radius where dlogrho/rlogr = -2
c   3: core radius as defined by BSI02

      CALL get_core_radius1
      CALL get_core_radius2
      CALL get_core_radius3

c---check to see if there is SMFP core (assumed to form BH)

      CALL get_BH_radius

c---compute density at characteristic radius

      CALL get_char_densities
      
      RETURN
      END
      
c***********************************************************************

      SUBROUTINE get_core_radius1
c-----------------------------------------------------------------------
c  Compute radius at which density falls to half the central density.
c  We shall refer to this as the core radius.
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   rhoc,rmed1,rmed2,rho1,rho2

c---

      rc1 = 0.0d0   
      rhoc = rho(1)/2.0d0
      
      DO i=1,Ngrid
        IF (rho(i).LT.rhoc) THEN
          rmed1 = (r(i-2) + r(i-1)) / 2.0d0
          rmed2 = (r(i-1) + r(i)) / 2.0d0
          rho1 = rho(i-1)
          rho2 = rho(i)
          rc1 = rmed1 + (rmed2-rmed1)*((rhoc-rho1)/(rho2-rho1))
          RETURN
        END IF  
      END DO

      rc1 = r(Ngrid)

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE get_core_radius2
c-----------------------------------------------------------------------
c  Compute radius at which logarithmic slope of the density profile
c  is equal to -2.0.
c  We shall refer to this as the core radius.
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   rho1,rho2,rmed1,rmed2,slope1,slope2

c---

      rc2 = 0.0d0   

      slope1 = 0.0d0      
      DO i=1,Ngrid
        rho1 = rho(i)
        rho2 = rho(i+1)
        rmed1 = (r(i) + r(i-1)) / 2.0d0
        rmed2 = (r(i+1) + r(i)) / 2.0d0
        slope2 = (DLOG10(rho2) - DLOG10(rho1)) / 
     &           (DLOG10(rmed2) - DLOG10(rmed1))
        IF (slope2.LT.-2.0d0.AND.rc2.EQ.0.0d0) THEN
          rc2 = rmed1 + (rmed2-rmed1)*((-2.0-slope1)/(slope2-slope1))
          RETURN  
        END IF  
        slope1 = slope2
      END DO

      rc2 = r(Ngrid)

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE get_core_radius3
c-----------------------------------------------------------------------
c  Following BSI02, we here the core radius is defined as the radius r
c  that is the root of
c              r = 3 SQRT(v^2(r) / 4 pi G rho(r))
c  In our dimensionless variables, this becomes
c              r = 3 SQRT(v^2/rho(r))
c  We solve this using interpolation
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   rmed1,rmed2,f1,f2

c---

      rc3 = 0.0d0   

      DO i=1,Ngrid
        rmed1 = (r(i) + r(i-1)) / 2.0d0
        rmed2 = (r(i+1) + r(i)) / 2.0d0
        f1 = rmed1 - 3.0d0 * DSQRT(v2(i)/rho(i))
        f2 = rmed2 - 3.0d0 * DSQRT(v2(i+1)/rho(i+1))
        IF ((f1*f2).LT.0.0d0) THEN
          rc3 = rmed1 + ((0.d0-f1)/(f2-f1)) * (rmed2 - rmed1)
          RETURN
        END IF        
      END DO

      rc3 = r(Ngrid)

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE get_BH_radius
c-----------------------------------------------------------------------
c  Compute radius at which dlogM(<r)/dlogr < 0.05
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   xlgM1,xlgM2,xlgr1,xlgr2,slope

c---

      DO i=2,Ngrid
        xlgM1 = DLOG10(M(i-1))
        xlgM2 = DLOG10(M(i))
        xlgr1 = DLOG10(r(i-1))
        xlgr2 = DLOG10(r(i))
        slope = (xlgM2 - xlgM1) / (xlgr2 - xlgr1)
        IF (slope.LT.0.05) THEN
          r_BH = (r(i) + r(i-1)) / 2.0d0
          M_BH = M(i)
          RETURN
        END IF            
      END DO

c---in case you make it to here, no BH radius was found

      r_BH = 0.0d0
      M_BH = 0.0d0

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE get_char_densities
c-----------------------------------------------------------------------
c  Compute density at a characteristic radius
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   rchar1,rchar2,rchar3,rchar4,rchar5,rchar6,rchar7,rchar8
      PARAMETER (rchar1=-1.25d0)
      PARAMETER (rchar2=-1.00d0)
      PARAMETER (rchar3=-0.75d0)
      PARAMETER (rchar4=-0.50d0)
      PARAMETER (rchar5=-0.25d0)
      PARAMETER (rchar6=+0.00d0)
      PARAMETER (rchar7=+0.25d0)
      PARAMETER (rchar8=+0.50d0)
      
      INTEGER  i
      REAL*8   rmed1,rmed2,rho1,rho2

c---

      rhochar1 = 0.0d0
      rhochar2 = 0.0d0
      rhochar3 = 0.0d0
      rhochar4 = 0.0d0
      rhochar5 = 0.0d0
      rhochar6 = 0.0d0
      rhochar7 = 0.0d0
      rhochar8 = 0.0d0
            
      DO i=1,Ngrid

        rmed2 = DLOG10((r(i-1) + r(i)) / 2.0d0)

        IF (rmed2.GT.rchar1.AND.rhochar1.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar1 = rho1 + (rchar1-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF
        
        IF (rmed2.GT.rchar2.AND.rhochar2.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar2 = rho1 + (rchar2-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF

        IF (rmed2.GT.rchar3.AND.rhochar3.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar3 = rho1 + (rchar3-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF

        IF (rmed2.GT.rchar4.AND.rhochar4.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar4 = rho1 + (rchar4-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF

        IF (rmed2.GT.rchar5.AND.rhochar5.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar5 = rho1 + (rchar5-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF
       
        IF (rmed2.GT.rchar6.AND.rhochar6.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar6 = rho1 + (rchar6-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF
       
        IF (rmed2.GT.rchar7.AND.rhochar7.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar7 = rho1 + (rchar7-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF

        IF (rmed2.GT.rchar8.AND.rhochar8.LT.1.0E-12) THEN
          rmed1 = DLOG10((r(i-2) + r(i-1)) / 2.0d0)
          rho1 = rho(i-1)
          rho2 = rho(i)
          rhochar8 = rho1 + (rchar8-rmed1)/(rmed2-rmed1) * (rho2-rho1)
        END IF

      END DO

      RETURN
      END
      
c***********************************************************************

      SUBROUTINE set_time_step
c-----------------------------------------------------------------------
c  Set time step to be used for integration step
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   trel,dt1,dt2
            
c---

c---compute the minimum relaxation time (in units of t0)

      trelmin = 9.9D+19
      DO i=1,Ngrid
        trel = 1.0d0 / (DSQRT(v2(i)) * rho(i))
        trelmin = MIN(trelmin,trel)
      END DO
       
c---make sure timesteps are very short during the phase in which the
c   cusp transforms into a core. We so do by exponentially decaying
c   the timestep control parameter.

      dt1 = eps_dt * trelmin
      dt2 = dt * 0.95d0 * (eps_du/dumax)
 
      dt = MIN(dt1,dt2)

      END
       
c***********************************************************************

      SUBROUTINE write_output(n)
c-----------------------------------------------------------------------
c  Output radial profiles
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i,n
      REAL*8   rmed,trel,Kn
            
      CHARACTER*4  ext
      CHARACTER*6  fmt  
      CHARACTER*60 outfile*60

      INTEGER  lblnk
      EXTERNAL lblnk

c---

      IF (n.LT.10) THEN
        fmt ='(I1)'
        WRITE(ext,fmt)n
      ELSE
        IF (n.LT.100) THEN
          fmt ='(I2)'
          WRITE(ext,fmt)n
        ELSE
          IF (n.LT.1000) THEN
            fmt ='(I3)'
            WRITE(ext,fmt)n
          ELSE
            IF (n.LT.10000) THEN
              fmt ='(I4)'
              WRITE(ext,fmt)n
            ELSE
              fmt ='(I5)'
              WRITE(ext,fmt)n
            END IF
          END IF
        END IF
      END IF

      outfile='/timestep_'//trim(ext)//'.dat'
      OPEN(11,file=moddir(1:lblnk(moddir))//outfile(1:lblnk(outfile)),
     &     status='UNKNOWN')

c---compute the minimum relaxation time (in units of t0)

      DO i=1,Ngrid
        rmed = (r(i) + r(i-1)) / 2.0d0
        trel = 1.0d0 / (DSQRT(v2(i)) * rho(i))
        Kn = (1.0d0/sigma_m) / DSQRT(P(i))
        WRITE(11,*)i,DLOG10(r(i)),DLOG10(rmed),M(i),rho(i),v2(i),trel,Kn
      END DO
      
      CLOSE(11)
      
      END
       
c***********************************************************************

      REAL*8 FUNCTION Mbaryon(rr)
c-----------------------------------------------------------------------
c  The enclosed baryonic mass within radius r
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   rr
      
c---

      Mbaryon = zeta_b * rr**0.6

      END
      
c***********************************************************************

      REAL*8 FUNCTION fNFW(x)
c-----------------------------------------------------------------------
c  The function f(x) used for computing enclosed mass of NFW profiles
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   x
      
c---

      fNFW = DLOG(1.0d0+x) - x/(1.0d0+x)

      END
      
c***********************************************************************

      REAL*8 FUNCTION gp(rr)
c-----------------------------------------------------------------------
c  The function g_p(r) used for computing the enclosed mass of a NFW
c  profile truncated with power-law (r/rt)^{-p}
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  ierr
      REAL*8   rr,rbuf,SS,SSerr
      
      REAL*8   toint2
      EXTERNAL toint2

c---

      IF (rr.LT.rtrunc) CALL Terminate(' r < r_t: not allowed')

      rbuf = rr
      CALL dqags(toint2,rtr,rbuf,1.0D-6,1.0D-4,SS,SSerr,
     &           Neval2,ierr,Nlimit,Nlenw,last2,iwork2,work2)

      gp = SS

      END
      
c***********************************************************************

      REAL*8 FUNCTION sigr_NFW(rr)
C-----------------------------------------------------------------------
C  Compute the radial velocity dispersion profile for an isotropic halo.
C  The input radius rr must be dimensionless (in units of r_s).
C  The method depends on the value of imode, which selects the halo model.
C-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  j,ierr,ierrb
      REAL*8   rr, rbuf, SS, SSerr, vesc, vmax, Pmax, sum
      REAL*8   v, ff, density, PT, SSb, SSberr
      
      REAL*8  psi
      COMMON /pot/psi

      REAL     ran3
      EXTERNAL ran3
      
      REAL*8   rhotilde2,rhotilde4,fevsq,potential,getQ,brent
      REAL*8   toint1,toint3,toint6,toint7,toint1b, toint1d
      EXTERNAL rhotilde2,rhotilde4,fevsq,potential,getQ,brent
      EXTERNAL toint1,toint3,toint6,toint7,toint1b, toint1d
      
C--- Store input radius

      rbuf = rr

C--- imode = 1 or 5: Untruncated NFW (standard or with central flat-top)
C      Use analytical integrands to compute dispersion integral
      IF (imode.EQ.1.OR.imode.EQ.5) THEN
        CALL dqagi(toint1,rbuf,1,1.0D-6,1.0D-4,SS,SSerr,Neval,
     &          ierr,Nlimit,Nlenw,last,iwork,work)
        CALL dqagi(toint1b,rbuf,1,1.0D-6,1.0D-4,SSb,SSberr,Neval,
     &          ierrb,Nlimit,Nlenw,last,iwork,work)
        sigr_NFW = rr * (1.0d0+rr)**2 * (SS + zeta_b * SSb)
      END IF
     
C--- imode = 2: Power-law truncation
C      Compute the integral of the distribution function from r to infinity
      IF (imode.EQ.2) THEN
        IF (rbuf.GT.rtr) THEN
          CALL dqagi(toint3,rbuf,1,1.0D-6,1.0D-4,SS,SSerr,Neval,
     &          ierr,Nlimit,Nlenw,last,iwork,work)
        ELSE 
          CALL dqags(toint3,rbuf,rtr,1.0D-6,1.0D-4,SS,SSerr,Neval,
     &          ierr,Nlimit,Nlenw,last,iwork,work)
          SS = SS + ftr
        END IF
        sigr_NFW = (1.0d0/rhotilde2(rr)) * SS
      END IF
      
C--- imode = 3: Energy-truncated NFW profile
C      Integrate to find dispersion, use DF-based density in outer region
      IF (imode.EQ.3) THEN
        IF (rr.GT.rcut) THEN
          sigr_NFW = 0.0d0
        ELSE
          CALL dqagi(toint6,rbuf,1,1.0D-6,1.0D-4,SS,SSerr,Neval,
     &          ierr,Nlimit,Nlenw,last,iwork,work)
          IF (rr.LT.rad(1)) THEN
            density = 1.0d0 / (rr * (1.0d0+rr)**2)    
          ELSE
            CALL locate(rad,Nmax,rr,j)
            PT = pot(j)+(rr-rad(j))/(rad(j+1)-rad(j))*(pot(j+1)-pot(j))
            density = getQ(PT) 
          END IF
          sigr_NFW = (1.0d0/density) * SS
        END IF
      END IF
       
C--- imode = 4: Exponential truncation
C      Use exponential cutoff profile to compute dispersion
      IF (imode.EQ.4) THEN
        CALL dqagi(toint7,rbuf,1,1.0D-6,1.0D-4,SS,SSerr,Neval,
     &          ierr,Nlimit,Nlenw,last,iwork,work)
        sigr_NFW = (1.0d0/rhotilde4(rr)) * SS
      END IF

C--- imode = 6: Alpha-Beta-Gamma profile
C      Use custom integrand to compute vdisp for ABG model
      IF (imode .EQ. 6) THEN
        CALL dqagi(toint1d, rbuf, 1, 1.0D-6, 1.0D-4, SS, SSerr, Neval,
     &            ierr, Nlimit, Nlenw, last, iwork, work)
        density = 1.0D0 / (rr**gamma * (1.0D0 + rr**(1.0D0 / alpha))
     &            **((beta - gamma) * alpha))
        sigr_NFW = (1.0D0 / density) * SS
      END IF
       
      END
      
c**********************************************************************

      REAL*8 FUNCTION Menc(rr)
c-----------------------------------------------------------------------
c  Compute the halo mass enclosed within radius rr (in units of r_s).
c  The result is returned in units of Mvir.
c  Different modes correspond to different profile models.
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  ierr
      REAL*8   rr,SS,SSerr,f1,f2
      
      REAL*8    fNFW,gp,toint4,toint4b,incomplete_gamma
      EXTERNAL  fNFW,gp,toint4,toint4b,incomplete_gamma

c--- imode = 1: Standard NFW profile (analytical)
      IF (imode.EQ.1) THEN
        Menc = fNFW(rr)
      END IF
      
c--- imode = 2: Truncated NFW profile
c      Use standard NFW within rtrunc; match power-law tail beyond
      IF (imode.EQ.2) THEN
        IF (rr.LT.rtrunc) THEN
          Menc = fNFW(rr)
        ELSE
          Menc = Mrt + gp(rr) * rtrunc**ppar
        END IF
      END IF
      
c--- imode = 3: Energy-truncated NFW profile
c      Mass profile from numerical integration of the DF over phase space
      IF (imode.EQ.3) THEN
        CALL dqags(toint4,0.0d0,rr,1.0D-5,1.0D-5,SS,SSerr,Neval2,
     &            ierr,Nlimit,Nlenw,last2,iwork2,work2)
        Menc = SS
      END IF

c--- imode = 4: NFW with exponential cutoff (e.g., Hernquist-like tail)
      IF (imode.EQ.4) THEN
        IF (rr.LE.rtrunc) THEN
          Menc = fNFW(rr)/fc
        ELSE
          f1 = incomplete_gamma(3.0d0+eps,1.0d0/eta)
          f2 = incomplete_gamma(3.0d0+eps,rr/(eta*rtrunc))
          Menc = (fNFW(rtrunc)/fc) + prefac * (f1-f2)
        END IF
        Menc = Menc * fc
      END IF

c--- imode = 5: Flat-topped core profile
c      Core has constant density = rho(1) inside r(1)
      IF (imode.EQ.5) THEN
        IF (rr.LE.r(1)) THEN
          Menc = (rho(1)/3.0d0) * rr**3
        ELSE
          Menc = fNFW(rr) - fNFW(r(1)) + (rho(1)/3.0d0) * r(1)**3
        END IF
      END IF

c--- imode = 6: Alpha-Beta-Gamma profile (numerical integration)
      IF (imode.EQ.6) THEN
        CALL dqags(toint4b, 0.0d0, rr, 1.0D-5, 1.0D-5, SS, SSerr, 
     &            Neval2, ierr, Nlimit, Nlenw, last2, iwork2, work2)
        Menc = SS
      END IF

      END
    
c**********************************************************************

      REAL*8 FUNCTION denss(rr)
c-----------------------------------------------------------------------
c  The normalized density at radius rr
c  NOTE: rr is in units of r_s
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  j
      REAL*8   rr,PT
 
      REAL*8   rhotilde,rhotilde2,rhotilde4,getQ      
      EXTERNAL rhotilde,rhotilde2,rhotilde4,getQ
      
c---

      IF (imode.EQ.1) denss = rhotilde(rr)

      IF (imode.EQ.2) denss = rhotilde2(rr)

      IF (imode.EQ.3) THEN
        IF (rr.LT.rad(1)) THEN
          denss = 1.0d0 / (rr * (1.0d0+rr)**2)    
        ELSE
          CALL locate(rad,Nmax,rr,j)
          PT = pot(j)+(rr-rad(j))/(rad(j+1)-rad(j))*(pot(j+1)-pot(j))
          denss = getQ(PT) 
        END IF
      END IF

      IF (imode.EQ.4) denss = rhotilde4(rr)
       
      IF (imode.EQ.5) denss = rhotilde(rr)

      END
      
c**********************************************************************

      REAL*8 FUNCTION rhotilde(rr)
c-----------------------------------------------------------------------
c  The normalized, flat-topped density at radius rad/r_s
c  NOTE: rad is in units of r_s
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   rr

c---

      IF (rr.LE.r(1)) THEN      
        rhotilde = rho(1)
      ELSE
        rhotilde = 1.0d0 / (rr * (1.0d0+rr)**2)
      END IF
      
      END
      
c**********************************************************************

      REAL*8 FUNCTION rhotilde2(rr)
c-----------------------------------------------------------------------
c  The normalized density at radius rad/r_s
c  NOTE: rad is in units of r_s
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   rr

c---

      IF (rr.GT.rtrunc) THEN      
        rhotilde2 = prefac * ((rr**pow) / ((1.0d0+rr)**2))
      ELSE
        rhotilde2 = 1.0d0 / (rr * (1.0d0+rr)**2)
      END IF
      
      END
      
c**********************************************************************

      REAL*8 FUNCTION rhotilde4(rr)
c-----------------------------------------------------------------------
c  The density at location r (in normalized units), computed from 
c  the initial enclosed mass profile.
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   rr,fac

c---

      IF (rr.LE.rtrunc) THEN
        rhotilde4 = 1.0d0/(rr * (1.0d0+rr)**2)
      ELSE
        fac = (rr/rtrunc)**eps * DEXP(-((rr/rtrunc) - 1.0d0)/eta) 
        rhotilde4 = fac/(rtrunc * (1.0d0+rtrunc)**2)
      END IF

      END
    
c***********************************************************************

      SUBROUTINE set_param
c-----------------------------------------------------------------------
c  Set parameters and convert units as needed
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  i
      REAL*8   x,rr,rtilde,rhotild,rhophys,rho0phys
      CHARACTER*60 outfile*60
      
      INTEGER  lblnk
      EXTERNAL lblnk
      
c---

c---compute characteristic quantities
c     rho_s   characteristic density  [Msun/Mpc^3]
c     v0      characteristic velocity [km/s]
c     t0      characteristic relaxation time 
c     sigma0  characteristic cross section  
c
      rho_s = Ms / (4.0d0 * pi * rs**3)    
      v0 = DSQRT(gee * Ms / rs)  
      sigma0 = 4.0d0 * pi * rs**2 / Ms
      
c---convert quantities to cgs and rho_s also to Msun/pc^3 and t0 to Gyr

      sigma0 = sigma0 * (Mpc_to_cm)**2 / (Msun_to_gram)
      rho_s_cgs = rho_s * (Msun_to_gram) / (Mpc_to_cm)**3
      rho_s_Msunpc3 = rho_s * 1.0E-18
      v0_cgs = v0 * 1.0D+5

      t0     = (1.0d0/(apar * sigma_m * v0_cgs * rho_s_cgs)) 
      t0_Gyr = t0 * sec_to_Gyr       

c---convert cross section per unit mass to dimensionless form

      sigma_m = sigma_m / sigma0

c---convert Mfly to characteristic mass units if flyby is on
      IF (flyby_on) THEN
         IF (imode .EQ. 6) THEN
            Mfly = Mfly * chi
         ELSE
            Mfly = Mfly * fc
         END IF
      END IF

c---write initial density profiles to file

      outfile='/initprof.dat'
      OPEN(88,file=moddir(1:lblnk(moddir))//outfile(1:lblnk(outfile)),
     &     status='UNKNOWN')

      rho0phys = Mvir/(4.0d0 * pi * rs**3 * fc)
      DO i=1,200
        x = -4.0d0 + (DBLE(i-1)/199.d0) * 4.0d0
        rr = 10.0d0**x * rvir
        rtilde = rr/rs
        rhotild = 1.0d0/(rtilde * (1.0d0+rtilde)**2)
        rhophys = rho0phys * rhotild * 1.0E-18
        WRITE(88,*)i,x,rtilde,DLOG10(rhotild),DLOG10(rhophys)
      END DO
      CLOSE(88)
            
      RETURN
      END
      
c***********************************************************************

      SUBROUTINE read_param
c-----------------------------------------------------------------------
c  Read parameters
c-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER  ianswer, ierr
      REAL*8   rmin,rmax,gc, SS, SSerr
       
      REAL*8   fNFW,xH,gp,Menc,df,Delta_crit, chiint
      EXTERNAL fNFW,xH,gp,Menc,df,Delta_crit, chiint
      
c---

c---give directory to which output should be written

      WRITE(*,*)' Give model directory '
      READ(*,'(A)')moddir
      WRITE(*,'(A)')moddir
      WRITE(*,*)' '

c---iseed needed for random number generator

      WRITE(*,*)' Give iseed (0=fiducial) '
      READ(*,*)iseed
      WRITE(*,*)' iseed = ',iseed
      WRITE(*,*)' '

c---iseed=0 corresponds to fiducial

      IF (iseed.EQ.0) iseed=-54991

c---give cosmology dependent quantities

      WRITE(*,*)' Give Omega_{m,0} '
      READ(*,*)omega_0
      WRITE(*,*)'   Omega_0 = ',omega_0
      WRITE(*,*)' '

      omega_lambda = 1.0d0 - omega_0
      
      WRITE(*,*)' Give hubble parameter (h=H_0/100 km/s/Mpc)'
      READ(*,*)xhubble
      WRITE(*,*)'         h = ',xhubble
      WRITE(*,*)' '

      xH_0 = 100.0d0 * xhubble

      WRITE(*,*)' Give Delta_vir (overdensity wrt critical) '
      READ(*,*)Delta_vir
      WRITE(*,*)'  Delta_vir = ',Delta_vir
      WRITE(*,*)' '

c---read halo properties

      WRITE(*,*)' Give halo virial mass (in Msun/h)'
      READ(*,*)Mvir
      WRITE(*,*)'  Mvir = ',Mvir
      WRITE(*,*)' '

      WRITE(*,*)' Give halo concentration '
      READ(*,*)cvir
      WRITE(*,*)'  cvir = ',cvir
      WRITE(*,*)' '

      WRITE(*,*)' Give redshift of halo in question'
      READ(*,*)zvir
      WRITE(*,*)'  zvir = ',zvir
      WRITE(*,*)' '

      WRITE(*,*)' Give baryonic mass fraction'
      READ(*,*)f_b
      WRITE(*,*)'  f_b = ',f_b
      WRITE(*,*)' '

c---specify what mode to use
c      1 = pure, untruncated NFW halo
c      2 = power-law truncation as in Nishikawa.etal.   
c      3 = energy truncation (as in Nbody simulations)
c      4 = exponential truncation (following Kazantzidis+03)

      WRITE(*,*) ' Give mode (1=NFW, 2=PL-trunc, 3=E-trunc,'//
     &          ' 4=exp-trunc, 6=abg)'
      READ(*,*)imode
      
      IF (imode.NE.1.AND.imode.NE.2.and.imode.NE.3.and.
     &    imode.NE.4.AND.imode.NE.5.AND.imode.NE.6) THEN
        CALL Terminate('invalid mode')
      END IF
      
      IF (f_b.GT.0.d0.AND.(imode.NE.1.AND.imode.NE.5)) THEN
        CALL Terminate(' Mode not supported with baryons')
      END IF

      WRITE(*,*)'  imode = ',imode
      WRITE(*,*)' '
      
c---Parameters used for PL and EXP-truncation (imode=2 or 4)

      WRITE(*,*)' Give truncation radius in units of rs (imode=2-4)'
      READ(*,*)rtrunc
      WRITE(*,*)'  r_trunc = ',rtrunc
      WRITE(*,*)' '
      IF (imode.EQ.1) rtrunc=0.0d0
      
c---Parameters used for PL-truncation (imode=2)

      WRITE(*,*)' Give truncation power law (imode=2)'
      READ(*,*)ppar
      WRITE(*,*)'  ppar = ',ppar
      WRITE(*,*)' '
      IF (imode.NE.2) ppar = 0.0d0
      
c---Parameters used for E-truncation (imode=3)

      WRITE(*,*)' Give truncation energy, Zt (imode=3)'
      READ(*,*)Zt
      WRITE(*,*)' Z_t = ',Zt
      WRITE(*,*)' '
      IF (imode.NE.3) Zt=0.0d0       
              
      WRITE(*,*)' Give deltaP (step size in P, fiducial=1.0D-5): '
      READ(*,*)deltaP
      deltaP = -DABS(deltaP)
      WRITE(*,*)' deltaP = ',deltaP
      WRITE(*,*)' '

c---Parameters used for Exp-truncation (imode=4)
      
      WRITE(*,*)' Give ratio rdecay/rtrunc (imode=4) '
      READ(*,*)eta
      WRITE(*,*)'  eta = ',eta
      WRITE(*,*)' '
      IF (imode.NE.4) eta=0.0d0

c---Parameters used to abg profile (imode=6)
      
      WRITE(*,*)' Give alpha (imode=6) '
      READ(*,*)alpha
      WRITE(*,*)'  alpha = ',alpha
      WRITE(*,*)' '
      IF (imode.NE.6) alpha=0.0d0

      WRITE(*,*)' Give beta (imode=6) '
      READ(*,*)beta
      WRITE(*,*)'  beta = ',beta
      WRITE(*,*)' '
      IF (imode.NE.6) beta=0.0d0

      WRITE(*,*)' Give gamma (imode=6) '
      READ(*,*)gamma
      WRITE(*,*)'  gamma = ',gamma
      WRITE(*,*)' '
      IF (imode.NE.6) gamma=0.0d0

c---compute chi parameter for abg profile
      If (imode.EQ.6) THEN
        CALL dqags(chiint, 0.0d0, 1.0d3, 1.0D-6, 1.0D-6, SS, SSerr,
     &             Neval2,ierr,Nlimit,Nlenw,last2,iwork2,work2)

        chi = SS
      ELSE
        chi = 0.0d0
      END IF

c---compute virial radius and scale radius in Mpc, Vvir in km/s
      
      Rvir = 0.169d0 * (Mvir/1.0D+12)**(1.d0/3.d0)
      Rvir = Rvir * (Delta_vir/178.d0)**(-1.d0/3.d0)
      Rvir = Rvir * (xH(zvir)/xH_0)**(-2.0/3.0)

      rvir = rvir / xhubble
      Mvir = Mvir / xhubble

      rs = rvir/cvir
      Vvir = DSQRT(gee * Mvir/rvir) 
      fc = fNFW(cvir)
      gc = gp(cvir) 

      If (imode.EQ.6) THEN
        Ms = Mvir / chi
      ELSE 
        Ms = Mvir / fc
      END IF

      zeta_b = (f_b/(1.0d0-f_b)) * fc/(cvir**(0.6d0))

c---mode-specific parameter definitions

      IF (imode.EQ.2) THEN
        Mrt = fNFW(rtrunc)
        rtr = SNGL(rtrunc)
        power = 1.0 - SNGL(ppar)
        pow = -ppar - 1.0d0
        prefac = rtrunc**ppar
      END IF
          
      IF (imode.EQ.3) THEN
        Ft = df(Zt)
      END IF
            
      IF (imode.EQ.4) THEN
        eps = (-1.0d0 - 3.0d0*rtrunc)/(1.0d0+rtrunc) + (1.0d0/eta)
        prefac = (rtrunc**2 / (fc*(1.0d0+rtrunc)**2)) * eta**(3.0d0+eps) 
        prefac = prefac * DEXP(1.0d0/eta)
      END IF

c---flyby parameters for impulsive heating

      WRITE(*,*)' flyby_on? '
      READ(*,*)flyby_on
      WRITE(*,*)'  flyby_on = ',flyby_on
      WRITE(*,*)' '

      WRITE(*,*)' tfly (in units of t0) '
      READ(*,*)tfly
      WRITE(*,*)'  tfly = ',tfly
      WRITE(*,*)' '
      IF (.NOT. flyby_on) tfly=0.0d0

      WRITE(*,*)' bfly (in units of rs) '
      READ(*,*)bfly
      WRITE(*,*)'  bfly = ',bfly
      WRITE(*,*)' '
      IF (.NOT. flyby_on) bfly=0.0d0

      WRITE(*,*)' Mfly (as fraction of Mtot or Mvir) '
      READ(*,*)Mfly
      WRITE(*,*)'  Mfly = ',Mfly
      WRITE(*,*)' '
      IF (.NOT. flyby_on) Mfly=0.0d0

      WRITE(*,*)' vfly (in units of v0) '
      READ(*,*)vfly
      WRITE(*,*)'  vfly = ',vfly
      WRITE(*,*)' '
      IF (.NOT. flyby_on) vfly=0.0d0
      
c----minimum and maximum radius of radial grid
      
      WRITE(*,*)' Give minimum normalized radius (r/rs)'
      READ(*,*)rmin
      WRITE(*,*)' rmin = ',rmin
      WRITE(*,*)' '   
            
      WRITE(*,*)' Give maximum normalized radius (r/rs)'
      READ(*,*)rmax
      WRITE(*,*)' rmax = ',rmax
      WRITE(*,*)' '      

      WRITE(*,*)' Give Nrad (number of radial grid cells)'
      READ(*,*)Ngrid
      WRITE(*,*)' Ngrid = ',Ngrid
      WRITE(*,*)' '      

      IF (Ngrid.GT.Ngrid_max) CALL Terminate('Ngrid too large')

      xlgrmin = DLOG10(rmin)
      xlgrmax = DLOG10(rmax)

c----read collisional cross section 

      WRITE(*,*)' Give cross section per unit mass [cm^2/g]'
      READ(*,*)sigma_m
      WRITE(*,*)'  sigma_m = ',sigma_m
      WRITE(*,*)' '

c---read characteristic parameters that specify collisions

      WRITE(*,*)' Give parameter a '
      READ(*,*)apar
      WRITE(*,*)'  a = ',apar
      WRITE(*,*)' '
      
      WRITE(*,*)' Give parameter b '
      READ(*,*)bpar
      WRITE(*,*)'  b = ',bpar
      WRITE(*,*)' '
      
      WRITE(*,*)' Give parameter c '
      READ(*,*)cpar
      WRITE(*,*)'  c = ',cpar
      WRITE(*,*)' '

c---read evaporation rate
            
      WRITE(*,*)' Give evaporation rate'
      READ(*,*)Gamma_evap
      WRITE(*,*)'  Gamma_evap = ',Gamma_evap
      WRITE(*,*)' '
                  
c---read characteristic parameters that specify collisions

      WRITE(*,*)' Give time step epsilon parameter '
      READ(*,*)eps_dt
      WRITE(*,*)'  eps_dt = ',eps_dt
      WRITE(*,*)' '

      WRITE(*,*)' Give convergence parameter for re-virialization'
      READ(*,*)eps_dr
      WRITE(*,*)'  eps_dr = ',eps_dr
      WRITE(*,*)' '

      WRITE(*,*)' Give convergence parameter for convection'
      READ(*,*)eps_du
      WRITE(*,*)'  eps_du = ',eps_du
      WRITE(*,*)' '

      WRITE(*,*)' Give t at which to stop integration (in units of t0)'
      READ(*,*)tstop
      WRITE(*,*)'  tstop = ',tstop
      WRITE(*,*)' '

c---pick method of revirialization
            
      WRITE(*,*)' Give revirialization method (1=Old, 2=New, 3=ZYY)'
      READ(*,*)iRevir
      
      IF (iRevir.NE.1.AND.iRevir.NE.2.AND.iRevir.NE.3) THEN
        CALL Terminate(' Invalid ')
      END IF
      
      RETURN
      END

c**********************************************************************

      SUBROUTINE integrate_potential

      INCLUDE 'paramfile.h'

      INTEGER  i,Nok,Nbad,N
      REAL*8   r0,r1,r2,dr,rr,step
      REAL*8   yy(2),yold,ynew,Mtot
             
      REAL*8   potential,getQ,Menc
      EXTERNAL potential,getQ,Menc

c---

      WRITE(*,*)' ...Computing potential...'
      WRITE(*,*)' '

      i = 1
      r0 = 0.01
      
c---compute potential P and logarithmic derivative dP/dlogr at rmin

      r1 = (1.0d0-epsilon) * r0
      r2 = (1.0d0+epsilon) * r0

      yy(1) = potential(r0)
      yy(2) =(potential(r2) - potential(r1)) / (r2-r1)
      
c---initialize starting step in radius

      r1 = r0
      dr = deltaP / yy(2)
      r2 = r1 + dr
      
c---store initial radius and density (normalized by rho0)

      pot(i) = yy(1) 
      rad(i) = r1
       
      Nok = 0
      Nbad= 0
      DO WHILE (yy(1).GT.0.0)
        step = (r2-r1)/FLOAT(Nstep)

        yold = yy(1)
        CALL odeint(yy,2,r1,r2,1.0D-5,step,0.0d0,Nok,Nbad)
        ynew = yy(1)

        i = i + 1
        rad(i) = r2
        pot(i) = ynew
         
        r1 = r2
        dr = deltaP/yy(2)
        IF (dr.LT.0.0) CALL Terminate(' dr negative ')

        dr = MIN(dr,0.01)
        r2 = r1 + dr
        rcut = r1
        N = i
      END DO 

c---set size of vector with potential and set boundary condition
       
      Nmax = N
      IF (Nmax.GT.Nmx) CALL Terminate(' Nmax too large ')

      pot(N) = 0.0d0

c---compute total mass of truncated profile and particle mass

      Mtot = Menc(rcut)
      
c---modify the maximum grid radius so not to exceed rcut

      xlgrmax = MIN(xlgrmax,DLOG10(0.99d0*rcut))
   
      WRITE(*,*)'                rcut : ',rcut
      WRITE(*,*)'           log[rmax] : ',xlgrmax
      WRITE(*,*)'      Mtot = M(rmax) : ',Mtot
      WRITE(*,*)'             M(rvir) : ',Menc(cvir)
      WRITE(*,*)' '
      
      RETURN
      END
      
c***********************************************************************

      REAL*8 FUNCTION fevsq(vv)
      
      INCLUDE 'paramfile.h'

      REAL*8  vv,ee,ff
     
      REAL*8  psi
      COMMON /pot/psi
      
      REAL*8   df_trunc
      EXTERNAL df_trunc
        
c---

      ee = psi - vv**2/2.0d0
      ff = df_trunc(ee)

      fevsq = 1.0d0/(ff * vv**2)
      
      END
      
c***********************************************************************

      REAL*8 FUNCTION df_trunc(e)

      INCLUDE 'paramfile.h'
      
      REAL*8    e
      
      REAL*8   df
      EXTERNAL df      

c---

      df_trunc = df(e+Zt) - Ft

      END
      
c***********************************************************************

      REAL*8 FUNCTION df(e)
      
      INCLUDE 'paramfile.h'
       
      REAL*8  e,e2,pp,fac1,fac2,fac3

c---

      e2 = 1.0d0 - e
      
      pp = p1 * e + p2 * e**2 + p3 * e**3 + p4 * e**4
      
      fac1 = e**(1.5) / e2**(2.5)
      fac2 = (-DLOG(e)/e2)**q
      fac3 = DEXP(pp)
      
      df = F0 * fac1 * fac2 * fac3
      
      END
      
c***********************************************************************

      REAL*8 FUNCTION getQ(PT)
      
      INCLUDE 'paramfile.h'
      
      INTEGER  ierr
      REAL*8   PT,SS,SSerr
      
      REAL*8   PP
      COMMON /Pbuf/ PP

      REAL*8   toint5
      EXTERNAL toint5

c---
      
      PP = PT
      IF (PT.EQ.0.0d0) THEN
        SS = 0.0d0
      ELSE
        CALL dqags(toint5,0.0d0,PT,1.0D-5,1.0D-5,SS,SSerr,Neval3,
     &            ierr,Nlimit,Nlenw,last3,iwork3,work3)
      END IF
       
      getQ = 4.0d0 * pi * SS

      END
            
c***********************************************************************

      REAL*8 FUNCTION potential(rr)
      
      INCLUDE 'paramfile.h'

      REAL*8 rr
      
c---
      
      potential = DLOG(1.0d0+rr)/rr
      potential = potential - Zt 
                  
      END
            
c***********************************************************************

      SUBROUTINE derivs(x,y,dydx)
      
      IMPLICIT NONE
      
      REAL*8    fac 
      PARAMETER (fac=0.434294482d0)
      
      REAL*8    r,x,y(2),dydx(2),y1,y2,Q
      
      REAL*8    getQ
      EXTERNAL  getQ
      
c---
      
      y1 = y(1)
      y2 = y(2)

      dydx(1) = y2
      Q = getQ(y1)

      r = x           
      dydx(2) = -Q - (2.0d0/r) * y2
      
      END
            
c***********************************************************************

      REAL*8 FUNCTION toint1(x)
      
      REAL*8     x,fac

c---

      fac = DLOG(1.0d0+x) - x/(1.0d0+x)
      toint1 = fac / (x**3 * (1.0d0+x)**2)
      
      END
      
c***********************************************************************

      REAL*8 FUNCTION toint1b(x)
      
      REAL*8     x,fac

c---

      toint1b = 1.0d0 / (x**(2.4d0) * (1.0d0+x)**2)
      
      END
      
c**********************************************************************

      REAL*8 FUNCTION toint1c(rr)
      
      INCLUDE 'paramfile.h'

      REAL*8   rr,rbuf
      
      REAL*8   rhotilde,Menc
      EXTERNAL rhotilde,Menc

c---

      rbuf = rr
      toint1c = rhotilde(rbuf) * Menc(rbuf) / rr**2

      END

c**********************************************************************

      REAL*8 FUNCTION toint1d(x)
C-----------------------------------------------------------------------
C  Integrand for radial velocity dispersion calculation using the
C  alpha-beta-gamma density profile:
C    rho(r) = C / [r^gamma * (1 + r^(1/alpha))^((beta - gamma) * alpha)]
C  This function returns rho(r) * M(r) / r^2, the integrand of the
C  Jeans equation under spherical, isotropic assumptions.
C-----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   x, C, density, integrand
      
      REAL*8   Menc
      EXTERNAL Menc

      C = 1.0D0  ! Placeholder normalization constant
      density = C / (x**gamma * (1.0D0 + x**(1.0D0 / alpha))
     &      **((beta - gamma) * alpha))
      integrand = density * Menc(x) / x**2
      toint1d = integrand

      RETURN
      END
      
c**********************************************************************

      REAL*8 FUNCTION toint2(rr)
      
      INCLUDE 'paramfile.h'

      REAL*8   rr

c---

      toint2 = (rr**power) / (1.0d0+rr)**2
      
      END
      
c**********************************************************************

      REAL*8 FUNCTION toint3(rr)
      
      INCLUDE 'paramfile.h'

      REAL*8   rr,rbuf
      
      REAL*8   rhotilde2,Menc
      EXTERNAL rhotilde2,Menc

c---

      rbuf = rr
      toint3 = rhotilde2(rbuf) * Menc(rbuf) / rr**2

      END
      
c***********************************************************************

      REAL*8 FUNCTION toint4(rr)
c-----------------------------------------------------------------------      
c  Compute integrand  rho(r) r^2  used for computing enclosed mass
c  If r<rad(1), we use analytical NFW profile
c  If r>rad(1), we find potential P(r) and compute rho(P)
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      INTEGER j
      REAL*8  rr,PT,density
      
      REAL*8   getQ     
      EXTERNAL getQ
 
c---
              
      IF (rr.LT.rad(1)) THEN
        density = 1.0d0 / (rr * (1.0d0+rr)**2)      
      ELSE
        CALL locate(rad,Nmax,rr,j)
        PT = pot(j) + (rr-rad(j))/(rad(j+1)-rad(j)) * (pot(j+1)-pot(j))
        density = getQ(PT) 
      END IF
      
      toint4 = density * rr**2      

      END

c***********************************************************************

      REAL*8 FUNCTION toint4b(rr)
C-----------------------------------------------------------------------      
C  Compute integrand  rho(r) r^2  used for computing enclosed mass
C  using an alpha-beta-gamma profile:
C  rho(r) = C / [r^gamma * (1 + r^(1/alpha))^((beta - gamma) * alpha)]
C  Parameters alpha, beta, gamma are provided by user input
C-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      REAL*8 rr, density, C

C---
      C = 1.0D0  ! Placeholder normalization constant

      density = C / (rr**gamma * (1.0D0 + rr**(1.0D0/alpha))**
     &          ((beta - gamma) * alpha))
      toint4b = density * rr**2

      RETURN
      END
      
c***********************************************************************

      REAL*8 FUNCTION toint5(Z)
      
      INCLUDE 'paramfile.h'

      REAL*8  Z

      REAL*8  PP
      COMMON /Pbuf/ PP

      REAL*8   df_trunc
      EXTERNAL df_trunc
        
c---
      
      IF (Z.GT.PP) THEN
        toint5 = 0.0d0 
      ELSE
        toint5 = df_trunc(Z) * DSQRT(2.0d0*(PP-Z))     
      END IF
      
      END
      
c***********************************************************************

      REAL*8 FUNCTION toint6(rr)
c-----------------------------------------------------------------------      
c  Compute integrand  M(<r) rho(r) / r^2  used for computing the
c  velocity dispersion as function of radius.
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      INTEGER j
      REAL*8  rr,rbuf,PT,density
      
      REAL*8   getQ,Menc
      EXTERNAL getQ,Menc
 
c---
              
      rbuf = rr

      IF (rr.LT.rad(1)) THEN
        density = 1.0d0 / (rr * (1.0d0+rr)**2)    
      ELSE
        CALL locate(rad,Nmax,rr,j)
        PT = pot(j) + (rr-rad(j))/(rad(j+1)-rad(j)) * (pot(j+1)-pot(j))
        density = getQ(PT) 
      END IF
            
      toint6 = density * Menc(rbuf) / rr**2

      END
      
c***********************************************************************

      REAL*8 FUNCTION toint7(rr)
c-----------------------------------------------------------------------      
c  Compute integrand  M(<r) rho(r) / r^2  used for computing the
c  velocity dispersion as function of radius.
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      INTEGER j
      REAL*8  rr,rbuf
      
      REAL*8   Menc,rhotilde4
      EXTERNAL Menc,rhotilde4
 
c---
              
      rbuf = rr

      toint7 = rhotilde4(rbuf) * Menc(rbuf) / rr**2

      END
      
c***********************************************************************

      REAL*8 FUNCTION toint8(rr)
c-----------------------------------------------------------------------      
c  Compute integrand  sig^2_r(r) rho(r) r^2  used for computing
c  the internal energy in central bin
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      INTEGER j
      REAL*8  rr,rbuf
      
      REAL*8    sigr_NFW,rhotilde
      EXTERNAL  sigr_NFW,rhotilde
 
c---
              
      rbuf = rr
      toint8 = sigr_NFW(rbuf) * rhotilde(rbuf) * rr**2

      END
      
c***********************************************************************

      REAL*8 FUNCTION chiint(x)
c-----------------------------------------------------------------------      
c  Integrand for the chi parameter in the abg profile
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      REAL*8  x, expo
 
c---

      expo = (beta - gamma) / alpha
      chiint = x**(2.d0 - gamma) / (1.d0 + x**alpha)**expo 

      END

c***********************************************************************

      REAL*8 FUNCTION findr(rr)
c-----------------------------------------------------------------------      
c  Compute integrand  sig^2_r(r) rho(r) r^2  used for computing
c  the internal energy in central bin
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      REAL*8  rr,rbuf
      
      REAL*8    Menc
      EXTERNAL  Menc
 
      REAL*8    Menclosed
      COMMON    Menclosed

c---
              
      rbuf = rr
      findr = Menc(rbuf) - Menclosed
      WRITE(*,*)rr, rbuf, Menc(rbuf), Menclosed, findr

      END
      
c**********************************************************************

      REAL*8 FUNCTION incomplete_gamma(a,x)
c----------------------------------------------------------------------
c
c  The functions gives the incomlete gamma function for parameters
c  a and x. We use the routine in Numerical Recipes. If a<0 then
c  we use the recurrence relations, using the fact that a is never
c  allowed to be smaller than -1.
c
c----------------------------------------------------------------------
  
      IMPLICIT NONE

      REAL*8   a,x,SS,a1

      REAL     expint
      REAL*8   incgam
      EXTERNAL expint,incgam
      
c---routine is only defined for x>=0 and a>=-1.0

      IF (a.LE.-1.0d0) CALL Terminate(' ERROR IncGam: a out of bound ')
      IF (x.LT.0.0d0)  CALL Terminate(' ERROR IncGam: x out of bound ')

c---if a=0, use exponential integral, otherwise use recursion and 
c   continued fraction developments
      
      IF (a.EQ.0.0d0) THEN
        SS = DBLE(expint(1,SNGL(x)))
      ELSE
        IF (a.GT.0.0d0) THEN
          SS = incgam(a,x)
        ELSE
          a1 = a+1.0d0 
          SS = (incgam(a1,x) - x**a * DEXP(-x)) / a
        END IF
      END IF

      incomplete_gamma = SS
      
      END

c**********************************************************************
       
      REAL*8 FUNCTION incgam(aa,xx)

      REAL*8 aa,xx,SS

      IF (xx.LT.aa+1.0d0) THEN
        CALL dgser2(SS,aa,xx)
      ELSE
        CALL dgcf2(SS,aa,xx)
      END IF
      
      incgam = SS

      END

c**********************************************************************

      REAL*8 FUNCTION xH(z)
c----------------------------------------------------------------------
c
c Calculate hubble constant (in physical units; km/s/Mpc) at redshift z.
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   z,z1,fac

c---

      z1 = 1.0d0 + z
    
      fac = omega_lambda + (1.0d0 - omega_lambda - omega_0) * z1**2 + 
     &      omega_0 * z1**3

      xH = xH_0 * DSQRT(fac)

      END

c**********************************************************************

      REAL*8 FUNCTION omega(z)
c----------------------------------------------------------------------
c
c Calculate the density parameter omega at redshift z. Note that omega 
c is only the mater contribution.
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8     z,xH
  
      EXTERNAL xH

c---

      omega = omega_0 * (1.0d0 + z)**3 / (xH(z)/xH_0)**2.d0

      END

c**********************************************************************

      REAL*8 FUNCTION Delta_crit(z)
c----------------------------------------------------------------------
c
c Calculate the virial density in terms of critical density of the 
c Universe. We use the fitting formulae of Bryan & Norman, 1998, ApJ,
c 495, 80. These are accurate to better than 1% for 0.1 < omega_0 < 1.0
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL*8   z,x,omega
 
      EXTERNAL omega

c---

      x = omega(z) - 1.0d0
      IF (DABS(1.0d0-omega_0-omega_lambda).LT.1.0D-4) THEN
        Delta_crit = 18.d0 * pi**2 + 82.d0*x - 39.d0 * x**2
        GOTO 3
      END IF
      
      IF (omega_0.LT.1.d0.AND.omega_lambda.EQ.0.d0) THEN
        Delta_crit = 18.d0 * pi**2 + 60.d0*x - 32.d0 * x**2
        GOTO 3
      END IF

      WRITE(*,*)' Delta_crit not defined for this cosmology'
      STOP

 3    CONTINUE

      END

c**********************************************************************

      SUBROUTINE OUTPUT_PARAMS
c-----------------------------------------------------------------------      
c  Write important parameters to screen & logfile
c-----------------------------------------------------------------------      

      INCLUDE 'paramfile.h'

      REAL*8   rmx,Mtot,rhochar
      
      REAL*8   Menc
      EXTERNAL Menc
      
c---

      rmx = 10.0**xlgrmax
      Mtot = Menc(rmx) * Ms
      WRITE(*,*)' >>>> ',xlgrmax, rmx, Mtot
      
c---compute rho_s in Msun/kpc^3

      rhochar = rho_s * 1.0D-9
            
      WRITE(*,*)' '
      WRITE(*,*)'====================================================='
      IF (imode.EQ.1) WRITE(*,*)'  Untruncated NFW profile'
      IF (imode.EQ.2) WRITE(*,*)'  Power-Law truncated NFW profile'
      IF (imode.EQ.3) WRITE(*,*)'  Energy truncated NFW profile'
      IF (imode.EQ.4) WRITE(*,*)'  Exponentially truncated NFW profile'
      IF (imode.EQ.6) WRITE(*,*)'  Alpha Beta Gamma profile'
      WRITE(*,*)'====================================================='
      WRITE(*,*)' '     
      WRITE(*,*)'          log[Mvir/Msun] = ',DLOG10(Mvir)
      WRITE(*,*)'          log[Mtot/Msun] = ',DLOG10(Mtot)
      WRITE(*,*)'             Vvir [km/s] = ',Vvir
      WRITE(*,*)'              v_0 [km/s] = ',v0
      WRITE(*,*)' log[rho_s/(Msun/kpc^3)] = ',DLOG10(rhochar)
      WRITE(*,*)'              r_s  [kpc] = ',rs*1.0D+3
      WRITE(*,*)'             r_trunc/r_s = ',rtrunc
      WRITE(*,*)'                    cvir = ',cvir
      WRITE(*,*)' '
      WRITE(*,*)'                   alpha = ',alpha
      WRITE(*,*)'                    beta = ',beta
      WRITE(*,*)'                   gamma = ',gamma
      WRITE(*,*)'               power-law = ',ppar
      WRITE(*,*)'                     Z_t = ',Zt
      WRITE(*,*)'         r_decay/r_trunc = ',eta
      WRITE(*,*)' '
      WRITE(*,*)'                flyby_on = ',flyby_on
      WRITE(*,*)'                    tfly = ',tfly
      WRITE(*,*)'                    bfly = ',bfly
      WRITE(*,*)'                    Mfly = ',Mfly
      WRITE(*,*)'                    vfly = ',vfly    
      WRITE(*,*)' '
      WRITE(*,*)'                 xlgrmin = ',xlgrmin
      WRITE(*,*)'                 xlgrmax = ',xlgrmax
      WRITE(*,*)'                   Ngrid = ',Ngrid
      WRITE(*,*)' '     
      WRITE(*,*)'          f_b [unitless] = ',f_b
      WRITE(*,*)'       zeta_b [unitless] = ',zeta_b
      WRITE(*,*)' '     
      WRITE(*,*)'      sigma_m [unitless] = ',sigma_m   
      WRITE(*,*)'   Gamma_evap [unitless] = ',Gamma_evap 
      WRITE(*,*)' '
      WRITE(*,*)'        sigma_m [cm^2/g] = ',sigma_m * sigma0  
      WRITE(*,*)'               t_0 [Gyr] = ',t0_Gyr
      WRITE(*,*)' '
      IF (iRevir.EQ.1) WRITE(*,*)'  Using Revirialization Method 1'
      IF (iRevir.EQ.2) WRITE(*,*)'  Using Revirialization Method 2'
      IF (iRevir.EQ.3) WRITE(*,*)'  Using Revirialization Method 3'
      WRITE(*,*)'====================================================='
      WRITE(*,*)' '

c---write same info to logfile

      WRITE(97,*)' '
      WRITE(97,*)'====================================================='
      IF (imode.EQ.1) WRITE(97,*)'  Untruncated NFW profile'
      IF (imode.EQ.2) WRITE(97,*)'  Power-Law truncated NFW profile'
      IF (imode.EQ.3) WRITE(97,*)'  Energy truncated NFW profile'
      IF (imode.EQ.4) WRITE(97,*)'  Exponentially truncated NFW profile'
      IF (imode.EQ.6) WRITE(97,*)'  Alpha Beta Gamma profile'
      WRITE(97,*)'====================================================='
      WRITE(97,*)' '     
      WRITE(97,*)'             Mvir [Msun] = ',Mvir
      WRITE(97,*)'          log[Mtot/Msun] = ',DLOG10(Mtot)
      WRITE(97,*)'             Vvir [km/s] = ',Vvir
      WRITE(97,*)'              v_0 [km/s] = ',v0
      WRITE(97,*)' log[rho_s/(Msun/kpc^3)] = ',DLOG10(rhochar)
      WRITE(97,*)'              r_s  [kpc] = ',rs*1.0D+3
      WRITE(97,*)'             r_trunc/r_s = ',rtrunc
      WRITE(97,*)'                    cvir = ',cvir
      WRITE(97,*)' '
      WRITE(97,*)'                   alpha = ',alpha
      WRITE(97,*)'                    beta = ',beta
      WRITE(97,*)'                   gamma = ',gamma
      WRITE(97,*)'               power-law = ',ppar
      WRITE(97,*)'                     Z_t = ',Zt
      WRITE(97,*)'         r_decay/r_trunc = ',eta
      WRITE(97,*)' '
      WRITE(97,*)'                flyby_on = ',flyby_on
      WRITE(97,*)'                    tfly = ',tfly
      WRITE(97,*)'                    bfly = ',bfly
      WRITE(97,*)'                    Mfly = ',Mfly
      WRITE(97,*)'                    vfly = ',vfly    
      WRITE(97,*)' '
      WRITE(97,*)'                 xlgrmin = ',xlgrmin
      WRITE(97,*)'                 xlgrmax = ',xlgrmax
      WRITE(97,*)'                   Ngrid = ',Ngrid
      WRITE(97,*)' '     
      WRITE(97,*)'          f_b [unitless] = ',f_b
      WRITE(97,*)'       zeta_b [unitless] = ',zeta_b
      WRITE(97,*)' '     
      WRITE(97,*)'      sigma_m [unitless] = ',sigma_m   
      WRITE(97,*)'   Gamma_evap [unitless] = ',Gamma_evap 
      WRITE(97,*)' '
      WRITE(97,*)'        sigma_m [cm^2/g] = ',sigma_m * sigma0  
      WRITE(97,*)'               t_0 [Gyr] = ',t0_Gyr
      WRITE(97,*)' '
      IF (iRevir.EQ.1) WRITE(97,*)'  Using Revirialization Method 1'
      IF (iRevir.EQ.2) WRITE(97,*)'  Using Revirialization Method 2'
      IF (iRevir.EQ.3) WRITE(97,*)'  Using Revirialization Method 3'
      WRITE(97,*)'====================================================='
      WRITE(97,*)' '

      RETURN
      END
            
c**********************************************************************
c  Numerical Recipes Routines
c**********************************************************************

      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      REAL*8 bet,gam(NMAX)
      if(b(1).eq.0.d0) CALL Terminate('tridag: rewrite equations')
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.d0) CALL Terminate('tridag failed')
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

c---

      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL*8 x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END

c---

      REAL*8 FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      REAL*8   ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.d0*tol1
        if(dabs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.d0) p=-p
          q=dabs(q)
          etemp=e
          e=d
          if(dabs(p).ge.dabs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-
     *x)) goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(dabs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      WRITE(*,*)'brent exceed maximum iterations'
      STOP 
3     xmin=x
      brent=fx
      return
      END

c---

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
      
c---

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs,bsstep
      PARAMETER (MAXSTP=10000,NMAX=2,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i)) + dabs(h*dydx(i)) + TINY
12      continue
        if(kmax.gt.0)then
          if(dabs(x-xsav).gt.dabs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
        call bsstep(y,dydx,nvar,x,h,eps,yscal,hdid,hnext)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.d0)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) THEN
          CALL Terminate('stepsize smaller than minimum in odeint')
        endif  
        h=hnext
16    continue
      CALL Terminate('too many steps in odeint')
      return
      END

c---

      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     *REDMAX=1.e-5,REDMIN=.7,TINY=1.e-30,SCALMX=.1)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.e29
        xnew=-1.e29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x) CALL Terminate('step size underflow in bsstep')
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,dabs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END

c---

      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout)
      INTEGER nstep,nvar,NMAX
      REAL*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL*8 h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END

c---

      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END

c---

      REAL FUNCTION expint(n,x)
      INTEGER n,MAXIT
      REAL    x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
      INTEGER i,ii,nm1
      REAL a,b,c,d,del,fact,h,psi
c---personal modification (I should only call expint for n=1)
      IF (n.NE.1) THEN
        WRITE(98,*)' >>> USING EXPINT : ',n,x
        WRITE(*,*)' >>> USING EXPINT : ',n,x
        STOP
      END IF  
      expint = 0.0  

      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        CALL Terminate('bad arguments in expint')
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0.)then
        expint=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
            expint=h*exp(-x)
            return
          endif
11      continue
        CALL Terminate('continued fraction failed in expint')
      else
        if(nm1.ne.0)then
          expint=1./nm1
        else
          expint=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*EPS) return
13      continue
        CALL Terminate('series failed in expint')
      endif
      return
      END

c---
       
      REAL*8 FUNCTION dgammln(xx)
      REAL*8  xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x+0.5d0) * DLOG(tmp) - tmp
      ser = 1.000000000190015d0
      do 11 j=1,6
        y = y + 1.d0
        ser = ser + cof(j)/y
11    continue
      dgammln = tmp + DLOG(stp*ser/x)
      return
      END

c---
              
      SUBROUTINE dgcf2(gammcf,a,x)

      INTEGER    ITMAX
      REAL*8     EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)

      INTEGER    i
      REAL*8     gammcf,a,x,an,b,c,d,del,h,xi,gln
        
      REAL*8     dgammln  
      EXTERNAL   dgammln 

      gln = dgammln(a)
      b = x + 1.0d0 - a
      c = 1.0d0 / FPMIN
      d = 1.0d0 / b
      h = d
      DO 11 i=1,ITMAX
        xi = DBLE(i)  
        an = -xi * (xi-a)
        b = b + 2.0d0
        d = an*d + b
        IF (DABS(d).LT.FPMIN) d = FPMIN
        c = b + an/c
        IF (DABS(c).LT.FPMIN) c = FPMIN
        d = 1.0d0/d
        del = d*c
        h = h*del
        IF (DABS(del-1.0d0).LT.EPS) GOTO 1
11    CONTINUE
      CALL Terminate('a too large, ITMAX too small in gcf')
1     gammcf = DEXP(-x + a*DLOG(x) - gln) * h
      gammcf = gammcf * DEXP(gln) 
      return
      END

c---
       
      SUBROUTINE dgser2(gamser,a,x)

      INTEGER ITMAX
      REAL*8  EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)

      INTEGER  n
      REAL*8   a,gamser,x,ap,del,sum,gln
      REAL*8   dgammln
      EXTERNAL dgammln
      
      gln = dgammln(a)
      ap = a
      sum = 1.0d0/a
      del = sum
      DO 11 n=1,ITMAX
        ap = ap + 1.0d0
        del = del * x/ap
        sum = sum + del
        IF (DABS(del).LT.DABS(sum)*EPS) GOTO 1
11    CONTINUE
      CALL Terminate('a too large, ITMAX too small in gser')
1     gamser = sum * DEXP(-x + a*DLOG(x) - gln)
      gamser = (1.0d0-gamser) * DEXP(gln) 
      return
      END

c---

c***********************************************************************
             
      
      