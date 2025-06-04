c=======================================================================
c
c                      INCLUDE FILE paramfile.h
c
c=======================================================================
c  Written by: Frank van den Bosch
c=======================================================================

      IMPLICIT NONE

c=======================================================================
c                        PARAMETERS
c=======================================================================

c---The gravitational constant in M_sun^-1 Mpc (km/s)^2

      REAL*8  gee
      PARAMETER (gee = 4.2994D-9)

c---The conversion of Mpc to cm

      REAL*8  Mpc_to_cm
      PARAMETER (Mpc_to_cm = 3.086D+24)

c---The conversion of Msun to grams

      REAL*8  Msun_to_gram
      PARAMETER (Msun_to_gram = 1.99D+33)

c---The conversion of Msun to grams

      REAL*8  sec_to_Gyr
      PARAMETER (sec_to_Gyr = 3.16881D-17)

c---the constant pi

      REAL*8  pi
      PARAMETER (pi=3.141592654d0)

c---small parameter, needed for odeint
      
      REAL*8   epsilon
      PARAMETER (epsilon=1.0D-2)      

c---parameters that specify DF of NFW profile

      REAL*8  F0,q,p1,p2,p3,p4
      PARAMETER (F0 = 9.1968E-2)
      PARAMETER (q = -2.7419d0)
      PARAMETER (p1 = 0.3620d0)
      PARAMETER (p2 =-0.5639d0)
      PARAMETER (p3 =-0.0859d0)
      PARAMETER (p4 =-0.4912d0)

c---number of radial grid points used

      INTEGER Ngrid_max
      PARAMETER (Ngrid_max=500)

c---number of velocities to draw for computing velocity dispersion

      INTEGER Nvel
      PARAMETER (Nvel=1000000)

c---parameters needed for odeint

      INTEGER  Nstep,Nmx
      PARAMETER (Nstep=10)
      PARAMETER (Nmx=1000000)

c---parameters needed for quadpack integration routines

      INTEGER Nlimit,Nlenw
      PARAMETER (Nlimit=10000)
      PARAMETER (Nlenw =41000)
      
c---the null character

      CHARACTER null*1
      PARAMETER (null='0')

c=======================================================================
c                        COMMON BLOCKS
c=======================================================================

c---debugging

      LOGICAL checkpoint
      COMMON /debug/ checkpoint

c---cosmology

      REAL*8   xhubble,xH_0,omega_0,omega_lambda
      COMMON /cospar/ xhubble,xH_0,omega_0,omega_lambda

c---radial grid

      INTEGER  Ngrid
      REAL*8   xlgrmin,xlgrmax,r(0:Ngrid_max),M(0:Ngrid_max)
      REAL*8   L(0:Ngrid_max)
      REAL*8   rho(Ngrid_max),P(Ngrid_max),u(Ngrid_max),v2(Ngrid_max)
      COMMON /radialgrid/ xlgrmin,xlgrmax,r,M,L,rho,P,u,v2,Ngrid
      
c---halo properties

      REAL*8   Delta_vir,Mvir,Ms,cvir,Vvir,rvir,rs,fc,Mtotal,zvir
      COMMON /haloprop/ Delta_vir,Mvir,Ms,cvir,Vvir,rvir,rs,fc,
     &                  Mtotal,zvir

c---baryon properties

      REAL*8   f_b, zeta_b
      COMMON /baryonprop/  f_b, zeta_b

c---parameters related to infalling perturber

      REAL*8 tsinki, tsinkf  ! start and end times of heat dump
      REAL*8 Rst             ! stalling radius in units of rs
      REAL*8 Mp              ! perturber mass
      REAL*8 Psink           ! power of heating
      REAL*8 Mrate           ! rate of point mass growth
      REAL*8 kecorei,eratio  ! initial KE of core
      INTEGER heat_func      ! option for location of heating
      LOGICAL infall_pert    ! whether to model an infalling perturber
      LOGICAL infall_triggered ! for the loop

      COMMON /infallparam/ tsinki, tsinkf, Mp, Rst, Psink, kecorei,
     &                    eratio, Mrate,heat_func, infall_pert,
     &                    infall_triggered

c---properties related to initial profile

      INTEGER  imode
      COMMON /mode/ imode

      REAL*8   alpha, beta, gamma, chi
      COMMON /abgpars/ alpha, beta, gamma, chi
      
      REAL*8   rtrunc,ppar,Mrt,prefac,pow,power,rtr,ftr
      REAL*8   Ft,Zt,rcut,eta,eps
      COMMON /trunc/ rtrunc,ppar,Mrt,prefac,pow,power,rtr,ftr,Ft,Zt,
     &               rcut,eta,eps

c---collision parameters

      REAL*8   sigma_m,apar,bpar,cpar,Gamma_evap
      COMMON /charprop/ sigma_m,apar,bpar,cpar,Gamma_evap
        
c---characteristic quantities

      REAL*8   v0,rho_s,t0,sigma0,maxvel,minKn,rc1,rc2,rc3,r_BH,M_BH
      REAL*8   rho_s_cgs,v0_cgs,t0_Gyr,rho_s_Msunpc3,v0_kms
      COMMON /charprop/ v0,rho_s,t0,sigma0,maxvel,minKn,rc1,rc2,rc3,
     &                  r_BH,M_BH,rho_s_cgs,v0_cgs,t0_Gyr,rho_s_Msunpc3,
     &                  v0_kms

c---characteristic densities

      REAL*8    rhochar1,rhochar2,rhochar3,rhochar4
      REAL*8    rhochar5,rhochar6,rhochar7,rhochar8
      COMMON /chardens/ rhochar1,rhochar2,rhochar3,rhochar4,
     &                  rhochar5,rhochar6,rhochar7,rhochar8

c---parameters related to time stepping and virialization convergence
        
      INTEGER  idrmax,idumax,itimestep
      REAL*8   t,eps_dt,eps_dr,eps_du,trelmin,dt,tstop,drmax,dumax,dtmax
      COMMON /timestepping/ t,eps_dt,eps_dr,eps_du,trelmin,dt,tstop,
     &                      drmax,dumax,dtmax,idrmax,idumax,itimestep
        
c---work space required for quadpack integration routines

      INTEGER Neval,Neval2,Neval3,last,last2,last3
      INTEGER iwork(Nlimit),iwork2(Nlimit),iwork3(Nlimit)
      REAL*8  work(Nlenw),work2(Nlenw),work3(Nlenw)
      COMMON /quadpackint/ work,work2,work3,iwork,iwork2,iwork3,
     &                     Neval,Neval2,Neval3,last,last2,last3

c---elapsed time

      INTEGER et_h,et_m,dt_h,dt_m
      REAL    e_time,d_time,et_s,dt_s,dt_tot
      COMMON /elapsedtime/ e_time,d_time,et_h,et_m,et_s,dt_h,dt_m,dt_s,
     &                     dt_tot

c---parameters needed for odeint of potential

      INTEGER  Nmax                   
      REAL*8   rad(Nmx),pot(Nmx)
      COMMON /profs/ rad,pot,Nmax

      REAL*8   deltaP
      COMMON /odeintegrate/ deltaP
      
c---random number seed
      
      INTEGER iseed
      COMMON /randomnr/ iseed
      
c---logicals

      LOGICAL  Failed
      COMMON /samplelogic/ Failed

c---Revirialization Method

      INTEGER  iRevir
      COMMON /methodology/ iRevir

c---Directories etc.

      CHARACTER moddir*60
      COMMON /modeldir/ moddir

c=======================================================================
c                             END
c=======================================================================





