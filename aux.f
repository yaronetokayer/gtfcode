c**********************************************************************
c          SOME USEFUL AUXILIARY FUNCTIONS & SUBROUTINES
c
c   - subroutine that reads time
c   - subroutine for linear interpolation
c   - ARCTANH
c   - TANH
c   - XEXP
c   - lblnk
c
c  Frank van den Bosch                                    Dec. 1999
c**********************************************************************

      SUBROUTINE get_time()
c----------------------------------------------------------------------
c
c  This subroutine computes time since start of execution, and
c  since last call to this subroutine. It uses two intrinsic functions:
c    - etime returns time (in sec) elapsed since start of execution
c    - dtime returns time (in sec) elapsed since last call to dtime
c  The argument ttime(2) on return yields the user time, ttime(1),
c  and the system time, ttime(2). The times (e1 and e2) returned by
c  the functions themselves are the sums of user and system time.
c
c See the following html for information
c   http://docs.oracle.com/cd/E19957-01/805-4942/6j4m3r8t4/index.html
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL    ttime1(2),ttime2(2)

c---
  
      e_time = etime(ttime1)
      d_time = dtime(ttime2)

      et_h = INT(e_time/3600.0)
      et_m = INT(MOD(e_time,3600.0)/60.0)
      et_s = MOD(MOD(e_time,3600.0),60.0)

      dt_tot = d_time
      dt_h = INT(d_time/3600.0)
      dt_m = INT(MOD(d_time,3600.0)/60.0)
      dt_s = MOD(MOD(d_time,3600.0),60.0)

      RETURN
      END      

c**********************************************************************

      SUBROUTINE write_time(ichoice)
c----------------------------------------------------------------------
c
c  IF (ichoice=1) write total elapsed time 
c  IF (ichoice=2) write elapsed time since last call
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      INTEGER ichoice

c---
  
      CALL get_time

      IF (ichoice.EQ.1) THEN
        WRITE(*,73)et_h,et_m,et_s
      ELSE
        WRITE(*,74)dt_h,dt_m,dt_s
      END IF

c---

 73   FORMAT('Total elapsed time: ',I2,'h',I2,'m',F5.2,'s')
 74   FORMAT('Time elapsed since last call: ',I2,'h',I2,'m',F5.2,'s')

      RETURN
      END      

c**********************************************************************

      REAL FUNCTION XEXP(x)
c--------------------------------------------------------------------
c
c Auxialiary function to compute EXP(x)
c
c--------------------------------------------------------------------

      REAL    x

c---

      IF (x.LT.-40.0) THEN
        XEXP = 0.0
      ELSE
        XEXP = EXP(x)
      END IF

      END

c**********************************************************************

      INTEGER FUNCTION lblnk(char)
c--------------------------------------------------------------------
c
c Function gives NUMBER of characters in a character variable `char'
c
c--------------------------------------------------------------------

      IMPLICIT NONE

      character char*(*)

c---

      lblnk=index(char,' ')-1

      RETURN
      END

c**********************************************************************

      SUBROUTINE Terminate(message)
c--------------------------------------------------------------------
c
c  Output error message and terminate program
c
c--------------------------------------------------------------------

      IMPLICIT NONE

      character  message*(*)

c---

      WRITE(*,'(A)')message

      STOP

      RETURN
      END

c**********************************************************************






