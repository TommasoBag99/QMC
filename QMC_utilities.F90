!======================================================================
      SUBROUTINE random_gauss(z,n,m)
!======================================================================
!     This subroutine is intended for generating random positions in a
!     Gaussian distribution, used to model the spatial distribution of 
!     electrons.      
!======================================================================

      IMPLICIT NONE

      INTEGER                       :: i, k
      INTEGER, INTENT(IN)           :: n, m

      DOUBLE PRECISION, INTENT(OUT) :: z(m,n)
      DOUBLE PRECISION              :: u(n+1)
      DOUBLE PRECISION, PARAMETER   :: two_pi = 2.d0*DACOS(-1.d0)

      DO k = 1, m
         CALL RANDOM_NUMBER(u)
         IF (IAND(n,1) == 0) THEN
            ! n is even
            DO i=1,n,2
               z(k,i)   = DSQRT(-2.d0*DLOG(u(i)))
               z(k,i+1) = z(k,i) * DSIN(two_pi*u(i+1))
               z(k,i)   = z(k,i) * DCOS(two_pi*u(i+1))
            END DO
         ELSE
            ! n is odd
            DO i=1,n-1,2
               z(k,i)   = DSQRT(-2.d0*DLOG(u(i)))
               z(k,i+1) = z(k,i) * DSIN(two_pi*u(i+1))
               z(k,i)   = z(k,i) * DCOS(two_pi*u(i+1))
            END DO
            z(k,n) = DSQRT(-2.d0*DLOG(u(n)))
            z(k,n) = z(k,n) * DCOS(two_pi*u(n+1))
         END IF
      END DO

      END SUBROUTINE random_gauss

!======================================================================
      SUBROUTINE ave_error(x,n,ave,err)
!======================================================================
!     Calculates mean (ave) and standard error (err) of array x.
!     Used mainly to compute mean energy and error after Monte Carlo
!     simulations. Also assesses acceptance rate and relative error,
!     crucial for verifying simulation convergence and accuracy.
!======================================================================

      IMPLICIT NONE
 
      INTEGER, INTENT(IN)           :: n

      DOUBLE PRECISION, INTENT(IN)  :: x(n)
      DOUBLE PRECISION, INTENT(OUT) :: ave, err
      DOUBLE PRECISION              :: variance

      IF (n < 1) THEN
         STOP 'n<1 in ave_error'
      ELSE IF (n == 1) THEN
         ave = x(1)
         err = 0.d0
      ELSE
         ave      = SUM(x(:)) / DBLE(n)
         variance = SUM((x(:) - ave)**2) / DBLE(n-1)
         err      = DSQRT(variance/DBLE(n))
      END IF

      END SUBROUTINE ave_error

!======================================================================
      SUBROUTINE output(N, r_N, q_tot, atom, E_ref,                   &
                        method, energy, accept, e_err, a_err)
!======================================================================
!     Outputs the results of the Monte Carlo calculation.
!     Displays the system and the atomic positions, total charge, and
!     details of the implemented method (variational or diffusion).
!     Shows reference and calculated energies and the acceptance rate.
!     All quantities are expressed in a.u.
!======================================================================

      IMPLICIT NONE
      
      INTEGER          :: N, q_tot, i
      CHARACTER*2      :: atom(N)
      CHARACTER*3      :: method
      DOUBLE PRECISION :: r_N(N,3), E_ref, energy, accept, e_err, a_err

      WRITE(6,*)
      WRITE(6,'(A)')   '+--------+-----------------------------------+'
      WRITE(6,'(A)')   '| SYSTEM |      X          Y          Z      |'
      WRITE(6,'(A)')   '+--------+-----------------------------------+'
      DO i = 1, N
         WRITE(6,1000) '|',atom(i),'|', r_N(i,:),                   '|'
      END DO
      WRITE(6,'(A)')   '+--------+-----------------------------------+'
      WRITE(6,*)
      WRITE(6,'(A)')   '+--------------------------------------------+'
      WRITE(6,2000)    '| CHARGE =', q_tot, '|'
      WRITE(6,'(A)')   '+--------------------------------------------+'
      WRITE(6,*)
      WRITE(6,'(A)')   '+--------------------------------------------+'
      IF (method == 'var') THEN
         WRITE(6,'(A)')'|      VARIATIONAL QUANTUM MONTE CARLO       |'
      ELSE IF (method == 'dif') THEN
         WRITE(6,'(A)')'|     PURE DIFFUSION QUANTUM MONTE CARLO     |'
      END IF
      WRITE(6,'(A)')   '+--------------------------------------------+'
      WRITE(6,4000)    '| ENERGY REF =', E_ref,                     '|'
      WRITE(6,3000)    '| ENERGY QMC =', energy, '+/-', e_err,      '|'
      WRITE(6,3000)    '| ACCEPT QMC =', accept, '+/-', a_err,      '|'
      WRITE(6,'(A)')   '+--------------------------------------------+'
      WRITE(6,*)

1000  FORMAT(A1,3X,A2,3X,A1,1X,3(F10.7,1X),1X,A1)
2000  FORMAT(A,1X,I2,32X,A)
3000  FORMAT(A,1X,F9.6,1X,A3,F9.6,8X,A)
4000  FORMAT(A,1X,F9.6,21X,A)

      END SUBROUTINE OUTPUT

!======================================================================
      SUBROUTINE data_out(nruns, X, name_sys, method)
!======================================================================
!     Generates an output file named based on 'name_sys' that records
!     each run's step number and corresponding energy value. The file's
!     naming is prefixed with 'energy_' and suffixed with '.out'.
!======================================================================

      IMPLICIT NONE

      INTEGER                      :: i
      INTEGER, INTENT(IN)          :: nruns

      CHARACTER*3                  :: method
      CHARACTER*4                  :: name_sys
      CHARACTER*7                  :: pre = "energy_"
      CHARACTER*1                  :: und = "_"
      CHARACTER*4                  :: suf = ".out"

      DOUBLE PRECISION, INTENT(IN) :: X(nruns)

      OPEN (30, FILE = pre//trim(name_sys)//und//trim(method)//suf)

      WRITE(30,1000) '#  i', 'energy'
      DO i = 1, nruns
         WRITE(30,2000) i, X(i)
      END DO

      CLOSE (30)

1000  FORMAT(A,2X,A)
2000  FORMAT(I4,1X,F9.6)

      END SUBROUTINE data_out
