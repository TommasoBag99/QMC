!======================================================================
      DOUBLE PRECISION FUNCTION psi(a, n_el, N, r_el, r_N)
!======================================================================
!     Computes the many-body wave function psi using the Jastrow factor
!     a. It iterates over electron-nucleus pairs, calculating their 
!     Euclidean distance and accumulating the exponential term for each 
!     pair. The final psi is the product of all psi_eN, representing 
!     the interaction of each electron with each nucleus.
!======================================================================
      
      IMPLICIT NONE

      INTEGER                      :: i, j
      INTEGER, INTENT(IN)          :: n_el, N

      DOUBLE PRECISION             :: d_eN, psi_eN
      DOUBLE PRECISION, INTENT(IN) :: a, r_el(n_el,3), r_N(N,3)

      psi = 1.d0

      DO i = 1, n_el
         psi_eN = 0.d0
         DO j = 1, N
            d_eN = SQRT((r_el(i,1)-r_N(j,1))**2.d0 +                  &
                        (r_el(i,2)-r_N(j,2))**2.d0 +                  &
                        (r_el(i,3)-r_N(j,3))**2.d0)
            psi_eN = psi_eN + DEXP(-a*d_eN)
         END DO
         psi = psi * psi_eN
      END DO

      END FUNCTION psi

!======================================================================
      SUBROUTINE drift(a, k, n_el, N, r_el, r_N, b)
!======================================================================
!     Calculating the drift vector for Monte Carlo sampling. Computes 
!     the first derivative of the wave function to guide electron moves
!     towards regions of higher probability, aiding in efficient 
!     exploration of configuration space in QMC simulations.
!======================================================================

      IMPLICIT NONE

      INTEGER                       :: i, j, k
      INTEGER, INTENT(IN)           :: n_el, N

      DOUBLE PRECISION              :: d_eN, psi_sum(3), d1_psi(3)
      DOUBLE PRECISION, INTENT(IN)  :: a, r_el(n_el,3), r_N(N,3)
      DOUBLE PRECISION, INTENT(OUT) :: b(3)
      DOUBLE PRECISION, EXTERNAL    :: psi

      d1_psi(:) = 1.d0

      DO i = 1, n_el
         psi_sum(:) = 0.d0
         IF (i /= k) THEN
            DO j = 1, N
               d_eN = SQRT((r_el(i,1)-r_N(j,1))**2.d0 +               &
                           (r_el(i,2)-r_N(j,2))**2.d0 +               &
                           (r_el(i,3)-r_N(j,3))**2.d0)
               psi_sum(:) = psi_sum(:) + DEXP(-a*d_eN)
            END DO
         ELSE
            DO j = 1, N
               d_eN = SQRT((r_el(i,1)-r_N(j,1))**2.d0 +               &
                           (r_el(i,2)-r_N(j,2))**2.d0 +               &
                           (r_el(i,3)-r_N(j,3))**2.d0)
               IF (d_eN < 1.d-8) THEN
                  WRITE(0,*)'warning drift: 1/d_eN diverges'
                  STOP
               END IF
               psi_sum(:) = psi_sum(:) -                              &
                            a*DEXP(-a*d_eN) *                         &
                            ((r_el(i,:)-r_N(j,:))/d_eN)
            END DO
         END IF
         d1_psi(:) = d1_psi(:) * psi_sum(:)
      END DO

      b(:) = d1_psi(:) / psi(a, n_el, N, r_el, r_N)

      END SUBROUTINE drift

!======================================================================
      DOUBLE PRECISION FUNCTION e_loc(a, n_el, N, r_el, r_N, V_NN, z)
!======================================================================
!     Calculates the local energy as the sum of kiinetic and potential
!     energies in Quantum Monte Carlo simulations.
!======================================================================

      IMPLICIT NONE

      INTEGER, INTENT(IN)          :: n_el, N, z(N)

      DOUBLE PRECISION, INTENT(IN) :: a, r_el(n_el,3), r_N(N,3), V_NN
      DOUBLE PRECISION, EXTERNAL   :: kinetic, potential

      e_loc = kinetic(a, n_el, N, r_el, r_N) +                        &
              potential(n_el, N, r_el, r_N, V_NN, z)

      END FUNCTION e_loc

!======================================================================
      DOUBLE PRECISION FUNCTION kinetic(a, n_el, N, r_el, r_N)
!======================================================================
!     Calculating the kinetic energy in QMC simulations. Computes the 
!     second derivative of the wave function and evaluates the kinetic 
!     energy as -0.5 times the ratio of the second derivative to the 
!     total wave function.
!======================================================================

      IMPLICIT NONE

      INTEGER                       :: i, j, m
      INTEGER, INTENT(IN)           :: n_el, N

      DOUBLE PRECISION              :: d_eN, k_sum, T, d2_psi
      DOUBLE PRECISION, INTENT(IN)  :: a, r_el(n_el,3), r_N(N,3)
      DOUBLE PRECISION, EXTERNAL    :: psi

      T = 0.d0

      DO m = 1, n_el
         d2_psi = 1.d0
         DO i = 1, n_el
            k_sum = 0.d0
            IF (i /= m) THEN
               DO j = 1, N
                  d_eN = SQRT((r_el(i,1)-r_N(j,1))**2.d0 +            &
                              (r_el(i,2)-r_N(j,2))**2.d0 +            &
                              (r_el(i,3)-r_N(j,3))**2.d0)
                  k_sum = k_sum + DEXP(-a * d_eN)
               END DO
            ELSE
               DO j = 1, N
                  d_eN = SQRT((r_el(i,1)-r_N(j,1))**2.d0 +            &
                              (r_el(i,2)-r_N(j,2))**2.d0 +            &
                              (r_el(i,3)-r_N(j,3))**2.d0)
                  IF (d_eN < 1.d-8) THEN
                     WRITE(0,*)'warning kinetic: 1/d_eN diverges'
                     STOP                  
                  END IF
                  k_sum = k_sum +                                     &
                          DEXP(-a*d_eN) *                             &
                          (a*a - ((2*a)/d_eN))
               END DO
            END IF
            d2_psi = d2_psi * k_sum
         END DO
         T = T + d2_psi
      END DO

      kinetic = - 0.5d0 * T / psi(a, n_el, N, r_el, r_N)

      END FUNCTION kinetic

!======================================================================
      DOUBLE PRECISION FUNCTION potential(n_el, N, r_el, r_N, V_NN, z)
!======================================================================
!     Computes the total potential energy by summing the contributions 
!     from electron-nucleus, electron-electron and nucleus-nucleus 
!     interaction term.     
!======================================================================

      IMPLICIT NONE

      INTEGER                      :: i, j
      INTEGER, INTENT(IN)          :: n_el, N, z(N)

      DOUBLE PRECISION             :: d_eN, d_ee, V_eN, V_ee
      DOUBLE PRECISION, INTENT(IN) :: r_el(n_el,3), r_N(N,3), V_NN

      ! electron-nucleus interaction

      V_eN = 0.d0

      DO i = 1, N
         DO j = 1, n_el
            d_eN = SQRT((r_el(j,1)-r_N(i,1))**2.d0 +                  &
                        (r_el(j,2)-r_N(i,2))**2.d0 +                  &
                        (r_el(j,3)-r_N(i,3))**2.d0)
            IF (d_eN < 1.d-8) THEN
               WRITE(0,*)'warning potential: potential-eN diverges'
               STOP
            END IF
            V_eN = V_eN - DBLE(z(i)) / d_eN
         END DO
      END DO

      ! electron-electron interaction

      V_ee = 0.d0

      DO i = 1, n_el-1
         DO j = i+1, n_el
            d_ee = SQRT((r_el(j,1)-r_el(i,1))**2.d0 +                 &
                        (r_el(j,2)-r_el(i,2))**2.d0 +                 &
                        (r_el(j,3)-r_el(i,3))**2.d0)
            IF (d_ee < 1.d-8) THEN
               WRITE(0,*)'warning potential: potential-ee diverges'
               STOP
            END IF
            V_ee = V_ee + 1.d0 / d_ee
         END DO
      END DO

      potential = V_NN + V_eN + V_ee

      END FUNCTION potential

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
      DOUBLE PRECISION FUNCTION potential_NN(N, r_N, z)
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      

      IMPLICIT NONE

      INTEGER                      :: i, j
      INTEGER, INTENT(IN)          :: N, z(N)

      DOUBLE PRECISION             :: V_NN, d_NN 
      DOUBLE PRECISION, INTENT(IN) :: r_N(N,3)

      ! nucleus-nucleus interaction

      V_NN = 0.d0

      DO i = 1, N-1
         DO j = i+1, N
            d_NN = SQRT((r_N(j,1)-r_N(i,1))**2.d0 + &
                        (r_N(j,2)-r_N(i,2))**2.d0 + &
                        (r_N(j,3)-r_N(i,3))**2.d0)
            IF (d_NN < 1.d-8) THEN
               WRITE(0,*)'warning potential: potential-NN diverges'
               STOP
            END IF
            V_NN = V_NN + DBLE(z(i)) * DBLE(z(j)) / d_NN
         END DO
      END DO
      
      potential_NN = V_NN

      END FUNCTION potential_NN
