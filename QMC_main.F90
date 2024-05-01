!======================================================================
      PROGRAM qmc
!======================================================================
!     This code implements the Quantum Monte Carlo (QMC) method,
!     allowing for both Pure Diffusion and Variational techniques.
!     It is designed to solve quantum mechanical systems by
!     statistically sampling over potential configurations.
!======================================================================

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     DECLARATION SECTION
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      IMPLICIT NONE

      CHARACTER*4                   :: name_sys
      CHARACTER*3                   :: method
      CHARACTER*8                   :: geom_file
      CHARACTER*2, ALLOCATABLE      :: atom(:)

      INTEGER                       :: irun, nruns, i, j
      INTEGER                       :: N, n_el, q_tot
      INTEGER*8                     :: nmax
      INTEGER, ALLOCATABLE          :: z(:)

      DOUBLE PRECISION              :: d_NN, V_NN
      DOUBLE PRECISION              :: a, dt, E_ref, tau, ave, err
      DOUBLE PRECISION              :: energy, e_err, accept, a_err
      DOUBLE PRECISION, ALLOCATABLE :: r_N(:,:), X(:), accep(:)
      DOUBLE PRECISION, EXTERNAL    :: potential_NN

      DOUBLE PRECISION, PARAMETER   :: A_to_Bohr = 1.8897259886d0

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     INPUT SECTION
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      OPEN(10, FILE = 'QMC_input', STATUS = 'old')

      READ(10,*) name_sys      ! name of the system
      READ(10,*) n_el          ! Number of electrons
      READ(10,*) method        ! var/dif 
      READ(10,*) geom_file     ! Geometry file .xyz
      READ(10,*) a             ! Jastrow factor
      READ(10,*)
      READ(10,*) dt            ! Time step 
      READ(10,*) nmax          ! Number of steps
      READ(10,*) nruns         ! Number of walkers
      READ(10,*)
      READ(10,*) E_ref         ! Reference Energy
      READ(10,*) tau           ! Projection time
      
      CLOSE (10)

      OPEN(20, FILE = geom_file, STATUS = 'old')

      READ(20,*) N             ! Number of atoms
      READ(20,*)

      ALLOCATE(atom(N), r_N(N,3), z(N), X(nruns), accep(nruns))
      
      ! Read atom types
      ! Read xyz positions
      ! Assign nuclear charges 

      DO i = 1, N
         READ(20,*) atom(i), r_N(i,:) 
         IF (atom(i) == 'H') THEN
            z(i) = 1.d0
         ELSE IF (atom(i) == 'He') THEN
            z(i) = 2.d0
         ELSE
            WRITE(0,*) 'INVALID ATOM: select "H" or "He" in .xyz file'
            STOP
         END IF
      END DO

      CLOSE (20)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     SIMULATION AND ANALYSIS SECTION
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      WRITE(6,*)
      WRITE(6,'(A)') 'WORKING...'

      r_N(:,:) = r_N(:,:) * A_to_Bohr      ! From Angstrom to Bohr
      q_tot    = SUM(z) - n_el             ! Global charge
      V_NN     = potential_NN(N, r_N, z)   ! Nucleus-nucleus potential

      ! Variational MC
      IF (method == 'var') THEN         
         DO irun = 1, nruns  
            CALL vmc(a, dt, nmax, X(irun), accep(irun),               &
                     n_el, N, r_N, V_NN, z)
         END DO
      ! Pure Diffusion MC
      ELSE IF (method == 'dif') THEN
         DO irun = 1, nruns
            CALL pdmc(a, dt, nmax, X(irun), accep(irun),              &
                      tau, E_ref, n_el, N, r_N, V_NN, z)
         END DO
      ELSE
         WRITE(0,*) 'INVALID METHOD: select "var" or "dif" in input'
         STOP
      END IF

      ! Acceptance rate and Energy calulations
      CALL ave_error(X,nruns,ave,err)
      energy = ave
      e_err  = err
      CALL ave_error(accep,nruns,ave,err)
      accept = ave
      a_err  = err

      ! Print results 
      CALL output(N, r_N, q_tot, atom, E_ref,                         &
                  method, energy, accept, e_err, a_err)
      CALL data_out(nruns, X, name_sys)

      END PROGRAM qmc

!======================================================================
      SUBROUTINE vmc(a, dt, nmax, energy, accep, n_el, N, r_N, V_NN, z)
!======================================================================
!     This subroutine implements the Variational Monte Carlo (VMC)
!     method for quantum systems. It performs Monte Carlo simulations
!     to estimate the ground state energy by sampling electron
!     configurations using a trial wave function and evaluating
!     their local energies. The subroutine includes the Metropolis
!     algorithm to accept or reject moves based on the ratio of
!     probabilities derived from the wave function.
!======================================================================

      IMPLICIT NONE

      INTEGER                       :: i
      INTEGER*8                     :: istep, n_accep
      INTEGER, INTENT(IN)           :: n_el, N, z(N)
      INTEGER*8, INTENT(IN)         :: nmax

      DOUBLE PRECISION              :: psi_old, psi_new
      DOUBLE PRECISION              :: sq_dt, u, q, argexpo
      DOUBLE PRECISION              :: prod(n_el), chi(n_el,3)
      DOUBLE PRECISION              :: r_old(n_el,3), r_new(n_el,3)
      DOUBLE PRECISION              :: d_old(n_el,3), d_new(n_el,3)
      DOUBLE PRECISION              :: d2_old(n_el), d2_new(n_el)
      DOUBLE PRECISION, INTENT(IN)  :: a, dt, r_N(N,3), V_NN
      DOUBLE PRECISION, INTENT(OUT) :: energy, accep
      DOUBLE PRECISION, EXTERNAL    :: e_loc, psi

      sq_dt = DSQRT(dt)

      energy  = 0.d0
      n_accep = 0_8

      DO i = 1, n_el
         CALL random_gauss(r_old(i,:),3)
      END DO

      DO i = 1, n_el
         CALL drift(a, i, n_el, N, r_old, r_N, d_old(i,:))
         d2_old(i)  = SUM(d_old(i,:)*d_old(i,:))
      END DO

      psi_old = psi(a, n_el, N, r_old, r_N)

      DO istep = 1,nmax
      
         energy = energy + e_loc(a, n_el, N, r_old, r_N, V_NN, z)

         DO i = 1, n_el
            CALL random_gauss(chi(i,:),3)
            r_new(i,:) = r_old(i,:) + dt*d_old(i,:) + chi(i,:)*sq_dt
         END DO

         DO i = 1, n_el
            CALL drift(a, i, n_el, N, r_new, r_N, d_new(i,:))
            d2_new(i) = SUM(d_new(i,:) * d_new(i,:))
         END DO
          
         psi_new = psi(a, n_el, N, r_new, r_N)
         
         ! Metropolis
         prod(:) = (d_new(:,1)+d_old(:,1))*(r_new(:,1)-r_old(:,1)) +  &
                   (d_new(:,2)+d_old(:,2))*(r_new(:,2)-r_old(:,2)) +  & 
                   (d_new(:,3)+d_old(:,3))*(r_new(:,3)-r_old(:,3))

         argexpo = 0.5d0 * SUM(d2_new(:) - d2_old(:))*dt + SUM(prod(:))

         q = psi_new / psi_old
         q = DEXP(-argexpo) * q*q

         CALL random_number(u)

         IF (u <= q) THEN
            n_accep    = n_accep + 1_8
            r_old(:,:) = r_new(:,:)
            d_old(:,:) = d_new(:,:)
            d2_old(:)  = d2_new(:)
            psi_old    = psi_new
         END IF

      END DO

      energy = energy / DBLE(nmax)
      accep  = DBLE(n_accep) / DBLE(nmax)

      END SUBROUTINE vmc

!======================================================================
      SUBROUTINE pdmc(a, dt, nmax, energy, accep,                     &
                      tau, E_ref, n_el, N, r_N, V_NN, z)
!======================================================================
!     This subroutine implements the Pure Diffusion Monte Carlo (PDMC)
!     method. It simulates the quantum evolution of a system using
!     imaginary time propagation to obtain ground state properties.
!     The simulation includes Metropolis sampling to ensure proper
!     statistical weighting and adjusts the system energy based on a
!     reference value.
!======================================================================

      IMPLICIT NONE

      INTEGER                       :: i
      INTEGER*8                     :: istep, n_accep
      INTEGER, INTENT(IN)           :: n_el, N, z(N)
      INTEGER*8, INTENT(IN)         :: nmax

      DOUBLE PRECISION              :: e, w
      DOUBLE PRECISION              :: tau_current
      DOUBLE PRECISION              :: normalization
      DOUBLE PRECISION              :: psi_old, psi_new
      DOUBLE PRECISION              :: sq_dt, u, q, argexpo
      DOUBLE PRECISION              :: prod(n_el), chi(n_el,3)
      DOUBLE PRECISION              :: r_old(n_el,3), r_new(n_el,3)
      DOUBLE PRECISION              :: d_old(n_el,3), d_new(n_el,3)
      DOUBLE PRECISION              :: d2_old(n_el), d2_new(n_el)
      DOUBLE PRECISION, INTENT(IN)  :: E_ref, tau, V_NN
      DOUBLE PRECISION, INTENT(IN)  :: a, dt, r_N(N,3)
      DOUBLE PRECISION, INTENT(OUT) :: energy, accep
      DOUBLE PRECISION, EXTERNAL    :: e_loc, psi

      sq_dt = DSQRT(dt)

      energy        = 0.d0
      n_accep       = 0_8
      normalization = 0.d0
      w             = 1.d0
      tau_current   = 0.d0

      DO i = 1, n_el
         CALL random_gauss(r_old(i,:),3)
      END DO

      DO i = 1, n_el
         CALL drift(a, i, n_el, N, r_old, r_N, d_old(i,:))
         d2_old(i)  = SUM(d_old(i,:)*d_old(i,:))
      END DO

      psi_old = psi(a, n_el, N, r_old, r_N)

      DO istep = 1,nmax

         e = e_loc(a, n_el, N, r_old, r_N, V_NN, z)
         w = w * DEXP(-dt*(e - E_ref))

         normalization = normalization + w
         energy = energy + w*e

         tau_current = tau_current + dt

         ! Reset when tau is reached
         IF (tau_current > tau) THEN
            w           = 1.d0
            tau_current = 0.d0
         ENDIF

         DO i = 1, n_el
            CALL random_gauss(chi(i,:),3)
            r_new(i,:) = r_old(i,:) + dt*d_old(i,:) + chi(i,:)*sq_dt
         END DO

         DO i = 1, n_el
            CALL drift(a, i, n_el, N, r_new, r_N, d_new(i,:))
            d2_new(i) = SUM(d_new(i,:)*d_new(i,:))
         END DO

         psi_new = psi(a, n_el, N, r_new, r_N)

         ! Metropolis
         prod(:) = (d_new(:,1)+d_old(:,1))*(r_new(:,1)-r_old(:,1)) +  &
                   (d_new(:,2)+d_old(:,2))*(r_new(:,2)-r_old(:,2)) +  &
                   (d_new(:,3)+d_old(:,3))*(r_new(:,3)-r_old(:,3))

         argexpo = 0.5d0 * SUM(d2_new(:) - d2_old(:))*dt + SUM(prod(:))

         q = psi_new / psi_old
         q = DEXP(-argexpo) * q*q

         CALL random_number(u)

         IF (u <= q) THEN
            n_accep    = n_accep + 1_8
            r_old(:,:) = r_new(:,:)
            d_old(:,:) = d_new(:,:)
            d2_old(:)  = d2_new(:)
            psi_old    = psi_new
         END IF

      END DO

      energy = energy / normalization
      accep  = DBLE(n_accep) / DBLE(nmax)

      END SUBROUTINE pdmc

