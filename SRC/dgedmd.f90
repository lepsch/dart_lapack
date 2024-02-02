      SUBROUTINE DGEDMD( JOBS, JOBZ, JOBR, JOBF,  WHTSVD,  &
                         M, N, X, LDX, Y, LDY, NRNK, TOL,  &
                         K, REIG,  IMEIG,   Z, LDZ,  RES,  &
                         B, LDB, W,  LDW,   S, LDS,        &
                         WORK, LWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine                                           --
!
!  -- LAPACK is a software package provided by University of          --
!  -- Tennessee, University of California Berkeley, University of     --
!  -- Colorado Denver and NAG Ltd..                                   --
!
!.....
      USE                   iso_fortran_env
      IMPLICIT NONE
      INTEGER, PARAMETER :: WP = real64
!
!     Scalar arguments
!     ~~~~~~~~~~~~~~~~
      CHARACTER, INTENT(IN)   :: JOBS,   JOBZ,  JOBR,  JOBF
      INTEGER,   INTENT(IN)   :: WHTSVD, M, N,   LDX,  LDY, &
                                 NRNK, LDZ, LDB, LDW,  LDS, &
                                 LWORK,  LIWORK
      INTEGER,   INTENT(OUT)  :: K, INFO
      REAL(KIND=WP), INTENT(IN)  :: TOL
!
!     Array arguments
!     ~~~~~~~~~~~~~~~
      REAL(KIND=WP), INTENT(INOUT) :: X(LDX,*), Y(LDY,*)
      REAL(KIND=WP), INTENT(OUT)   :: Z(LDZ,*), B(LDB,*), &
                                      W(LDW,*), S(LDS,*)
      REAL(KIND=WP), INTENT(OUT)   :: REIG(*),  IMEIG(*), &
                                      RES(*)
      REAL(KIND=WP), INTENT(OUT)   :: WORK(*)
      INTEGER,       INTENT(OUT)   :: IWORK(*)
!
!     Parameters
!     ~~~~~~~~~~
      REAL(KIND=WP), PARAMETER ::  ONE = 1.0_WP
      REAL(KIND=WP), PARAMETER :: ZERO = 0.0_WP
!
!     Local scalars
!     ~~~~~~~~~~~~~
      REAL(KIND=WP) :: OFL,    ROOTSC, SCALE,  SMALL,  &
                       SSUM,   XSCL1,  XSCL2
      INTEGER       :: i,   j, IMINWR,  INFO1, INFO2,  &
                       LWRKEV, LWRSDD, LWRSVD,         &
                       LWRSVQ, MLWORK, MWRKEV, MWRSDD, &
                       MWRSVD, MWRSVJ, MWRSVQ, NUMRNK, &
                       OLWORK
      LOGICAL       :: BADXY,  LQUERY, SCCOLX, SCCOLY, &
                       WNTEX,  WNTREF, WNTRES, WNTVEC
      CHARACTER     :: JOBZL,  T_OR_N
      CHARACTER     :: JSVOPT
!
!     Local arrays
!     ~~~~~~~~~~~~
      REAL(KIND=WP) :: AB(2,2), RDUMMY(2), RDUMMY2(2)
!
!     External functions (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~
      REAL(KIND=WP) DLANGE, DLAMCH, DNRM2
      EXTERNAL      DLANGE, DLAMCH, DNRM2, IDAMAX
      INTEGER       IDAMAX
      LOGICAL       DISNAN, LSAME
      EXTERNAL      DISNAN, LSAME
!
!     External subroutines (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~~~~
      EXTERNAL      DAXPY,  DGEMM,  DSCAL
      EXTERNAL      DGEEV,  DGEJSV, DGESDD, DGESVD, DGESVDQ, &
                    DLACPY, DLASCL, DLASSQ, XERBLA
!
!     Intrinsic functions
!     ~~~~~~~~~~~~~~~~~~~
      INTRINSIC     DBLE, INT, MAX, SQRT
!............................................................
!
!    Test the input arguments
!
      WNTRES = LSAME(JOBR,'R')
      SCCOLX = LSAME(JOBS,'S') .OR. LSAME(JOBS,'C')
      SCCOLY = LSAME(JOBS,'Y')
      WNTVEC = LSAME(JOBZ,'V')
      WNTREF = LSAME(JOBF,'R')
      WNTEX  = LSAME(JOBF,'E')
      INFO   = 0
      LQUERY = ( ( LWORK == -1 ) .OR. ( LIWORK == -1 ) )
!
      IF ( .NOT. (SCCOLX .OR. SCCOLY .OR. &
                                  LSAME(JOBS,'N')) )   THEN
          INFO = -1
      ELSE IF ( .NOT. (WNTVEC .OR. LSAME(JOBZ,'N')        &
                              .OR. LSAME(JOBZ,'F')) )  THEN
          INFO = -2
      ELSE IF ( .NOT. (WNTRES .OR. LSAME(JOBR,'N')) .OR.  &
                ( WNTRES .AND. (.NOT.WNTVEC) ) )       THEN
          INFO = -3
      ELSE IF ( .NOT. (WNTREF .OR. WNTEX .OR.             &
                LSAME(JOBF,'N') ) )                    THEN
          INFO = -4
      ELSE IF ( .NOT.((WHTSVD == 1) .OR. (WHTSVD == 2) .OR.  &
                      (WHTSVD == 3) .OR. (WHTSVD == 4) )) THEN
          INFO = -5
      ELSE IF ( M < 0 )   THEN
          INFO = -6
      ELSE IF ( ( N < 0 ) .OR. ( N > M ) ) THEN
          INFO = -7
      ELSE IF ( LDX < M ) THEN
          INFO = -9
      ELSE IF ( LDY < M ) THEN
          INFO = -11
      ELSE IF ( .NOT. (( NRNK == -2).OR.(NRNK == -1).OR. &
                ((NRNK >= 1).AND.(NRNK <=N ))) )      THEN
          INFO = -12
      ELSE IF ( ( TOL < ZERO ) .OR. ( TOL >= ONE ) )  THEN
          INFO = -13
      ELSE IF ( LDZ < M ) THEN
          INFO = -18
      ELSE IF ( (WNTREF .OR. WNTEX ) .AND. ( LDB < M ) ) THEN
          INFO = -21
      ELSE IF ( LDW < N ) THEN
          INFO = -23
      ELSE IF ( LDS < N ) THEN
          INFO = -25
      END IF
!
      IF ( INFO == 0 ) THEN
          ! Compute the minimal and the optimal workspace
          ! requirements. Simulate running the code and
          ! determine minimal and optimal sizes of the
          ! workspace at any moment of the run.
         IF ( N == 0 ) THEN
             ! Quick return. All output except K is void.
             ! INFO=1 signals the void input.
             ! In case of a workspace query, the default
             ! minimal workspace lengths are returned.
            IF ( LQUERY ) THEN
                IWORK(1) = 1
                WORK(1)  = 2
                WORK(2)  = 2
            ELSE
               K = 0
            END IF
            INFO = 1
            RETURN
         END IF
         MLWORK = MAX(2,N)
         OLWORK = MAX(2,N)
         IMINWR = 1
         SELECT CASE ( WHTSVD )
         CASE (1)
             ! The following is specified as the minimal
             ! length of WORK in the definition of DGESVD:
             ! MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
             MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
             MLWORK = MAX(MLWORK,N + MWRSVD)
             IF ( LQUERY ) THEN
                CALL DGESVD( 'O', 'S', M, N, X, LDX, WORK, &
                           B, LDB, W, LDW, RDUMMY, -1, INFO1 )
                LWRSVD = MAX( MWRSVD, INT( RDUMMY(1) ) )
                OLWORK = MAX(OLWORK,N + LWRSVD)
             END IF
         CASE (2)
             ! The following is specified as the minimal
             ! length of WORK in the definition of DGESDD:
             ! MWRSDD = 3*MIN(M,N)*MIN(M,N) +
             ! MAX( MAX(M,N),5*MIN(M,N)*MIN(M,N)+4*MIN(M,N) )
             ! IMINWR = 8*MIN(M,N)
             MWRSDD = 3*MIN(M,N)*MIN(M,N) +                &
              MAX( MAX(M,N),5*MIN(M,N)*MIN(M,N)+4*MIN(M,N) )
             MLWORK = MAX(MLWORK,N + MWRSDD)
             IMINWR = 8*MIN(M,N)
             IF ( LQUERY ) THEN
                CALL DGESDD( 'O', M, N, X, LDX, WORK, B,     &
                     LDB, W, LDW, RDUMMY, -1, IWORK, INFO1 )
                LWRSDD = MAX( MWRSDD, INT( RDUMMY(1) ) )
                OLWORK = MAX(OLWORK,N + LWRSDD)
             END IF
         CASE (3)
             !LWQP3 = 3*N+1
             !LWORQ = MAX(N, 1)
             !MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
             !MWRSVQ = N + MAX( LWQP3, MWRSVD, LWORQ ) + MAX(M,2)
             !MLWORK = N +  MWRSVQ
             !IMINWR = M+N-1
             CALL DGESVDQ( 'H', 'P', 'N', 'R', 'R', M, N, &
                             X, LDX, WORK, Z, LDZ, W, LDW,   &
                             NUMRNK, IWORK, LIWORK, RDUMMY,  &
                             -1, RDUMMY2, -1, INFO1 )
             IMINWR = IWORK(1)
             MWRSVQ = INT(RDUMMY(2))
             MLWORK = MAX(MLWORK,N+MWRSVQ+INT(RDUMMY2(1)))
             IF ( LQUERY ) THEN
                LWRSVQ = MAX( MWRSVQ, INT(RDUMMY(1)) )
                OLWORK = MAX(OLWORK,N+LWRSVQ+INT(RDUMMY2(1)))
             END IF
         CASE (4)
             JSVOPT = 'J'
             !MWRSVJ = MAX( 7, 2*M+N, 6*N+2*N*N ) ! for JSVOPT='V'
             MWRSVJ = MAX( 7, 2*M+N, 4*N+N*N, 2*N+N*N+6 )
             MLWORK = MAX(MLWORK,N+MWRSVJ)
             IMINWR = MAX( 3, M+3*N )
             IF ( LQUERY ) THEN
                OLWORK =  MAX(OLWORK,N+MWRSVJ)
             END IF
         END SELECT
         IF ( WNTVEC .OR. WNTEX .OR. LSAME(JOBZ,'F') ) THEN
             JOBZL = 'V'
         ELSE
             JOBZL = 'N'
         END IF
         ! Workspace calculation to the DGEEV call
         IF ( LSAME(JOBZL,'V') ) THEN
             MWRKEV = MAX( 1, 4*N )
         ELSE
             MWRKEV = MAX( 1, 3*N )
         END IF
         MLWORK = MAX(MLWORK,N+MWRKEV)
         IF ( LQUERY ) THEN
                CALL DGEEV( 'N', JOBZL, N, S, LDS, REIG, &
                    IMEIG, W, LDW, W, LDW, RDUMMY, -1, INFO1 )
                LWRKEV = MAX( MWRKEV, INT(RDUMMY(1)) )
                OLWORK = MAX( OLWORK, N+LWRKEV )
         END IF
!
         IF ( LIWORK < IMINWR .AND. (.NOT.LQUERY) ) INFO = -29
         IF (  LWORK < MLWORK .AND. (.NOT.LQUERY) ) INFO = -27
      END IF
!
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'DGEDMD', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
!     Return minimal and optimal workspace sizes
          IWORK(1) = IMINWR
          WORK(1)  = MLWORK
          WORK(2)  = OLWORK
          RETURN
      END IF
!............................................................
!
      OFL   = DLAMCH('O')
      SMALL = DLAMCH('S')
      BADXY = .FALSE.
!
!     <1> Optional scaling of the snapshots (columns of X, Y)
!     ==========================================================
      IF ( SCCOLX ) THEN
          ! The columns of X will be normalized.
          ! To prevent overflows, the column norms of X are
          ! carefully computed using DLASSQ.
          K = 0
          DO i = 1, N
            !WORK(i) = DNRM2( M, X(1,i), 1 )
            SSUM  = ONE
            SCALE = ZERO
            CALL DLASSQ( M, X(1,i), 1, SCALE, SSUM )
            IF ( DISNAN(SCALE) .OR. DISNAN(SSUM) ) THEN
                K    =  0
                INFO = -8
                CALL XERBLA('DGEDMD',-INFO)
            END IF
            IF ( (SCALE /= ZERO) .AND. (SSUM /= ZERO) ) THEN
               ROOTSC = SQRT(SSUM)
               IF ( SCALE .GE. (OFL / ROOTSC) ) THEN
!                 Norm of X(:,i) overflows. First, X(:,i)
!                 is scaled by
!                 ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2.
!                 Next, the norm of X(:,i) is stored without
!                 overflow as WORK(i) = - SCALE * (ROOTSC/M),
!                 the minus sign indicating the 1/M factor.
!                 Scaling is performed without overflow, and
!                 underflow may occur in the smallest entries
!                 of X(:,i). The relative backward and forward
!                 errors are small in the ell_2 norm.
                  CALL DLASCL( 'G', 0, 0, SCALE, ONE/ROOTSC, &
                               M, 1, X(1,i), M, INFO2 )
                  WORK(i) = - SCALE * ( ROOTSC / DBLE(M) )
               ELSE
!                 X(:,i) will be scaled to unit 2-norm
                  WORK(i) =   SCALE * ROOTSC
                  CALL DLASCL( 'G',0, 0, WORK(i), ONE, M, 1, &
                               X(1,i), M, INFO2 )              ! LAPACK CALL
!                 X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)          ! INTRINSIC
               END IF
            ELSE
               WORK(i) = ZERO
               K = K + 1
            END IF
          END DO
          IF ( K == N ) THEN
          ! All columns of X are zero. Return error code -8.
          ! (the 8th input variable had an illegal value)
          K = 0
          INFO = -8
          CALL XERBLA('DGEDMD',-INFO)
          RETURN
          END IF
          DO i = 1, N
!           Now, apply the same scaling to the columns of Y.
            IF ( WORK(i) >  ZERO ) THEN
                CALL DSCAL( M, ONE/WORK(i), Y(1,i), 1 )  ! BLAS CALL
!               Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)      ! INTRINSIC
            ELSE IF ( WORK(i) < ZERO ) THEN
                CALL DLASCL( 'G', 0, 0, -WORK(i),          &
                     ONE/DBLE(M), M, 1, Y(1,i), M, INFO2 ) ! LAPACK CALL
            ELSE IF ( Y(IDAMAX(M, Y(1,i),1),i )  &
                                            /= ZERO ) THEN
!               X(:,i) is zero vector. For consistency,
!               Y(:,i) should also be zero. If Y(:,i) is not
!               zero, then the data might be inconsistent or
!               corrupted. If JOBS == 'C', Y(:,i) is set to
!               zero and a warning flag is raised.
!               The computation continues but the
!               situation will be reported in the output.
                BADXY = .TRUE.
                IF ( LSAME(JOBS,'C')) &
                CALL DSCAL( M, ZERO, Y(1,i), 1 )  ! BLAS CALL
            END IF
          END DO
      END IF
  !
      IF ( SCCOLY ) THEN
          ! The columns of Y will be normalized.
          ! To prevent overflows, the column norms of Y are
          ! carefully computed using DLASSQ.
          DO i = 1, N
            !WORK(i) = DNRM2( M, Y(1,i), 1 )
            SSUM  = ONE
            SCALE = ZERO
            CALL DLASSQ( M, Y(1,i), 1, SCALE, SSUM )
            IF ( DISNAN(SCALE) .OR. DISNAN(SSUM) ) THEN
                K    =  0
                INFO = -10
                CALL XERBLA('DGEDMD',-INFO)
            END IF
            IF ( SCALE /= ZERO  .AND. (SSUM /= ZERO) ) THEN
               ROOTSC = SQRT(SSUM)
               IF ( SCALE .GE. (OFL / ROOTSC) ) THEN
!                 Norm of Y(:,i) overflows. First, Y(:,i)
!                 is scaled by
!                 ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2.
!                 Next, the norm of Y(:,i) is stored without
!                 overflow as WORK(i) = - SCALE * (ROOTSC/M),
!                 the minus sign indicating the 1/M factor.
!                 Scaling is performed without overflow, and
!                 underflow may occur in the smallest entries
!                 of Y(:,i). The relative backward and forward
!                 errors are small in the ell_2 norm.
                  CALL DLASCL( 'G', 0, 0, SCALE, ONE/ROOTSC, &
                               M, 1, Y(1,i), M, INFO2 )
                  WORK(i) = - SCALE * ( ROOTSC / DBLE(M) )
               ELSE
!                 X(:,i) will be scaled to unit 2-norm
                  WORK(i) =   SCALE * ROOTSC
                  CALL DLASCL( 'G',0, 0, WORK(i), ONE, M, 1, &
                               Y(1,i), M, INFO2 )              ! LAPACK CALL
!                 Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)          ! INTRINSIC
               END IF
            ELSE
               WORK(i) = ZERO
            END IF
         END DO
         DO i = 1, N
!           Now, apply the same scaling to the columns of X.
            IF ( WORK(i) >  ZERO ) THEN
                CALL DSCAL( M, ONE/WORK(i), X(1,i), 1 )  ! BLAS CALL
!               X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)      ! INTRINSIC
            ELSE IF ( WORK(i) < ZERO ) THEN
                CALL DLASCL( 'G', 0, 0, -WORK(i),          &
                     ONE/DBLE(M), M, 1, X(1,i), M, INFO2 ) ! LAPACK CALL
            ELSE IF ( X(IDAMAX(M, X(1,i),1),i )  &
                                           /= ZERO ) THEN
!               Y(:,i) is zero vector.  If X(:,i) is not
!               zero, then a warning flag is raised.
!               The computation continues but the
!               situation will be reported in the output.
                BADXY = .TRUE.
            END IF
         END DO
       END IF
!
!     <2> SVD of the data snapshot matrix X.
!     =====================================
!     The left singular vectors are stored in the array X.
!     The right singular vectors are in the array W.
!     The array W will later on contain the eigenvectors
!     of a Rayleigh quotient.
      NUMRNK = N
      SELECT CASE ( WHTSVD )
         CASE (1)
             CALL DGESVD( 'O', 'S', M, N, X, LDX, WORK, B, &
                  LDB, W, LDW, WORK(N+1), LWORK-N, INFO1 ) ! LAPACK CALL
             T_OR_N = 'T'
         CASE (2)
            CALL DGESDD( 'O', M, N, X, LDX, WORK, B, LDB, W, &
                 LDW, WORK(N+1), LWORK-N, IWORK, INFO1 )   ! LAPACK CALL
            T_OR_N = 'T'
         CASE (3)
              CALL DGESVDQ( 'H', 'P', 'N', 'R', 'R', M, N, &
                   X, LDX, WORK, Z, LDZ, W, LDW, &
                   NUMRNK, IWORK, LIWORK, WORK(N+MAX(2,M)+1),&
                   LWORK-N-MAX(2,M), WORK(N+1), MAX(2,M), INFO1)     ! LAPACK CALL
              CALL DLACPY( 'A', M, NUMRNK, Z, LDZ, X, LDX )   ! LAPACK CALL
         T_OR_N = 'T'
         CASE (4)
              CALL DGEJSV( 'F', 'U', JSVOPT, 'N', 'N', 'P', M, &
                   N, X, LDX, WORK, Z, LDZ, W, LDW, &
                   WORK(N+1), LWORK-N, IWORK, INFO1 )    ! LAPACK CALL
              CALL DLACPY( 'A', M, N, Z, LDZ, X, LDX )   ! LAPACK CALL
              T_OR_N = 'N'
              XSCL1 = WORK(N+1)
              XSCL2 = WORK(N+2)
              IF ( XSCL1 /=  XSCL2 ) THEN
                 ! This is an exceptional situation. If the
                 ! data matrices are not scaled and the
                 ! largest singular value of X overflows.
                 ! In that case DGEJSV can return the SVD
                 ! in scaled form. The scaling factor can be used
                 ! to rescale the data (X and Y).
                 CALL DLASCL( 'G', 0, 0, XSCL1, XSCL2, M, N, Y, LDY, INFO2  )
              END IF
      END SELECT
!
      IF ( INFO1 > 0 ) THEN
         ! The SVD selected subroutine did not converge.
         ! Return with an error code.
         INFO = 2
         RETURN
      END IF
!
      IF ( WORK(1) == ZERO ) THEN
          ! The largest computed singular value of (scaled)
          ! X is zero. Return error code -8
          ! (the 8th input variable had an illegal value).
          K = 0
          INFO = -8
          CALL XERBLA('DGEDMD',-INFO)
          RETURN
      END IF
!
      !<3> Determine the numerical rank of the data
      !    snapshots matrix X. This depends on the
      !    parameters NRNK and TOL.

      SELECT CASE ( NRNK )
          CASE ( -1 )
               K = 1
               DO i = 2, NUMRNK
                 IF ( ( WORK(i) <= WORK(1)*TOL ) .OR. &
                      ( WORK(i) <= SMALL ) ) EXIT
                 K = K + 1
               END DO
          CASE ( -2 )
               K = 1
               DO i = 1, NUMRNK-1
                 IF ( ( WORK(i+1) <= WORK(i)*TOL  ) .OR. &
                      ( WORK(i) <= SMALL ) ) EXIT
                 K = K + 1
               END DO
          CASE DEFAULT
               K = 1
               DO i = 2, NRNK
                  IF ( WORK(i) <= SMALL ) EXIT
                  K = K + 1
               END DO
          END SELECT
      !   Now, U = X(1:M,1:K) is the SVD/POD basis for the
      !   snapshot data in the input matrix X.

      !<4> Compute the Rayleigh quotient S = U^T * A * U.
      !    Depending on the requested outputs, the computation
      !    is organized to compute additional auxiliary
      !    matrices (for the residuals and refinements).
      !
      !    In all formulas below, we need V_k*Sigma_k^(-1)
      !    where either V_k is in W(1:N,1:K), or V_k^T is in
      !    W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
      IF ( LSAME(T_OR_N, 'N') ) THEN
          DO i = 1, K
           CALL DSCAL( N, ONE/WORK(i), W(1,i), 1 )    ! BLAS CALL
           ! W(1:N,i) = (ONE/WORK(i)) * W(1:N,i)      ! INTRINSIC
          END DO
      ELSE
          ! This non-unit stride access is due to the fact
          ! that DGESVD, DGESVDQ and DGESDD return the
          ! transposed matrix of the right singular vectors.
          !DO i = 1, K
          ! CALL DSCAL( N, ONE/WORK(i), W(i,1), LDW )    ! BLAS CALL
          ! ! W(i,1:N) = (ONE/WORK(i)) * W(i,1:N)      ! INTRINSIC
          !END DO
          DO i = 1, K
              WORK(N+i) = ONE/WORK(i)
          END DO
          DO j = 1, N
             DO i = 1, K
                 W(i,j) = (WORK(N+i))*W(i,j)
             END DO
          END DO
      END IF
!
      IF ( WNTREF ) THEN
         !
         ! Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K)))
         ! for computing the refined Ritz vectors
         ! (optionally, outside DGEDMD).
          CALL DGEMM( 'N', T_OR_N, M, K, N, ONE, Y, LDY, W, &
                      LDW, ZERO, Z, LDZ )                        ! BLAS CALL
          ! Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(W(1:K,1:N)))  ! INTRINSIC, for T_OR_N=='T'
          ! Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))             ! INTRINSIC, for T_OR_N=='N'
          !
          ! At this point Z contains
          ! A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and
          ! this is needed for computing the residuals.
          ! This matrix is  returned in the array B and
          ! it can be used to compute refined Ritz vectors.
          CALL DLACPY( 'A', M, K, Z, LDZ, B, LDB )   ! BLAS CALL
          ! B(1:M,1:K) = Z(1:M,1:K)                  ! INTRINSIC

          CALL DGEMM( 'T', 'N', K, K, M, ONE, X, LDX, Z, &
                      LDZ, ZERO, S, LDS )                        ! BLAS CALL
          ! S(1:K,1:K) = MATMUL(TANSPOSE(X(1:M,1:K)),Z(1:M,1:K)) ! INTRINSIC
          ! At this point S = U^T * A * U is the Rayleigh quotient.
      ELSE
        ! A * U(:,1:K) is not explicitly needed and the
        ! computation is organized differently. The Rayleigh
        ! quotient is computed more efficiently.
        CALL DGEMM( 'T', 'N', K, N, M, ONE, X, LDX, Y, LDY, &
                   ZERO, Z, LDZ )                                   ! BLAS CALL
        ! Z(1:K,1:N) = MATMUL( TRANSPOSE(X(1:M,1:K)), Y(1:M,1:N) )  ! INTRINSIC
        ! In the two DGEMM calls here, can use K for LDZ.
        CALL DGEMM( 'N', T_OR_N, K, K, N, ONE, Z, LDZ, W, &
                    LDW, ZERO, S, LDS )                         ! BLAS CALL
        ! S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRINSIC, for T_OR_N=='T'
        ! S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))          ! INTRINSIC, for T_OR_N=='N'
        ! At this point S = U^T * A * U is the Rayleigh quotient.
        ! If the residuals are requested, save scaled V_k into Z.
        ! Recall that V_k or V_k^T is stored in W.
        IF ( WNTRES .OR. WNTEX ) THEN
          IF ( LSAME(T_OR_N, 'N') ) THEN
              CALL DLACPY( 'A', N, K, W, LDW, Z, LDZ )
          ELSE
              CALL DLACPY( 'A', K, N, W, LDW, Z, LDZ )
          END IF
        END IF
      END IF
!
      !<5> Compute the Ritz values and (if requested) the
      !   right eigenvectors of the Rayleigh quotient.
      !
      CALL DGEEV( 'N', JOBZL, K, S, LDS, REIG, IMEIG, W, &
                  LDW, W, LDW, WORK(N+1), LWORK-N, INFO1 )   ! LAPACK CALL
      !
      ! W(1:K,1:K) contains the eigenvectors of the Rayleigh
      ! quotient. Even in the case of complex spectrum, all
      ! computation is done in real arithmetic. REIG and
      ! IMEIG are the real and the imaginary parts of the
      ! eigenvalues, so that the spectrum is given as
      ! REIG(:) + sqrt(-1)*IMEIG(:). Complex conjugate pairs
      ! are listed at consecutive positions. For such a
      ! complex conjugate pair of the eigenvalues, the
      ! corresponding eigenvectors are also a complex
      ! conjugate pair with the real and imaginary parts
      ! stored column-wise in W at the corresponding
      ! consecutive column indices. See the description of Z.
      ! Also, see the description of DGEEV.
      IF ( INFO1 > 0 ) THEN
         ! DGEEV failed to compute the eigenvalues and
         ! eigenvectors of the Rayleigh quotient.
         INFO = 3
         RETURN
      END IF
!
      ! <6> Compute the eigenvectors (if requested) and,
      ! the residuals (if requested).
      !
      IF ( WNTVEC .OR. WNTEX ) THEN
      IF ( WNTRES ) THEN
          IF ( WNTREF ) THEN
            ! Here, if the refinement is requested, we have
            ! A*U(:,1:K) already computed and stored in Z.
            ! For the residuals, need Y = A * U(:,1;K) * W.
            CALL DGEMM( 'N', 'N', M, K, K, ONE, Z, LDZ, W, &
                       LDW, ZERO, Y, LDY )               ! BLAS CALL
            ! Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)       ! INTRINSIC
            ! This frees Z; Y contains A * U(:,1:K) * W.
          ELSE
            ! Compute S = V_k * Sigma_k^(-1) * W, where
            ! V_k * Sigma_k^(-1) is stored in Z
            CALL DGEMM( T_OR_N, 'N', N, K, K, ONE, Z, LDZ, &
                       W, LDW, ZERO, S, LDS)
            ! Then, compute Z = Y * S =
            ! = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            ! = A * U(:,1:K) * W(1:K,1:K)
            CALL DGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
                       LDS, ZERO, Z, LDZ)
            ! Save a copy of Z into Y and free Z for holding
            ! the Ritz vectors.
            CALL DLACPY( 'A', M, K, Z, LDZ, Y, LDY )
            IF ( WNTEX ) CALL DLACPY( 'A', M, K, Z, LDZ, B, LDB )
          END IF
      ELSE IF ( WNTEX ) THEN
          ! Compute S = V_k * Sigma_k^(-1) * W, where
            ! V_k * Sigma_k^(-1) is stored in Z
            CALL DGEMM( T_OR_N, 'N', N, K, K, ONE, Z, LDZ, &
                       W, LDW, ZERO, S, LDS )
            ! Then, compute Z = Y * S =
            ! = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            ! = A * U(:,1:K) * W(1:K,1:K)
            CALL DGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
                       LDS, ZERO, B, LDB )
            ! The above call replaces the following two calls
            ! that were used in the developing-testing phase.
            ! CALL DGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
            !           LDS, ZERO, Z, LDZ)
            ! Save a copy of Z into B and free Z for holding
            ! the Ritz vectors.
            ! CALL DLACPY( 'A', M, K, Z, LDZ, B, LDB )
      END IF
!
      ! Compute the real form of the Ritz vectors
      IF ( WNTVEC ) CALL DGEMM( 'N', 'N', M, K, K, ONE, X, LDX, W, LDW, &
                   ZERO, Z, LDZ )                           ! BLAS CALL
      ! Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))         ! INTRINSIC
!
      IF ( WNTRES ) THEN
         i = 1
         DO WHILE ( i <= K )
            IF ( IMEIG(i) == ZERO ) THEN
                ! have a real eigenvalue with real eigenvector
                CALL DAXPY( M, -REIG(i), Z(1,i), 1, Y(1,i), 1 )       ! BLAS CALL
                ! Y(1:M,i) = Y(1:M,i) - REIG(i) * Z(1:M,i)            ! INTRINSIC
                RES(i) = DNRM2( M, Y(1,i), 1)                         ! BLAS CALL
                i = i + 1
            ELSE
               ! Have a complex conjugate pair
               ! REIG(i) +- sqrt(-1)*IMEIG(i).
               ! Since all computation is done in real
               ! arithmetic, the formula for the residual
               ! is recast for real representation of the
               ! complex conjugate eigenpair. See the
               ! description of RES.
               AB(1,1) =  REIG(i)
               AB(2,1) = -IMEIG(i)
               AB(1,2) =  IMEIG(i)
               AB(2,2) =  REIG(i)
               CALL DGEMM( 'N', 'N', M, 2, 2, -ONE, Z(1,i), &
                           LDZ, AB, 2, ONE, Y(1,i), LDY )          ! BLAS CALL
               ! Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB   ! INTRINSIC
               RES(i)   = DLANGE( 'F', M, 2, Y(1,i), LDY, &
                                  WORK(N+1) )                      ! LAPACK CALL
               RES(i+1) = RES(i)
               i = i + 2
            END IF
         END DO
      END IF
      END IF
!
      IF ( WHTSVD == 4 ) THEN
          WORK(N+1) = XSCL1
          WORK(N+2) = XSCL2
      END IF
!
!     Successful exit.
      IF ( .NOT. BADXY ) THEN
         INFO = 0
      ELSE
         ! A warning on possible data inconsistency.
         ! This should be a rare event.
         INFO = 4
      END IF
!............................................................
      RETURN
!     ......
      END SUBROUTINE DGEDMD