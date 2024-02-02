      SUBROUTINE ZGEDMD( JOBS, JOBZ, JOBR, JOBF,  WHTSVD,   &
                         M, N, X, LDX, Y, LDY, NRNK, TOL,   &
                         K, EIGS, Z, LDZ, RES, B,    LDB,   &
                         W, LDW,  S, LDS, ZWORK,  LZWORK,   &
                         RWORK, LRWORK, IWORK, LIWORK, INFO )
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
                                 LIWORK, LRWORK, LZWORK
      INTEGER,       INTENT(OUT)  :: K, INFO
      REAL(KIND=WP), INTENT(IN)   ::    TOL
!
!     Array arguments
!     ~~~~~~~~~~~~~~~
      COMPLEX(KIND=WP), INTENT(INOUT) :: X(LDX,*), Y(LDY,*)
      COMPLEX(KIND=WP), INTENT(OUT)   :: Z(LDZ,*), B(LDB,*), &
                                         W(LDW,*), S(LDS,*)
      COMPLEX(KIND=WP), INTENT(OUT)   :: EIGS(*)
      COMPLEX(KIND=WP), INTENT(OUT)   :: ZWORK(*)
      REAL(KIND=WP),    INTENT(OUT)   :: RES(*)
      REAL(KIND=WP),    INTENT(OUT)   :: RWORK(*)
      INTEGER,          INTENT(OUT)   :: IWORK(*)
!
!     Parameters
!     ~~~~~~~~~~
      REAL(KIND=WP),    PARAMETER ::  ONE = 1.0_WP
      REAL(KIND=WP),    PARAMETER :: ZERO = 0.0_WP
      COMPLEX(KIND=WP), PARAMETER ::  ZONE = ( 1.0_WP, 0.0_WP )
      COMPLEX(KIND=WP), PARAMETER :: ZZERO = ( 0.0_WP, 0.0_WP )
!
!     Local scalars
!     ~~~~~~~~~~~~~
      REAL(KIND=WP) :: OFL,   ROOTSC, SCALE,  SMALL,    &
                       SSUM,  XSCL1,  XSCL2
      INTEGER       ::  i,  j,  IMINWR,  INFO1, INFO2,  &
                        LWRKEV, LWRSDD, LWRSVD, LWRSVJ, &
                        LWRSVQ, MLWORK, MWRKEV, MWRSDD, &
                        MWRSVD, MWRSVJ, MWRSVQ, NUMRNK, &
                        OLWORK, MLRWRK
      LOGICAL       ::  BADXY, LQUERY, SCCOLX, SCCOLY,  &
                        WNTEX, WNTREF, WNTRES, WNTVEC
      CHARACTER     ::  JOBZL, T_OR_N
      CHARACTER     ::  JSVOPT
!
!     Local arrays
!     ~~~~~~~~~~~~
      REAL(KIND=WP) :: RDUMMY(2)
!
!     External functions (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~
      REAL(KIND=WP) ZLANGE, DLAMCH, DZNRM2
      EXTERNAL      ZLANGE, DLAMCH, DZNRM2, IZAMAX
      INTEGER                               IZAMAX
      LOGICAL       DISNAN, LSAME
      EXTERNAL      DISNAN, LSAME
!
!     External subroutines (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~~~~
      EXTERNAL      ZAXPY,  ZGEMM,  ZDSCAL
      EXTERNAL      ZGEEV,  ZGEJSV, ZGESDD, ZGESVD, ZGESVDQ, &
                    ZLACPY, ZLASCL, ZLASSQ, XERBLA
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
      LQUERY = ( ( LZWORK == -1 ) .OR. ( LIWORK == -1 ) &
                                  .OR. ( LRWORK == -1 ) )
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
          INFO = -17
      ELSE IF ( (WNTREF .OR. WNTEX ) .AND. ( LDB < M ) ) THEN
          INFO = -20
      ELSE IF ( LDW < N ) THEN
          INFO = -22
      ELSE IF ( LDS < N ) THEN
          INFO = -24
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
                RWORK(1) = 1
                ZWORK(1) = 2
                ZWORK(2) = 2
            ELSE
               K   =  0
            END IF
            INFO = 1
            RETURN
         END IF

         IMINWR = 1
         MLRWRK = MAX(1,N)
         MLWORK = 2
         OLWORK = 2
         SELECT CASE ( WHTSVD )
         CASE (1)
             ! The following is specified as the minimal
             ! length of WORK in the definition of ZGESVD:
             ! MWRSVD = MAX(1,2*MIN(M,N)+MAX(M,N))
             MWRSVD = MAX(1,2*MIN(M,N)+MAX(M,N))
             MLWORK = MAX(MLWORK,MWRSVD)
             MLRWRK = MAX(MLRWRK,N + 5*MIN(M,N))
             IF ( LQUERY ) THEN
                CALL ZGESVD( 'O', 'S', M, N, X, LDX, RWORK, &
                     B, LDB, W, LDW, ZWORK, -1, RDUMMY, INFO1 )
                LWRSVD = INT( ZWORK(1) )
                OLWORK = MAX(OLWORK,LWRSVD)
             END IF
         CASE (2)
             ! The following is specified as the minimal
             ! length of WORK in the definition of ZGESDD:
             ! MWRSDD = 2*min(M,N)*min(M,N)+2*min(M,N)+max(M,N).
             ! RWORK length: 5*MIN(M,N)*MIN(M,N)+7*MIN(M,N)
             ! In LAPACK 3.10.1 RWORK is defined differently.
             ! Below we take max over the two versions.
             ! IMINWR = 8*MIN(M,N)
             MWRSDD = 2*MIN(M,N)*MIN(M,N)+2*MIN(M,N)+MAX(M,N)
             MLWORK = MAX(MLWORK,MWRSDD)
             IMINWR = 8*MIN(M,N)
             MLRWRK = MAX( MLRWRK,  N +                    &
                      MAX( 5*MIN(M,N)*MIN(M,N)+7*MIN(M,N), &
                           5*MIN(M,N)*MIN(M,N)+5*MIN(M,N), &
                           2*MAX(M,N)*MIN(M,N)+            &
                           2*MIN(M,N)*MIN(M,N)+MIN(M,N) ) )
             IF ( LQUERY ) THEN
                CALL ZGESDD( 'O', M, N, X, LDX, RWORK, B,LDB,&
                     W, LDW, ZWORK, -1, RDUMMY, IWORK, INFO1 )
                LWRSDD = MAX( MWRSDD,INT( ZWORK(1) ))
                ! Possible bug in ZGESDD optimal workspace size.
                OLWORK = MAX(OLWORK,LWRSDD)
             END IF
         CASE (3)
             CALL ZGESVDQ( 'H', 'P', 'N', 'R', 'R', M, N, &
                  X, LDX, RWORK, Z, LDZ, W, LDW, NUMRNK,  &
                  IWORK, -1, ZWORK, -1, RDUMMY, -1, INFO1 )
             IMINWR = IWORK(1)
             MWRSVQ = INT(ZWORK(2))
             MLWORK = MAX(MLWORK,MWRSVQ)
             MLRWRK  = MAX(MLRWRK,N + INT(RDUMMY(1)))
             IF ( LQUERY ) THEN
                LWRSVQ = INT(ZWORK(1))
                OLWORK = MAX(OLWORK,LWRSVQ)
             END IF
         CASE (4)
             JSVOPT = 'J'
             CALL ZGEJSV( 'F', 'U', JSVOPT, 'R', 'N', 'P', M, &
                   N, X, LDX, RWORK, Z, LDZ, W, LDW,       &
                   ZWORK, -1, RDUMMY, -1, IWORK, INFO1 )
             IMINWR = IWORK(1)
             MWRSVJ = INT(ZWORK(2))
             MLWORK = MAX(MLWORK,MWRSVJ)
             MLRWRK = MAX(MLRWRK,N + MAX(7,INT(RDUMMY(1))))
             IF ( LQUERY ) THEN
                LWRSVJ = INT(ZWORK(1))
                OLWORK = MAX(OLWORK,LWRSVJ)
             END IF
         END SELECT
         IF ( WNTVEC .OR. WNTEX .OR. LSAME(JOBZ,'F') ) THEN
             JOBZL = 'V'
         ELSE
             JOBZL = 'N'
         END IF
         ! Workspace calculation to the ZGEEV call
         MWRKEV = MAX( 1, 2*N )
         MLWORK = MAX(MLWORK,MWRKEV)
         MLRWRK = MAX(MLRWRK,N+2*N)
         IF ( LQUERY ) THEN
             CALL ZGEEV( 'N', JOBZL, N, S, LDS, EIGS, &
              W, LDW, W, LDW, ZWORK, -1, RWORK, INFO1 )
                LWRKEV = INT(ZWORK(1))
                OLWORK = MAX( OLWORK, LWRKEV )
         END IF
!
         IF ( LIWORK < IMINWR .AND. (.NOT.LQUERY) ) INFO = -30
         IF ( LRWORK < MLRWRK .AND. (.NOT.LQUERY) ) INFO = -28
         IF ( LZWORK < MLWORK .AND. (.NOT.LQUERY) ) INFO = -26

      END IF
!
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'ZGEDMD', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
!     Return minimal and optimal workspace sizes
          IWORK(1) = IMINWR
          RWORK(1) = MLRWRK
          ZWORK(1) = MLWORK
          ZWORK(2) = OLWORK
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
          ! carefully computed using ZLASSQ.
          K = 0
          DO i = 1, N
            !WORK(i) = DZNRM2( M, X(1,i), 1 )
            SSUM  = ONE
            SCALE = ZERO
            CALL ZLASSQ( M, X(1,i), 1, SCALE, SSUM )
            IF ( DISNAN(SCALE) .OR. DISNAN(SSUM) ) THEN
                K    =  0
                INFO = -8
                CALL XERBLA('ZGEDMD',-INFO)
            END IF
            IF ( (SCALE /= ZERO) .AND. (SSUM /= ZERO) ) THEN
               ROOTSC = SQRT(SSUM)
               IF ( SCALE .GE. (OFL / ROOTSC) ) THEN
!                 Norm of X(:,i) overflows. First, X(:,i)
!                 is scaled by
!                 ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2.
!                 Next, the norm of X(:,i) is stored without
!                 overflow as RWORK(i) = - SCALE * (ROOTSC/M),
!                 the minus sign indicating the 1/M factor.
!                 Scaling is performed without overflow, and
!                 underflow may occur in the smallest entries
!                 of X(:,i). The relative backward and forward
!                 errors are small in the ell_2 norm.
                  CALL ZLASCL( 'G', 0, 0, SCALE, ONE/ROOTSC, &
                               M, 1, X(1,i), LDX, INFO2 )
                  RWORK(i) = - SCALE * ( ROOTSC / DBLE(M) )
               ELSE
!                 X(:,i) will be scaled to unit 2-norm
                  RWORK(i) =   SCALE * ROOTSC
                  CALL ZLASCL( 'G',0, 0, RWORK(i), ONE, M, 1, &
                               X(1,i), LDX, INFO2 )      ! LAPACK CALL
!                 X(1:M,i) = (ONE/RWORK(i)) * X(1:M,i)   ! INTRINSIC
               END IF
            ELSE
               RWORK(i) = ZERO
               K = K + 1
            END IF
          END DO
          IF ( K == N ) THEN
             ! All columns of X are zero. Return error code -8.
             ! (the 8th input variable had an illegal value)
             K = 0
             INFO = -8
             CALL XERBLA('ZGEDMD',-INFO)
             RETURN
          END IF
          DO i = 1, N
!           Now, apply the same scaling to the columns of Y.
            IF ( RWORK(i) >  ZERO ) THEN
                CALL ZDSCAL( M, ONE/RWORK(i), Y(1,i), 1 )  ! BLAS CALL
!               Y(1:M,i) = (ONE/RWORK(i)) * Y(1:M,i)       ! INTRINSIC
            ELSE IF ( RWORK(i) < ZERO ) THEN
                CALL ZLASCL( 'G', 0, 0, -RWORK(i),          &
                     ONE/DBLE(M), M, 1, Y(1,i), LDY, INFO2 ) ! LAPACK CALL
            ELSE IF ( ABS(Y(IZAMAX(M, Y(1,i),1),i ))  &
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
                CALL ZDSCAL( M, ZERO, Y(1,i), 1 )  ! BLAS CALL
            END IF
          END DO
      END IF
  !
      IF ( SCCOLY ) THEN
          ! The columns of Y will be normalized.
          ! To prevent overflows, the column norms of Y are
          ! carefully computed using ZLASSQ.
          DO i = 1, N
            !RWORK(i) = DZNRM2( M, Y(1,i), 1 )
            SSUM  = ONE
            SCALE = ZERO
            CALL ZLASSQ( M, Y(1,i), 1, SCALE, SSUM )
            IF ( DISNAN(SCALE) .OR. DISNAN(SSUM) ) THEN
                K    =  0
                INFO = -10
                CALL XERBLA('ZGEDMD',-INFO)
            END IF
            IF ( SCALE /= ZERO  .AND. (SSUM /= ZERO) ) THEN
               ROOTSC = SQRT(SSUM)
               IF ( SCALE .GE. (OFL / ROOTSC) ) THEN
!                 Norm of Y(:,i) overflows. First, Y(:,i)
!                 is scaled by
!                 ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2.
!                 Next, the norm of Y(:,i) is stored without
!                 overflow as RWORK(i) = - SCALE * (ROOTSC/M),
!                 the minus sign indicating the 1/M factor.
!                 Scaling is performed without overflow, and
!                 underflow may occur in the smallest entries
!                 of Y(:,i). The relative backward and forward
!                 errors are small in the ell_2 norm.
                  CALL ZLASCL( 'G', 0, 0, SCALE, ONE/ROOTSC, &
                               M, 1, Y(1,i), LDY, INFO2 )
                  RWORK(i) = - SCALE * ( ROOTSC / DBLE(M) )
               ELSE
!                 Y(:,i) will be scaled to unit 2-norm
                  RWORK(i) =   SCALE * ROOTSC
                  CALL ZLASCL( 'G',0, 0, RWORK(i), ONE, M, 1, &
                               Y(1,i), LDY, INFO2 )             ! LAPACK CALL
!                 Y(1:M,i) = (ONE/RWORK(i)) * Y(1:M,i)          ! INTRINSIC
               END IF
            ELSE
               RWORK(i) = ZERO
            END IF
         END DO
         DO i = 1, N
!           Now, apply the same scaling to the columns of X.
            IF ( RWORK(i) >  ZERO ) THEN
                CALL ZDSCAL( M, ONE/RWORK(i), X(1,i), 1 ) ! BLAS CALL
!               X(1:M,i) = (ONE/RWORK(i)) * X(1:M,i)      ! INTRINSIC
            ELSE IF ( RWORK(i) < ZERO ) THEN
                CALL ZLASCL( 'G', 0, 0, -RWORK(i),          &
                     ONE/DBLE(M), M, 1, X(1,i), LDX, INFO2 ) ! LAPACK CALL
            ELSE IF ( ABS(X(IZAMAX(M, X(1,i),1),i ))  &
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
             CALL ZGESVD( 'O', 'S', M, N, X, LDX, RWORK, B, &
                  LDB, W, LDW, ZWORK, LZWORK,  RWORK(N+1), INFO1 ) ! LAPACK CALL
             T_OR_N = 'C'
         CASE (2)
            CALL ZGESDD( 'O', M, N, X, LDX, RWORK, B, LDB, W, &
                 LDW, ZWORK, LZWORK, RWORK(N+1), IWORK, INFO1 )   ! LAPACK CALL
            T_OR_N = 'C'
         CASE (3)
              CALL ZGESVDQ( 'H', 'P', 'N', 'R', 'R', M, N, &
                   X, LDX, RWORK, Z, LDZ, W, LDW, &
                   NUMRNK, IWORK, LIWORK, ZWORK,     &
                   LZWORK, RWORK(N+1), LRWORK-N, INFO1)     ! LAPACK CALL
              CALL ZLACPY( 'A', M, NUMRNK, Z, LDZ, X, LDX )   ! LAPACK CALL
         T_OR_N = 'C'
         CASE (4)
              CALL ZGEJSV( 'F', 'U', JSVOPT, 'R', 'N', 'P', M, &
                   N, X, LDX, RWORK, Z, LDZ, W, LDW, &
                   ZWORK, LZWORK, RWORK(N+1), LRWORK-N, IWORK, INFO1 )    ! LAPACK CALL
              CALL ZLACPY( 'A', M, N, Z, LDZ, X, LDX )   ! LAPACK CALL
              T_OR_N = 'N'
              XSCL1 = RWORK(N+1)
              XSCL2 = RWORK(N+2)
              IF ( XSCL1 /=  XSCL2 ) THEN
                 ! This is an exceptional situation. If the
                 ! data matrices are not scaled and the
                 ! largest singular value of X overflows.
                 ! In that case ZGEJSV can return the SVD
                 ! in scaled form. The scaling factor can be used
                 ! to rescale the data (X and Y).
                 CALL ZLASCL( 'G', 0, 0, XSCL1, XSCL2, M, N, Y, LDY, INFO2  )
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
      IF ( RWORK(1) == ZERO ) THEN
          ! The largest computed singular value of (scaled)
          ! X is zero. Return error code -8
          ! (the 8th input variable had an illegal value).
          K = 0
          INFO = -8
          CALL XERBLA('ZGEDMD',-INFO)
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
                 IF ( ( RWORK(i) <= RWORK(1)*TOL ) .OR. &
                      ( RWORK(i) <= SMALL ) ) EXIT
                 K = K + 1
               END DO
          CASE ( -2 )
               K = 1
               DO i = 1, NUMRNK-1
                 IF ( ( RWORK(i+1) <= RWORK(i)*TOL  ) .OR. &
                      ( RWORK(i) <= SMALL ) ) EXIT
                 K = K + 1
               END DO
          CASE DEFAULT
               K = 1
               DO i = 2, NRNK
                  IF ( RWORK(i) <= SMALL ) EXIT
                  K = K + 1
               END DO
          END SELECT
      !   Now, U = X(1:M,1:K) is the SVD/POD basis for the
      !   snapshot data in the input matrix X.

      !<4> Compute the Rayleigh quotient S = U^H * A * U.
      !    Depending on the requested outputs, the computation
      !    is organized to compute additional auxiliary
      !    matrices (for the residuals and refinements).
      !
      !    In all formulas below, we need V_k*Sigma_k^(-1)
      !    where either V_k is in W(1:N,1:K), or V_k^H is in
      !    W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
      IF ( LSAME(T_OR_N, 'N') ) THEN
          DO i = 1, K
           CALL ZDSCAL( N, ONE/RWORK(i), W(1,i), 1 )    ! BLAS CALL
           ! W(1:N,i) = (ONE/RWORK(i)) * W(1:N,i)      ! INTRINSIC
          END DO
      ELSE
          ! This non-unit stride access is due to the fact
          ! that ZGESVD, ZGESVDQ and ZGESDD return the
          ! adjoint matrix of the right singular vectors.
          !DO i = 1, K
          ! CALL ZDSCAL( N, ONE/RWORK(i), W(i,1), LDW )    ! BLAS CALL
          ! ! W(i,1:N) = (ONE/RWORK(i)) * W(i,1:N)      ! INTRINSIC
          !END DO
          DO i = 1, K
              RWORK(N+i) = ONE/RWORK(i)
          END DO
          DO j = 1, N
             DO i = 1, K
                 W(i,j) = CMPLX(RWORK(N+i),ZERO,KIND=WP)*W(i,j)
             END DO
          END DO
      END IF
!
      IF ( WNTREF ) THEN
         !
         ! Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K)))
         ! for computing the refined Ritz vectors
         ! (optionally, outside ZGEDMD).
          CALL ZGEMM( 'N', T_OR_N, M, K, N, ZONE, Y, LDY, W, &
                      LDW, ZZERO, Z, LDZ )                       ! BLAS CALL
          ! Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(CONJG(W(1:K,1:N)))) ! INTRINSIC, for T_OR_N=='C'
          ! Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))                   ! INTRINSIC, for T_OR_N=='N'
          !
          ! At this point Z contains
          ! A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and
          ! this is needed for computing the residuals.
          ! This matrix is  returned in the array B and
          ! it can be used to compute refined Ritz vectors.
          CALL ZLACPY( 'A', M, K, Z, LDZ, B, LDB )   ! BLAS CALL
          ! B(1:M,1:K) = Z(1:M,1:K)                  ! INTRINSIC

          CALL ZGEMM( 'C', 'N', K, K, M, ZONE, X, LDX, Z, &
                      LDZ, ZZERO, S, LDS )                        ! BLAS CALL
          ! S(1:K,1:K) = MATMUL(TRANSPOSE(CONJG(X(1:M,1:K))),Z(1:M,1:K)) ! INTRINSIC
          ! At this point S = U^H * A * U is the Rayleigh quotient.
      ELSE
        ! A * U(:,1:K) is not explicitly needed and the
        ! computation is organized differently. The Rayleigh
        ! quotient is computed more efficiently.
        CALL ZGEMM( 'C', 'N', K, N, M, ZONE, X, LDX, Y, LDY, &
                   ZZERO, Z, LDZ )                                         ! BLAS CALL
        ! Z(1:K,1:N) = MATMUL( TRANSPOSE(CONJG(X(1:M,1:K))), Y(1:M,1:N) )  ! INTRINSIC
        !
        CALL ZGEMM( 'N', T_OR_N, K, K, N, ZONE, Z, LDZ, W, &
                    LDW, ZZERO, S, LDS )                         ! BLAS CALL
        ! S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(CONJG(W(1:K,1:N)))) ! INTRINSIC, for T_OR_N=='T'
        ! S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))                 ! INTRINSIC, for T_OR_N=='N'
        ! At this point S = U^H * A * U is the Rayleigh quotient.
        ! If the residuals are requested, save scaled V_k into Z.
        ! Recall that V_k or V_k^H is stored in W.
        IF ( WNTRES .OR. WNTEX ) THEN
          IF ( LSAME(T_OR_N, 'N') ) THEN
              CALL ZLACPY( 'A', N, K, W, LDW, Z, LDZ )
          ELSE
              CALL ZLACPY( 'A', K, N, W, LDW, Z, LDZ )
          END IF
        END IF
      END IF
!
      !<5> Compute the Ritz values and (if requested) the
      !   right eigenvectors of the Rayleigh quotient.
      !
      CALL ZGEEV( 'N', JOBZL, K, S, LDS, EIGS, W, LDW, &
            W, LDW, ZWORK, LZWORK, RWORK(N+1), INFO1 )  ! LAPACK CALL
      !
      ! W(1:K,1:K) contains the eigenvectors of the Rayleigh
      ! quotient.  See the description of Z.
      ! Also, see the description of ZGEEV.
      IF ( INFO1 > 0 ) THEN
         ! ZGEEV failed to compute the eigenvalues and
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
            CALL ZGEMM( 'N', 'N', M, K, K, ZONE, Z, LDZ, W, &
                       LDW, ZZERO, Y, LDY )               ! BLAS CALL
            ! Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)        ! INTRINSIC
            ! This frees Z; Y contains A * U(:,1:K) * W.
          ELSE
            ! Compute S = V_k * Sigma_k^(-1) * W, where
            ! V_k * Sigma_k^(-1) (or its adjoint) is stored in Z
            CALL ZGEMM( T_OR_N, 'N', N, K, K, ZONE, Z, LDZ, &
                       W, LDW, ZZERO, S, LDS )
            ! Then, compute Z = Y * S =
            ! = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            ! = A * U(:,1:K) * W(1:K,1:K)
            CALL ZGEMM( 'N', 'N', M, K, N, ZONE, Y, LDY, S, &
                       LDS, ZZERO, Z, LDZ )
            ! Save a copy of Z into Y and free Z for holding
            ! the Ritz vectors.
            CALL ZLACPY( 'A', M, K, Z, LDZ, Y, LDY )
            IF ( WNTEX ) CALL ZLACPY( 'A', M, K, Z, LDZ, B, LDB )
          END IF
      ELSE IF ( WNTEX ) THEN
          ! Compute S = V_k * Sigma_k^(-1) * W, where
            ! V_k * Sigma_k^(-1) is stored in Z
            CALL ZGEMM( T_OR_N, 'N', N, K, K, ZONE, Z, LDZ, &
                       W, LDW, ZZERO, S, LDS )
            ! Then, compute Z = Y * S =
            ! = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            ! = A * U(:,1:K) * W(1:K,1:K)
            CALL ZGEMM( 'N', 'N', M, K, N, ZONE, Y, LDY, S, &
                       LDS, ZZERO, B, LDB )
            ! The above call replaces the following two calls
            ! that were used in the developing-testing phase.
            ! CALL ZGEMM( 'N', 'N', M, K, N, ZONE, Y, LDY, S, &
            !           LDS, ZZERO, Z, LDZ)
            ! Save a copy of Z into B and free Z for holding
            ! the Ritz vectors.
            ! CALL ZLACPY( 'A', M, K, Z, LDZ, B, LDB )
      END IF
!
      ! Compute the Ritz vectors
      IF ( WNTVEC ) CALL ZGEMM( 'N', 'N', M, K, K, ZONE, X, LDX, W, LDW, &
                   ZZERO, Z, LDZ )                          ! BLAS CALL
      ! Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))         ! INTRINSIC
!
      IF ( WNTRES ) THEN
         DO i = 1, K
            CALL ZAXPY( M, -EIGS(i), Z(1,i), 1, Y(1,i), 1 )       ! BLAS CALL
            ! Y(1:M,i) = Y(1:M,i) - EIGS(i) * Z(1:M,i)            ! INTRINSIC
            RES(i) = DZNRM2( M, Y(1,i), 1 )                       ! BLAS CALL
         END DO
      END IF
      END IF
!
      IF ( WHTSVD == 4 ) THEN
          RWORK(N+1) = XSCL1
          RWORK(N+2) = XSCL2
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
      END SUBROUTINE ZGEDMD