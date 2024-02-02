                    WHTSVD,   M, N, F, LDF,  X, LDX,  Y,   &
                    LDY,   NRNK,  TOL,   K,  REIG, IMEIG,  &
                    Z, LDZ, RES,  B,     LDB,   V, LDV,    &
                    S, LDS, WORK, LWORK, IWORK, LIWORK, INFO )
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
      INTEGER, PARAMETER :: WP = real32
!
!     Scalar arguments
!     ~~~~~~~~~~~~~~~~
      CHARACTER, INTENT(IN)  :: JOBS, JOBZ, JOBR, JOBQ,    &
                                JOBT, JOBF
      INTEGER,   INTENT(IN)  :: WHTSVD, M, N,   LDF, LDX,  &
                                LDY, NRNK, LDZ, LDB, LDV,  &
                                LDS, LWORK,  LIWORK
      INTEGER,   INTENT(OUT) :: INFO,   K
      REAL(KIND=WP), INTENT(IN)    ::   TOL
!
!     Array arguments
!     ~~~~~~~~~~~~~~~~
      REAL(KIND=WP), INTENT(INOUT) :: F(LDF,*)
      REAL(KIND=WP), INTENT(OUT)   :: X(LDX,*), Y(LDY,*),  &
                                      Z(LDZ,*), B(LDB,*),  &
                                      V(LDV,*), S(LDS,*)
      REAL(KIND=WP), INTENT(OUT)   :: REIG(*),  IMEIG(*),  &
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
      INTEGER           :: IMINWR, INFO1,  MLWDMD, MLWGQR, &
                           MLWMQR, MLWORK, MLWQR,  MINMN,  &
                           OLWDMD, OLWGQR, OLWMQR, OLWORK, &
                           OLWQR
      LOGICAL           :: LQUERY, SCCOLX, SCCOLY, WANTQ,  &
                           WNTTRF, WNTRES, WNTVEC, WNTVCF, &
                           WNTVCQ, WNTREF, WNTEX
      CHARACTER(LEN=1)  :: JOBVL
!
!     Local array
!     ~~~~~~~~~~~
      REAL(KIND=WP) :: RDUMMY(2)
!
!     External functions (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~
      LOGICAL       LSAME
      EXTERNAL      LSAME
!
!     External subroutines (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~~~~
      EXTERNAL      SGEMM
      EXTERNAL      SGEDMD, SGEQRF, SLACPY, SLASET, SORGQR, &
                    SORMQR, XERBLA
!
!     Intrinsic functions
!     ~~~~~~~~~~~~~~~~~~~
      INTRINSIC      MAX, MIN, INT
!..........................................................
!
!     Test the input arguments
      WNTRES = LSAME(JOBR,'R')
      SCCOLX = LSAME(JOBS,'S') .OR. LSAME( JOBS, 'C' )
      SCCOLY = LSAME(JOBS,'Y')
      WNTVEC = LSAME(JOBZ,'V')
      WNTVCF = LSAME(JOBZ,'F')
      WNTVCQ = LSAME(JOBZ,'Q')
      WNTREF = LSAME(JOBF,'R')
      WNTEX  = LSAME(JOBF,'E')
      WANTQ  = LSAME(JOBQ,'Q')
      WNTTRF = LSAME(JOBT,'R')
      MINMN  = MIN(M,N)
      INFO = 0
      LQUERY = ( ( LWORK == -1 ) .OR. ( LIWORK == -1 ) )
!
      IF ( .NOT. (SCCOLX .OR. SCCOLY .OR. LSAME(JOBS,'N')) )  THEN
          INFO = -1
      ELSE IF ( .NOT. (WNTVEC .OR. WNTVCF .OR. WNTVCQ       &
                              .OR. LSAME(JOBZ,'N')) ) THEN
          INFO = -2
      ELSE IF ( .NOT. (WNTRES .OR. LSAME(JOBR,'N')) .OR.    &
          ( WNTRES .AND. LSAME(JOBZ,'N') ) ) THEN
          INFO = -3
      ELSE IF ( .NOT. (WANTQ .OR. LSAME(JOBQ,'N')) ) THEN
          INFO = -4
      ELSE IF ( .NOT. ( WNTTRF .OR. LSAME(JOBT,'N') ) )  THEN
          INFO = -5
      ELSE IF ( .NOT. (WNTREF .OR. WNTEX .OR.             &
                LSAME(JOBF,'N') ) )                    THEN
          INFO = -6
      ELSE IF ( .NOT. ((WHTSVD == 1).OR.(WHTSVD == 2).OR.   &
                       (WHTSVD == 3).OR.(WHTSVD == 4)) ) THEN
          INFO = -7
      ELSE IF ( M < 0 ) THEN
          INFO = -8
      ELSE IF ( ( N < 0 ) .OR. ( N > M+1 ) ) THEN
          INFO = -9
      ELSE IF ( LDF < M ) THEN
          INFO = -11
      ELSE IF ( LDX < MINMN ) THEN
          INFO = -13
      ELSE IF ( LDY < MINMN ) THEN
          INFO = -15
      ELSE IF ( .NOT. (( NRNK == -2).OR.(NRNK == -1).OR.    &
                       ((NRNK >= 1).AND.(NRNK <=N ))) )  THEN
          INFO = -16
      ELSE IF ( ( TOL < ZERO ) .OR. ( TOL >= ONE ) ) THEN
          INFO = -17
      ELSE IF ( LDZ < M ) THEN
          INFO = -22
      ELSE IF ( (WNTREF.OR.WNTEX ).AND.( LDB < MINMN ) ) THEN
          INFO = -25
      ELSE IF ( LDV < N-1 ) THEN
          INFO = -27
      ELSE IF ( LDS < N-1 ) THEN
          INFO = -29
      END IF
!
      IF ( WNTVEC .OR. WNTVCF ) THEN
          JOBVL = 'V'
      ELSE
          JOBVL = 'N'
      END IF
      IF ( INFO == 0 ) THEN
          ! Compute the minimal and the optimal workspace
          ! requirements. Simulate running the code and
          ! determine minimal and optimal sizes of the
          ! workspace at any moment of the run.
         IF ( ( N == 0 ) .OR. ( N == 1 ) ) THEN
             ! All output except K is void. INFO=1 signals
             ! the void input. In case of a workspace query,
             ! the minimal workspace lengths are returned.
            IF ( LQUERY ) THEN
               IWORK(1) = 1
                WORK(1) = 2
                WORK(2) = 2
            ELSE
               K = 0
            END IF
            INFO = 1
            RETURN
         END IF
         MLWQR  = MAX(1,N)  ! Minimal workspace length for SGEQRF.
         MLWORK = MIN(M,N) + MLWQR
         IF ( LQUERY ) THEN
             CALL SGEQRF( M, N, F, LDF, WORK, RDUMMY, -1, &
                          INFO1 )
             OLWQR = INT(RDUMMY(1))
             OLWORK = MIN(M,N) + OLWQR
         END IF
         CALL SGEDMD( JOBS, JOBVL, JOBR, JOBF, WHTSVD, MINMN,&
                      N-1, X, LDX, Y, LDY, NRNK, TOL, K,     &
                      REIG, IMEIG, Z, LDZ, RES,  B, LDB,     &
                      V, LDV, S, LDS, WORK, -1, IWORK,       &
                      LIWORK, INFO1 )
         MLWDMD = INT(WORK(1))
         MLWORK = MAX(MLWORK, MINMN + MLWDMD)
         IMINWR = IWORK(1)
         IF ( LQUERY ) THEN
             OLWDMD = INT(WORK(2))
             OLWORK = MAX(OLWORK, MINMN+OLWDMD)
         END IF
         IF ( WNTVEC .OR. WNTVCF ) THEN
            MLWMQR = MAX(1,N)
            MLWORK = MAX(MLWORK,MINMN+N-1+MLWMQR)
            IF ( LQUERY ) THEN
               CALL SORMQR( 'L','N', M, N, MINMN, F, LDF,  &
                            WORK, Z, LDZ, WORK, -1, INFO1 )
               OLWMQR = INT(WORK(1))
               OLWORK = MAX(OLWORK,MINMN+N-1+OLWMQR)
            END IF
         END IF
         IF ( WANTQ ) THEN
            MLWGQR = N
            MLWORK = MAX(MLWORK,MINMN+N-1+MLWGQR)
            IF ( LQUERY ) THEN
                CALL SORGQR( M, MINMN, MINMN, F, LDF, WORK, &
                             WORK, -1, INFO1 )
                OLWGQR = INT(WORK(1))
                OLWORK = MAX(OLWORK,MINMN+N-1+OLWGQR)
            END IF
         END IF
         IMINWR = MAX( 1, IMINWR )
         MLWORK = MAX( 2, MLWORK )
         IF (  LWORK < MLWORK .AND. (.NOT.LQUERY) ) INFO = -31
         IF ( LIWORK < IMINWR .AND. (.NOT.LQUERY) ) INFO = -33
      END IF
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'SGEDMDQ', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
!     Return minimal and optimal workspace sizes
          IWORK(1) = IMINWR
          WORK(1)  = MLWORK
          WORK(2)  = OLWORK
          RETURN
      END IF
!.....
!     Initial QR factorization that is used to represent the
!     snapshots as elements of lower dimensional subspace.
!     For large scale computation with M >>N , at this place
!     one can use an out of core QRF.
!
      CALL SGEQRF( M, N, F, LDF, WORK,         &
                   WORK(MINMN+1), LWORK-MINMN, INFO1 )
!
!     Define X and Y as the snapshots representations in the
!     orthogonal basis computed in the QR factorization.
!     X corresponds to the leading N-1 and Y to the trailing
!     N-1 snapshots.
      CALL SLASET( 'L', MINMN, N-1, ZERO,  ZERO, X, LDX )
      CALL SLACPY( 'U', MINMN, N-1, F,      LDF, X, LDX )
      CALL SLACPY( 'A', MINMN, N-1, F(1,2), LDF, Y, LDY )
      IF ( M >= 3 ) THEN
          CALL SLASET( 'L', MINMN-2, N-2, ZERO,  ZERO, &
                       Y(3,1), LDY )
      END IF
!
!     Compute the DMD of the projected snapshot pairs (X,Y)
      CALL SGEDMD( JOBS, JOBVL, JOBR, JOBF, WHTSVD, MINMN,  &
                   N-1, X, LDX, Y, LDY, NRNK,   TOL, K,     &
                   REIG, IMEIG, Z, LDZ, RES, B, LDB, V,     &
                   LDV, S, LDS, WORK(MINMN+1), LWORK-MINMN, IWORK,  &
                   LIWORK, INFO1 )
      IF ( INFO1 == 2 .OR. INFO1 == 3 ) THEN
          ! Return with error code.
          INFO = INFO1
          RETURN
      ELSE
          INFO = INFO1
      END IF
!
!     The Ritz vectors (Koopman modes) can be explicitly
!     formed or returned in factored form.
      IF ( WNTVEC ) THEN
        ! Compute the eigenvectors explicitly.
        IF ( M > MINMN ) CALL SLASET( 'A', M-MINMN, K, ZERO, &
                                     ZERO, Z(MINMN+1,1), LDZ )
        CALL SORMQR( 'L','N', M, K, MINMN, F, LDF, WORK, Z,  &
                     LDZ, WORK(MINMN+N), LWORK-(MINMN+N-1), INFO1 )
      ELSE IF ( WNTVCF ) THEN
        !   Return the Ritz vectors (eigenvectors) in factored
        !   form Z*V, where Z contains orthonormal matrix (the
        !   product of Q from the initial QR factorization and
        !   the SVD/POD_basis returned by SGEDMD in X) and the
        !   second factor (the eigenvectors of the Rayleigh
        !   quotient) is in the array V, as returned by SGEDMD.
        CALL SLACPY( 'A', N, K, X, LDX, Z, LDZ )
        IF ( M > N ) CALL SLASET( 'A', M-N, K, ZERO, ZERO,   &
                                  Z(N+1,1), LDZ )
        CALL SORMQR( 'L','N', M, K, MINMN, F, LDF, WORK, Z,  &
             LDZ, WORK(MINMN+N), LWORK-(MINMN+N-1), INFO1 )
      END IF
!
!     Some optional output variables:
!
!     The upper triangular factor in the initial QR
!     factorization is optionally returned in the array Y.
!     This is useful if this call to SGEDMDQ is to be
!     followed by a streaming DMD that is implemented in a
!     QR compressed form.
      IF ( WNTTRF ) THEN ! Return the upper triangular R in Y
         CALL SLASET( 'A', MINMN, N, ZERO,  ZERO, Y, LDY )
         CALL SLACPY( 'U', MINMN, N, F, LDF,      Y, LDY )
      END IF
!
!     The orthonormal/orthogonal factor in the initial QR
!     factorization is optionally returned in the array F.
!     Same as with the triangular factor above, this is
!     useful in a streaming DMD.
      IF ( WANTQ ) THEN  ! Q overwrites F
         CALL SORGQR( M, MINMN, MINMN, F, LDF, WORK, &
              WORK(MINMN+N), LWORK-(MINMN+N-1), INFO1 )
      END IF
!
      RETURN
!
      END SUBROUTINE SGEDMDQ