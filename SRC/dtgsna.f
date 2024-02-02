      SUBROUTINE DTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,
     $                   IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), DIF( * ), S( * ),
     $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            DIFDRI
      PARAMETER          ( DIFDRI = 3 )
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                   FOUR = 4.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, PAIR, SOMCON, WANTBH, WANTDF, WANTS
      INTEGER            I, IERR, IFST, ILST, IZ, K, KS, LWMIN, N1, N2
      DOUBLE PRECISION   ALPHAI, ALPHAR, ALPRQT, BETA, C1, C2, COND,
     $                   EPS, LNRM, RNRM, ROOT1, ROOT2, SCALE, SMLNUM,
     $                   TMPII, TMPIR, TMPRI, TMPRR, UHAV, UHAVI, UHBV,
     $                   UHBVI
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUMMY( 1 ), DUMMY1( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLAPY2, DNRM2
      EXTERNAL           LSAME, DDOT, DLAMCH, DLAPY2, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DLACPY, DLAG2, DTGEXC, DTGSYL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTDF = LSAME( JOB, 'V' ) .OR. WANTBH
*
      SOMCON = LSAME( HOWMNY, 'S' )
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.WANTS .AND. .NOT.WANTDF ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( WANTS .AND. LDVL.LT.N ) THEN
         INFO = -10
      ELSE IF( WANTS .AND. LDVR.LT.N ) THEN
         INFO = -12
      ELSE
*
*        Set M to the number of eigenpairs for which condition numbers
*        are required, and test MM.
*
         IF( SOMCON ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 K = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
               ELSE
                  IF( K.LT.N ) THEN
                     IF( A( K+1, K ).EQ.ZERO ) THEN
                        IF( SELECT( K ) )
     $                     M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( K ) .OR. SELECT( K+1 ) )
     $                     M = M + 2
                     END IF
                  ELSE
                     IF( SELECT( N ) )
     $                  M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
*
         IF( N.EQ.0 ) THEN
            LWMIN = 1
         ELSE IF( LSAME( JOB, 'V' ) .OR. LSAME( JOB, 'B' ) ) THEN
            LWMIN = 2*N*( N + 2 ) + 16
         ELSE
            LWMIN = N
         END IF
         WORK( 1 ) = LWMIN
*
         IF( MM.LT.M ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTGSNA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      KS = 0
      PAIR = .FALSE.
*
      DO 20 K = 1, N
*
*        Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block.
*
         IF( PAIR ) THEN
            PAIR = .FALSE.
            GO TO 20
         ELSE
            IF( K.LT.N )
     $         PAIR = A( K+1, K ).NE.ZERO
         END IF
*
*        Determine whether condition numbers are required for the k-th
*        eigenpair.
*
         IF( SOMCON ) THEN
            IF( PAIR ) THEN
               IF( .NOT.SELECT( K ) .AND. .NOT.SELECT( K+1 ) )
     $            GO TO 20
            ELSE
               IF( .NOT.SELECT( K ) )
     $            GO TO 20
            END IF
         END IF
*
         KS = KS + 1
*
         IF( WANTS ) THEN
*
*           Compute the reciprocal condition number of the k-th
*           eigenvalue.
*
            IF( PAIR ) THEN
*
*              Complex eigenvalue pair.
*
               RNRM = DLAPY2( DNRM2( N, VR( 1, KS ), 1 ),
     $                DNRM2( N, VR( 1, KS+1 ), 1 ) )
               LNRM = DLAPY2( DNRM2( N, VL( 1, KS ), 1 ),
     $                DNRM2( N, VL( 1, KS+1 ), 1 ) )
               CALL DGEMV( 'N', N, N, ONE, A, LDA, VR( 1, KS ), 1, ZERO,
     $                     WORK, 1 )
               TMPRR = DDOT( N, WORK, 1, VL( 1, KS ), 1 )
               TMPRI = DDOT( N, WORK, 1, VL( 1, KS+1 ), 1 )
               CALL DGEMV( 'N', N, N, ONE, A, LDA, VR( 1, KS+1 ), 1,
     $                     ZERO, WORK, 1 )
               TMPII = DDOT( N, WORK, 1, VL( 1, KS+1 ), 1 )
               TMPIR = DDOT( N, WORK, 1, VL( 1, KS ), 1 )
               UHAV = TMPRR + TMPII
               UHAVI = TMPIR - TMPRI
               CALL DGEMV( 'N', N, N, ONE, B, LDB, VR( 1, KS ), 1, ZERO,
     $                     WORK, 1 )
               TMPRR = DDOT( N, WORK, 1, VL( 1, KS ), 1 )
               TMPRI = DDOT( N, WORK, 1, VL( 1, KS+1 ), 1 )
               CALL DGEMV( 'N', N, N, ONE, B, LDB, VR( 1, KS+1 ), 1,
     $                     ZERO, WORK, 1 )
               TMPII = DDOT( N, WORK, 1, VL( 1, KS+1 ), 1 )
               TMPIR = DDOT( N, WORK, 1, VL( 1, KS ), 1 )
               UHBV = TMPRR + TMPII
               UHBVI = TMPIR - TMPRI
               UHAV = DLAPY2( UHAV, UHAVI )
               UHBV = DLAPY2( UHBV, UHBVI )
               COND = DLAPY2( UHAV, UHBV )
               S( KS ) = COND / ( RNRM*LNRM )
               S( KS+1 ) = S( KS )
*
            ELSE
*
*              Real eigenvalue.
*
               RNRM = DNRM2( N, VR( 1, KS ), 1 )
               LNRM = DNRM2( N, VL( 1, KS ), 1 )
               CALL DGEMV( 'N', N, N, ONE, A, LDA, VR( 1, KS ), 1, ZERO,
     $                     WORK, 1 )
               UHAV = DDOT( N, WORK, 1, VL( 1, KS ), 1 )
               CALL DGEMV( 'N', N, N, ONE, B, LDB, VR( 1, KS ), 1, ZERO,
     $                     WORK, 1 )
               UHBV = DDOT( N, WORK, 1, VL( 1, KS ), 1 )
               COND = DLAPY2( UHAV, UHBV )
               IF( COND.EQ.ZERO ) THEN
                  S( KS ) = -ONE
               ELSE
                  S( KS ) = COND / ( RNRM*LNRM )
               END IF
            END IF
         END IF
*
         IF( WANTDF ) THEN
            IF( N.EQ.1 ) THEN
               DIF( KS ) = DLAPY2( A( 1, 1 ), B( 1, 1 ) )
               GO TO 20
            END IF
*
*           Estimate the reciprocal condition number of the k-th
*           eigenvectors.
            IF( PAIR ) THEN
*
*              Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)).
*              Compute the eigenvalue(s) at position K.
*
               WORK( 1 ) = A( K, K )
               WORK( 2 ) = A( K+1, K )
               WORK( 3 ) = A( K, K+1 )
               WORK( 4 ) = A( K+1, K+1 )
               WORK( 5 ) = B( K, K )
               WORK( 6 ) = B( K+1, K )
               WORK( 7 ) = B( K, K+1 )
               WORK( 8 ) = B( K+1, K+1 )
               CALL DLAG2( WORK, 2, WORK( 5 ), 2, SMLNUM*EPS, BETA,
     $                     DUMMY1( 1 ), ALPHAR, DUMMY( 1 ), ALPHAI )
               ALPRQT = ONE
               C1 = TWO*( ALPHAR*ALPHAR+ALPHAI*ALPHAI+BETA*BETA )
               C2 = FOUR*BETA*BETA*ALPHAI*ALPHAI
               ROOT1 = C1 + SQRT( C1*C1-4.0D0*C2 )
               ROOT1 = ROOT1 / TWO
               ROOT2 = C2 / ROOT1
               COND = MIN( SQRT( ROOT1 ), SQRT( ROOT2 ) )
            END IF
*
*           Copy the matrix (A, B) to the array WORK and swap the
*           diagonal block beginning at A(k,k) to the (1,1) position.
*
            CALL DLACPY( 'Full', N, N, A, LDA, WORK, N )
            CALL DLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
            IFST = K
            ILST = 1
*
            CALL DTGEXC( .FALSE., .FALSE., N, WORK, N, WORK( N*N+1 ), N,
     $                   DUMMY, 1, DUMMY1, 1, IFST, ILST,
     $                   WORK( N*N*2+1 ), LWORK-2*N*N, IERR )
*
            IF( IERR.GT.0 ) THEN
*
*              Ill-conditioned problem - swap rejected.
*
               DIF( KS ) = ZERO
            ELSE
*
*              Reordering successful, solve generalized Sylvester
*              equation for R and L,
*                         A22 * R - L * A11 = A12
*                         B22 * R - L * B11 = B12,
*              and compute estimate of Difl((A11,B11), (A22, B22)).
*
               N1 = 1
               IF( WORK( 2 ).NE.ZERO )
     $            N1 = 2
               N2 = N - N1
               IF( N2.EQ.0 ) THEN
                  DIF( KS ) = COND
               ELSE
                  I = N*N + 1
                  IZ = 2*N*N + 1
                  CALL DTGSYL( 'N', DIFDRI, N2, N1, WORK( N*N1+N1+1 ),
     $                         N, WORK, N, WORK( N1+1 ), N,
     $                         WORK( N*N1+N1+I ), N, WORK( I ), N,
     $                         WORK( N1+I ), N, SCALE, DIF( KS ),
     $                         WORK( IZ+1 ), LWORK-2*N*N, IWORK, IERR )
*
                  IF( PAIR )
     $               DIF( KS ) = MIN( MAX( ONE, ALPRQT )*DIF( KS ),
     $                           COND )
               END IF
            END IF
            IF( PAIR )
     $         DIF( KS+1 ) = DIF( KS )
         END IF
         IF( PAIR )
     $      KS = KS + 1
*
   20 CONTINUE
      WORK( 1 ) = LWMIN
      RETURN
*
*     End of DTGSNA
*
      END