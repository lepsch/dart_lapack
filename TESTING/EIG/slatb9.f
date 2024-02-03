      SUBROUTINE SLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                IMAT, KLA, KUA, KLB, KUB, M, P, MODEA, MODEB, N;
      REAL               ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               SHRINK, TENTH
      PARAMETER          ( SHRINK = 0.25E0, TENTH = 0.1E+0 )
      REAL               ONE, TEN
      PARAMETER          ( ONE = 1.0E+0, TEN = 1.0E+1 )
*     ..
*     .. Local Scalars ..
      bool               FIRST;
      REAL               BADC1, BADC2, EPS, LARGE, SMALL
*     ..
*     .. External Functions ..
      bool               LSAMEN;
      REAL               SLAMCH
      EXTERNAL           LSAMEN, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Save statement ..
      SAVE               EPS, SMALL, LARGE, BADC1, BADC2, FIRST
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     Set some constants for use in the subroutine.
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         EPS = SLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = SLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      END IF
*
*     Set some parameters we don't plan to change.
*
      TYPE = 'N'
      DISTA = 'S'
      DISTB = 'S'
      MODEA = 3
      MODEB = 4
*
*     Set the lower and upper bandwidths.
*
      IF( LSAMEN( 3, PATH, 'GRQ') .OR. LSAMEN( 3, PATH, 'LSE') .OR. LSAMEN( 3, PATH, 'GSV') )THEN
*
*        A: M by N, B: P by N
*
         IF( IMAT.EQ.1 ) THEN
*
*           A: diagonal, B: upper triangular
*
            KLA = 0
            KUA = 0
            KLB = 0
            KUB = MAX( N-1,0 )
*
         ELSE IF( IMAT.EQ.2 ) THEN
*
*           A: upper triangular, B: upper triangular
*
            KLA = 0
            KUA = MAX( N-1, 0 )
            KLB = 0
            KUB = MAX( N-1, 0 )
*
         ELSE IF( IMAT.EQ.3 ) THEN
*
*           A: lower triangular, B: upper triangular
*
            KLA = MAX( M-1, 0 )
            KUA = 0
            KLB = 0
            KUB = MAX( N-1, 0 )
*
         ELSE
*
*           A: general dense, B: general dense
*
            KLA = MAX( M-1, 0 )
            KUA = MAX( N-1, 0 )
            KLB = MAX( P-1, 0 )
            KUB = MAX( N-1, 0 )
*
         END IF
*
      ELSE IF( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GLM') )THEN
*
*        A: N by M, B: N by P
*
         IF( IMAT.EQ.1 ) THEN
*
*           A: diagonal, B: lower triangular
*
            KLA = 0
            KUA = 0
            KLB = MAX( N-1,0 )
            KUB = 0
         ELSE IF( IMAT.EQ.2 ) THEN
*
*           A: lower triangular, B: diagonal
*
            KLA = MAX( N-1, 0 )
            KUA = 0
            KLB = 0
            KUB = 0
*
         ELSE IF( IMAT.EQ.3 ) THEN
*
*           A: lower triangular, B: upper triangular
*
            KLA = MAX( N-1, 0 )
            KUA = 0
            KLB = 0
            KUB = MAX( P-1, 0 )
*
         ELSE
*
*           A: general dense, B: general dense
*
            KLA = MAX( N-1, 0 )
            KUA = MAX( M-1, 0 )
            KLB = MAX( N-1, 0 )
            KUB = MAX( P-1, 0 )
         END IF
*
      END IF
*
*     Set the condition number and norm.
*
      CNDNMA = TEN*TEN
      CNDNMB = TEN
      IF( LSAMEN( 3, PATH, 'GQR') .OR. LSAMEN( 3, PATH, 'GRQ') .OR. LSAMEN( 3, PATH, 'GSV') )THEN
         IF( IMAT.EQ.5 ) THEN
            CNDNMA = BADC1
            CNDNMB = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNMA = BADC2
            CNDNMB = BADC2
         ELSE IF( IMAT.EQ.7 ) THEN
            CNDNMA = BADC1
            CNDNMB = BADC2
         ELSE IF( IMAT.EQ.8 ) THEN
            CNDNMA = BADC2
            CNDNMB = BADC1
         END IF
      END IF
*
      ANORM = TEN
      BNORM = TEN*TEN*TEN
      IF( LSAMEN( 3, PATH, 'GQR') .OR. LSAMEN( 3, PATH, 'GRQ') )THEN
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
            BNORM = LARGE
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
            BNORM = SMALL
         END IF
      END IF
*
      IF( N.LE.1 )THEN
         CNDNMA = ONE
         CNDNMB = ONE
      END IF
*
      RETURN
*
*     End of SLATB9
*
      END
