      SUBROUTINE DPTT05( N, NRHS, D, E, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDB, LDX, LDXACT, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), E( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IMAX, J, K, NZ
      DOUBLE PRECISION   AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IDAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      NZ = 4
*
*     Test 1:  Compute the maximum of
*        norm(X - XACT) / ( norm(X) * FERR )
*     over all the vectors X and XACT using the infinity-norm.
*
      ERRBND = ZERO
      DO 30 J = 1, NRHS
         IMAX = IDAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( ABS( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         DO 10 I = 1, N
            DIFF = MAX( DIFF, ABS( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE
*
         IF( XNORM.GT.ONE ) THEN
            GO TO 20
         ELSE IF( DIFF.LE.OVFL*XNORM ) THEN
            GO TO 20
         ELSE
            ERRBND = ONE / EPS
            GO TO 30
         END IF
*
   20    CONTINUE
         IF( DIFF / XNORM.LE.FERR( J ) ) THEN
            ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
         ELSE
            ERRBND = ONE / EPS
         END IF
   30 CONTINUE
      RESLTS( 1 ) = ERRBND
*
*     Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
*     (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
*
      DO 50 K = 1, NRHS
         IF( N.EQ.1 ) THEN
            AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) )
         ELSE
            AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) ) + ABS( E( 1 )*X( 2, K ) )
            DO 40 I = 2, N - 1
               TMP = ABS( B( I, K ) ) + ABS( E( I-1 )*X( I-1, K ) ) + ABS( D( I )*X( I, K ) ) + ABS( E( I )*X( I+1, K ) )
               AXBI = MIN( AXBI, TMP )
   40       CONTINUE
            TMP = ABS( B( N, K ) ) + ABS( E( N-1 )*X( N-1, K ) ) + ABS( D( N )*X( N, K ) )
            AXBI = MIN( AXBI, TMP )
         END IF
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / MAX( AXBI, NZ*UNFL ) )
         IF( K.EQ.1 ) THEN
            RESLTS( 2 ) = TMP
         ELSE
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         END IF
   50 CONTINUE
*
      RETURN
*
*     End of DPTT05
*
      END
