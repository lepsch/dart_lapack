      SUBROUTINE CTPT05( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, LDXACT, N, NRHS;
*     ..
*     .. Array Arguments ..
      REAL               BERR( * ), FERR( * ), RESLTS( * )
      COMPLEX            AP( * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, JC, K;
      REAL               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
      COMPLEX            ZDUM
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      EXTERNAL           LSAME, ICAMAX, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
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
      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      UNIT = LSAME( DIAG, 'U' )
*
*     Test 1:  Compute the maximum of
*        norm(X - XACT) / ( norm(X) * FERR )
*     over all the vectors X and XACT using the infinity-norm.
*
      ERRBND = ZERO
      DO 30 J = 1, NRHS
         IMAX = ICAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( CABS1( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         DO 10 I = 1, N
            DIFF = MAX( DIFF, CABS1( X( I, J )-XACT( I, J ) ) )
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
*     Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
*     (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )
*
      IFU = 0
      IF( UNIT ) IFU = 1
      DO 90 K = 1, NRHS
         DO 80 I = 1, N
            TMP = CABS1( B( I, K ) )
            IF( UPPER ) THEN
               JC = ( ( I-1 )*I ) / 2
               IF( .NOT.NOTRAN ) THEN
                  DO 40 J = 1, I - IFU
                     TMP = TMP + CABS1( AP( JC+J ) )*CABS1( X( J, K ) )
   40             CONTINUE
                  IF( UNIT ) TMP = TMP + CABS1( X( I, K ) )
               ELSE
                  JC = JC + I
                  IF( UNIT ) THEN
                     TMP = TMP + CABS1( X( I, K ) )
                     JC = JC + I
                  END IF
                  DO 50 J = I + IFU, N
                     TMP = TMP + CABS1( AP( JC ) )*CABS1( X( J, K ) )
                     JC = JC + J
   50             CONTINUE
               END IF
            ELSE
               IF( NOTRAN ) THEN
                  JC = I
                  DO 60 J = 1, I - IFU
                     TMP = TMP + CABS1( AP( JC ) )*CABS1( X( J, K ) )
                     JC = JC + N - J
   60             CONTINUE
                  IF( UNIT ) TMP = TMP + CABS1( X( I, K ) )
               ELSE
                  JC = ( I-1 )*( N-I ) + ( I*( I+1 ) ) / 2
                  IF( UNIT ) TMP = TMP + CABS1( X( I, K ) )
                  DO 70 J = I + IFU, N
                     TMP = TMP + CABS1( AP( JC+J-I ) )* CABS1( X( J, K ) )
   70             CONTINUE
               END IF
            END IF
            IF( I.EQ.1 ) THEN
               AXBI = TMP
            ELSE
               AXBI = MIN( AXBI, TMP )
            END IF
   80    CONTINUE
         TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / MAX( AXBI, ( N+1 )*UNFL ) )
         IF( K.EQ.1 ) THEN
            RESLTS( 2 ) = TMP
         ELSE
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         END IF
   90 CONTINUE
*
      RETURN
*
*     End of CTPT05
*
      END
