      SUBROUTINE DGTT05( TRANS, N, NRHS, DL, D, DU, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, LDXACT, N, NRHS;
*     ..
*     .. Array Arguments ..
      double             B( LDB, * ), BERR( * ), D( * ), DL( * ), DU( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOTRAN;
      int                I, IMAX, J, K, NZ;
      double             AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
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
      NOTRAN = LSAME( TRANS, 'N' )
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
*     (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
*
      DO 60 K = 1, NRHS
         IF( NOTRAN ) THEN
            IF( N.EQ.1 ) THEN
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) )
            ELSE
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) ) + ABS( DU( 1 )*X( 2, K ) )
               DO 40 I = 2, N - 1
                  TMP = ABS( B( I, K ) ) + ABS( DL( I-1 )*X( I-1, K ) ) + ABS( D( I )*X( I, K ) ) + ABS( DU( I )*X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
   40          CONTINUE
               TMP = ABS( B( N, K ) ) + ABS( DL( N-1 )*X( N-1, K ) ) + ABS( D( N )*X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            END IF
         ELSE
            IF( N.EQ.1 ) THEN
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) )
            ELSE
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) ) + ABS( DL( 1 )*X( 2, K ) )
               DO 50 I = 2, N - 1
                  TMP = ABS( B( I, K ) ) + ABS( DU( I-1 )*X( I-1, K ) ) + ABS( D( I )*X( I, K ) ) + ABS( DL( I )*X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
   50          CONTINUE
               TMP = ABS( B( N, K ) ) + ABS( DU( N-1 )*X( N-1, K ) ) + ABS( D( N )*X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            END IF
         END IF
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / MAX( AXBI, NZ*UNFL ) )
         IF( K.EQ.1 ) THEN
            RESLTS( 2 ) = TMP
         ELSE
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         END IF
   60 CONTINUE
*
      RETURN
*
*     End of DGTT05
*
      END
