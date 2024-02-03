      SUBROUTINE CGTT05( TRANS, N, NRHS, DL, D, DU, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               BERR( * ), FERR( * ), RESLTS( * )
      COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * ), XACT( LDXACT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IMAX, J, K, NZ;
      REAL               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
      COMPLEX            ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      END IF

      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      NOTRAN = LSAME( TRANS, 'N' )
      NZ = 4

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO
      DO 30 J = 1, NRHS
         IMAX = ICAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( CABS1( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         DO 10 I = 1, N
            DIFF = MAX( DIFF, CABS1( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE

         IF( XNORM.GT.ONE ) THEN
            GO TO 20
         ELSE IF( DIFF.LE.OVFL*XNORM ) THEN
            GO TO 20
         ELSE
            ERRBND = ONE / EPS
            GO TO 30
         END IF

   20    CONTINUE
         IF( DIFF / XNORM.LE.FERR( J ) ) THEN
            ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
         ELSE
            ERRBND = ONE / EPS
         END IF
   30 CONTINUE
      RESLTS( 1 ) = ERRBND

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

      DO 60 K = 1, NRHS
         IF( NOTRAN ) THEN
            IF( N.EQ.1 ) THEN
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) )
            ELSE
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) ) + CABS1( DU( 1 ) )*CABS1( X( 2, K ) )
               DO 40 I = 2, N - 1
                  TMP = CABS1( B( I, K ) ) + CABS1( DL( I-1 ) )*CABS1( X( I-1, K ) ) + CABS1( D( I ) )*CABS1( X( I, K ) ) + CABS1( DU( I ) )*CABS1( X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
   40          CONTINUE
               TMP = CABS1( B( N, K ) ) + CABS1( DL( N-1 ) )* CABS1( X( N-1, K ) ) + CABS1( D( N ) )* CABS1( X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            END IF
         ELSE
            IF( N.EQ.1 ) THEN
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) )
            ELSE
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) ) + CABS1( DL( 1 ) )*CABS1( X( 2, K ) )
               DO 50 I = 2, N - 1
                  TMP = CABS1( B( I, K ) ) + CABS1( DU( I-1 ) )*CABS1( X( I-1, K ) ) + CABS1( D( I ) )*CABS1( X( I, K ) ) + CABS1( DL( I ) )*CABS1( X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
   50          CONTINUE
               TMP = CABS1( B( N, K ) ) + CABS1( DU( N-1 ) )* CABS1( X( N-1, K ) ) + CABS1( D( N ) )* CABS1( X( N, K ) )
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

      RETURN

      // End of CGTT05

      }
