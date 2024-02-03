      SUBROUTINE DTBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), B( LDB, * ), BERR( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, K, NZ;
      double             AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      }

      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      UNIT = LSAME( DIAG, 'U' )
      NZ = MIN( KD, N-1 ) + 1

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = IDAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( ABS( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         for (I = 1; I <= N; I++) { // 10
            DIFF = MAX( DIFF, ABS( X( I, J )-XACT( I, J ) ) )
         } // 10

         if ( XNORM.GT.ONE ) {
            GO TO 20
         } else if ( DIFF.LE.OVFL*XNORM ) {
            GO TO 20
         } else {
            ERRBND = ONE / EPS
            GO TO 30
         }

         } // 20
         if ( DIFF / XNORM.LE.FERR( J ) ) {
            ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
         } else {
            ERRBND = ONE / EPS
         }
      } // 30
      RESLTS( 1 ) = ERRBND

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      IFU = 0
      IF( UNIT ) IFU = 1
      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = ABS( B( I, K ) )
            if ( UPPER ) {
               if ( .NOT.NOTRAN ) {
                  DO 40 J = MAX( I-KD, 1 ), I - IFU
                     TMP = TMP + ABS( AB( KD+1-I+J, I ) )* ABS( X( J, K ) )
                  } // 40
                  IF( UNIT ) TMP = TMP + ABS( X( I, K ) )
               } else {
                  IF( UNIT ) TMP = TMP + ABS( X( I, K ) )
                  DO 50 J = I + IFU, MIN( I+KD, N )
                     TMP = TMP + ABS( AB( KD+1+I-J, J ) )* ABS( X( J, K ) )
                  } // 50
               }
            } else {
               if ( NOTRAN ) {
                  DO 60 J = MAX( I-KD, 1 ), I - IFU
                     TMP = TMP + ABS( AB( 1+I-J, J ) )*ABS( X( J, K ) )
                  } // 60
                  IF( UNIT ) TMP = TMP + ABS( X( I, K ) )
               } else {
                  IF( UNIT ) TMP = TMP + ABS( X( I, K ) )
                  DO 70 J = I + IFU, MIN( I+KD, N )
                     TMP = TMP + ABS( AB( 1+J-I, I ) )*ABS( X( J, K ) )
                  } // 70
               }
            }
            if ( I.EQ.1 ) {
               AXBI = TMP
            } else {
               AXBI = MIN( AXBI, TMP )
            }
         } // 80
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / MAX( AXBI, NZ*UNFL ) )
         if ( K.EQ.1 ) {
            RESLTS( 2 ) = TMP
         } else {
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         }
      } // 90

      RETURN

      // End of DTBT05

      }
