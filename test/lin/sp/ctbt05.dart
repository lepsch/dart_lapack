      void ctbt05(final int UPLO, final int TRANS, final int DIAG, final int N, final int KD, final int NRHS, final Matrix<double> AB, final int LDAB, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, final Matrix<double> XACT, final int LDXACT, final int FERR, final int BERR, final int RESLTS) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, LDXACT, N, NRHS;
      double               BERR( * ), FERR( * ), RESLTS( * );
      Complex            AB( LDAB, * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, K, NZ;
      double               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      Complex            ZDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS[1] = ZERO;
         RESLTS[2] = ZERO;
         return;
      }

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      UPPER = lsame( UPLO, 'U' );
      NOTRAN = lsame( TRANS, 'N' );
      UNIT = lsame( DIAG, 'U' );
      NZ = min( KD, N-1 ) + 1;

      // Test 1:  Compute the maximum of
      //    norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = ICAMAX( N, X( 1, J ), 1 );
         XNORM = max( CABS1( X( IMAX, J ) ), UNFL );
         DIFF = ZERO;
         for (I = 1; I <= N; I++) { // 10
            DIFF = max( DIFF, CABS1( X( I, J )-XACT( I, J ) ) );
         } // 10

         if ( XNORM > ONE ) {
            GO TO 20;
         } else if ( DIFF <= OVFL*XNORM ) {
            GO TO 20;
         } else {
            ERRBND = ONE / EPS;
            GO TO 30;
         }

         } // 20
         if ( DIFF / XNORM <= FERR( J ) ) {
            ERRBND = max( ERRBND, ( DIFF / XNORM ) / FERR( J ) );
         } else {
            ERRBND = ONE / EPS;
         }
      } // 30
      RESLTS[1] = ERRBND;

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      IFU = 0;
      if (UNIT) IFU = 1;
      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = CABS1( B( I, K ) );
            if ( UPPER ) {
               if ( !NOTRAN ) {
                  for (J = max( I-KD, 1 ); J <= I - IFU; J++) { // 40
                     TMP = TMP + CABS1( AB( KD+1-I+J, I ) )* CABS1( X( J, K ) );
                  } // 40
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
               } else {
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
                  for (J = I + IFU; J <= min( I+KD, N ); J++) { // 50
                     TMP = TMP + CABS1( AB( KD+1+I-J, J ) )* CABS1( X( J, K ) );
                  } // 50
               }
            } else {
               if ( NOTRAN ) {
                  for (J = max( I-KD, 1 ); J <= I - IFU; J++) { // 60
                     TMP = TMP + CABS1( AB( 1+I-J, J ) )* CABS1( X( J, K ) );
                  } // 60
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
               } else {
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
                  for (J = I + IFU; J <= min( I+KD, N ); J++) { // 70
                     TMP = TMP + CABS1( AB( 1+J-I, I ) )* CABS1( X( J, K ) );
                  } // 70
               }
            }
            if ( I == 1 ) {
               AXBI = TMP;
            } else {
               AXBI = min( AXBI, TMP );
            }
         } // 80
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / max( AXBI, NZ*UNFL ) );
         if ( K == 1 ) {
            RESLTS[2] = TMP;
         } else {
            RESLTS[2] = max( RESLTS( 2 ), TMP );
         }
      } // 90

      }
