      void ctrt05(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Matrix<double> X_, final int LDX, final Matrix<double> XACT_, final int LDXACT, final int FERR, final int BERR, final int RESLTS,) {
  final A = A_.dim();
  final B = B_.dim();
  final X = X_.dim();
  final XACT = XACT_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, LDXACT, N, NRHS;
      double               BERR( * ), FERR( * ), RESLTS( * );
      Complex            A( LDA, * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, K;
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

      // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
      // (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      IFU = 0;
      if (UNIT) IFU = 1;
      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = CABS1( B( I, K ) );
            if ( UPPER ) {
               if ( !NOTRAN ) {
                  for (J = 1; J <= I - IFU; J++) { // 40
                     TMP = TMP + CABS1( A( J, I ) )*CABS1( X( J, K ) );
                  } // 40
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
               } else {
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
                  for (J = I + IFU; J <= N; J++) { // 50
                     TMP = TMP + CABS1( A( I, J ) )*CABS1( X( J, K ) );
                  } // 50
               }
            } else {
               if ( NOTRAN ) {
                  for (J = 1; J <= I - IFU; J++) { // 60
                     TMP = TMP + CABS1( A( I, J ) )*CABS1( X( J, K ) );
                  } // 60
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
               } else {
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
                  for (J = I + IFU; J <= N; J++) { // 70
                     TMP = TMP + CABS1( A( J, I ) )*CABS1( X( J, K ) );
                  } // 70
               }
            }
            if ( I == 1 ) {
               AXBI = TMP;
            } else {
               AXBI = min( AXBI, TMP );
            }
         } // 80
         TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / max( AXBI, ( N+1 )*UNFL ) );
         if ( K == 1 ) {
            RESLTS[2] = TMP;
         } else {
            RESLTS[2] = max( RESLTS( 2 ), TMP );
         }
      } // 90

      }
