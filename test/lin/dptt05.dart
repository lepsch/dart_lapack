      void dptt05(final int N, final int NRHS, final int D, final int E, final Matrix<double> B_, final int LDB, final Matrix<double> X_, final int LDX, final Matrix<double> XACT_, final int LDXACT, final int FERR, final int BERR, final int RESLTS,) {
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDB, LDX, LDXACT, N, NRHS;
      double             B( LDB, * ), BERR( * ), D( * ), E( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, IMAX, J, K, NZ;
      double             AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      // ..
      // .. External Functions ..
      //- int                idamax;
      //- double             DLAMCH;
      // EXTERNAL idamax, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS[1] = ZERO;
         RESLTS[2] = ZERO;
         return;
      }

      EPS = dlamch( 'Epsilon' );
      UNFL = dlamch( 'Safe minimum' );
      OVFL = ONE / UNFL;
      NZ = 4;

      // Test 1:  Compute the maximum of
      //    norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = idamax( N, X( 1, J ), 1 );
         XNORM = max( ( X( IMAX, J ) ).abs(), UNFL );
         DIFF = ZERO;
         for (I = 1; I <= N; I++) { // 10
            DIFF = max( DIFF, ABS( X( I, J )-XACT( I, J ) ) );
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

      for (K = 1; K <= NRHS; K++) { // 50
         if ( N == 1 ) {
            AXBI = ( B( 1, K ) ).abs() + ABS( D( 1 )*X( 1, K ) );
         } else {
            AXBI = ( B( 1, K ) ).abs() + ABS( D( 1 )*X( 1, K ) ) + ABS( E( 1 )*X( 2, K ) );
            for (I = 2; I <= N - 1; I++) { // 40
               TMP = ( B( I, K ) ).abs() + ABS( E( I-1 )*X( I-1, K ) ) + ABS( D( I )*X( I, K ) ) + ABS( E( I )*X( I+1, K ) );
               AXBI = min( AXBI, TMP );
            } // 40
            TMP = ( B( N, K ) ).abs() + ABS( E( N-1 )*X( N-1, K ) ) + ABS( D( N )*X( N, K ) );
            AXBI = min( AXBI, TMP );
         }
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / max( AXBI, NZ*UNFL ) );
         if ( K == 1 ) {
            RESLTS[2] = TMP;
         } else {
            RESLTS[2] = max( RESLTS( 2 ), TMP );
         }
      } // 50

      }
