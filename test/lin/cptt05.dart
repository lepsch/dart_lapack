      void cptt05(N, NRHS, D, E, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      double               BERR( * ), D( * ), FERR( * ), RESLTS( * );
      Complex            B( LDB, * ), E( * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IMAX, J, K, NZ;
      double               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      Complex            ZDUM;
      // ..
      // .. External Functions ..
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS[1] = ZERO;
         RESLTS[2] = ZERO;
         return;
      }

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      NZ = 4;

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
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

      for (K = 1; K <= NRHS; K++) { // 50
         if ( N == 1 ) {
            AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 )*X( 1, K ) );
         } else {
            AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 )*X( 1, K ) ) + CABS1( E( 1 ) )*CABS1( X( 2, K ) );
            for (I = 2; I <= N - 1; I++) { // 40
               TMP = CABS1( B( I, K ) ) + CABS1( E( I-1 ) )* CABS1( X( I-1, K ) ) + CABS1( D( I )*X( I, K ) ) + CABS1( E( I ) )*CABS1( X( I+1, K ) );
               AXBI = min( AXBI, TMP );
            } // 40
            TMP = CABS1( B( N, K ) ) + CABS1( E( N-1 ) )* CABS1( X( N-1, K ) ) + CABS1( D( N )*X( N, K ) );
            AXBI = min( AXBI, TMP );
         }
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / max( AXBI, NZ*UNFL ) );
         if ( K == 1 ) {
            RESLTS[2] = TMP;
         } else {
            RESLTS[2] = max( RESLTS( 2 ), TMP );
         }
      } // 50

      return;
      }