      void dppt05(UPLO, N, NRHS, AP, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             AP( * ), B( LDB, * ), BERR( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IMAX, J, JC, K;
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

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS( 1 ) = ZERO;
         RESLTS( 2 ) = ZERO;
         return;
      }

      EPS = DLAMCH( 'Epsilon' );
      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      UPPER = LSAME( UPLO, 'U' );

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = IDAMAX( N, X( 1, J ), 1 );
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
      RESLTS( 1 ) = ERRBND;

      // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
      // (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = ( B( I, K ) ).abs();
            if ( UPPER ) {
               JC = ( ( I-1 )*I ) / 2;
               for (J = 1; J <= I; J++) { // 40
                  TMP = TMP + ( AP( JC+J ) ).abs()*( X( J, K ) ).abs();
               } // 40
               JC = JC + I;
               for (J = I + 1; J <= N; J++) { // 50
                  TMP = TMP + ( AP( JC ) ).abs()*( X( J, K ) ).abs();
                  JC = JC + J;
               } // 50
            } else {
               JC = I;
               for (J = 1; J <= I - 1; J++) { // 60
                  TMP = TMP + ( AP( JC ) ).abs()*( X( J, K ) ).abs();
                  JC = JC + N - J;
               } // 60
               for (J = I; J <= N; J++) { // 70
                  TMP = TMP + ( AP( JC+J-I ) ).abs()*( X( J, K ) ).abs();
               } // 70
            }
            if ( I == 1 ) {
               AXBI = TMP;
            } else {
               AXBI = min( AXBI, TMP );
            }
         } // 80
         TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / max( AXBI, ( N+1 )*UNFL ) );
         if ( K == 1 ) {
            RESLTS( 2 ) = TMP;
         } else {
            RESLTS( 2 ) = max( RESLTS( 2 ), TMP );
         }
      } // 90

      return;
      }
