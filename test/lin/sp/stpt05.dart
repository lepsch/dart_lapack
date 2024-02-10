      void stpt05(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final int AP, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, final Matrix<double> XACT, final int LDXACT, final int FERR, final int BERR, final int RESLTS) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, LDXACT, N, NRHS;
      double               AP( * ), B( LDB, * ), BERR( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, JC, K;
      double               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ISAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

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
         IMAX = ISAMAX( N, X( 1, J ), 1 );
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

      // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
      // (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      IFU = 0;
      if (UNIT) IFU = 1;
      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = ( B( I, K ) ).abs();
            if ( UPPER ) {
               JC = ( ( I-1 )*I ) / 2;
               if ( !NOTRAN ) {
                  for (J = 1; J <= I - IFU; J++) { // 40
                     TMP = TMP + ( AP( JC+J ) ).abs()*( X( J, K ) ).abs();
                  } // 40
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
               } else {
                  JC = JC + I;
                  if ( UNIT ) {
                     TMP = TMP + ( X( I, K ) ).abs();
                     JC = JC + I;
                  }
                  for (J = I + IFU; J <= N; J++) { // 50
                     TMP = TMP + ( AP( JC ) ).abs()*( X( J, K ) ).abs();
                     JC = JC + J;
                  } // 50
               }
            } else {
               if ( NOTRAN ) {
                  JC = I;
                  for (J = 1; J <= I - IFU; J++) { // 60
                     TMP = TMP + ( AP( JC ) ).abs()*( X( J, K ) ).abs();
                     JC = JC + N - J;
                  } // 60
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
               } else {
                  JC = ( I-1 )*( N-I ) + ( I*( I+1 ) ) / 2;
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
                  for (J = I + IFU; J <= N; J++) { // 70
                     TMP = TMP + ( AP( JC+J-I ) ).abs()*( X( J, K ) ).abs();
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
