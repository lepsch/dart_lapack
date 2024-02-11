      void zget07(final int TRANS, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, final Matrix<double> XACT, final int LDXACT, final int FERR, final int CHKFERR, final int BERR, final int RESLTS,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      bool               CHKFERR;
      int                LDA, LDB, LDX, LDXACT, N, NRHS;
      double             BERR( * ), FERR( * ), RESLTS( * );
      Complex         A( LDA, * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               NOTRAN;
      int                I, IMAX, J, K;
      double             AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      Complex         ZDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                IZAMAX;
      //- double             DLAMCH;
      // EXTERNAL lsame, IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS[1] = ZERO;
         RESLTS[2] = ZERO;
         return;
      }

      EPS = dlamch( 'Epsilon' );
      UNFL = dlamch( 'Safe minimum' );
      OVFL = ONE / UNFL;
      NOTRAN = lsame( TRANS, 'N' );

      // Test 1:  Compute the maximum of
      //    norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO;
      if ( CHKFERR ) {
         for (J = 1; J <= NRHS; J++) { // 30
            IMAX = IZAMAX( N, X( 1, J ), 1 );
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
      }
      RESLTS[1] = ERRBND;

      // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
      // (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

      for (K = 1; K <= NRHS; K++) { // 70
         for (I = 1; I <= N; I++) { // 60
            TMP = CABS1( B( I, K ) );
            if ( NOTRAN ) {
               for (J = 1; J <= N; J++) { // 40
                  TMP = TMP + CABS1( A( I, J ) )*CABS1( X( J, K ) );
               } // 40
            } else {
               for (J = 1; J <= N; J++) { // 50
                  TMP = TMP + CABS1( A( J, I ) )*CABS1( X( J, K ) );
               } // 50
            }
            if ( I == 1 ) {
               AXBI = TMP;
            } else {
               AXBI = min( AXBI, TMP );
            }
         } // 60
         TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / max( AXBI, ( N+1 )*UNFL ) );
         if ( K == 1 ) {
            RESLTS[2] = TMP;
         } else {
            RESLTS[2] = max( RESLTS( 2 ), TMP );
         }
      } // 70

      }
