      void zgbt05(final int TRANS, final int N, final int KL, final int KU, final int NRHS, final Matrix<double> AB_, final int LDAB, final Matrix<double> B_, final int LDB, final Matrix<double> X_, final int LDX, final Matrix<double> XACT_, final int LDXACT, final int FERR, final int BERR, final int RESLTS,) {
  final AB = AB_.dim();
  final B = B_.dim();
  final X = X_.dim();
  final XACT = XACT_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                KL, KU, LDAB, LDB, LDX, LDXACT, N, NRHS;
      double             BERR( * ), FERR( * ), RESLTS( * );
      Complex         AB( LDAB, * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               NOTRAN;
      int                I, IMAX, J, K, NZ;
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
      NZ = min( KL+KU+2, N+1 );

      // Test 1:  Compute the maximum of
      //    norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO;
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
      RESLTS[1] = ERRBND;

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

      for (K = 1; K <= NRHS; K++) { // 70
         for (I = 1; I <= N; I++) { // 60
            TMP = CABS1( B( I, K ) );
            if ( NOTRAN ) {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) { // 40
                  TMP = TMP + CABS1( AB( KU+1+I-J, J ) )* CABS1( X( J, K ) );
               } // 40
            } else {
               for (J = max( I-KU, 1 ); J <= min( I+KL, N ); J++) { // 50
                  TMP = TMP + CABS1( AB( KU+1+J-I, I ) )* CABS1( X( J, K ) );
               } // 50
            }
            if ( I == 1 ) {
               AXBI = TMP;
            } else {
               AXBI = min( AXBI, TMP );
            }
         } // 60
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / max( AXBI, NZ*UNFL ) );
         if ( K == 1 ) {
            RESLTS[2] = TMP;
         } else {
            RESLTS[2] = max( RESLTS( 2 ), TMP );
         }
      } // 70

      }
