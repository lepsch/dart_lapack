      void cgbt05(TRANS, N, KL, KU, NRHS, AB, LDAB, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                KL, KU, LDAB, LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               BERR( * ), FERR( * ), RESLTS( * );
      COMPLEX            AB( LDAB, * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IMAX, J, K, NZ;
      REAL               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      COMPLEX            ZDUM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ( REAL( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS( 1 ) = ZERO;
         RESLTS( 2 ) = ZERO;
         return;
      }

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      NOTRAN = LSAME( TRANS, 'N' );
      NZ = min( KL+KU+2, N+1 );

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
      RESLTS( 1 ) = ERRBND;

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
            RESLTS( 2 ) = TMP;
         } else {
            RESLTS( 2 ) = max( RESLTS( 2 ), TMP );
         }
      } // 70

      return;
      }
