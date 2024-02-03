      void dtbt05(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), B( LDB, * ), BERR( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, K, NZ;
      double             AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                IDAMAX;
      //- double             DLAMCH;
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
      NOTRAN = LSAME( TRANS, 'N' );
      UNIT = LSAME( DIAG, 'U' );
      NZ = min( KD, N-1 ) + 1;

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

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      IFU = 0;
      if (UNIT) IFU = 1;
      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = ( B( I, K ) ).abs();
            if ( UPPER ) {
               if ( !NOTRAN ) {
                  DO 40 J = max( I-KD, 1 ), I - IFU;
                     TMP = TMP + ( AB( KD+1-I+J, I ) ).abs()* ( X( J, K ) ).abs();
                  } // 40
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
               } else {
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
                  DO 50 J = I + IFU, min( I+KD, N );
                     TMP = TMP + ( AB( KD+1+I-J, J ) ).abs()* ( X( J, K ) ).abs();
                  } // 50
               }
            } else {
               if ( NOTRAN ) {
                  DO 60 J = max( I-KD, 1 ), I - IFU;
                     TMP = TMP + ( AB( 1+I-J, J ) ).abs()*( X( J, K ) ).abs();
                  } // 60
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
               } else {
                  if (UNIT) TMP = TMP + ( X( I, K ) ).abs();
                  DO 70 J = I + IFU, min( I+KD, N );
                     TMP = TMP + ( AB( 1+J-I, I ) ).abs()*( X( J, K ) ).abs();
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
            RESLTS( 2 ) = TMP;
         } else {
            RESLTS( 2 ) = max( RESLTS( 2 ), TMP );
         }
      } // 90

      return;
      }
