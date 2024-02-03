      REAL ssxt1(IJOB, D1, N1, D2, N2, ABSTOL, ULP, UNFL ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IJOB, N1, N2;
      REAL               ABSTOL, ULP, UNFL;
      // ..
      // .. Array Arguments ..
      REAL               D1( * ), D2( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               TEMP1, TEMP2;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      TEMP1 = ZERO;

      J = 1;
      for (I = 1; I <= N1; I++) { // 20
         } // 10
         if ( D2( J ) < D1( I ) && J < N2 ) {
            J = J + 1;
            GO TO 10;
         }
         if ( J == 1 ) {
            TEMP2 = ABS( D2( J )-D1( I ) );
            if (IJOB == 2) TEMP2 = TEMP2 / max( UNFL, ABSTOL+ULP*ABS( D1( I ) ) );
         } else {
            TEMP2 = min( ABS( D2( J )-D1( I ) ), ABS( D1( I )-D2( J-1 ) ) )             IF( IJOB == 2 ) TEMP2 = TEMP2 / max( UNFL, ABSTOL+ULP*ABS( D1( I ) ) );
         }
         TEMP1 = max( TEMP1, TEMP2 );
      } // 20

      SSXT1 = TEMP1;
      return;
      }
