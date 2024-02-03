      double zlangb(NORM, N, KL, KU, AB, LDAB, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                KL, KU, LDAB, N;
      // ..
      // .. Array Arguments ..
      double             WORK( * );
      Complex         AB( LDAB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K, L;
      double             SCALE, SUM, VALUE, TEMP;
      // ..
      // .. External Functions ..
      //- bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 20
            DO 10 I = max( KU+2-J, 1 ), min( N+KU+1-J, KL+KU+1 );
               TEMP = ( AB( I, J ) ).abs();
               if( VALUE < TEMP || DISNAN( TEMP ) ) VALUE = TEMP;
            } // 10
         } // 20
      } else if ( ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO;
            DO 30 I = max( KU+2-J, 1 ), min( N+KU+1-J, KL+KU+1 );
               SUM = SUM + ( AB( I, J ) ).abs();
            } // 30
            if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
         } // 40
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= N; I++) { // 50
            WORK( I ) = ZERO;
         } // 50
         for (J = 1; J <= N; J++) { // 70
            K = KU + 1 - J;
            DO 60 I = max( 1, J-KU ), min( N, J+KL );
               WORK( I ) = WORK( I ) + ( AB( K+I, J ) ).abs();
            } // 60
         } // 70
         VALUE = ZERO;
         for (I = 1; I <= N; I++) { // 80
            TEMP = WORK( I );
            if( VALUE < TEMP || DISNAN( TEMP ) ) VALUE = TEMP;
         } // 80
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         for (J = 1; J <= N; J++) { // 90
            L = max( 1, J-KU );
            K = KU + 1 - J + L;
            zlassq(min( N, J+KL )-L+1, AB( K, J ), 1, SCALE, SUM );
         } // 90
         VALUE = SCALE*sqrt( SUM );
      }

      ZLANGB = VALUE;
      return;
      }
