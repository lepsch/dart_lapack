      double clanhb(NORM, UPLO, N, K, AB, LDAB, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                K, LDAB, N;
      // ..
      // .. Array Arguments ..
      double               WORK( * );
      Complex            AB( LDAB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, L;
      double               ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      // EXTERNAL lsame, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = max( K+2-J, 1 ); I <= K; I++) { // 10
                  SUM = ( AB( I, J ) ).abs();
                  if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               } // 10
               SUM = ABS( double( AB( K+1, J ) ) );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( double( AB( 1, J ) ) );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               for (I = 2; I <= min( N+1-J, K+1 ); I++) { // 30
                  SUM = ( AB( I, J ) ).abs();
                  if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               } // 30
            } // 40
         }
      } else if ( ( lsame( NORM, 'I' ) ) || ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               L = K + 1 - J;
               for (I = max( 1, J-K ); I <= J - 1; I++) { // 50
                  ABSA = ( AB( L+I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 50
               WORK[J] = SUM + ABS( double( AB( K+1, J ) ) );
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK[I] = ZERO;
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( double( AB( 1, J ) ) );
               L = 1 - J;
               for (I = J + 1; I <= min( N, J+K ); I++) { // 90
                  ABSA = ( AB( L+I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 90
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 100
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( K > 0 ) {
            if ( lsame( UPLO, 'U' ) ) {
               for (J = 2; J <= N; J++) { // 110
                  classq(min( J-1, K ), AB( max( K+2-J, 1 ), J ), 1, SCALE, SUM );
               } // 110
               L = K + 1;
            } else {
               for (J = 1; J <= N - 1; J++) { // 120
                  classq(min( N-J, K ), AB( 2, J ), 1, SCALE, SUM );
               } // 120
               L = 1;
            }
            SUM = 2*SUM;
         } else {
            L = 1;
         }
         for (J = 1; J <= N; J++) { // 130
            if ( double( AB( L, J ) ) != ZERO ) {
               ABSA = ABS( double( AB( L, J ) ) );
               if ( SCALE < ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2;
                  SCALE = ABSA;
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2;
               }
            }
         } // 130
         VALUE = SCALE*sqrt( SUM );
      }

      CLANHB = VALUE;
      return;
      }