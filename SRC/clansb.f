      REAL clansb(NORM, UPLO, N, K, AB, LDAB, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                K, LDAB, N;
      // ..
      // .. Array Arguments ..
      REAL               WORK( * );
      COMPLEX            AB( LDAB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, L;
      REAL               ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      //- bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
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
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               DO 10 I = max( K+2-J, 1 ), K + 1;
                  SUM = ( AB( I, J ) ).abs();
                  if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               } // 10
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               DO 30 I = 1, min( N+1-J, K+1 );
                  SUM = ( AB( I, J ) ).abs();
                  if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               } // 30
            } // 40
         }
      } else if ( ( LSAME( NORM, 'I' ) ) || ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is symmetric).

         VALUE = ZERO;
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               L = K + 1 - J;
               DO 50 I = max( 1, J-K ), J - 1;
                  ABSA = ( AB( L+I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK( I ) = WORK( I ) + ABSA;
               } // 50
               WORK( J ) = SUM + ( AB( K+1, J ) ).abs();
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK( I ) = ZERO;
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ( AB( 1, J ) ).abs();
               L = 1 - J;
               DO 90 I = J + 1, min( N, J+K );
                  ABSA = ( AB( L+I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK( I ) = WORK( I ) + ABSA;
               } // 90
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 100
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( K > 0 ) {
            if ( LSAME( UPLO, 'U' ) ) {
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
         classq(N, AB( L, 1 ), LDAB, SCALE, SUM );
         VALUE = SCALE*sqrt( SUM );
      }

      CLANSB = VALUE;
      return;
      }
