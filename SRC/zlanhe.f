      double zlanhe(NORM, UPLO, N, A, LDA, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                LDA, N;
      // ..
      // .. Array Arguments ..
      double             WORK( * );
      COMPLEX*16         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  SUM = ABS( A( I, J ) );
                  if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
               } // 10
               SUM = ABS( DBLE( A( J, J ) ) );
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( DBLE( A( J, J ) ) );
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
               for (I = J + 1; I <= N; I++) { // 30
                  SUM = ABS( A( I, J ) );
                  if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
               } // 30
            } // 40
         }
      } else if ( ( LSAME( NORM, 'I' ) ) || ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO;
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ABS( A( I, J ) );
                  SUM = SUM + ABSA;
                  WORK( I ) = WORK( I ) + ABSA;
               } // 50
               WORK( J ) = SUM + ABS( DBLE( A( J, J ) ) );
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I );
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK( I ) = ZERO;
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( DBLE( A( J, J ) ) );
               for (I = J + 1; I <= N; I++) { // 90
                  ABSA = ABS( A( I, J ) );
                  SUM = SUM + ABSA;
                  WORK( I ) = WORK( I ) + ABSA;
               } // 90
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 100
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               zlassq(J-1, A( 1, J ), 1, SCALE, SUM );
            } // 110
         } else {
            for (J = 1; J <= N - 1; J++) { // 120
               zlassq(N-J, A( J+1, J ), 1, SCALE, SUM );
            } // 120
         }
         SUM = 2*SUM;
         for (I = 1; I <= N; I++) { // 130
            if ( DBLE( A( I, I ) ) != ZERO ) {
               ABSA = ABS( DBLE( A( I, I ) ) );
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

      ZLANHE = VALUE;
      return;

      // End of ZLANHE

      }
