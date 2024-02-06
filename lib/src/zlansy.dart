      double zlansy(NORM, UPLO, N, A, LDA, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM, UPLO;
      int                LDA, N;
      double             WORK( * );
      Complex         A( LDA, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J;
      double             ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      // EXTERNAL lsame, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J; I++) { // 10
                  SUM = ( A( I, J ) ).abs();
                  if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               } // 10
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               for (I = J; I <= N; I++) { // 30
                  SUM = ( A( I, J ) ).abs();
                  if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               } // 30
            } // 40
         }
      } else if ( ( lsame( NORM, 'I' ) ) || ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is symmetric).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ( A( I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 50
               WORK[J] = SUM + ( A( J, J ) ).abs();
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I );
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK[I] = ZERO;
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ( A( J, J ) ).abs();
               for (I = J + 1; I <= N; I++) { // 90
                  ABSA = ( A( I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 90
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
            } // 100
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               zlassq(J-1, A( 1, J ), 1, SCALE, SUM );
            } // 110
         } else {
            for (J = 1; J <= N - 1; J++) { // 120
               zlassq(N-J, A( J+1, J ), 1, SCALE, SUM );
            } // 120
         }
         SUM = 2*SUM;
         zlassq(N, A, LDA+1, SCALE, SUM );
         VALUE = SCALE*sqrt( SUM );
      }

      ZLANSY = VALUE;
      return;
      }
