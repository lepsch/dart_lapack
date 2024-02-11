      double zlanhe(final int NORM, final int UPLO, final int N, final Matrix<double> A, final int LDA, final Array<double> WORK,) {

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
      // INTRINSIC ABS, DBLE, SQRT

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  SUM = ( A( I, J ) ).abs();
                  if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               } // 10
               SUM = ABS( (A( J, J )).toDouble() );
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( (A( J, J )).toDouble() );
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               for (I = J + 1; I <= N; I++) { // 30
                  SUM = ( A( I, J ) ).abs();
                  if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               } // 30
            } // 40
         }
      } else if ( ( lsame( NORM, 'I' ) ) || ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ( A( I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 50
               WORK[J] = SUM + ABS( (A( J, J )).toDouble() );
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
               SUM = WORK( J ) + ABS( (A( J, J )).toDouble() );
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
         for (I = 1; I <= N; I++) { // 130
            if ( (A( I, I )).toDouble() != ZERO ) {
               ABSA = ABS( (A( I, I )).toDouble() );
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
      }
