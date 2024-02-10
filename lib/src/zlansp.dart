      double zlansp(NORM, UPLO, N, AP, final Array<double> WORK) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM, UPLO;
      int                N;
      double             WORK( * );
      Complex         AP( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J, K;
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
      // INTRINSIC ABS, DBLE, DIMAG, SQRT

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            K = 1;
            for (J = 1; J <= N; J++) { // 20
               for (I = K; I <= K + J - 1; I++) { // 10
                  SUM = ( AP( I ) ).abs();
                  if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               } // 10
               K = K + J;
            } // 20
         } else {
            K = 1;
            for (J = 1; J <= N; J++) { // 40
               for (I = K; I <= K + N - J; I++) { // 30
                  SUM = ( AP( I ) ).abs();
                  if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
               } // 30
               K = K + N - J + 1;
            } // 40
         }
      } else if ( ( lsame( NORM, 'I' ) ) || ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is symmetric).

         VALUE = ZERO;
         K = 1;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ( AP( K ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
                  K = K + 1;
               } // 50
               WORK[J] = SUM + ( AP( K ) ).abs();
               K = K + 1;
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
               SUM = WORK( J ) + ( AP( K ) ).abs();
               K = K + 1;
               for (I = J + 1; I <= N; I++) { // 90
                  ABSA = ( AP( K ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
                  K = K + 1;
               } // 90
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
            } // 100
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         K = 2;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               zlassq(J-1, AP( K ), 1, SCALE, SUM );
               K = K + J;
            } // 110
         } else {
            for (J = 1; J <= N - 1; J++) { // 120
               zlassq(N-J, AP( K ), 1, SCALE, SUM );
               K = K + N - J + 1;
            } // 120
         }
         SUM = 2*SUM;
         K = 1;
         for (I = 1; I <= N; I++) { // 130
            if ( (AP( K )).toDouble() != ZERO ) {
               ABSA = ABS( (AP( K )).toDouble() );
               if ( SCALE < ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2;
                  SCALE = ABSA;
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2;
               }
            }
            if ( DIMAG( AP( K ) ) != ZERO ) {
               ABSA = ABS( DIMAG( AP( K ) ) );
               if ( SCALE < ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2;
                  SCALE = ABSA;
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2;
               }
            }
            if ( lsame( UPLO, 'U' ) ) {
               K = K + I + 1;
            } else {
               K = K + N - I + 1;
            }
         } // 130
         VALUE = SCALE*sqrt( SUM );
      }

      ZLANSP = VALUE;
      }
