      double zlantp(final int NORM, final int UPLO, final int DIAG, final int N, final int AP, final Array<double> WORK,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, NORM, UPLO;
      int                N;
      double             WORK( * );
      Complex         AP( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UDIAG;
      int                I, J, K;
      double             SCALE, SUM, VALUE;
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

         K = 1;
         if ( lsame( DIAG, 'U' ) ) {
            VALUE = ONE;
            if ( lsame( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 20
                  for (I = K; I <= K + J - 2; I++) { // 10
                     SUM = ( AP( I ) ).abs();
                     if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
                  } // 10
                  K = K + J;
               } // 20
            } else {
               for (J = 1; J <= N; J++) { // 40
                  for (I = K + 1; I <= K + N - J; I++) { // 30
                     SUM = ( AP( I ) ).abs();
                     if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
                  } // 30
                  K = K + N - J + 1;
               } // 40
            }
         } else {
            VALUE = ZERO;
            if ( lsame( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 60
                  for (I = K; I <= K + J - 1; I++) { // 50
                     SUM = ( AP( I ) ).abs();
                     if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
                  } // 50
                  K = K + J;
               } // 60
            } else {
               for (J = 1; J <= N; J++) { // 80
                  for (I = K; I <= K + N - J; I++) { // 70
                     SUM = ( AP( I ) ).abs();
                     if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
                  } // 70
                  K = K + N - J + 1;
               } // 80
            }
         }
      } else if ( ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO;
         K = 1;
         UDIAG = lsame( DIAG, 'U' );
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 110
               if ( UDIAG ) {
                  SUM = ONE;
                  for (I = K; I <= K + J - 2; I++) { // 90
                     SUM = SUM + ( AP( I ) ).abs();
                  } // 90
               } else {
                  SUM = ZERO;
                  for (I = K; I <= K + J - 1; I++) { // 100
                     SUM = SUM + ( AP( I ) ).abs();
                  } // 100
               }
               K = K + J;
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
            } // 110
         } else {
            for (J = 1; J <= N; J++) { // 140
               if ( UDIAG ) {
                  SUM = ONE;
                  for (I = K + 1; I <= K + N - J; I++) { // 120
                     SUM = SUM + ( AP( I ) ).abs();
                  } // 120
               } else {
                  SUM = ZERO;
                  for (I = K; I <= K + N - J; I++) { // 130
                     SUM = SUM + ( AP( I ) ).abs();
                  } // 130
               }
               K = K + N - J + 1;
               if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
            } // 140
         }
      } else if ( lsame( NORM, 'I' ) ) {

         // Find normI(A).

         K = 1;
         if ( lsame( UPLO, 'U' ) ) {
            if ( lsame( DIAG, 'U' ) ) {
               for (I = 1; I <= N; I++) { // 150
                  WORK[I] = ONE;
               } // 150
               for (J = 1; J <= N; J++) { // 170
                  for (I = 1; I <= J - 1; I++) { // 160
                     WORK[I] = WORK( I ) + ( AP( K ) ).abs();
                     K = K + 1;
                  } // 160
                  K = K + 1;
               } // 170
            } else {
               for (I = 1; I <= N; I++) { // 180
                  WORK[I] = ZERO;
               } // 180
               for (J = 1; J <= N; J++) { // 200
                  for (I = 1; I <= J; I++) { // 190
                     WORK[I] = WORK( I ) + ( AP( K ) ).abs();
                     K = K + 1;
                  } // 190
               } // 200
            }
         } else {
            if ( lsame( DIAG, 'U' ) ) {
               for (I = 1; I <= N; I++) { // 210
                  WORK[I] = ONE;
               } // 210
               for (J = 1; J <= N; J++) { // 230
                  K = K + 1;
                  for (I = J + 1; I <= N; I++) { // 220
                     WORK[I] = WORK( I ) + ( AP( K ) ).abs();
                     K = K + 1;
                  } // 220
               } // 230
            } else {
               for (I = 1; I <= N; I++) { // 240
                  WORK[I] = ZERO;
               } // 240
               for (J = 1; J <= N; J++) { // 260
                  for (I = J; I <= N; I++) { // 250
                     WORK[I] = WORK( I ) + ( AP( K ) ).abs();
                     K = K + 1;
                  } // 250
               } // 260
            }
         }
         VALUE = ZERO;
         for (I = 1; I <= N; I++) { // 270
            SUM = WORK( I );
            if( VALUE < SUM || disnan( SUM ) ) VALUE = SUM;
         } // 270
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( lsame( UPLO, 'U' ) ) {
            if ( lsame( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = N;
               K = 2;
               for (J = 2; J <= N; J++) { // 280
                  zlassq(J-1, AP( K ), 1, SCALE, SUM );
                  K = K + J;
               } // 280
            } else {
               SCALE = ZERO;
               SUM = ONE;
               K = 1;
               for (J = 1; J <= N; J++) { // 290
                  zlassq(J, AP( K ), 1, SCALE, SUM );
                  K = K + J;
               } // 290
            }
         } else {
            if ( lsame( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = N;
               K = 2;
               for (J = 1; J <= N - 1; J++) { // 300
                  zlassq(N-J, AP( K ), 1, SCALE, SUM );
                  K = K + N - J + 1;
               } // 300
            } else {
               SCALE = ZERO;
               SUM = ONE;
               K = 1;
               for (J = 1; J <= N; J++) { // 310
                  zlassq(N-J+1, AP( K ), 1, SCALE, SUM );
                  K = K + N - J + 1;
               } // 310
            }
         }
         VALUE = SCALE*sqrt( SUM );
      }

      ZLANTP = VALUE;
      }
