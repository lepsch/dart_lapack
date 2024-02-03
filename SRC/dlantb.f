      double dlantb(NORM, UPLO, DIAG, N, K, AB, LDAB, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                K, LDAB, N;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J, L;
      double             SCALE, SUM, VALUE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ
      // ..
      // .. External Functions ..
      //- bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         if ( LSAME( DIAG, 'U' ) ) {
            VALUE = ONE;
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 20
                  DO 10 I = max( K+2-J, 1 ), K;
                     SUM = ( AB( I, J ) ).abs();
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 10
               } // 20
            } else {
               for (J = 1; J <= N; J++) { // 40
                  for (I = 2; I <= min( N+1-J, K+1 ); I++) { // 30
                     SUM = ( AB( I, J ) ).abs();
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 30
               } // 40
            }
         } else {
            VALUE = ZERO;
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 60
                  DO 50 I = max( K+2-J, 1 ), K + 1;
                     SUM = ( AB( I, J ) ).abs();
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 50
               } // 60
            } else {
               for (J = 1; J <= N; J++) { // 80
                  for (I = 1; I <= min( N+1-J, K+1 ); I++) { // 70
                     SUM = ( AB( I, J ) ).abs();
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 70
               } // 80
            }
         }
      } else if ( ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO;
         UDIAG = LSAME( DIAG, 'U' );
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 110
               if ( UDIAG ) {
                  SUM = ONE;
                  DO 90 I = max( K+2-J, 1 ), K;
                     SUM = SUM + ( AB( I, J ) ).abs();
                  } // 90
               } else {
                  SUM = ZERO;
                  DO 100 I = max( K+2-J, 1 ), K + 1;
                     SUM = SUM + ( AB( I, J ) ).abs();
                  } // 100
               }
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 110
         } else {
            for (J = 1; J <= N; J++) { // 140
               if ( UDIAG ) {
                  SUM = ONE;
                  for (I = 2; I <= min( N+1-J, K+1 ); I++) { // 120
                     SUM = SUM + ( AB( I, J ) ).abs();
                  } // 120
               } else {
                  SUM = ZERO;
                  for (I = 1; I <= min( N+1-J, K+1 ); I++) { // 130
                     SUM = SUM + ( AB( I, J ) ).abs();
                  } // 130
               }
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 140
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         VALUE = ZERO;
         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               for (I = 1; I <= N; I++) { // 150
                  WORK( I ) = ONE;
               } // 150
               for (J = 1; J <= N; J++) { // 170
                  L = K + 1 - J;
                  DO 160 I = max( 1, J-K ), J - 1;
                     WORK( I ) = WORK( I ) + ( AB( L+I, J ) ).abs();
                  } // 160
               } // 170
            } else {
               for (I = 1; I <= N; I++) { // 180
                  WORK( I ) = ZERO;
               } // 180
               for (J = 1; J <= N; J++) { // 200
                  L = K + 1 - J;
                  DO 190 I = max( 1, J-K ), J;
                     WORK( I ) = WORK( I ) + ( AB( L+I, J ) ).abs();
                  } // 190
               } // 200
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               for (I = 1; I <= N; I++) { // 210
                  WORK( I ) = ONE;
               } // 210
               for (J = 1; J <= N; J++) { // 230
                  L = 1 - J;
                  DO 220 I = J + 1, min( N, J+K );
                     WORK( I ) = WORK( I ) + ( AB( L+I, J ) ).abs();
                  } // 220
               } // 230
            } else {
               for (I = 1; I <= N; I++) { // 240
                  WORK( I ) = ZERO;
               } // 240
               for (J = 1; J <= N; J++) { // 260
                  L = 1 - J;
                  for (I = J; I <= min( N, J+K ); I++) { // 250
                     WORK( I ) = WORK( I ) + ( AB( L+I, J ) ).abs();
                  } // 250
               } // 260
            }
         }
         for (I = 1; I <= N; I++) { // 270
            SUM = WORK( I );
            if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
         } // 270
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = N;
               if ( K > 0 ) {
                  for (J = 2; J <= N; J++) { // 280
                     dlassq(min( J-1, K ), AB( max( K+2-J, 1 ), J ), 1, SCALE, SUM );
                  } // 280
               }
            } else {
               SCALE = ZERO;
               SUM = ONE;
               for (J = 1; J <= N; J++) { // 290
                  dlassq(min( J, K+1 ), AB( max( K+2-J, 1 ), J ), 1, SCALE, SUM );
               } // 290
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = N;
               if ( K > 0 ) {
                  for (J = 1; J <= N - 1; J++) { // 300
                     dlassq(min( N-J, K ), AB( 2, J ), 1, SCALE, SUM );
                  } // 300
               }
            } else {
               SCALE = ZERO;
               SUM = ONE;
               for (J = 1; J <= N; J++) { // 310
                  dlassq(min( N-J+1, K+1 ), AB( 1, J ), 1, SCALE, SUM );
               } // 310
            }
         }
         VALUE = SCALE*sqrt( SUM );
      }

      DLANTB = VALUE;
      return;
      }
