      double dlantr(NORM, UPLO, DIAG, M, N, A, LDA, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J;
      double             SCALE, SUM, VALUE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( min( M, N ) == 0 ) {
         VALUE = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         if ( LSAME( DIAG, 'U' ) ) {
            VALUE = ONE;
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 20
                  DO 10 I = 1, min( M, J-1 );
                     SUM = ABS( A( I, J ) );
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 10
               } // 20
            } else {
               for (J = 1; J <= N; J++) { // 40
                  for (I = J + 1; I <= M; I++) { // 30
                     SUM = ABS( A( I, J ) );
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 30
               } // 40
            }
         } else {
            VALUE = ZERO;
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 60
                  DO 50 I = 1, min( M, J );
                     SUM = ABS( A( I, J ) );
                     if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
                  } // 50
               } // 60
            } else {
               for (J = 1; J <= N; J++) { // 80
                  for (I = J; I <= M; I++) { // 70
                     SUM = ABS( A( I, J ) );
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
               if ( ( UDIAG ) && ( J <= M ) ) {
                  SUM = ONE;
                  for (I = 1; I <= J - 1; I++) { // 90
                     SUM = SUM + ABS( A( I, J ) );
                  } // 90
               } else {
                  SUM = ZERO;
                  DO 100 I = 1, min( M, J );
                     SUM = SUM + ABS( A( I, J ) );
                  } // 100
               }
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 110
         } else {
            for (J = 1; J <= N; J++) { // 140
               if ( UDIAG ) {
                  SUM = ONE;
                  for (I = J + 1; I <= M; I++) { // 120
                     SUM = SUM + ABS( A( I, J ) );
                  } // 120
               } else {
                  SUM = ZERO;
                  for (I = J; I <= M; I++) { // 130
                     SUM = SUM + ABS( A( I, J ) );
                  } // 130
               }
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 140
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               for (I = 1; I <= M; I++) { // 150
                  WORK( I ) = ONE;
               } // 150
               for (J = 1; J <= N; J++) { // 170
                  DO 160 I = 1, min( M, J-1 );
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) );
                  } // 160
               } // 170
            } else {
               for (I = 1; I <= M; I++) { // 180
                  WORK( I ) = ZERO;
               } // 180
               for (J = 1; J <= N; J++) { // 200
                  DO 190 I = 1, min( M, J );
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) );
                  } // 190
               } // 200
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               DO 210 I = 1, min( M, N );
                  WORK( I ) = ONE;
               } // 210
               for (I = N + 1; I <= M; I++) { // 220
                  WORK( I ) = ZERO;
               } // 220
               for (J = 1; J <= N; J++) { // 240
                  for (I = J + 1; I <= M; I++) { // 230
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) );
                  } // 230
               } // 240
            } else {
               for (I = 1; I <= M; I++) { // 250
                  WORK( I ) = ZERO;
               } // 250
               for (J = 1; J <= N; J++) { // 270
                  for (I = J; I <= M; I++) { // 260
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) );
                  } // 260
               } // 270
            }
         }
         VALUE = ZERO;
         for (I = 1; I <= M; I++) { // 280
            SUM = WORK( I );
            if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
         } // 280
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = min( M, N );
               for (J = 2; J <= N; J++) { // 290
                  dlassq(min( M, J-1 ), A( 1, J ), 1, SCALE, SUM );
               } // 290
            } else {
               SCALE = ZERO;
               SUM = ONE;
               for (J = 1; J <= N; J++) { // 300
                  dlassq(min( M, J ), A( 1, J ), 1, SCALE, SUM );
               } // 300
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = min( M, N );
               for (J = 1; J <= N; J++) { // 310
                  dlassq(M-J, A( min( M, J+1 ), J ), 1, SCALE, SUM );
               } // 310
            } else {
               SCALE = ZERO;
               SUM = ONE;
               for (J = 1; J <= N; J++) { // 320
                  dlassq(M-J+1, A( J, J ), 1, SCALE, SUM );
               } // 320
            }
         }
         VALUE = SCALE*sqrt( SUM );
      }

      DLANTR = VALUE;
      return;
      }
