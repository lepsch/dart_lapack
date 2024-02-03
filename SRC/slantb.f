      REAL             FUNCTION SLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                K, LDAB, N;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), WORK( * );
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J, L;
      REAL               SCALE, SUM, VALUE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
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
                  DO 10 I = MAX( K+2-J, 1 ), K;
                     SUM = ABS( AB( I, J ) );
                     IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
                  } // 10
               } // 20
            } else {
               for (J = 1; J <= N; J++) { // 40
                  DO 30 I = 2, MIN( N+1-J, K+1 );
                     SUM = ABS( AB( I, J ) );
                     IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
                  } // 30
               } // 40
            }
         } else {
            VALUE = ZERO;
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 60
                  DO 50 I = MAX( K+2-J, 1 ), K + 1;
                     SUM = ABS( AB( I, J ) );
                     IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
                  } // 50
               } // 60
            } else {
               for (J = 1; J <= N; J++) { // 80
                  DO 70 I = 1, MIN( N+1-J, K+1 );
                     SUM = ABS( AB( I, J ) );
                     IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
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
                  DO 90 I = MAX( K+2-J, 1 ), K;
                     SUM = SUM + ABS( AB( I, J ) );
                  } // 90
               } else {
                  SUM = ZERO;
                  DO 100 I = MAX( K+2-J, 1 ), K + 1;
                     SUM = SUM + ABS( AB( I, J ) );
                  } // 100
               }
               IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 110
         } else {
            for (J = 1; J <= N; J++) { // 140
               if ( UDIAG ) {
                  SUM = ONE;
                  DO 120 I = 2, MIN( N+1-J, K+1 );
                     SUM = SUM + ABS( AB( I, J ) );
                  } // 120
               } else {
                  SUM = ZERO;
                  DO 130 I = 1, MIN( N+1-J, K+1 );
                     SUM = SUM + ABS( AB( I, J ) );
                  } // 130
               }
               IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
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
                  DO 160 I = MAX( 1, J-K ), J - 1;
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) );
                  } // 160
               } // 170
            } else {
               for (I = 1; I <= N; I++) { // 180
                  WORK( I ) = ZERO;
               } // 180
               for (J = 1; J <= N; J++) { // 200
                  L = K + 1 - J;
                  DO 190 I = MAX( 1, J-K ), J;
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) );
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
                  DO 220 I = J + 1, MIN( N, J+K );
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) );
                  } // 220
               } // 230
            } else {
               for (I = 1; I <= N; I++) { // 240
                  WORK( I ) = ZERO;
               } // 240
               for (J = 1; J <= N; J++) { // 260
                  L = 1 - J;
                  DO 250 I = J, MIN( N, J+K );
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) );
                  } // 250
               } // 260
            }
         }
         for (I = 1; I <= N; I++) { // 270
            SUM = WORK( I );
            IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
         } // 270
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = N;
               if ( K > 0 ) {
                  for (J = 2; J <= N; J++) { // 280
                     slassq(MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SUM );
                  } // 280
               }
            } else {
               SCALE = ZERO;
               SUM = ONE;
               for (J = 1; J <= N; J++) { // 290
                  slassq(MIN( J, K+1 ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SUM );
               } // 290
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE;
               SUM = N;
               if ( K > 0 ) {
                  for (J = 1; J <= N - 1; J++) { // 300
                     slassq(MIN( N-J, K ), AB( 2, J ), 1, SCALE, SUM );
                  } // 300
               }
            } else {
               SCALE = ZERO;
               SUM = ONE;
               for (J = 1; J <= N; J++) { // 310
                  slassq(MIN( N-J+1, K+1 ), AB( 1, J ), 1, SCALE, SUM );
               } // 310
            }
         }
         VALUE = SCALE*SQRT( SUM );
      }

      SLANTB = VALUE;
      RETURN;

      // End of SLANTB

      }
