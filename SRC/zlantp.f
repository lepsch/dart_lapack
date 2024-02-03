      double           FUNCTION ZLANTP( NORM, UPLO, DIAG, N, AP, WORK );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      double             WORK( * );
      COMPLEX*16         AP( * )
      // ..

* =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J, K;
      double             SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         K = 1
         if ( LSAME( DIAG, 'U' ) ) {
            VALUE = ONE
            if ( LSAME( UPLO, 'U' ) ) {
               DO 20 J = 1, N
                  DO 10 I = K, K + J - 2
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10             CONTINUE
                  K = K + J
   20          CONTINUE
            } else {
               DO 40 J = 1, N
                  DO 30 I = K + 1, K + N - J
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30             CONTINUE
                  K = K + N - J + 1
   40          CONTINUE
            }
         } else {
            VALUE = ZERO
            if ( LSAME( UPLO, 'U' ) ) {
               DO 60 J = 1, N
                  DO 50 I = K, K + J - 1
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   50             CONTINUE
                  K = K + J
   60          CONTINUE
            } else {
               DO 80 J = 1, N
                  DO 70 I = K, K + N - J
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70             CONTINUE
                  K = K + N - J + 1
   80          CONTINUE
            }
         }
      } else if ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find norm1(A).

         VALUE = ZERO
         K = 1
         UDIAG = LSAME( DIAG, 'U' )
         if ( LSAME( UPLO, 'U' ) ) {
            DO 110 J = 1, N
               if ( UDIAG ) {
                  SUM = ONE
                  DO 90 I = K, K + J - 2
                     SUM = SUM + ABS( AP( I ) )
   90             CONTINUE
               } else {
                  SUM = ZERO
                  DO 100 I = K, K + J - 1
                     SUM = SUM + ABS( AP( I ) )
  100             CONTINUE
               }
               K = K + J
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  110       CONTINUE
         } else {
            DO 140 J = 1, N
               if ( UDIAG ) {
                  SUM = ONE
                  DO 120 I = K + 1, K + N - J
                     SUM = SUM + ABS( AP( I ) )
  120             CONTINUE
               } else {
                  SUM = ZERO
                  DO 130 I = K, K + N - J
                     SUM = SUM + ABS( AP( I ) )
  130             CONTINUE
               }
               K = K + N - J + 1
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  140       CONTINUE
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         K = 1
         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               DO 150 I = 1, N
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  DO 160 I = 1, J - 1
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  160             CONTINUE
                  K = K + 1
  170          CONTINUE
            } else {
               DO 180 I = 1, N
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  DO 190 I = 1, J
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  190             CONTINUE
  200          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               DO 210 I = 1, N
                  WORK( I ) = ONE
  210          CONTINUE
               DO 230 J = 1, N
                  K = K + 1
                  DO 220 I = J + 1, N
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  220             CONTINUE
  230          CONTINUE
            } else {
               DO 240 I = 1, N
                  WORK( I ) = ZERO
  240          CONTINUE
               DO 260 J = 1, N
                  DO 250 I = J, N
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  250             CONTINUE
  260          CONTINUE
            }
         }
         VALUE = ZERO
         DO 270 I = 1, N
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  270    CONTINUE
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = N
               K = 2
               DO 280 J = 2, N
                  zlassq(J-1, AP( K ), 1, SCALE, SUM );
                  K = K + J
  280          CONTINUE
            } else {
               SCALE = ZERO
               SUM = ONE
               K = 1
               DO 290 J = 1, N
                  zlassq(J, AP( K ), 1, SCALE, SUM );
                  K = K + J
  290          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = N
               K = 2
               DO 300 J = 1, N - 1
                  zlassq(N-J, AP( K ), 1, SCALE, SUM );
                  K = K + N - J + 1
  300          CONTINUE
            } else {
               SCALE = ZERO
               SUM = ONE
               K = 1
               DO 310 J = 1, N
                  zlassq(N-J+1, AP( K ), 1, SCALE, SUM );
                  K = K + N - J + 1
  310          CONTINUE
            }
         }
         VALUE = SCALE*SQRT( SUM )
      }

      ZLANTP = VALUE
      RETURN

      // End of ZLANTP

      }
