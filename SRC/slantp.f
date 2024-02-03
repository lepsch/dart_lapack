      REAL             FUNCTION SLANTP( NORM, UPLO, DIAG, N, AP, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), WORK( * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J, K;
      REAL               SCALE, SUM, VALUE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
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
               for (J = 1; J <= N; J++) { // 20
                  DO 10 I = K, K + J - 2
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   10             CONTINUE
                  K = K + J
   20          CONTINUE
            } else {
               for (J = 1; J <= N; J++) { // 40
                  DO 30 I = K + 1, K + N - J
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   30             CONTINUE
                  K = K + N - J + 1
   40          CONTINUE
            }
         } else {
            VALUE = ZERO
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 60
                  DO 50 I = K, K + J - 1
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   50             CONTINUE
                  K = K + J
   60          CONTINUE
            } else {
               for (J = 1; J <= N; J++) { // 80
                  DO 70 I = K, K + N - J
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
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
            for (J = 1; J <= N; J++) { // 110
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
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  110       CONTINUE
         } else {
            for (J = 1; J <= N; J++) { // 140
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
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  140       CONTINUE
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         K = 1
         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               for (I = 1; I <= N; I++) { // 150
                  WORK( I ) = ONE
  150          CONTINUE
               for (J = 1; J <= N; J++) { // 170
                  DO 160 I = 1, J - 1
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  160             CONTINUE
                  K = K + 1
  170          CONTINUE
            } else {
               for (I = 1; I <= N; I++) { // 180
                  WORK( I ) = ZERO
  180          CONTINUE
               for (J = 1; J <= N; J++) { // 200
                  for (I = 1; I <= J; I++) { // 190
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  190             CONTINUE
  200          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               for (I = 1; I <= N; I++) { // 210
                  WORK( I ) = ONE
  210          CONTINUE
               for (J = 1; J <= N; J++) { // 230
                  K = K + 1
                  DO 220 I = J + 1, N
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  220             CONTINUE
  230          CONTINUE
            } else {
               for (I = 1; I <= N; I++) { // 240
                  WORK( I ) = ZERO
  240          CONTINUE
               for (J = 1; J <= N; J++) { // 260
                  for (I = J; I <= N; I++) { // 250
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  250             CONTINUE
  260          CONTINUE
            }
         }
         VALUE = ZERO
         for (I = 1; I <= N; I++) { // 270
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  270    CONTINUE
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = N
               K = 2
               for (J = 2; J <= N; J++) { // 280
                  slassq(J-1, AP( K ), 1, SCALE, SUM );
                  K = K + J
  280          CONTINUE
            } else {
               SCALE = ZERO
               SUM = ONE
               K = 1
               for (J = 1; J <= N; J++) { // 290
                  slassq(J, AP( K ), 1, SCALE, SUM );
                  K = K + J
  290          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = N
               K = 2
               DO 300 J = 1, N - 1
                  slassq(N-J, AP( K ), 1, SCALE, SUM );
                  K = K + N - J + 1
  300          CONTINUE
            } else {
               SCALE = ZERO
               SUM = ONE
               K = 1
               for (J = 1; J <= N; J++) { // 310
                  slassq(N-J+1, AP( K ), 1, SCALE, SUM );
                  K = K + N - J + 1
  310          CONTINUE
            }
         }
         VALUE = SCALE*SQRT( SUM )
      }

      SLANTP = VALUE
      RETURN

      // End of SLANTP

      }
