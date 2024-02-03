      REAL             FUNCTION SLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J;
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
      // INTRINSIC ABS, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( MIN( M, N ).EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         if ( LSAME( DIAG, 'U' ) ) {
            VALUE = ONE
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 20
                  DO 10 I = 1, MIN( M, J-1 )
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   10             CONTINUE
   20          CONTINUE
            } else {
               for (J = 1; J <= N; J++) { // 40
                  DO 30 I = J + 1, M
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   30             CONTINUE
   40          CONTINUE
            }
         } else {
            VALUE = ZERO
            if ( LSAME( UPLO, 'U' ) ) {
               for (J = 1; J <= N; J++) { // 60
                  DO 50 I = 1, MIN( M, J )
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   50             CONTINUE
   60          CONTINUE
            } else {
               for (J = 1; J <= N; J++) { // 80
                  for (I = J; I <= M; I++) { // 70
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   70             CONTINUE
   80          CONTINUE
            }
         }
      } else if ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find norm1(A).

         VALUE = ZERO
         UDIAG = LSAME( DIAG, 'U' )
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 110
               if ( ( UDIAG ) .AND. ( J.LE.M ) ) {
                  SUM = ONE
                  DO 90 I = 1, J - 1
                     SUM = SUM + ABS( A( I, J ) )
   90             CONTINUE
               } else {
                  SUM = ZERO
                  DO 100 I = 1, MIN( M, J )
                     SUM = SUM + ABS( A( I, J ) )
  100             CONTINUE
               }
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  110       CONTINUE
         } else {
            for (J = 1; J <= N; J++) { // 140
               if ( UDIAG ) {
                  SUM = ONE
                  DO 120 I = J + 1, M
                     SUM = SUM + ABS( A( I, J ) )
  120             CONTINUE
               } else {
                  SUM = ZERO
                  for (I = J; I <= M; I++) { // 130
                     SUM = SUM + ABS( A( I, J ) )
  130             CONTINUE
               }
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  140       CONTINUE
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               for (I = 1; I <= M; I++) { // 150
                  WORK( I ) = ONE
  150          CONTINUE
               for (J = 1; J <= N; J++) { // 170
                  DO 160 I = 1, MIN( M, J-1 )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  160             CONTINUE
  170          CONTINUE
            } else {
               for (I = 1; I <= M; I++) { // 180
                  WORK( I ) = ZERO
  180          CONTINUE
               for (J = 1; J <= N; J++) { // 200
                  DO 190 I = 1, MIN( M, J )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  190             CONTINUE
  200          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               DO 210 I = 1, MIN( M, N )
                  WORK( I ) = ONE
  210          CONTINUE
               DO 220 I = N + 1, M
                  WORK( I ) = ZERO
  220          CONTINUE
               for (J = 1; J <= N; J++) { // 240
                  DO 230 I = J + 1, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  230             CONTINUE
  240          CONTINUE
            } else {
               for (I = 1; I <= M; I++) { // 250
                  WORK( I ) = ZERO
  250          CONTINUE
               for (J = 1; J <= N; J++) { // 270
                  for (I = J; I <= M; I++) { // 260
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  260             CONTINUE
  270          CONTINUE
            }
         }
         VALUE = ZERO
         for (I = 1; I <= M; I++) { // 280
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  280    CONTINUE
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = MIN( M, N )
               for (J = 2; J <= N; J++) { // 290
                  slassq(MIN( M, J-1 ), A( 1, J ), 1, SCALE, SUM );
  290          CONTINUE
            } else {
               SCALE = ZERO
               SUM = ONE
               for (J = 1; J <= N; J++) { // 300
                  slassq(MIN( M, J ), A( 1, J ), 1, SCALE, SUM );
  300          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = MIN( M, N )
               for (J = 1; J <= N; J++) { // 310
                  slassq(M-J, A( MIN( M, J+1 ), J ), 1, SCALE, SUM );
  310          CONTINUE
            } else {
               SCALE = ZERO
               SUM = ONE
               for (J = 1; J <= N; J++) { // 320
                  slassq(M-J+1, A( J, J ), 1, SCALE, SUM );
  320          CONTINUE
            }
         }
         VALUE = SCALE*SQRT( SUM )
      }

      SLANTR = VALUE
      RETURN

      // End of SLANTR

      }
