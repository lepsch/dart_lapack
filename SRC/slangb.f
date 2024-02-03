      REAL             FUNCTION SLANGB( NORM, N, KL, KU, AB, LDAB, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                KL, KU, LDAB, N;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), WORK( * )
      // ..

* =====================================================================


      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K, L;
      REAL               SCALE, SUM, VALUE, TEMP
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

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         for (J = 1; J <= N; J++) { // 20
            DO 10 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )
               TEMP = ABS( AB( I, J ) )
               IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
   10       CONTINUE
   20    CONTINUE
      } else if ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find norm1(A).

         VALUE = ZERO
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO
            DO 30 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )
               SUM = SUM + ABS( AB( I, J ) )
   30       CONTINUE
            IF( VALUE.LT.SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   40    CONTINUE
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= N; I++) { // 50
            WORK( I ) = ZERO
   50    CONTINUE
         for (J = 1; J <= N; J++) { // 70
            K = KU + 1 - J
            DO 60 I = MAX( 1, J-KU ), MIN( N, J+KL )
               WORK( I ) = WORK( I ) + ABS( AB( K+I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         for (I = 1; I <= N; I++) { // 80
            TEMP = WORK( I )
            IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
   80    CONTINUE
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         for (J = 1; J <= N; J++) { // 90
            L = MAX( 1, J-KU )
            K = KU + 1 - J + L
            slassq(MIN( N, J+KL )-L+1, AB( K, J ), 1, SCALE, SUM );
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      }

      SLANGB = VALUE
      RETURN

      // End of SLANGB

      }
