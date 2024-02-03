      double           FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
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
      int                I, J, K;
      double             ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            K = 0
            for (J = 1; J <= N; J++) { // 20
               DO 10 I = K + 1, K + J - 1
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               } // 10
               K = K + J
               SUM = ABS( DBLE( AP( K ) ) )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            } // 20
         } else {
            K = 1
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( DBLE( AP( K ) ) )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               DO 30 I = K + 1, K + N - J
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               } // 30
               K = K + N - J + 1
            } // 40
         }
      } else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO
         K = 1
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
               } // 50
               WORK( J ) = SUM + ABS( DBLE( AP( K ) ) )
               K = K + 1
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK( I ) = ZERO
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( DBLE( AP( K ) ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
               } // 90
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
            } // 100
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         K = 2
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               zlassq(J-1, AP( K ), 1, SCALE, SUM );
               K = K + J
            } // 110
         } else {
            DO 120 J = 1, N - 1
               zlassq(N-J, AP( K ), 1, SCALE, SUM );
               K = K + N - J + 1
            } // 120
         }
         SUM = 2*SUM
         K = 1
         for (I = 1; I <= N; I++) { // 130
            if ( DBLE( AP( K ) ).NE.ZERO ) {
               ABSA = ABS( DBLE( AP( K ) ) )
               if ( SCALE.LT.ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2
               }
            }
            if ( LSAME( UPLO, 'U' ) ) {
               K = K + I + 1
            } else {
               K = K + N - I + 1
            }
         } // 130
         VALUE = SCALE*SQRT( SUM )
      }

      ZLANHP = VALUE
      RETURN

      // End of ZLANHP

      }
