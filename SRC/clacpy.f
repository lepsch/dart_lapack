      SUBROUTINE CLACPY( UPLO, M, N, A, LDA, B, LDB )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDB, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE

      } else if ( LSAME( UPLO, 'L' ) ) {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= M; I++) { // 30
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE

      } else {
         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= M; I++) { // 50
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      }

      RETURN

      // End of CLACPY

      }
