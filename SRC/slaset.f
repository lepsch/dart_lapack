      SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, M, N;
      REAL               ALPHA, BETA
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..

* =====================================================================

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

         // Set the strictly upper triangular or trapezoidal part of the
         // array to ALPHA.

         for (J = 2; J <= N; J++) { // 20
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE

      } else if ( LSAME( UPLO, 'L' ) ) {

         // Set the strictly lower triangular or trapezoidal part of the
         // array to ALPHA.

         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE

      } else {

         // Set the leading m-by-n submatrix to ALPHA.

         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= M; I++) { // 50
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      }

      // Set the first min(M,N) diagonal elements to BETA.

      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE

      RETURN

      // End of SLASET

      }
