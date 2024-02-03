      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, M, N;
      COMPLEX*16         ALPHA, BETA
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
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

      IF( LSAME( UPLO, 'U' ) ) THEN

         // Set the diagonal to BETA and the strictly upper triangular
         // part of the array to ALPHA.

         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE

      ELSE IF( LSAME( UPLO, 'L' ) ) THEN

         // Set the diagonal to BETA and the strictly lower triangular
         // part of the array to ALPHA.

         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE

      ELSE

         // Set the array to BETA on the diagonal and ALPHA on the
         // offdiagonal.

         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF

      RETURN

      // End of ZLASET

      END
