      SUBROUTINE SPOTRI( UPLO, N, A, LDA, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * );
      // ..

*  =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAUUM, STRTRI, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SPOTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Invert the triangular Cholesky factor U or L.

      strtri(UPLO, 'Non-unit', N, A, LDA, INFO );
      if (INFO > 0) RETURN;

      // Form inv(U) * inv(U)**T or inv(L)**T * inv(L).

      slauum(UPLO, N, A, LDA, INFO );

      return;

      // End of SPOTRI

      }
