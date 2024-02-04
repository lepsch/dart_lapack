      void sgeqrs(M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORMQR, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || N > M ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -8;
      } else if ( LWORK < 1 || LWORK < NRHS && M > 0 && N > 0 ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('SGEQRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) return;

      // B := Q' * B

      sormqr('Left', 'Transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      // Solve R*X = B(1:n,:)

      strsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

      return;
      }
