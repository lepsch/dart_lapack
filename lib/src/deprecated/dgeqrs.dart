      void dgeqrs(M, N, NRHS, final Matrix<double> A, final int LDA, TAU, final Matrix<double> B, final int LDB, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      double             A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DORMQR, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

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
         xerbla('DGEQRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) return;

      // B := Q' * B

      dormqr('Left', 'Transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      // Solve R*X = B(1:n,:)

      dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

      }
