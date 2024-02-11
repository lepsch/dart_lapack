      void sgelqs(final int M, final int N, final int NRHS, final Matrix<double> A, final int LDA, final int TAU, final Matrix<double> B, final int LDB, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      double               A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SORMLQ, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || M > N ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < 1 || LWORK < NRHS && M > 0 && N > 0 ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('SGELQS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) return;

      // Solve L*X = B(1:m,:)

      strsm('Left', 'Lower', 'No transpose', 'Non-unit', M, NRHS, ONE, A, LDA, B, LDB );

      // Set B(m+1:n,:) to zero

      if (M < N) slaset( 'Full', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB );

      // B := Q' * B

      sormlq('Left', 'Transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      }
