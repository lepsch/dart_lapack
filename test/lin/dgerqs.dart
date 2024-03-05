      void dgerqs(final int M, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final int TAU, final Matrix<double> B_, final int LDB, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.having();
  final B = B_.having();
  final WORK = WORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      double             A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DORMRQ, DTRSM, XERBLA
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
         xerbla('DGERQS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) return;

      // Solve R*X = B(n-m+1:n,:)

      dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', M, NRHS, ONE, A( 1, N-M+1 ), LDA, B( N-M+1, 1 ), LDB );

      // Set B(1:n-m,:) to zero

      dlaset('Full', N-M, NRHS, ZERO, ZERO, B, LDB );

      // B := Q' * B

      dormrq('Left', 'Transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      }
