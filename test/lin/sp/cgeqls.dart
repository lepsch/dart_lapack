      void cgeqls(final int M, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final int TAU, final Matrix<double> B_, final int LDB, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      Complex            A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRSM, CUNMQL, XERBLA
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
         xerbla('CGEQLS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) return;

      // B := Q' * B

      cunmql('Left', 'Conjugate transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      // Solve L*X = B(m-n+1:m,:)

      ctrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A( M-N+1, 1 ), LDA, B( M-N+1, 1 ), LDB );

      }