      void cgerqs(M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CTRSM, CUNMRQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

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
         xerbla('CGERQS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) return;

      // Solve R*X = B(n-m+1:n,:)

      ctrsm('Left', 'Upper', 'No transpose', 'Non-unit', M, NRHS, CONE, A( 1, N-M+1 ), LDA, B( N-M+1, 1 ), LDB );

      // Set B(1:n-m,:) to zero

      claset('Full', N-M, NRHS, CZERO, CZERO, B, LDB );

      // B := Q' * B

      cunmrq('Left', 'Conjugate transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      return;
      }
