      SUBROUTINE CGEQLS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRSM, CUNMQL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 .OR. N.GT.M ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -8
      } else if ( LWORK.LT.1 .OR. LWORK.LT.NRHS && M.GT.0 && N.GT.0 ) {
         INFO = -10
      }
      if ( INFO != 0 ) {
         xerbla('CGEQLS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 .OR. NRHS == 0 .OR. M == 0) RETURN;

      // B := Q' * B

      cunmql('Left', 'Conjugate transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      // Solve L*X = B(m-n+1:m,:)

      ctrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A( M-N+1, 1 ), LDA, B( M-N+1, 1 ), LDB );

      RETURN

      // End of CGEQLS

      }
