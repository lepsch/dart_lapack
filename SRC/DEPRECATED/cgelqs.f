      SUBROUTINE CGELQS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )

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
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CTRSM, CUNMLQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < 0 || M.GT.N ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8
      } else if ( LWORK < 1 || LWORK < NRHS && M.GT.0 && N.GT.0 ) {
         INFO = -10
      }
      if ( INFO != 0 ) {
         xerbla('CGELQS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) RETURN;

      // Solve L*X = B(1:m,:)

      ctrsm('Left', 'Lower', 'No transpose', 'Non-unit', M, NRHS, CONE, A, LDA, B, LDB );

      // Set B(m+1:n,:) to zero

      if (M < N) CALL CLASET( 'Full', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB );

      // B := Q' * B

      cunmlq('Left', 'Conjugate transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      RETURN

      // End of CGELQS

      }
