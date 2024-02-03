      SUBROUTINE SGERQS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SORMRQ, STRSM, XERBLA
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
         xerbla('SGERQS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0 || M == 0) RETURN;

      // Solve R*X = B(n-m+1:n,:)

      strsm('Left', 'Upper', 'No transpose', 'Non-unit', M, NRHS, ONE, A( 1, N-M+1 ), LDA, B( N-M+1, 1 ), LDB );

      // Set B(1:n-m,:) to zero

      slaset('Full', N-M, NRHS, ZERO, ZERO, B, LDB );

      // B := Q' * B

      sormrq('Left', 'Transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      RETURN

      // End of SGERQS

      }
