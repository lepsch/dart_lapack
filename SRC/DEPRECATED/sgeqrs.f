      SUBROUTINE SGEQRS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )

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
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORMQR, STRSM, XERBLA
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
      } else if ( LWORK.LT.1 .OR. LWORK.LT.NRHS .AND. M.GT.0 .AND. N.GT.0 ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         xerbla('SGEQRS', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 .OR. M.EQ.0 ) RETURN

      // B := Q' * B

      sormqr('Left', 'Transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK, LWORK, INFO );

      // Solve R*X = B(1:n,:)

      strsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

      RETURN

      // End of SGEQRS

      }