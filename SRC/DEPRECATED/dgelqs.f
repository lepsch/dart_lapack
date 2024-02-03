      SUBROUTINE DGELQS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DORMLQ, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. M.GT.N ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.1 .OR. LWORK.LT.NRHS .AND. M.GT.0 .AND. N.GT.0 ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELQS', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 .OR. M.EQ.0 ) RETURN
*
      // Solve L*X = B(1:m,:)
*
      CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', M, NRHS, ONE, A, LDA, B, LDB )
*
      // Set B(m+1:n,:) to zero
*
      IF( M.LT.N ) CALL DLASET( 'Full', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
*
      // B := Q' * B
*
      CALL DORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )
*
      RETURN
*
      // End of DGELQS
*
      END
