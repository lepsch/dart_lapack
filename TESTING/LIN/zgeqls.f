      SUBROUTINE ZGEQLS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), TAU( * ), WORK( LWORK )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTRSM, ZUNMQL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.1 .OR. LWORK.LT.NRHS .AND. M.GT.0 .AND. N.GT.0 ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQLS', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 .OR. M.EQ.0 ) RETURN
*
      // B := Q' * B
*
      CALL ZUNMQL( 'Left', 'Conjugate transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK, LWORK, INFO )
*
      // Solve L*X = B(m-n+1:m,:)
*
      CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A( M-N+1, 1 ), LDA, B( M-N+1, 1 ), LDB )
*
      RETURN
*
      // End of ZGEQLS
*
      END
