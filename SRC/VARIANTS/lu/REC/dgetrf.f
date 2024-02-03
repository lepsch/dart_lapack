      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.X) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO, NEGONE;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      PARAMETER          ( NEGONE = -1.0D+0 )
      // ..
      // .. Local Scalars ..
      double             SFMIN, TMP;
      int                I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
      int                KSTART, IPIVSTART, JPIVSTART, KCOLS;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      int                IDAMAX;
      bool               DISNAN;
      // EXTERNAL DLAMCH, IDAMAX, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTRSM, DSCAL, XERBLA, DLASWP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, IAND
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      // Compute machine safe minimum
*
      SFMIN = DLAMCH( 'S' )
*
      NSTEP = MIN( M, N )
      DO J = 1, NSTEP
         KAHEAD = IAND( J, -J )
         KSTART = J + 1 - KAHEAD
         KCOLS = MIN( KAHEAD, M-J )
*
         // Find pivot.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP

         // Permute just this column.
         IF (JP .NE. J) THEN
            TMP = A( J, J )
            A( J, J ) = A( JP, J )
            A( JP, J ) = TMP
         END IF

         // Apply pending permutations to L
         NTOPIV = 1
         IPIVSTART = J
         JPIVSTART = J - NTOPIV
         DO WHILE ( NTOPIV .LT. KAHEAD )
            CALL DLASWP( NTOPIV, A( 1, JPIVSTART ), LDA, IPIVSTART, J, IPIV, 1 )
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
         END DO

         // Permute U block to match L
         CALL DLASWP( KCOLS, A( 1,J+1 ), LDA, KSTART, J, IPIV, 1 )

         // Factor the current column
         IF( A( J, J ).NE.ZERO .AND. .NOT.DISNAN( A( J, J ) ) ) THEN
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                 DO I = 1, M-J
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
                 END DO
               END IF
         ELSE IF( A( J,J ) .EQ. ZERO .AND. INFO .EQ. 0 ) THEN
            INFO = J
         END IF

         // Solve for U block.
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', KAHEAD, KCOLS, ONE, A( KSTART, KSTART ), LDA, A( KSTART, J+1 ), LDA )
         // Schur complement.
         CALL DGEMM( 'No transpose', 'No transpose', M-J, KCOLS, KAHEAD, NEGONE, A( J+1, KSTART ), LDA, A( KSTART, J+1 ), LDA, ONE, A( J+1, J+1 ), LDA )
      END DO

      // Handle pivot permutations on the way out of the recursion
      NPIVED = IAND( NSTEP, -NSTEP )
      J = NSTEP - NPIVED
      DO WHILE ( J .GT. 0 )
         NTOPIV = IAND( J, -J )
         CALL DLASWP( NTOPIV, A( 1, J-NTOPIV+1 ), LDA, J+1, NSTEP, IPIV, 1 )
         J = J - NTOPIV
      END DO

      // If short and wide, handle the rest of the columns.
      IF ( M .LT. N ) THEN
         CALL DLASWP( N-M, A( 1, M+KCOLS+1 ), LDA, 1, M, IPIV, 1 )
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', M, N-M, ONE, A, LDA, A( 1,M+KCOLS+1 ), LDA )
      END IF

      RETURN
*
      // End of DGETRF
*
      END
