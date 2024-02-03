      SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine (version 3.X) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, NEGONE
      REAL               ZERO
      PARAMETER          ( ONE = (1.0E+0, 0.0E+0) )
      PARAMETER          ( NEGONE = (-1.0E+0, 0.0E+0) )
      PARAMETER          ( ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      REAL               SFMIN, PIVMAG
      COMPLEX            TMP
      int                I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
      int                KSTART, IPIVSTART, JPIVSTART, KCOLS;
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      int                ICAMAX;
      bool               SISNAN;
      // EXTERNAL SLAMCH, ICAMAX, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRSM, CSCAL, XERBLA, CLASWP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, IAND, ABS
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Compute machine safe minimum

      SFMIN = SLAMCH( 'S' )

      NSTEP = MIN( M, N )
      DO J = 1, NSTEP
         KAHEAD = IAND( J, -J )
         KSTART = J + 1 - KAHEAD
         KCOLS = MIN( KAHEAD, M-J )

         // Find pivot.

         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 )
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
            CALL CLASWP( NTOPIV, A( 1, JPIVSTART ), LDA, IPIVSTART, J, IPIV, 1 )
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
         END DO

         // Permute U block to match L
         CALL CLASWP( KCOLS, A( 1,J+1 ), LDA, KSTART, J, IPIV, 1 )

         // Factor the current column
         PIVMAG = ABS( A( J, J ) )
         IF( PIVMAG.NE.ZERO .AND. .NOT.SISNAN( PIVMAG ) ) THEN
               IF( PIVMAG .GE. SFMIN ) THEN
                  CALL CSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                 DO I = 1, M-J
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
                 END DO
               END IF
         ELSE IF( PIVMAG .EQ. ZERO .AND. INFO .EQ. 0 ) THEN
            INFO = J
         END IF

         // Solve for U block.
         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', KAHEAD, KCOLS, ONE, A( KSTART, KSTART ), LDA, A( KSTART, J+1 ), LDA )
         // Schur complement.
         CALL CGEMM( 'No transpose', 'No transpose', M-J, KCOLS, KAHEAD, NEGONE, A( J+1, KSTART ), LDA, A( KSTART, J+1 ), LDA, ONE, A( J+1, J+1 ), LDA )
      END DO

      // Handle pivot permutations on the way out of the recursion
      NPIVED = IAND( NSTEP, -NSTEP )
      J = NSTEP - NPIVED
      DO WHILE ( J .GT. 0 )
         NTOPIV = IAND( J, -J )
         CALL CLASWP( NTOPIV, A( 1, J-NTOPIV+1 ), LDA, J+1, NSTEP, IPIV, 1 )
         J = J - NTOPIV
      END DO

      // If short and wide, handle the rest of the columns.
      IF ( M .LT. N ) THEN
         CALL CLASWP( N-M, A( 1, M+KCOLS+1 ), LDA, 1, M, IPIV, 1 )
         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', M, N-M, ONE, A, LDA, A( 1,M+KCOLS+1 ), LDA )
      END IF

      RETURN

      // End of CGETRF

      END
