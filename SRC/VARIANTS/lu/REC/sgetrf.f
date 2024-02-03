      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine (version 3.X) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO, NEGONE
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      const              NEGONE = -1.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               SFMIN, TMP
      int                I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
      int                KSTART, IPIVSTART, JPIVSTART, KCOLS;
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      int                ISAMAX;
      bool               SISNAN;
      // EXTERNAL SLAMCH, ISAMAX, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL STRSM, SSCAL, XERBLA, SLASWP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, IAND
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('SGETRF', -INFO );
         RETURN
      }

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

         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP

         // Permute just this column.
         if (JP .NE. J) {
            TMP = A( J, J )
            A( J, J ) = A( JP, J )
            A( JP, J ) = TMP
         }

         // Apply pending permutations to L
         NTOPIV = 1
         IPIVSTART = J
         JPIVSTART = J - NTOPIV
         DO WHILE ( NTOPIV .LT. KAHEAD )
            slaswp(NTOPIV, A( 1, JPIVSTART ), LDA, IPIVSTART, J, IPIV, 1 );
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
         END DO

         // Permute U block to match L
         slaswp(KCOLS, A( 1,J+1 ), LDA, KSTART, J, IPIV, 1 );

         // Factor the current column
         if ( A( J, J ).NE.ZERO .AND. .NOT.SISNAN( A( J, J ) ) ) {
               if ( ABS(A( J, J )) .GE. SFMIN ) {
                  sscal(M-J, ONE / A( J, J ), A( J+1, J ), 1 );
               } else {
                 DO I = 1, M-J
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
                 END DO
               }
         } else if ( A( J,J ) .EQ. ZERO .AND. INFO .EQ. 0 ) {
            INFO = J
         }

         // Solve for U block.
         strsm('Left', 'Lower', 'No transpose', 'Unit', KAHEAD, KCOLS, ONE, A( KSTART, KSTART ), LDA, A( KSTART, J+1 ), LDA );
         // Schur complement.
         sgemm('No transpose', 'No transpose', M-J, KCOLS, KAHEAD, NEGONE, A( J+1, KSTART ), LDA, A( KSTART, J+1 ), LDA, ONE, A( J+1, J+1 ), LDA );
      END DO

      // Handle pivot permutations on the way out of the recursion
      NPIVED = IAND( NSTEP, -NSTEP )
      J = NSTEP - NPIVED
      DO WHILE ( J .GT. 0 )
         NTOPIV = IAND( J, -J )
         slaswp(NTOPIV, A( 1, J-NTOPIV+1 ), LDA, J+1, NSTEP, IPIV, 1 );
         J = J - NTOPIV
      END DO

      // If short and wide, handle the rest of the columns.
      if ( M .LT. N ) {
         slaswp(N-M, A( 1, M+KCOLS+1 ), LDA, 1, M, IPIV, 1 );
         strsm('Left', 'Lower', 'No transpose', 'Unit', M, N-M, ONE, A, LDA, A( 1,M+KCOLS+1 ), LDA );
      }

      RETURN

      // End of SGETRF

      }
