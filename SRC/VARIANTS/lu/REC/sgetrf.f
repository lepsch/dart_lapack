      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO );

// -- LAPACK computational routine (version 3.X) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO, NEGONE;
      const              ONE = 1.0, ZERO = 0.0 ;
      const              NEGONE = -1.0 ;
      // ..
      // .. Local Scalars ..
      REAL               SFMIN, TMP;
      int                I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
      int                KSTART, IPIVSTART, JPIVSTART, KCOLS;
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
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

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SGETRF', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Compute machine safe minimum

      SFMIN = SLAMCH( 'S' );

      NSTEP = MIN( M, N );
      for (J = 1; J <= NSTEP; J++) {
         KAHEAD = IAND( J, -J );
         KSTART = J + 1 - KAHEAD;
         KCOLS = MIN( KAHEAD, M-J );

         // Find pivot.

         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 );
         IPIV( J ) = JP;

         // Permute just this column.
         if (JP != J) {
            TMP = A( J, J );
            A( J, J ) = A( JP, J );
            A( JP, J ) = TMP;
         }

         // Apply pending permutations to L
         NTOPIV = 1;
         IPIVSTART = J;
         JPIVSTART = J - NTOPIV;
         DO WHILE ( NTOPIV < KAHEAD );
            slaswp(NTOPIV, A( 1, JPIVSTART ), LDA, IPIVSTART, J, IPIV, 1 );
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
         }

         // Permute U block to match L
         slaswp(KCOLS, A( 1,J+1 ), LDA, KSTART, J, IPIV, 1 );

         // Factor the current column
         if ( A( J, J ) != ZERO && !SISNAN( A( J, J ) ) ) {
               if ( ABS(A( J, J )) >= SFMIN ) {
                  sscal(M-J, ONE / A( J, J ), A( J+1, J ), 1 );
               } else {
                 for (I = 1; I <= M-J; I++) {
                    A( J+I, J ) = A( J+I, J ) / A( J, J );
                 }
               }
         } else if ( A( J,J ) == ZERO && INFO == 0 ) {
            INFO = J;
         }

         // Solve for U block.
         strsm('Left', 'Lower', 'No transpose', 'Unit', KAHEAD, KCOLS, ONE, A( KSTART, KSTART ), LDA, A( KSTART, J+1 ), LDA );
         // Schur complement.
         sgemm('No transpose', 'No transpose', M-J, KCOLS, KAHEAD, NEGONE, A( J+1, KSTART ), LDA, A( KSTART, J+1 ), LDA, ONE, A( J+1, J+1 ), LDA );
      }

      // Handle pivot permutations on the way out of the recursion
      NPIVED = IAND( NSTEP, -NSTEP );
      J = NSTEP - NPIVED;
      DO WHILE ( J > 0 );
         NTOPIV = IAND( J, -J );
         slaswp(NTOPIV, A( 1, J-NTOPIV+1 ), LDA, J+1, NSTEP, IPIV, 1 );
         J = J - NTOPIV;
      }

      // If short and wide, handle the rest of the columns.
      if ( M < N ) {
         slaswp(N-M, A( 1, M+KCOLS+1 ), LDA, 1, M, IPIV, 1 );
         strsm('Left', 'Lower', 'No transpose', 'Unit', M, N-M, ONE, A, LDA, A( 1,M+KCOLS+1 ), LDA );
      }

      return;

      // End of SGETRF

      }
