      void zgetrf(M, N, final Matrix<double> A, final int LDA, IPIV, final Box<int> INFO ) {

// -- LAPACK computational routine (version 3.X) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      int                IPIV( * );
      Complex         A( LDA, * );
      // ..

      Complex         ONE, NEGONE;
      double             ZERO;
      const              ONE = (1.0, 0.0) ;
      const              NEGONE = (-1.0, 0.0) ;
      const              ZERO = 0.0 ;
      double             SFMIN, PIVMAG;
      Complex         TMP;
      int                I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
      int                KSTART, IPIVSTART, JPIVSTART, KCOLS;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      //- int                IZAMAX;
      //- bool               DISNAN;
      // EXTERNAL DLAMCH, IZAMAX, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZTRSM, ZSCAL, XERBLA, ZLASWP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, IAND, ABS

      // Test the input parameters.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZGETRF', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Compute machine safe minimum

      SFMIN = dlamch( 'S' );

      NSTEP = min( M, N );
      for (J = 1; J <= NSTEP; J++) {
         KAHEAD = IAND( J, -J );
         KSTART = J + 1 - KAHEAD;
         KCOLS = min( KAHEAD, M-J );

         // Find pivot.

         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 );
         IPIV[J] = JP;

         // Permute just this column.
         if (JP != J) {
            TMP = A( J, J );
            A[J][J] = A( JP, J );
            A[JP][J] = TMP;
         }

         // Apply pending permutations to L
         NTOPIV = 1;
         IPIVSTART = J;
         JPIVSTART = J - NTOPIV;
         while (NTOPIV < KAHEAD) {
            zlaswp(NTOPIV, A( 1, JPIVSTART ), LDA, IPIVSTART, J, IPIV, 1 );
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
         }

         // Permute U block to match L
         zlaswp(KCOLS, A( 1,J+1 ), LDA, KSTART, J, IPIV, 1 );

         // Factor the current column
         PIVMAG = ( A( J, J ) ).abs();
         if ( PIVMAG != ZERO && !disnan( PIVMAG ) ) {
               if ( PIVMAG >= SFMIN ) {
                  zscal(M-J, ONE / A( J, J ), A( J+1, J ), 1 );
               } else {
                 for (I = 1; I <= M-J; I++) {
                    A[J+I][J] = A( J+I, J ) / A( J, J );
                 }
               }
         } else if ( PIVMAG == ZERO && INFO == 0 ) {
            INFO = J;
         }

         // Solve for U block.
         ztrsm('Left', 'Lower', 'No transpose', 'Unit', KAHEAD, KCOLS, ONE, A( KSTART, KSTART ), LDA, A( KSTART, J+1 ), LDA );
         // Schur complement.
         zgemm('No transpose', 'No transpose', M-J, KCOLS, KAHEAD, NEGONE, A( J+1, KSTART ), LDA, A( KSTART, J+1 ), LDA, ONE, A( J+1, J+1 ), LDA );
      }

      // Handle pivot permutations on the way out of the recursion
      NPIVED = IAND( NSTEP, -NSTEP );
      J = NSTEP - NPIVED;
      while (J > 0) {
         NTOPIV = IAND( J, -J );
         zlaswp(NTOPIV, A( 1, J-NTOPIV+1 ), LDA, J+1, NSTEP, IPIV, 1 );
         J = J - NTOPIV;
      }

      // If short and wide, handle the rest of the columns.
      if ( M < N ) {
         zlaswp(N-M, A( 1, M+KCOLS+1 ), LDA, 1, M, IPIV, 1 );
         ztrsm('Left', 'Lower', 'No transpose', 'Unit', M, N-M, ONE, A, LDA, A( 1,M+KCOLS+1 ), LDA );
      }

      }
