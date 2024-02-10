      void sgetrf(M, N, final Matrix<double> A, final int LDA, final Array<int> IPIV, INFO) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      int                IPIV( * );
      double               A( LDA, * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      int                I, IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGETF2, SLASWP, STRSM, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

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
         xerbla('SGETRF', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'SGETRF', ' ', M, N, -1, -1 );
      if ( NB <= 1 || NB >= min( M, N ) ) {

         // Use unblocked code.

         sgetf2(M, N, A, LDA, IPIV, INFO );
      } else {

         // Use blocked code.

         for (J = 1; NB < 0 ? J >= min( M, N ) : J <= min( M, N ); J += NB) { // 20
            JB = min( min( M, N )-J+1, NB );

            // Update current block.

            sgemm('No transpose', 'No transpose', M-J+1, JB, J-1, -ONE, A( J, 1 ), LDA, A( 1, J ), LDA, ONE, A( J, J ), LDA );


            // Factor diagonal and subdiagonal blocks and test for exact
            // singularity.

            sgetf2(M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO );

            // Adjust INFO and the pivot indices.

            if (INFO == 0 && IINFO > 0) INFO = IINFO + J - 1;
            for (I = J; I <= min( M, J+JB-1 ); I++) { // 10
               IPIV[I] = J - 1 + IPIV( I );
            } // 10

            // Apply interchanges to column 1:J-1

            slaswp(J-1, A, LDA, J, J+JB-1, IPIV, 1 );

            if ( J+JB <= N ) {

               // Apply interchanges to column J+JB:N

               slaswp(N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, IPIV, 1 );

               sgemm('No transpose', 'No transpose', JB, N-J-JB+1, J-1, -ONE, A( J, 1 ), LDA, A( 1, J+JB ), LDA, ONE, A( J, J+JB ), LDA );

               // Compute block row of U.

               strsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA );
            }

         } // 20

      }
      }
