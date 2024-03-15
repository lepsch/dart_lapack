      void cgetrf(final int M, final int N, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      int                IPIV( * );
      Complex            A( LDA, * );
      // ..

      Complex            ONE;
      const              ONE = (1.0, 0.0) ;
      int                I, IINFO, J, JB, K, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
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
         xerbla('CGETRF', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'CGETRF', ' ', M, N, -1, -1 );
      if ( NB <= 1 || NB >= min( M, N ) ) {

         // Use unblocked code.

         cgetf2(M, N, A, LDA, IPIV, INFO );

      } else {

         // Use blocked code.

         for (J = 1; NB < 0 ? J >= min( M, N ) : J <= min( M, N ); J += NB) { // 20
            JB = min( min( M, N )-J+1, NB );


            // Update before factoring the current panel

            for (K = 1; NB < 0 ? K >= J-NB : K <= J-NB; K += NB) { // 30

               // Apply interchanges to rows K:K+NB-1.

               claswp(JB, A(1, J), LDA, K, K+NB-1, IPIV, 1 );

               // Compute block row of U.

               ctrsm('Left', 'Lower', 'No transpose', 'Unit', NB, JB, ONE, A( K, K ), LDA, A( K, J ), LDA );

               // Update trailing submatrix.

               cgemm('No transpose', 'No transpose', M-K-NB+1, JB, NB, -ONE, A( K+NB, K ), LDA, A( K, J ), LDA, ONE, A( K+NB, J ), LDA );
            } // 30

            // Factor diagonal and subdiagonal blocks and test for exact
            // singularity.

            cgetf2(M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO );

            // Adjust INFO and the pivot indices.

            if (INFO == 0 && IINFO > 0) INFO = IINFO + J - 1;
            for (I = J; I <= min( M, J+JB-1 ); I++) { // 10
               IPIV[I] = J - 1 + IPIV( I );
            } // 10

         } // 20


         // Apply interchanges to the left-overs

         for (K = 1; NB < 0 ? K >= min( M, N ) : K <= min( M, N ); K += NB) { // 40
            claswp(K-1, A( 1, 1 ), LDA, K, min(K+NB-1, min( M, N )), IPIV, 1 );
         } // 40

         // Apply update to the M+1:N columns when N > M

         if ( N > M ) {

            claswp(N-M, A(1, M+1), LDA, 1, M, IPIV, 1 );

            for (K = 1; NB < 0 ? K >= M : K <= M; K += NB) { // 50

               JB = min( M-K+1, NB );

               ctrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-M, ONE, A( K, K ), LDA, A( K, M+1 ), LDA );


               if ( K+NB <= M ) {
                    cgemm('No transpose', 'No transpose', M-K-NB+1, N-M, NB, -ONE, A( K+NB, K ), LDA, A( K, M+1 ), LDA, ONE, A( K+NB, M+1 ), LDA );
               }
            } // 50
         }

      }
      }