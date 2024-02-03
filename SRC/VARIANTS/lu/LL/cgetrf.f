      void cgetrf(M, N, A, LDA, IPIV, INFO) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE;
      const              ONE = (1.0, 0.0) ;
      // ..
      // .. Local Scalars ..
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
      // ..
      // .. Executable Statements ..

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

      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 );
      if ( NB <= 1 || NB >= min( M, N ) ) {

         // Use unblocked code.

         cgetf2(M, N, A, LDA, IPIV, INFO );

      } else {

         // Use blocked code.

         DO 20 J = 1, min( M, N ), NB;
            JB = min( min( M, N )-J+1, NB );


            // Update before factoring the current panel

            DO 30 K = 1, J-NB, NB;

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
            DO 10 I = J, min( M, J+JB-1 );
               IPIV( I ) = J - 1 + IPIV( I );
            } // 10

         } // 20


         // Apply interchanges to the left-overs

         DO 40 K = 1, min( M, N ), NB;
            claswp(K-1, A( 1, 1 ), LDA, K, min(K+NB-1, min( M, N )), IPIV, 1 );
         } // 40

         // Apply update to the M+1:N columns when N > M

         if ( N > M ) {

            claswp(N-M, A(1, M+1), LDA, 1, M, IPIV, 1 );

            DO 50 K = 1, M, NB;

               JB = min( M-K+1, NB );

               ctrsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-M, ONE, A( K, K ), LDA, A( K, M+1 ), LDA );


               if ( K+NB <= M ) {
                    cgemm('No transpose', 'No transpose', M-K-NB+1, N-M, NB, -ONE, A( K+NB, K ), LDA, A( K, M+1 ), LDA, ONE, A( K+NB, M+1 ), LDA );
               }
            } // 50
         }

      }
      return;
      }
