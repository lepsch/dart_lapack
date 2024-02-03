      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine --
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
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGETRF2, SLASWP, STRSM, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
      if ( INFO != 0 ) {
         xerbla('SGETRF', -INFO );
         RETURN
      }

      // Quick return if possible

      if (M == 0 .OR. N == 0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) {

         // Use unblocked code.

         sgetrf2(M, N, A, LDA, IPIV, INFO );
      } else {

         // Use blocked code.

         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

            // Factor diagonal and subdiagonal blocks and test for exact
            // singularity.

            sgetrf2(M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO );

            // Adjust INFO and the pivot indices.

            if (INFO == 0 && IINFO.GT.0) INFO = IINFO + J - 1;
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
            } // 10

            // Apply interchanges to columns 1:J-1.

            slaswp(J-1, A, LDA, J, J+JB-1, IPIV, 1 );

            if ( J+JB.LE.N ) {

               // Apply interchanges to columns J+JB:N.

               slaswp(N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, IPIV, 1 );

               // Compute block row of U.

               strsm('Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA );
               if ( J+JB.LE.M ) {

                  // Update trailing submatrix.

                  sgemm('No transpose', 'No transpose', M-J-JB+1, N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
            }
         } // 20
      }
      RETURN

      // End of SGETRF

      }
