      SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )

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
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               SFMIN
      int                I, J, JP;
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      int                ISAMAX;
      // EXTERNAL SLAMCH, ISAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGER, SSCAL, SSWAP, XERBLA
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
      if ( INFO.NE.0 ) {
         xerbla('SGETF2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (M.EQ.0 .OR. N.EQ.0) RETURN;

      // Compute machine safe minimum

      SFMIN = SLAMCH('S')

      DO 10 J = 1, MIN( M, N )

         // Find pivot and test for singularity.

         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         if ( A( JP, J ).NE.ZERO ) {

            // Apply the interchange to columns 1:N.

            if (JP.NE.J) CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA );

            // Compute elements J+1:M of J-th column.

            if ( J.LT.M ) {
               if ( ABS(A( J, J )) .GE. SFMIN ) {
                  sscal(M-J, ONE / A( J, J ), A( J+1, J ), 1 );
               } else {
                 for (I = 1; I <= M-J; I++) { // 20
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
                 } // 20
               }
            }

         } else if ( INFO.EQ.0 ) {

            INFO = J
         }

         if ( J.LT.MIN( M, N ) ) {

            // Update trailing submatrix.

            sger(M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA );
         }
      } // 10
      RETURN

      // End of SGETF2

      }
