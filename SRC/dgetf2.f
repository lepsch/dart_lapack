      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             SFMIN;
      int                I, J, JP;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      int                IDAMAX;
      // EXTERNAL DLAMCH, IDAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGER, DSCAL, DSWAP, XERBLA
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
         xerbla('DGETF2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Compute machine safe minimum

      SFMIN = DLAMCH('S')

      DO 10 J = 1, MIN( M, N )

         // Find pivot and test for singularity.

         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         if ( A( JP, J ).NE.ZERO ) {

            // Apply the interchange to columns 1:N.

            IF( JP.NE.J ) CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )

            // Compute elements J+1:M of J-th column.

            if ( J.LT.M ) {
               if ( ABS(A( J, J )) .GE. SFMIN ) {
                  dscal(M-J, ONE / A( J, J ), A( J+1, J ), 1 );
               } else {
                 DO 20 I = 1, M-J
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
                 } // 20
               }
            }

         } else if ( INFO.EQ.0 ) {

            INFO = J
         }

         if ( J.LT.MIN( M, N ) ) {

            // Update trailing submatrix.

            dger(M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA );
         }
      } // 10
      RETURN

      // End of DGETF2

      }
