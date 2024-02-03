      SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO );

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
      COMPLEX            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, JP;
      // ..
      // .. External Functions ..
      int                ICAMAX;
      // EXTERNAL ICAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGERU, CRSCL, CSWAP, XERBLA
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
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CGETF2', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) RETURN;

      DO 10 J = 1, MIN( M, N );

         // Find pivot and test for singularity.

         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 );
         IPIV( J ) = JP;
         if ( A( JP, J ) != ZERO ) {

            // Apply the interchange to columns 1:N.

            if (JP != J) CALL CSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA );

            // Compute elements J+1:M of J-th column.

            if (J < M) CALL CRSCL( M-J, A( J, J ), A( J+1, J ), 1 );

         } else if ( INFO == 0 ) {

            INFO = J;
         }

         if ( J < MIN( M, N ) ) {

            // Update trailing submatrix.

            cgeru(M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA );
         }
      } // 10
      return;

      // End of CGETF2

      }
