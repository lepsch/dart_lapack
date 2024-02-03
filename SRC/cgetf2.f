      SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) )
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

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETF2', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      DO 10 J = 1, MIN( M, N )

         // Find pivot and test for singularity.

         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN

            // Apply the interchange to columns 1:N.

            IF( JP.NE.J ) CALL CSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )

            // Compute elements J+1:M of J-th column.

            IF( J.LT.M ) CALL CRSCL( M-J, A( J, J ), A( J+1, J ), 1 )

         ELSE IF( INFO.EQ.0 ) THEN

            INFO = J
         END IF

         IF( J.LT.MIN( M, N ) ) THEN

            // Update trailing submatrix.

            CALL CGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN

      // End of CGETF2

      }
