      SUBROUTINE CGETRF ( M, N, A, LDA, IPIV, INFO)

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
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, IINFO, J, JB, NB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
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
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN

         // Use unblocked code.

         CALL CGETF2( M, N, A, LDA, IPIV, INFO )
      } else {

         // Use blocked code.

         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

            // Update current block.

            CALL CGEMM( 'No transpose', 'No transpose', M-J+1, JB, J-1, -ONE, A( J, 1 ), LDA, A( 1, J ), LDA, ONE, A( J, J ), LDA )


            // Factor diagonal and subdiagonal blocks and test for exact
            // singularity.

            CALL CGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )

            // Adjust INFO and the pivot indices.

            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE

            // Apply interchanges to column 1:J-1

            CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )

            IF ( J+JB.LE.N ) THEN

               // Apply interchanges to column J+JB:N

               CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, IPIV, 1 )

               CALL CGEMM( 'No transpose', 'No transpose', JB, N-J-JB+1, J-1, -ONE, A( J, 1 ), LDA, A( 1, J+JB ), LDA, ONE, A( J, J+JB ), LDA )

               // Compute block row of U.

               CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )
            END IF

   20    CONTINUE

      END IF
      RETURN

      // End of CGETRF

      }
