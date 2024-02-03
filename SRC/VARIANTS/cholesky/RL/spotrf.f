      SUBROUTINE SPOTRF ( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SPOTRF2, SSYRK, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('SPOTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'SPOTRF', UPLO, N, -1, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code.

         spotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            DO 10 J = 1, N, NB

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = MIN( NB, N-J+1 )

               spotrf2('Upper', JB, A( J, J ), LDA, INFO );
                IF( INFO.NE.0 ) GO TO 30

               if ( J+JB.LE.N ) {

                  // Updating the trailing submatrix.

                  strsm('Left', 'Upper', 'Transpose', 'Non-unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )                   CALL SSYRK( 'Upper', 'Transpose', N-J-JB+1, JB, -ONE, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
   10       CONTINUE

         } else {

            // Compute the Cholesky factorization A = L*L'.

            DO 20 J = 1, N, NB

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = MIN( NB, N-J+1 )

               spotrf2('Lower', JB, A( J, J ), LDA, INFO );
                IF( INFO.NE.0 ) GO TO 30

               if ( J+JB.LE.N ) {

                 // Updating the trailing submatrix.

                 strsm('Right', 'Lower', 'Transpose', 'Non-unit', N-J-JB+1, JB, ONE, A( J, J ), LDA, A( J+JB, J ), LDA )                   CALL SSYRK( 'Lower', 'No Transpose', N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
   20       CONTINUE
         }
      }
      GO TO 40

   30 CONTINUE
      INFO = INFO + J - 1

   40 CONTINUE
      RETURN

      // End of SPOTRF

      }
