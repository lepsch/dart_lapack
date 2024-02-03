      SUBROUTINE ZPOTRF ( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      COMPLEX*16         CONE
      const              ONE = 1.0D+0, CONE = ( 1.0D+0, 0.0D+0 ) ;
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
      // EXTERNAL XERBLA, ZGEMM, ZHERK, ZPOTRF2, ZTRSM
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
         xerbla('ZPOTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code.

         zpotrf2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            DO 10 J = 1, N, NB

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = MIN( NB, N-J+1 )

               zpotrf2('Upper', JB, A( J, J ), LDA, INFO );
                IF( INFO.NE.0 ) GO TO 30

               if ( J+JB.LE.N ) {

                  // Updating the trailing submatrix.

                  ztrsm('Left', 'Upper', 'Conjugate Transpose', 'Non-unit', JB, N-J-JB+1, CONE, A( J, J ), LDA, A( J, J+JB ), LDA )                   CALL ZHERK( 'Upper', 'Conjugate transpose', N-J-JB+1, JB, -ONE, A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
   10       CONTINUE

         } else {

            // Compute the Cholesky factorization A = L*L'.

            DO 20 J = 1, N, NB

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               JB = MIN( NB, N-J+1 )

               zpotrf2('Lower', JB, A( J, J ), LDA, INFO );
                IF( INFO.NE.0 ) GO TO 30

               if ( J+JB.LE.N ) {

                 // Updating the trailing submatrix.

                 ztrsm('Right', 'Lower', 'Conjugate Transpose', 'Non-unit', N-J-JB+1, JB, CONE, A( J, J ), LDA, A( J+JB, J ), LDA )                   CALL ZHERK( 'Lower', 'No Transpose', N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, ONE, A( J+JB, J+JB ), LDA );
               }
   20       CONTINUE
         }
      }
      GO TO 40

   30 CONTINUE
      INFO = INFO + J - 1

   40 CONTINUE
      RETURN

      // End of ZPOTRF

      }
