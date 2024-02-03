      SUBROUTINE ZPOTRF ( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
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
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOTRF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN

         // Use unblocked code.

         CALL ZPOTRF2( UPLO, N, A, LDA, INFO )
      } else {

         // Use blocked code.

         IF( UPPER ) THEN

            // Compute the Cholesky factorization A = U'*U.

            DO 10 J = 1, N, NB

               JB = MIN( NB, N-J+1 )

               // Compute the current block.

               CALL ZTRSM( 'Left', 'Upper', 'Conjugate Transpose', 'Non-unit', J-1, JB, CONE, A( 1, 1 ), LDA, A( 1, J ), LDA )
                CALL ZHERK( 'Upper', 'Conjugate Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA )

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               CALL ZPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30

   10       CONTINUE

         } else {

            // Compute the Cholesky factorization A = L*L'.

            DO 20 J = 1, N, NB

               JB = MIN( NB, N-J+1 )

               // Compute the current block.

               CALL ZTRSM( 'Right', 'Lower', 'Conjugate Transpose', 'Non-unit', JB, J-1, CONE, A( 1, 1 ), LDA, A( J, 1 ), LDA )                 CALL ZHERK( 'Lower', 'No Transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA )

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               CALL ZPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30

   20       CONTINUE
         END IF
      END IF
      GO TO 40

   30 CONTINUE
      INFO = INFO + J - 1

   40 CONTINUE
      RETURN

      // End of ZPOTRF

      }
