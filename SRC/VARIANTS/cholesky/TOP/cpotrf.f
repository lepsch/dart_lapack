      SUBROUTINE CPOTRF ( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      COMPLEX            CONE
      const              ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ) ;
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
      // EXTERNAL CGEMM, CHERK, CPOTRF2, CTRSM, XERBLA
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
         CALL XERBLA( 'CPOTRF', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'CPOTRF', UPLO, N, -1, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code.

         CALL CPOTRF2( UPLO, N, A, LDA, INFO )
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            DO 10 J = 1, N, NB

               JB = MIN( NB, N-J+1 )

               // Compute the current block.

               CALL CTRSM( 'Left', 'Upper', 'Conjugate Transpose', 'Non-unit', J-1, JB, CONE, A( 1, 1 ), LDA, A( 1, J ), LDA )
                CALL CHERK( 'Upper', 'Conjugate Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA )

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               CALL CPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30

   10       CONTINUE

         } else {

            // Compute the Cholesky factorization A = L*L'.

            DO 20 J = 1, N, NB

               JB = MIN( NB, N-J+1 )

               // Compute the current block.

               CALL CTRSM( 'Right', 'Lower', 'Conjugate Transpose', 'Non-unit', JB, J-1, CONE, A( 1, 1 ), LDA, A( J, 1 ), LDA )                 CALL CHERK( 'Lower', 'No Transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA )

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               CALL CPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30

   20       CONTINUE
         }
      }
      GO TO 40

   30 CONTINUE
      INFO = INFO + J - 1

   40 CONTINUE
      RETURN

      // End of CPOTRF

      }
