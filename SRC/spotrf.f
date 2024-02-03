      SUBROUTINE SPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
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
*
      // Test the input parameters.
*
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
         CALL XERBLA( 'SPOTRF', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // Determine the block size for this environment.
*
      NB = ILAENV( 1, 'SPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
         // Use unblocked code.
*
         CALL SPOTRF2( UPLO, N, A, LDA, INFO )
      ELSE
*
         // Use blocked code.
*
         IF( UPPER ) THEN
*
            // Compute the Cholesky factorization A = U**T*U.
*
            DO 10 J = 1, N, NB
*
               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL SSYRK( 'Upper', 'Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL SPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30
               IF( J+JB.LE.N ) THEN
*
                  // Compute the current block row.
*
                  CALL SGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1, J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ), LDA, ONE, A( J, J+JB ), LDA )                   CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', JB, N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
*
         ELSE
*
            // Compute the Cholesky factorization A = L*L**T.
*
            DO 20 J = 1, N, NB
*
               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL SSYRK( 'Lower', 'No transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL SPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30
               IF( J+JB.LE.N ) THEN
*
                  // Compute the current block column.
*
                  CALL SGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+JB, J ), LDA )                   CALL STRSM( 'Right', 'Lower', 'Transpose', 'Non-unit', N-J-JB+1, JB, ONE, A( J, J ), LDA, A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
      // End of SPOTRF
*
      END
