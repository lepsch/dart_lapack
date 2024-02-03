      SUBROUTINE SLAUUM( UPLO, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
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
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IB, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLAUU2, SSYRK, STRMM, XERBLA
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
         CALL XERBLA( 'SLAUUM', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'SLAUUM', UPLO, N, -1, -1, -1 )

      IF( NB.LE.1 .OR. NB.GE.N ) THEN

         // Use unblocked code

         CALL SLAUU2( UPLO, N, A, LDA, INFO )
      ELSE

         // Use blocked code

         IF( UPPER ) THEN

            // Compute the product U * U**T.

            DO 10 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Non-unit', I-1, IB, ONE, A( I, I ), LDA, A( 1, I ), LDA )
               CALL SLAUU2( 'Upper', IB, A( I, I ), LDA, INFO )
               IF( I+IB.LE.N ) THEN
                  CALL SGEMM( 'No transpose', 'Transpose', I-1, IB, N-I-IB+1, ONE, A( 1, I+IB ), LDA, A( I, I+IB ), LDA, ONE, A( 1, I ), LDA )                   CALL SSYRK( 'Upper', 'No transpose', IB, N-I-IB+1, ONE, A( I, I+IB ), LDA, ONE, A( I, I ), LDA )
               END IF
   10       CONTINUE
         ELSE

            // Compute the product L**T * L.

            DO 20 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               CALL STRMM( 'Left', 'Lower', 'Transpose', 'Non-unit', IB, I-1, ONE, A( I, I ), LDA, A( I, 1 ), LDA )
               CALL SLAUU2( 'Lower', IB, A( I, I ), LDA, INFO )
               IF( I+IB.LE.N ) THEN
                  CALL SGEMM( 'Transpose', 'No transpose', IB, I-1, N-I-IB+1, ONE, A( I+IB, I ), LDA, A( I+IB, 1 ), LDA, ONE, A( I, 1 ), LDA )
                  CALL SSYRK( 'Lower', 'Transpose', IB, N-I-IB+1, ONE, A( I+IB, I ), LDA, ONE, A( I, I ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF

      RETURN

      // End of SLAUUM

      END
