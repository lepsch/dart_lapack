      SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      IF( UPPER ) THEN

         // Solve A*X = B where A = U**T *U.

         // Solve U**T *X = B, overwriting B with X.

         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )

         // Solve U*X = B, overwriting B with X.

         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
      } else {

         // Solve A*X = B where A = L*L**T.

         // Solve L*X = B, overwriting B with X.

         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )

         // Solve L**T *X = B, overwriting B with X.

         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
      END IF

      RETURN

      // End of DPOTRS

      }
