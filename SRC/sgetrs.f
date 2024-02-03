      SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASWP, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      IF( NOTRAN ) THEN

         // Solve A * X = B.

         // Apply row interchanges to the right hand sides.

         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )

         // Solve L*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB )

         // Solve U*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
      ELSE

         // Solve A**T * X = B.

         // Solve U**T *X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )

         // Solve L**T *X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB )

         // Apply row interchanges to the solution vectors.

         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF

      RETURN

      // End of SGETRS

      }
