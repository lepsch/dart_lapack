      SUBROUTINE SPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
      // ..
*
*  =====================================================================
*
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPOTRF, SPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
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
         CALL XERBLA( 'SPOSV ', -INFO )
         RETURN
      END IF
*
      // Compute the Cholesky factorization A = U**T*U or A = L*L**T.
*
      CALL SPOTRF( UPLO, N, A, LDA, INFO )
      IF( INFO.EQ.0 ) THEN
*
         // Solve the system A*X = B, overwriting B with X.
*
         CALL SPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
      END IF
      RETURN
*
      // End of SPOSV
*
      END
