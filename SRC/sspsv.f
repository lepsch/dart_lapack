      SUBROUTINE SSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               AP( * ), B( LDB, * )
      // ..
*
*  =====================================================================
*
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSPTRF, SSPTRS, XERBLA
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
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSPSV ', -INFO )
         RETURN
      END IF
*
      // Compute the factorization A = U*D*U**T or A = L*D*L**T.
*
      CALL SSPTRF( UPLO, N, AP, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
         // Solve the system A*X = B, overwriting B with X.
*
         CALL SSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
*
      END IF
      RETURN
*
      // End of SSPSV
*
      END
