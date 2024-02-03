      SUBROUTINE CPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            AB( LDAB, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CPBTRF, CPBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPBSV ', -INFO )
         RETURN
      END IF

      // Compute the Cholesky factorization A = U**H*U or A = L*L**H.

      CALL CPBTRF( UPLO, N, KD, AB, LDAB, INFO )
      IF( INFO.EQ.0 ) THEN

         // Solve the system A*X = B, overwriting B with X.

         CALL CPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

      END IF
      RETURN

      // End of CPBSV

      END
