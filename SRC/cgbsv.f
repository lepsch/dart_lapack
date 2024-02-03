      SUBROUTINE CGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            AB( LDAB, * ), B( LDB, * )
      // ..
*
*  =====================================================================
*
      // .. External Subroutines ..
      // EXTERNAL CGBTRF, CGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( KL.LT.0 ) THEN
         INFO = -2
      ELSE IF( KU.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGBSV ', -INFO )
         RETURN
      END IF
*
      // Compute the LU factorization of the band matrix A.
*
      CALL CGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
         // Solve the system A*X = B, overwriting B with X.
*
         CALL CGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
      END IF
      RETURN
*
      // End of CGBSV
*
      END
