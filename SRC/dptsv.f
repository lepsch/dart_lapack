      SUBROUTINE DPTSV( N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), E( * );
      // ..
*
*  =====================================================================
*
      // .. External Subroutines ..
      // EXTERNAL DPTTRF, DPTTRS, XERBLA
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
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPTSV ', -INFO )
         RETURN
      END IF
*
      // Compute the L*D*L**T (or U**T*D*U) factorization of A.
*
      CALL DPTTRF( N, D, E, INFO )
      IF( INFO.EQ.0 ) THEN
*
         // Solve the system A*X = B, overwriting B with X.
*
         CALL DPTTRS( N, NRHS, D, E, B, LDB, INFO )
      END IF
      RETURN
*
      // End of DPTSV
*
      END
