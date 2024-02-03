      SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. External Subroutines ..
      // EXTERNAL SGETRF, SGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( NRHS.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SGESV ', -INFO )
         RETURN
      }

      // Compute the LU factorization of A.

      CALL SGETRF( N, N, A, LDA, IPIV, INFO )
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         CALL SGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      }
      RETURN

      // End of SGESV

      }
