      SUBROUTINE SPTSV( N, NRHS, D, E, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), E( * )
      // ..

*  =====================================================================

      // .. External Subroutines ..
      // EXTERNAL SPTTRF, SPTTRS, XERBLA
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
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('SPTSV ', -INFO );
         RETURN
      }

      // Compute the L*D*L**T (or U**T*D*U) factorization of A.

      spttrf(N, D, E, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         spttrs(N, NRHS, D, E, B, LDB, INFO );
      }
      RETURN

      // End of SPTSV

      }
