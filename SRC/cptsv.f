      SUBROUTINE CPTSV( N, NRHS, D, E, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               D( * )
      COMPLEX            B( LDB, * ), E( * )
      // ..

*  =====================================================================

      // .. External Subroutines ..
      // EXTERNAL CPTTRF, CPTTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( N < 0 ) {
         INFO = -1
      } else if ( NRHS < 0 ) {
         INFO = -2
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('CPTSV ', -INFO );
         RETURN
      }

      // Compute the L*D*L**H (or U**H*D*U) factorization of A.

      cpttrf(N, D, E, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         cpttrs('Lower', N, NRHS, D, E, B, LDB, INFO );
      }
      RETURN

      // End of CPTSV

      }
