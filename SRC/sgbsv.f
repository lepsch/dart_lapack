      SUBROUTINE SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               AB( LDAB, * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. External Subroutines ..
      // EXTERNAL SGBTRF, SGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( KL < 0 ) {
         INFO = -2;
      } else if ( KU < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDAB < 2*KL+KU+1 ) {
         INFO = -6;
      } else if ( LDB < MAX( N, 1 ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SGBSV ', -INFO );
         RETURN;
      }

      // Compute the LU factorization of the band matrix A.

      sgbtrf(N, N, KL, KU, AB, LDAB, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         sgbtrs('No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO );
      }
      RETURN;

      // End of SGBSV

      }
