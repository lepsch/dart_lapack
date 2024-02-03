      SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             AB( LDAB, * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. External Subroutines ..
      // EXTERNAL DGBTRF, DGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( KL.LT.0 ) {
         INFO = -2
      } else if ( KU.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.2*KL+KU+1 ) {
         INFO = -6
      } else if ( LDB.LT.MAX( N, 1 ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('DGBSV ', -INFO );
         RETURN
      }

      // Compute the LU factorization of the band matrix A.

      dgbtrf(N, N, KL, KU, AB, LDAB, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         dgbtrs('No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO );
      }
      RETURN

      // End of DGBSV

      }
