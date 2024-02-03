      void zpbsv(UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AB( LDAB, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZPBTRF, ZPBTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDAB < KD+1 ) {
         INFO = -6;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZPBSV ', -INFO );
         return;
      }

      // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

      zpbtrf(UPLO, N, KD, AB, LDAB, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         zpbtrs(UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO );

      }
      return;

      // End of ZPBSV

      }
