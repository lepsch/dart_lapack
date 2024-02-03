      void zposv(UPLO, N, NRHS, A, LDA, B, LDB, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZPOTRF, ZPOTRS
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
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('ZPOSV ', -INFO );
         return;
      }

      // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

      zpotrf(UPLO, N, A, LDA, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         zpotrs(UPLO, N, NRHS, A, LDA, B, LDB, INFO );

      }
      return;
      }
