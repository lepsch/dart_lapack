      void sposv(UPLO, N, NRHS, A, LDA, B, LDB, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      double               A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPOTRF, SPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
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
         xerbla('SPOSV ', -INFO );
         return;
      }

      // Compute the Cholesky factorization A = U**T*U or A = L*L**T.

      spotrf(UPLO, N, A, LDA, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         spotrs(UPLO, N, NRHS, A, LDA, B, LDB, INFO );

      }
      return;
      }
