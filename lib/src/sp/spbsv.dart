      void spbsv(final int UPLO, final int N, final int KD, final int NRHS, final Matrix<double> AB_, final int LDAB, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final AB = AB_.dim();
  final B = B_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      double               AB( LDAB, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPBTRF, SPBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
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
         xerbla('SPBSV ', -INFO );
         return;
      }

      // Compute the Cholesky factorization A = U**T*U or A = L*L**T.

      spbtrf(UPLO, N, KD, AB, LDAB, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         spbtrs(UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO );

      }
      }
