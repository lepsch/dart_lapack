      int ilauplo(UPLO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      // ..

// =====================================================================

      // .. Parameters ..
      int     BLAS_UPPER, BLAS_LOWER;
      const     BLAS_UPPER = 121, BLAS_LOWER = 122 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..
      if ( LSAME( UPLO, 'U' ) ) {
         ILAUPLO = BLAS_UPPER;
      } else if ( LSAME( UPLO, 'L' ) ) {
         ILAUPLO = BLAS_LOWER;
      } else {
         ILAUPLO = -1;
      }
      return;

      // End of ILAUPLO

      }
