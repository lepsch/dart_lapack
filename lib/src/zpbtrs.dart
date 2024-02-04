      void zpbtrs(UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      Complex         AB( LDAB, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTBSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
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
         xerbla('ZPBTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( UPPER ) {

         // Solve A*X = B where A = U**H *U.

         for (J = 1; J <= NRHS; J++) { // 10

            // Solve U**H *X = B, overwriting B with X.

            ztbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve U*X = B, overwriting B with X.

            ztbsv('Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
         } // 10
      } else {

         // Solve A*X = B where A = L*L**H.

         for (J = 1; J <= NRHS; J++) { // 20

            // Solve L*X = B, overwriting B with X.

            ztbsv('Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve L**H *X = B, overwriting B with X.

            ztbsv('Lower', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
         } // 20
      }

      return;
      }