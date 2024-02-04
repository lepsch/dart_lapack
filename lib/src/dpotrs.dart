      void dpotrs(UPLO, N, NRHS, A, LDA, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTRSM, XERBLA
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
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DPOTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( UPPER ) {

         // Solve A*X = B where A = U**T *U.

         // Solve U**T *X = B, overwriting B with X.

         dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve U*X = B, overwriting B with X.

         dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      } else {

         // Solve A*X = B where A = L*L**T.

         // Solve L*X = B, overwriting B with X.

         dtrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve L**T *X = B, overwriting B with X.

         dtrsm('Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      }

      return;
      }