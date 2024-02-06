      void zpotrs(UPLO, N, NRHS, A, LDA, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      Complex         A( LDA, * ), B( LDB, * );
      // ..

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

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
         xerbla('ZPOTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( UPPER ) {

         // Solve A*X = B where A = U**H *U.

         // Solve U**H *X = B, overwriting B with X.

         ztrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve U*X = B, overwriting B with X.

         ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      } else {

         // Solve A*X = B where A = L*L**H.

         // Solve L*X = B, overwriting B with X.

         ztrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve L**H *X = B, overwriting B with X.

         ztrsm('Left', 'Lower', 'Conjugate transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      }

      return;
      }
