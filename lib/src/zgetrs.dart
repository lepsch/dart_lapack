      void zgetrs(TRANS, N, NRHS, final Matrix<double> A, final int LDA, IPIV, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                INFO, LDA, LDB, N, NRHS;
      int                IPIV( * );
      Complex         A( LDA, * ), B( LDB, * );
      // ..

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               NOTRAN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLASWP, ZTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZGETRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( NOTRAN ) {

         // Solve A * X = B.

         // Apply row interchanges to the right hand sides.

         zlaswp(NRHS, B, LDB, 1, N, IPIV, 1 );

         // Solve L*X = B, overwriting B with X.

         ztrsm('Left', 'Lower', 'No transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve U*X = B, overwriting B with X.

         ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      } else {

         // Solve A**T * X = B  or A**H * X = B.

         // Solve U**T *X = B or U**H *X = B, overwriting B with X.

         ztrsm('Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve L**T *X = B, or L**H *X = B overwriting B with X.

         ztrsm('Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Apply row interchanges to the solution vectors.

         zlaswp(NRHS, B, LDB, 1, N, IPIV, -1 );
      }

      }
