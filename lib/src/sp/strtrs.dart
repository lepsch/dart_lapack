      void strtrs(UPLO, TRANS, DIAG, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      double               A( LDA, * ), B( LDB, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               NOUNIT;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NOUNIT = lsame( DIAG, 'N' );
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !lsame( TRANS, 'N' ) && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('STRTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check for singularity.

      if ( NOUNIT ) {
         for (INFO = 1; INFO <= N; INFO++) { // 10
            if( A( INFO, INFO ) == ZERO ) return;
         } // 10
      }
      INFO = 0;

      // Solve A * x = b  or  A**T * x = b.

      strsm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      }
