      void ztptrs(UPLO, TRANS, DIAG, N, NRHS, AP, final Matrix<double> B, final int LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
      Complex         AP( * ), B( LDB, * );
      // ..

      Complex         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      bool               NOUNIT, UPPER;
      int                J, JC;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOUNIT = lsame( DIAG, 'N' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !lsame( TRANS, 'N' ) && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZTPTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check for singularity.

      if ( NOUNIT ) {
         if ( UPPER ) {
            JC = 1;
            for (INFO = 1; INFO <= N; INFO++) { // 10
               if( AP( JC+INFO-1 ) == ZERO ) return;
               JC = JC + INFO;
            } // 10
         } else {
            JC = 1;
            for (INFO = 1; INFO <= N; INFO++) { // 20
               if( AP( JC ) == ZERO ) return;
               JC = JC + N - INFO + 1;
            } // 20
         }
      }
      INFO = 0;

      // Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.

      for (J = 1; J <= NRHS; J++) { // 30
         ztpsv(UPLO, TRANS, DIAG, N, AP, B( 1, J ), 1 );
      } // 30

      }
