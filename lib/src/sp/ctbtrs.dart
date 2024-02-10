      void ctbtrs(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      Complex            AB( LDAB, * ), B( LDB, * );
      // ..

      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      bool               NOUNIT, UPPER;
      int                J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NOUNIT = lsame( DIAG, 'N' );
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !lsame( TRANS, 'N' ) && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( KD < 0 ) {
         INFO = -5;
      } else if ( NRHS < 0 ) {
         INFO = -6;
      } else if ( LDAB < KD+1 ) {
         INFO = -8;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CTBTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check for singularity.

      if ( NOUNIT ) {
         if ( UPPER ) {
            for (INFO = 1; INFO <= N; INFO++) { // 10
               if( AB( KD+1, INFO ) == ZERO ) return;
            } // 10
         } else {
            for (INFO = 1; INFO <= N; INFO++) { // 20
               if( AB( 1, INFO ) == ZERO ) return;
            } // 20
         }
      }
      INFO = 0;

      // Solve A * X = B,  A**T * X = B,  or  A**H * X = B.

      for (J = 1; J <= NRHS; J++) { // 30
         ctbsv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, B( 1, J ), 1 );
      } // 30

      }
