      void cpftrs(TRANSR, UPLO, N, NRHS, A, B, LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANSR, UPLO;
      int                INFO, LDB, N, NRHS;
      Complex            A( 0: * ), B( LDB, * );
      // ..

      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               LOWER, NORMALTRANSR;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CTFSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LOWER = lsame( UPLO, 'L' );
      if ( !NORMALTRANSR && !lsame( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CPFTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // start execution: there are two triangular solves

      if ( LOWER ) {
         ctfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B, LDB );
         ctfsm(TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B, LDB );
      } else {
         ctfsm(TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B, LDB );
         ctfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B, LDB );
      }

      }
