      SUBROUTINE DPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             A( 0: * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, NORMALTRANSR;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DTFSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      if ( !NORMALTRANSR && !LSAME( TRANSR, 'T' ) ) {
         INFO = -1
      } else if ( !LOWER && !LSAME( UPLO, 'U' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 0 ) {
         INFO = -4
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('DPFTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      // start execution: there are two triangular solves

      if ( LOWER ) {
         dtfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB );
         dtfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB );
      } else {
         dtfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB );
         dtfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB );
      }

      RETURN

      // End of DPFTRS

      }
