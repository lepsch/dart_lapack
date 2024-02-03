      SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NOUNIT = LSAME( DIAG, 'N' );
      if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !LSAME( TRANS, 'N' ) && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('ZTRTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Check for singularity.

      if ( NOUNIT ) {
         for (INFO = 1; INFO <= N; INFO++) { // 10
            if( A( INFO, INFO ) == ZERO ) RETURN;
         } // 10
      }
      INFO = 0;

      // Solve A * x = b,  A**T * x = b,  or  A**H * x = b.

      ztrsm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      return;

      // End of ZTRTRS

      }
