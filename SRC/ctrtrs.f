      SUBROUTINE CTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE
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
      // EXTERNAL CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) && .NOT. LSAME( TRANS, 'T' ) && .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT && .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( NRHS < 0 ) {
         INFO = -5
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('CTRTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Check for singularity.

      if ( NOUNIT ) {
         for (INFO = 1; INFO <= N; INFO++) { // 10
            IF( A( INFO, INFO ) == ZERO ) RETURN
         } // 10
      }
      INFO = 0

      // Solve A * x = b,  A**T * x = b,  or  A**H * x = b.

      ctrsm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      RETURN

      // End of CTRTRS

      }
