      SUBROUTINE ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JC;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      NOUNIT = LSAME( DIAG, 'N' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !LSAME( TRANS, 'N' ) && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZTPTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Check for singularity.

      if ( NOUNIT ) {
         if ( UPPER ) {
            JC = 1;
            for (INFO = 1; INFO <= N; INFO++) { // 10
               if( AP( JC+INFO-1 ) == ZERO ) RETURN;
               JC = JC + INFO;
            } // 10
         } else {
            JC = 1;
            for (INFO = 1; INFO <= N; INFO++) { // 20
               if( AP( JC ) == ZERO ) RETURN;
               JC = JC + N - INFO + 1;
            } // 20
         }
      }
      INFO = 0;

      // Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.

      for (J = 1; J <= NRHS; J++) { // 30
         ztpsv(UPLO, TRANS, DIAG, N, AP, B( 1, J ), 1 );
      } // 30

      return;

      // End of ZTPTRS

      }
