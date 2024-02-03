      SUBROUTINE STBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL STBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) && .NOT. LSAME( TRANS, 'T' ) && .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT && .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( KD.LT.0 ) {
         INFO = -5
      } else if ( NRHS.LT.0 ) {
         INFO = -6
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -10
      }
      if ( INFO != 0 ) {
         xerbla('STBTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Check for singularity.

      if ( NOUNIT ) {
         if ( UPPER ) {
            for (INFO = 1; INFO <= N; INFO++) { // 10
               IF( AB( KD+1, INFO ) == ZERO ) RETURN
            } // 10
         } else {
            for (INFO = 1; INFO <= N; INFO++) { // 20
               IF( AB( 1, INFO ) == ZERO ) RETURN
            } // 20
         }
      }
      INFO = 0

      // Solve A * X = B  or  A**T * X = B.

      for (J = 1; J <= NRHS; J++) { // 30
         stbsv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, B( 1, J ), 1 );
      } // 30

      RETURN

      // End of STBTRS

      }
