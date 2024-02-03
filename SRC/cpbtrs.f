      SUBROUTINE CPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            AB( LDAB, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('CPBTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 .OR. NRHS == 0) RETURN;

      if ( UPPER ) {

         // Solve A*X = B where A = U**H *U.

         for (J = 1; J <= NRHS; J++) { // 10

            // Solve U**H *X = B, overwriting B with X.

            ctbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve U*X = B, overwriting B with X.

            ctbsv('Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
         } // 10
      } else {

         // Solve A*X = B where A = L*L**H.

         for (J = 1; J <= NRHS; J++) { // 20

            // Solve L*X = B, overwriting B with X.

            ctbsv('Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve L**H *X = B, overwriting B with X.

            ctbsv('Lower', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
         } // 20
      }

      RETURN

      // End of CPBTRS

      }
