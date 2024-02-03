      SUBROUTINE SPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), B( LDB, * )
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
      // EXTERNAL STBSV, XERBLA
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
         xerbla('SPBTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      if ( UPPER ) {

         // Solve A*X = B where A = U**T *U.

         DO 10 J = 1, NRHS

            // Solve U**T *X = B, overwriting B with X.

            stbsv('Upper', 'Transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve U*X = B, overwriting B with X.

            stbsv('Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
   10    CONTINUE
      } else {

         // Solve A*X = B where A = L*L**T.

         DO 20 J = 1, NRHS

            // Solve L*X = B, overwriting B with X.

            stbsv('Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve L**T *X = B, overwriting B with X.

            stbsv('Lower', 'Transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
   20    CONTINUE
      }

      RETURN

      // End of SPBTRS

      }
