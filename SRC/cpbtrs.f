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
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPBTRS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      IF( UPPER ) THEN

         // Solve A*X = B where A = U**H *U.

         DO 10 J = 1, NRHS

            // Solve U**H *X = B, overwriting B with X.

            CALL CTBSV( 'Upper', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 )

            // Solve U*X = B, overwriting B with X.

            CALL CTBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 )
   10    CONTINUE
      ELSE

         // Solve A*X = B where A = L*L**H.

         DO 20 J = 1, NRHS

            // Solve L*X = B, overwriting B with X.

            CALL CTBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 )

            // Solve L**H *X = B, overwriting B with X.

            CALL CTBSV( 'Lower', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 )
   20    CONTINUE
      END IF

      RETURN

      // End of CPBTRS

      END
