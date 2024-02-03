      SUBROUTINE SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               AB( LDAB, * ), B( LDB, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               LNOTI, NOTRAN;
      int                I, J, KD, L, LM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SGER, SSWAP, STBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRS', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
         // Solve  A*X = B.
*
         // Solve L*X = B, overwriting B with X.
*
         // L is represented as a product of permutations and unit lower
        t // riangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
         // where each transformation L(i) is a rank-one modification of
        t // he identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )                CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
*
         DO 20 I = 1, NRHS
*
            // Solve U*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
*
      ELSE
*
         // Solve A**T*X = B.
*
         DO 30 I = 1, NRHS
*
            // Solve U**T*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 )
   30    CONTINUE
*
         // Solve L**T*X = B, overwriting B with X.
*
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
*
      // End of SGBTRS
*
      END
