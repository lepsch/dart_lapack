      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      int                INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      int                I, J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
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
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )                CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
*
         DO 20 I = 1, NRHS
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
*
      ELSE
*
*        Solve A**T*X = B.
*
         DO 30 I = 1, NRHS
*
*           Solve U**T*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 )
   30    CONTINUE
*
*        Solve L**T*X = B, overwriting B with X.
*
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGBTRS
*
      END
