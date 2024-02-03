      SUBROUTINE DORT03( RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK, RESULT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      List<String>       RC;
      int                INFO, K, LDU, LDV, LWORK, MU, MV, N;
      double             RESULT;
*     ..
*     .. Array Arguments ..
      double             U( LDU, * ), V( LDV, * ), WORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      int                I, IRC, J, LMX;
      double             RES1, RES2, S, ULP;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORT01, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Check inputs
*
      INFO = 0
      IF( LSAME( RC, 'R' ) ) THEN
         IRC = 0
      ELSE IF( LSAME( RC, 'C' ) ) THEN
         IRC = 1
      ELSE
         IRC = -1
      END IF
      IF( IRC.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( MU.LT.0 ) THEN
         INFO = -2
      ELSE IF( MV.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.MAX( MU, MV ) ) THEN
         INFO = -5
      ELSE IF( ( IRC.EQ.0 .AND. LDU.LT.MAX( 1, MU ) ) .OR. ( IRC.EQ.1 .AND. LDU.LT.MAX( 1, N ) ) ) THEN
         INFO = -7
      ELSE IF( ( IRC.EQ.0 .AND. LDV.LT.MAX( 1, MV ) ) .OR. ( IRC.EQ.1 .AND. LDV.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORT03', -INFO )
         RETURN
      END IF
*
*     Initialize result
*
      RESULT = ZERO
      IF( MU.EQ.0 .OR. MV.EQ.0 .OR. N.EQ.0 ) RETURN
*
*     Machine constants
*
      ULP = DLAMCH( 'Precision' )
*
      IF( IRC.EQ.0 ) THEN
*
*        Compare rows
*
         RES1 = ZERO
         DO 20 I = 1, K
            LMX = IDAMAX( N, U( I, 1 ), LDU )
            S = SIGN( ONE, U( I, LMX ) )*SIGN( ONE, V( I, LMX ) )
            DO 10 J = 1, N
               RES1 = MAX( RES1, ABS( U( I, J )-S*V( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
         RES1 = RES1 / ( DBLE( N )*ULP )
*
*        Compute orthogonality of rows of V.
*
         CALL DORT01( 'Rows', MV, N, V, LDV, WORK, LWORK, RES2 )
*
      ELSE
*
*        Compare columns
*
         RES1 = ZERO
         DO 40 I = 1, K
            LMX = IDAMAX( N, U( 1, I ), 1 )
            S = SIGN( ONE, U( LMX, I ) )*SIGN( ONE, V( LMX, I ) )
            DO 30 J = 1, N
               RES1 = MAX( RES1, ABS( U( J, I )-S*V( J, I ) ) )
   30       CONTINUE
   40    CONTINUE
         RES1 = RES1 / ( DBLE( N )*ULP )
*
*        Compute orthogonality of columns of V.
*
         CALL DORT01( 'Columns', N, MV, V, LDV, WORK, LWORK, RES2 )
      END IF
*
      RESULT = MIN( MAX( RES1, RES2 ), ONE / ULP )
      RETURN
*
*     End of DORT03
*
      END
