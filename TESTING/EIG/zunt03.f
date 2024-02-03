      SUBROUTINE ZUNT03( RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK, RWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       RC;
      int                INFO, K, LDU, LDV, LWORK, MU, MV, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      int                I, IRC, J, LMX;
      double             RES1, RES2, ULP;
      COMPLEX*16         S, SU, SV
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZUNT01
      // ..
      // .. Executable Statements ..

      // Check inputs

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
         CALL XERBLA( 'ZUNT03', -INFO )
         RETURN
      END IF

      // Initialize result

      RESULT = ZERO
      IF( MU.EQ.0 .OR. MV.EQ.0 .OR. N.EQ.0 ) RETURN

      // Machine constants

      ULP = DLAMCH( 'Precision' )

      IF( IRC.EQ.0 ) THEN

         // Compare rows

         RES1 = ZERO
         DO 20 I = 1, K
            LMX = IZAMAX( N, U( I, 1 ), LDU )
            IF( V( I, LMX ).EQ.DCMPLX( ZERO ) ) THEN
               SV = ONE
            ELSE
               SV = ABS( V( I, LMX ) ) / V( I, LMX )
            END IF
            IF( U( I, LMX ).EQ.DCMPLX( ZERO ) ) THEN
               SU = ONE
            ELSE
               SU = ABS( U( I, LMX ) ) / U( I, LMX )
            END IF
            S = SV / SU
            DO 10 J = 1, N
               RES1 = MAX( RES1, ABS( U( I, J )-S*V( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
         RES1 = RES1 / ( DBLE( N )*ULP )

         // Compute orthogonality of rows of V.

         CALL ZUNT01( 'Rows', MV, N, V, LDV, WORK, LWORK, RWORK, RES2 )

      ELSE

         // Compare columns

         RES1 = ZERO
         DO 40 I = 1, K
            LMX = IZAMAX( N, U( 1, I ), 1 )
            IF( V( LMX, I ).EQ.DCMPLX( ZERO ) ) THEN
               SV = ONE
            ELSE
               SV = ABS( V( LMX, I ) ) / V( LMX, I )
            END IF
            IF( U( LMX, I ).EQ.DCMPLX( ZERO ) ) THEN
               SU = ONE
            ELSE
               SU = ABS( U( LMX, I ) ) / U( LMX, I )
            END IF
            S = SV / SU
            DO 30 J = 1, N
               RES1 = MAX( RES1, ABS( U( J, I )-S*V( J, I ) ) )
   30       CONTINUE
   40    CONTINUE
         RES1 = RES1 / ( DBLE( N )*ULP )

         // Compute orthogonality of columns of V.

         CALL ZUNT01( 'Columns', N, MV, V, LDV, WORK, LWORK, RWORK, RES2 )
      END IF

      RESULT = MIN( MAX( RES1, RES2 ), ONE / ULP )
      RETURN

      // End of ZUNT03

      END
