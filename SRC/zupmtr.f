      SUBROUTINE ZUPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             SIDE, TRANS, UPLO;
      int                INFO, LDC, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * ), C( LDC, * ), TAU( * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      bool               FORWRD, LEFT, NOTRAN, UPPER;
      int                I, I1, I2, I3, IC, II, JC, MI, NI, NQ;
      COMPLEX*16         AII, TAUI
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARF
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      UPPER = LSAME( UPLO, 'U' )
*
      // NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUPMTR', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
         // Q was determined by a call to ZHPTRD with UPLO = 'U'
*
         FORWRD = ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN )
*
         IF( FORWRD ) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         END IF
*
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
*
         DO 10 I = I1, I2, I3
            IF( LEFT ) THEN
*
               // H(i) or H(i)**H is applied to C(1:i,1:n)
*
               MI = I
            ELSE
*
               // H(i) or H(i)**H is applied to C(1:m,1:i)
*
               NI = I
            END IF
*
            // Apply H(i) or H(i)**H
*
            IF( NOTRAN ) THEN
               TAUI = TAU( I )
            ELSE
               TAUI = DCONJG( TAU( I ) )
            END IF
            AII = AP( II )
            AP( II ) = ONE
            CALL ZLARF( SIDE, MI, NI, AP( II-I+1 ), 1, TAUI, C, LDC, WORK )
            AP( II ) = AII
*
            IF( FORWRD ) THEN
               II = II + I + 2
            ELSE
               II = II - I - 1
            END IF
   10    CONTINUE
      ELSE
*
         // Q was determined by a call to ZHPTRD with UPLO = 'L'.
*
         FORWRD = ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN )
*
         IF( FORWRD ) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 20 I = I1, I2, I3
            AII = AP( II )
            AP( II ) = ONE
            IF( LEFT ) THEN
*
               // H(i) or H(i)**H is applied to C(i+1:m,1:n)
*
               MI = M - I
               IC = I + 1
            ELSE
*
               // H(i) or H(i)**H is applied to C(1:m,i+1:n)
*
               NI = N - I
               JC = I + 1
            END IF
*
            // Apply H(i) or H(i)**H
*
            IF( NOTRAN ) THEN
               TAUI = TAU( I )
            ELSE
               TAUI = DCONJG( TAU( I ) )
            END IF
            CALL ZLARF( SIDE, MI, NI, AP( II ), 1, TAUI, C( IC, JC ), LDC, WORK )
            AP( II ) = AII
*
            IF( FORWRD ) THEN
               II = II + NQ - I + 1
            ELSE
               II = II - NQ + I - 2
            END IF
   20    CONTINUE
      END IF
      RETURN
*
      // End of ZUPMTR
*
      END
