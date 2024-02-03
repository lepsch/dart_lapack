      SUBROUTINE SOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS, UPLO;
      int                INFO, LDC, M, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), C( LDC, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               FORWRD, LEFT, NOTRAN, UPPER;
      int                I, I1, I2, I3, IC, II, JC, MI, NI, NQ;
      REAL               AII
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      UPPER = LSAME( UPLO, 'U' )

      // NQ is the order of Q

      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SOPMTR', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      IF( UPPER ) THEN

         // Q was determined by a call to SSPTRD with UPLO = 'U'

         FORWRD = ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN )

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

         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF

         DO 10 I = I1, I2, I3
            IF( LEFT ) THEN

               // H(i) is applied to C(1:i,1:n)

               MI = I
            ELSE

               // H(i) is applied to C(1:m,1:i)

               NI = I
            END IF

            // Apply H(i)

            AII = AP( II )
            AP( II ) = ONE
            CALL SLARF( SIDE, MI, NI, AP( II-I+1 ), 1, TAU( I ), C, LDC, WORK )
            AP( II ) = AII

            IF( FORWRD ) THEN
               II = II + I + 2
            ELSE
               II = II - I - 1
            END IF
   10    CONTINUE
      ELSE

         // Q was determined by a call to SSPTRD with UPLO = 'L'.

         FORWRD = ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN )

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

         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF

         DO 20 I = I1, I2, I3
            AII = AP( II )
            AP( II ) = ONE
            IF( LEFT ) THEN

               // H(i) is applied to C(i+1:m,1:n)

               MI = M - I
               IC = I + 1
            ELSE

               // H(i) is applied to C(1:m,i+1:n)

               NI = N - I
               JC = I + 1
            END IF

            // Apply H(i)

            CALL SLARF( SIDE, MI, NI, AP( II ), 1, TAU( I ), C( IC, JC ), LDC, WORK )
            AP( II ) = AII

            IF( FORWRD ) THEN
               II = II + NQ - I + 1
            ELSE
               II = II - NQ + I - 2
            END IF
   20    CONTINUE
      END IF
      RETURN

      // End of SOPMTR

      END
