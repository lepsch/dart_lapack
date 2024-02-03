      double           FUNCTION DSXT1( IJOB, D1, N1, D2, N2, ABSTOL, ULP, UNFL );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IJOB, N1, N2;
      double             ABSTOL, ULP, UNFL;
      // ..
      // .. Array Arguments ..
      double             D1( * ), D2( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             TEMP1, TEMP2;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      TEMP1 = ZERO

      J = 1
      DO 20 I = 1, N1
   10    CONTINUE
         IF( D2( J ).LT.D1( I ) .AND. J.LT.N2 ) THEN
            J = J + 1
            GO TO 10
         END IF
         IF( J.EQ.1 ) THEN
            TEMP2 = ABS( D2( J )-D1( I ) )
            IF( IJOB.EQ.2 ) TEMP2 = TEMP2 / MAX( UNFL, ABSTOL+ULP*ABS( D1( I ) ) )
         ELSE
            TEMP2 = MIN( ABS( D2( J )-D1( I ) ), ABS( D1( I )-D2( J-1 ) ) )             IF( IJOB.EQ.2 ) TEMP2 = TEMP2 / MAX( UNFL, ABSTOL+ULP*ABS( D1( I ) ) )
         END IF
         TEMP1 = MAX( TEMP1, TEMP2 )
   20 CONTINUE

      DSXT1 = TEMP1
      RETURN

      // End of DSXT1

      }
