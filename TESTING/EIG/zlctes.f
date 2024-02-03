      bool             FUNCTION ZLCTES( Z, D );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      COMPLEX*16         D, Z
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
*
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      double             ZMAX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, SIGN
      // ..
      // .. Executable Statements ..
*
      IF( D.EQ.CZERO ) THEN
         ZLCTES = ( DBLE( Z ).LT.ZERO )
      ELSE
         IF( DBLE( Z ).EQ.ZERO .OR. DBLE( D ).EQ.ZERO ) THEN
            ZLCTES = ( SIGN( ONE, DIMAG( Z ) ).NE. SIGN( ONE, DIMAG( D ) ) )
         ELSE IF( DIMAG( Z ).EQ.ZERO .OR. DIMAG( D ).EQ.ZERO ) THEN
            ZLCTES = ( SIGN( ONE, DBLE( Z ) ).NE. SIGN( ONE, DBLE( D ) ) )
         ELSE
            ZMAX = MAX( ABS( DBLE( Z ) ), ABS( DIMAG( Z ) ) )
            ZLCTES = ( ( DBLE( Z ) / ZMAX )*DBLE( D )+ ( DIMAG( Z ) / ZMAX )*DIMAG( D ).LT.ZERO )
         END IF
      END IF
*
      RETURN
*
      // End of ZLCTES
*
      END
