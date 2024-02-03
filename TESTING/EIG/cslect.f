      bool             FUNCTION CSLECT( Z );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            Z
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               RMIN, X
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      REAL               SELWI( 20 ), SELWR( 20 )
      // ..
      // .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, REAL
      // ..
      // .. Executable Statements ..

      IF( SELOPT.EQ.0 ) THEN
         CSLECT = ( REAL( Z ).LT.ZERO )
      ELSE
         RMIN = ABS( Z-CMPLX( SELWR( 1 ), SELWI( 1 ) ) )
         CSLECT = SELVAL( 1 )
         DO 10 I = 2, SELDIM
            X = ABS( Z-CMPLX( SELWR( I ), SELWI( I ) ) )
            IF( X.LE.RMIN ) THEN
               RMIN = X
               CSLECT = SELVAL( I )
            END IF
   10    CONTINUE
      END IF
      RETURN

      // End of CSLECT

      }
