      bool             FUNCTION ZSLECT( Z );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16         Z
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             RMIN, X;
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      IF( SELOPT.EQ.0 ) THEN
         ZSLECT = ( DBLE( Z ).LT.ZERO )
      ELSE
         RMIN = ABS( Z-DCMPLX( SELWR( 1 ), SELWI( 1 ) ) )
         ZSLECT = SELVAL( 1 )
         DO 10 I = 2, SELDIM
            X = ABS( Z-DCMPLX( SELWR( I ), SELWI( I ) ) )
            IF( X.LE.RMIN ) THEN
               RMIN = X
               ZSLECT = SELVAL( I )
            END IF
   10    CONTINUE
      END IF
      RETURN

      // End of ZSLECT

      }
