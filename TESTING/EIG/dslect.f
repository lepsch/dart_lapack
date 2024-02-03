      bool             FUNCTION DSLECT( ZR, ZI );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      double             ZI, ZR;
*     ..
*
*  =====================================================================
*
*     .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
*     ..
*     .. Scalars in Common ..
      int                SELDIM, SELOPT
*     ..
*     .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
*     ..
*     .. Local Scalars ..
      int                I
      double             RMIN, X;
*     ..
*     .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. External Functions ..
      double             DLAPY2;
      EXTERNAL           DLAPY2
*     ..
*     .. Executable Statements ..
*
      IF( SELOPT.EQ.0 ) THEN
         DSLECT = ( ZR.LT.ZERO )
      ELSE
         RMIN = DLAPY2( ZR-SELWR( 1 ), ZI-SELWI( 1 ) )
         DSLECT = SELVAL( 1 )
         DO 10 I = 2, SELDIM
            X = DLAPY2( ZR-SELWR( I ), ZI-SELWI( I ) )
            IF( X.LE.RMIN ) THEN
               RMIN = X
               DSLECT = SELVAL( I )
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DSLECT
*
      END
