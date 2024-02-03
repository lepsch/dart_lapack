      bool             FUNCTION ZLCTSX( ALPHA, BETA );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      // ..

*  =====================================================================

      // .. Parameters ..
      // DOUBLE PRECISION               ZERO
      // PARAMETER          ( ZERO = 0.0E+0 )
      // COMPLEX*16            CZERO
      // PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      // ..
      // .. Scalars in Common ..
      bool               FS;
      int                I, M, MPLUSN, N;
      // ..
      // .. Common blocks ..
      COMMON             / MN / M, N, MPLUSN, I, FS
      // ..
      // .. Save statement ..
      SAVE
      // ..
      // .. Executable Statements ..

      IF( FS ) THEN
         I = I + 1
         IF( I.LE.M ) THEN
            ZLCTSX = .FALSE.
         ELSE
            ZLCTSX = .TRUE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .FALSE.
            I = 0
         END IF
      ELSE
         I = I + 1
         IF( I.LE.N ) THEN
            ZLCTSX = .TRUE.
         ELSE
            ZLCTSX = .FALSE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .TRUE.
            I = 0
         END IF
      END IF

       // IF( BETA.EQ.CZERO ) THEN
          // ZLCTSX = ( DBLE( ALPHA ).GT.ZERO )
       // ELSE
          // ZLCTSX = ( DBLE( ALPHA/BETA ).GT.ZERO )
       // END IF

      RETURN

      // End of ZLCTSX

      }
