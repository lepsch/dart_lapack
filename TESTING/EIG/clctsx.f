      bool             FUNCTION CLCTSX( ALPHA, BETA );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*     REAL               ZERO
*     PARAMETER          ( ZERO = 0.0E+0 )
*     COMPLEX            CZERO
*     PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Scalars in Common ..
      bool               FS;
      int                I, M, MPLUSN, N;
*     ..
*     .. Common blocks ..
      COMMON             / MN / M, N, MPLUSN, I, FS
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
      IF( FS ) THEN
         I = I + 1
         IF( I.LE.M ) THEN
            CLCTSX = .FALSE.
         ELSE
            CLCTSX = .TRUE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .FALSE.
            I = 0
         END IF
      ELSE
         I = I + 1
         IF( I.LE.N ) THEN
            CLCTSX = .TRUE.
         ELSE
            CLCTSX = .FALSE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .TRUE.
            I = 0
         END IF
      END IF
*
*      IF( BETA.EQ.CZERO ) THEN
*         CLCTSX = ( REAL( ALPHA ).GT.ZERO )
*      ELSE
*         CLCTSX = ( REAL( ALPHA/BETA ).GT.ZERO )
*      END IF
*
      RETURN
*
*     End of CLCTSX
*
      END
