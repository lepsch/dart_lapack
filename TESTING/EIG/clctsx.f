      bool             FUNCTION CLCTSX( ALPHA, BETA );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      // ..

*  =====================================================================

      // .. Parameters ..
      // REAL               ZERO
      // PARAMETER          ( ZERO = 0.0E+0 )
      // COMPLEX            CZERO
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

      if ( FS ) {
         I = I + 1
         if ( I.LE.M ) {
            CLCTSX = .FALSE.
         } else {
            CLCTSX = .TRUE.
         }
         if ( I.EQ.MPLUSN ) {
            FS = .FALSE.
            I = 0
         }
      } else {
         I = I + 1
         if ( I.LE.N ) {
            CLCTSX = .TRUE.
         } else {
            CLCTSX = .FALSE.
         }
         if ( I.EQ.MPLUSN ) {
            FS = .TRUE.
            I = 0
         }
      }

       // IF( BETA.EQ.CZERO ) THEN
          // CLCTSX = ( REAL( ALPHA ).GT.ZERO )
       // ELSE
          // CLCTSX = ( REAL( ALPHA/BETA ).GT.ZERO )
       // END IF

      RETURN

      // End of CLCTSX

      }
