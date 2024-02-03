      bool             FUNCTION SLCTSX( AR, AI, BETA );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               AI, AR, BETA
      // ..

*  =====================================================================

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
            SLCTSX = .FALSE.
         } else {
            SLCTSX = .TRUE.
         }
         if ( I.EQ.MPLUSN ) {
            FS = .FALSE.
            I = 0
         }
      } else {
         I = I + 1
         if ( I.LE.N ) {
            SLCTSX = .TRUE.
         } else {
            SLCTSX = .FALSE.
         }
         if ( I.EQ.MPLUSN ) {
            FS = .TRUE.
            I = 0
         }
      }

        // IF( AR/BETA.GT.0.0 )THEN
           // SLCTSX = .TRUE.
        // ELSE
           // SLCTSX = .FALSE.
        // END IF

      RETURN

      // End of SLCTSX

      }
