      bool             FUNCTION SLCTSX( AR, AI, BETA );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               AI, AR, BETA;
      // ..

*  =====================================================================

      // .. Scalars in Common ..
      bool               FS;
      int                I, M, MPLUSN, N;
      // ..
      // .. Common blocks ..
      // COMMON / MN / M, N, MPLUSN, I, FS
      // ..
      // .. Save statement ..
      SAVE;
      // ..
      // .. Executable Statements ..

      if ( FS ) {
         I = I + 1;
         if ( I <= M ) {
            SLCTSX = false;
         } else {
            SLCTSX = true;
         }
         if ( I == MPLUSN ) {
            FS = false;
            I = 0;
         }
      } else {
         I = I + 1;
         if ( I <= N ) {
            SLCTSX = true;
         } else {
            SLCTSX = false;
         }
         if ( I == MPLUSN ) {
            FS = true;
            I = 0;
         }
      }

        // IF( AR/BETA > 0.0 )THEN
           // SLCTSX = true;
        // ELSE
           // SLCTSX = false;
        // END IF

      return;

      // End of SLCTSX

      }
