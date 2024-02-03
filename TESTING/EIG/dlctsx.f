      bool             FUNCTION DLCTSX( AR, AI, BETA );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             AI, AR, BETA;
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
      SAVE
      // ..
      // .. Executable Statements ..

      if ( FS ) {
         I = I + 1
         if ( I <= M ) {
            DLCTSX = false;
         } else {
            DLCTSX = true;
         }
         if ( I == MPLUSN ) {
            FS = false;
            I = 0
         }
      } else {
         I = I + 1
         if ( I <= N ) {
            DLCTSX = true;
         } else {
            DLCTSX = false;
         }
         if ( I == MPLUSN ) {
            FS = true;
            I = 0
         }
      }

        // IF( AR/BETA > 0.0 )THEN
           // DLCTSX = true;
        // ELSE
           // DLCTSX = false;
        // END IF

      RETURN

      // End of DLCTSX

      }
