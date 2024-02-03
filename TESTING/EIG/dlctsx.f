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
      COMMON             / MN / M, N, MPLUSN, I, FS
      // ..
      // .. Save statement ..
      SAVE
      // ..
      // .. Executable Statements ..

      if ( FS ) {
         I = I + 1
         if ( I.LE.M ) {
            DLCTSX = .FALSE.
         } else {
            DLCTSX = .TRUE.
         }
         if ( I.EQ.MPLUSN ) {
            FS = .FALSE.
            I = 0
         }
      } else {
         I = I + 1
         if ( I.LE.N ) {
            DLCTSX = .TRUE.
         } else {
            DLCTSX = .FALSE.
         }
         if ( I.EQ.MPLUSN ) {
            FS = .TRUE.
            I = 0
         }
      }

        // IF( AR/BETA.GT.0.0 )THEN
           // DLCTSX = .TRUE.
        // ELSE
           // DLCTSX = .FALSE.
        // END IF

      RETURN

      // End of DLCTSX

      }
