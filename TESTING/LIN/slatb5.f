      SUBROUTINE SLATB5( PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               ANORM, CNDNUM
      int                IMAT, KL, KU, MODE, N;
      String             DIST, TYPE;
      String             PATH;
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               SHRINK, TENTH
      const              SHRINK = 0.25E0, TENTH = 0.1E+0 ;
      REAL               ONE
      const              ONE = 1.0E+0 ;
      REAL               TWO
      const              TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               BADC1, BADC2, EPS, LARGE, SMALL
      bool               FIRST;
      String             C2;
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Save statement ..
      SAVE               EPS, SMALL, LARGE, BADC1, BADC2, FIRST
      // ..
      // .. Data statements ..
      DATA               FIRST / .TRUE. /
      // ..
      // .. Executable Statements ..

      // Set some constants for use in the subroutine.

      if ( FIRST ) {
         FIRST = .FALSE.
         EPS = SLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = SLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      }

      C2 = PATH( 2: 3 )

      // Set some parameters

      DIST = 'S'
      MODE = 3

      // Set TYPE, the type of matrix to be generated.

      TYPE = C2( 1: 1 )

      // Set the lower and upper bandwidths.

      if ( IMAT.EQ.1 ) {
         KL = 0
      } else {
         KL = MAX( N-1, 0 )
      }
      KU = KL

      // Set the condition number and norm.etc

      if ( IMAT.EQ.3 ) {
         CNDNUM = 1.0E4
         MODE = 2
      } else if ( IMAT.EQ.4 ) {
         CNDNUM = 1.0E4
         MODE = 1
      } else if ( IMAT.EQ.5 ) {
         CNDNUM = 1.0E4
         MODE = 3
      } else if ( IMAT.EQ.6 ) {
         CNDNUM = BADC1
      } else if ( IMAT.EQ.7 ) {
         CNDNUM = BADC2
      } else {
         CNDNUM = TWO
      }

      if ( IMAT.EQ.8 ) {
         ANORM = SMALL
      } else if ( IMAT.EQ.9 ) {
         ANORM = LARGE
      } else {
         ANORM = ONE
      }

      IF( N.LE.1 ) CNDNUM = ONE

      RETURN

      // End of SLATB5

      }
