      SUBROUTINE DLATB5( PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ANORM, CNDNUM;
      int                IMAT, KL, KU, MODE, N;
      String             DIST, TYPE;
      String             PATH;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             SHRINK, TENTH;
      const              SHRINK = 0.25D0, TENTH = 0.1D+0 ;
      double             ONE;
      const              ONE = 1.0D+0 ;
      double             TWO;
      const              TWO = 2.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             BADC1, BADC2, EPS, LARGE, SMALL;
      bool               FIRST;
      String             C2;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Save statement ..
      SAVE               EPS, SMALL, LARGE, BADC1, BADC2, FIRST
      // ..
      // .. Data statements ..
      DATA               FIRST / true /
      // ..
      // .. Executable Statements ..

      // Set some constants for use in the subroutine.

      if ( FIRST ) {
         FIRST = false;
         EPS = DLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = DLAMCH( 'Safe minimum' )
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

      if ( IMAT == 1 ) {
         KL = 0
      } else {
         KL = MAX( N-1, 0 )
      }
      KU = KL

      // Set the condition number and norm.etc

      if ( IMAT == 3 ) {
         CNDNUM = 1.0D12
         MODE = 2
      } else if ( IMAT == 4 ) {
         CNDNUM = 1.0D12
         MODE = 1
      } else if ( IMAT == 5 ) {
         CNDNUM = 1.0D12
         MODE = 3
      } else if ( IMAT == 6 ) {
         CNDNUM = BADC1
      } else if ( IMAT == 7 ) {
         CNDNUM = BADC2
      } else {
         CNDNUM = TWO
      }

      if ( IMAT == 8 ) {
         ANORM = SMALL
      } else if ( IMAT == 9 ) {
         ANORM = LARGE
      } else {
         ANORM = ONE
      }

      if (N <= 1) CNDNUM = ONE;

      RETURN

      // End of DLATB5

      }
