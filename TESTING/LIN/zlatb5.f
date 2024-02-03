      SUBROUTINE ZLATB5( PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

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
      PARAMETER          ( SHRINK = 0.25D0, TENTH = 0.1D+0 )
      double             ONE;
      PARAMETER          ( ONE = 1.0D+0 )
      double             TWO;
      PARAMETER          ( TWO = 2.0D+0 )
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
      DATA               FIRST / .TRUE. /
      // ..
      // .. Executable Statements ..

      // Set some constants for use in the subroutine.

      IF( FIRST ) THEN
         FIRST = .FALSE.
         EPS = DLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = DLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      END IF

      C2 = PATH( 2: 3 )

      // Set some parameters

      DIST = 'S'
      MODE = 3

      // Set TYPE, the type of matrix to be generated.

      TYPE = C2( 1: 1 )

      // Set the lower and upper bandwidths.

      IF( IMAT.EQ.1 ) THEN
         KL = 0
      ELSE
         KL = MAX( N-1, 0 )
      END IF
      KU = KL

      // Set the condition number and norm.etc

      IF( IMAT.EQ.3 ) THEN
         CNDNUM = 1.0D12
         MODE = 2
      ELSE IF( IMAT.EQ.4 ) THEN
         CNDNUM = 1.0D12
         MODE = 1
      ELSE IF( IMAT.EQ.5 ) THEN
         CNDNUM = 1.0D12
         MODE = 3
      ELSE IF( IMAT.EQ.6 ) THEN
         CNDNUM = BADC1
      ELSE IF( IMAT.EQ.7 ) THEN
         CNDNUM = BADC2
      ELSE
         CNDNUM = TWO
      END IF

      IF( IMAT.EQ.8 ) THEN
         ANORM = SMALL
      ELSE IF( IMAT.EQ.9 ) THEN
         ANORM = LARGE
      ELSE
         ANORM = ONE
      END IF

      IF( N.LE.1 ) CNDNUM = ONE

      RETURN

      // End of ZLATB5

      }
