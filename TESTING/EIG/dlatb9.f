      SUBROUTINE DLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                IMAT, KLA, KLB, KUA, KUB, M, MODEA, MODEB, N, P;
      double             ANORM, BNORM, CNDNMA, CNDNMB;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             SHRINK, TENTH;
      const              SHRINK = 0.25D0, TENTH = 0.1D+0 ;
      double             ONE, TEN;
      const              ONE = 1.0D+0, TEN = 1.0D+1 ;
      // ..
      // .. Local Scalars ..
      bool               FIRST;
      double             BADC1, BADC2, EPS, LARGE, SMALL;
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      double             DLAMCH;
      // EXTERNAL LSAMEN, DLAMCH
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

      // Set some parameters we don't plan to change.

      TYPE = 'N'
      DISTA = 'S'
      DISTB = 'S'
      MODEA = 3
      MODEB = 4

      // Set the lower and upper bandwidths.

      if ( LSAMEN( 3, PATH, 'GRQ' ) .OR. LSAMEN( 3, PATH, 'LSE' ) .OR. LSAMEN( 3, PATH, 'GSV' ) ) {

         // A: M by N, B: P by N

         if ( IMAT.EQ.1 ) {

            // A: diagonal, B: upper triangular

            KLA = 0
            KUA = 0
            KLB = 0
            KUB = MAX( N-1, 0 )

         } else if ( IMAT.EQ.2 ) {

            // A: upper triangular, B: upper triangular

            KLA = 0
            KUA = MAX( N-1, 0 )
            KLB = 0
            KUB = MAX( N-1, 0 )

         } else if ( IMAT.EQ.3 ) {

            // A: lower triangular, B: upper triangular

            KLA = MAX( M-1, 0 )
            KUA = 0
            KLB = 0
            KUB = MAX( N-1, 0 )

         } else {

            // A: general dense, B: general dense

            KLA = MAX( M-1, 0 )
            KUA = MAX( N-1, 0 )
            KLB = MAX( P-1, 0 )
            KUB = MAX( N-1, 0 )

         }

      } else if ( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GLM' ) ) {

         // A: N by M, B: N by P

         if ( IMAT.EQ.1 ) {

            // A: diagonal, B: lower triangular

            KLA = 0
            KUA = 0
            KLB = MAX( N-1, 0 )
            KUB = 0
         } else if ( IMAT.EQ.2 ) {

            // A: lower triangular, B: diagonal

            KLA = MAX( N-1, 0 )
            KUA = 0
            KLB = 0
            KUB = 0

         } else if ( IMAT.EQ.3 ) {

            // A: lower triangular, B: upper triangular

            KLA = MAX( N-1, 0 )
            KUA = 0
            KLB = 0
            KUB = MAX( P-1, 0 )

         } else {

            // A: general dense, B: general dense

            KLA = MAX( N-1, 0 )
            KUA = MAX( M-1, 0 )
            KLB = MAX( N-1, 0 )
            KUB = MAX( P-1, 0 )
         }

      }

      // Set the condition number and norm.

      CNDNMA = TEN*TEN
      CNDNMB = TEN
      if ( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' ) .OR. LSAMEN( 3, PATH, 'GSV' ) ) {
         if ( IMAT.EQ.5 ) {
            CNDNMA = BADC1
            CNDNMB = BADC1
         } else if ( IMAT.EQ.6 ) {
            CNDNMA = BADC2
            CNDNMB = BADC2
         } else if ( IMAT.EQ.7 ) {
            CNDNMA = BADC1
            CNDNMB = BADC2
         } else if ( IMAT.EQ.8 ) {
            CNDNMA = BADC2
            CNDNMB = BADC1
         }
      }

      ANORM = TEN
      BNORM = TEN*TEN*TEN
      if ( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' ) ) {
         if ( IMAT.EQ.7 ) {
            ANORM = SMALL
            BNORM = LARGE
         } else if ( IMAT.EQ.8 ) {
            ANORM = LARGE
            BNORM = SMALL
         }
      }

      if ( N.LE.1 ) {
         CNDNMA = ONE
         CNDNMB = ONE
      }

      RETURN

      // End of DLATB9

      }
