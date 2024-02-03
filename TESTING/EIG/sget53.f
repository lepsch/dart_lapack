      void sget53(A, LDA, B, LDB, SCALE, WR, WI, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB;
      REAL               RESULT, SCALE, WI, WR;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      REAL               ABSW, ANORM, BNORM, CI11, CI12, CI22, CNORM, CR11, CR12, CR21, CR22, CSCALE, DETI, DETR, S1, SAFMIN, SCALES, SIGMIN, TEMP, ULP, WIS, WRS;
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Initialize

      INFO = 0;
      RESULT = ZERO;
      SCALES = SCALE;
      WRS = WR;
      WIS = WI;

      // Machine constants and norms

      SAFMIN = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
      ABSW = ABS( WRS ) + ABS( WIS );
      ANORM = max( ABS( A( 1, 1 ) )+ABS( A( 2, 1 ) ), ABS( A( 1, 2 ) )+ABS( A( 2, 2 ) ), SAFMIN )       BNORM = max( ABS( B( 1, 1 ) ), ABS( B( 1, 2 ) )+ABS( B( 2, 2 ) ), SAFMIN );

      // Check for possible overflow.

      TEMP = ( SAFMIN*BNORM )*ABSW + ( SAFMIN*ANORM )*SCALES;
      if ( TEMP >= ONE ) {

         // Scale down to avoid overflow

         INFO = 1;
         TEMP = ONE / TEMP;
         SCALES = SCALES*TEMP;
         WRS = WRS*TEMP;
         WIS = WIS*TEMP;
         ABSW = ABS( WRS ) + ABS( WIS );
      }
      S1 = max( ULP*max( SCALES*ANORM, ABSW*BNORM ), SAFMIN*max( SCALES, ABSW ) );

      // Check for W and SCALE essentially zero.

      if ( S1 < SAFMIN ) {
         INFO = 2;
         if ( SCALES < SAFMIN && ABSW < SAFMIN ) {
            INFO = 3;
            RESULT = ONE / ULP;
            return;
         }

         // Scale up to avoid underflow

         TEMP = ONE / max( SCALES*ANORM+ABSW*BNORM, SAFMIN );
         SCALES = SCALES*TEMP;
         WRS = WRS*TEMP;
         WIS = WIS*TEMP;
         ABSW = ABS( WRS ) + ABS( WIS );
         S1 = max( ULP*max( SCALES*ANORM, ABSW*BNORM ), SAFMIN*max( SCALES, ABSW ) );
         if ( S1 < SAFMIN ) {
            INFO = 3;
            RESULT = ONE / ULP;
            return;
         }
      }

      // Compute C = s A - w B

      CR11 = SCALES*A( 1, 1 ) - WRS*B( 1, 1 );
      CI11 = -WIS*B( 1, 1 );
      CR21 = SCALES*A( 2, 1 );
      CR12 = SCALES*A( 1, 2 ) - WRS*B( 1, 2 );
      CI12 = -WIS*B( 1, 2 );
      CR22 = SCALES*A( 2, 2 ) - WRS*B( 2, 2 );
      CI22 = -WIS*B( 2, 2 );

      // Compute the smallest singular value of s A - w B:

                  // |det( s A - w B )|
      // sigma_min = ------------------
                  // norm( s A - w B )

      CNORM = max( ABS( CR11 )+ABS( CI11 )+ABS( CR21 ), ABS( CR12 )+ABS( CI12 )+ABS( CR22 )+ABS( CI22 ), SAFMIN );
      CSCALE = ONE / sqrt( CNORM );
      DETR = ( CSCALE*CR11 )*( CSCALE*CR22 ) - ( CSCALE*CI11 )*( CSCALE*CI22 ) - ( CSCALE*CR12 )*( CSCALE*CR21 )       DETI = ( CSCALE*CR11 )*( CSCALE*CI22 ) + ( CSCALE*CI11 )*( CSCALE*CR22 ) - ( CSCALE*CI12 )*( CSCALE*CR21 );
      SIGMIN = ABS( DETR ) + ABS( DETI );
      RESULT = SIGMIN / S1;
      return;
      }
