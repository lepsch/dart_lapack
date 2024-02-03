      SUBROUTINE DGET53( A, LDA, B, LDB, SCALE, WR, WI, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB;
      double             RESULT, SCALE, WI, WR;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      double             ABSW, ANORM, BNORM, CI11, CI12, CI22, CNORM, CR11, CR12, CR21, CR22, CSCALE, DETI, DETR, S1, SAFMIN, SCALES, SIGMIN, TEMP, ULP, WIS, WRS;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Initialize

      INFO = 0
      RESULT = ZERO
      SCALES = SCALE
      WRS = WR
      WIS = WI

      // Machine constants and norms

      SAFMIN = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ABSW = ABS( WRS ) + ABS( WIS )
      ANORM = MAX( ABS( A( 1, 1 ) )+ABS( A( 2, 1 ) ), ABS( A( 1, 2 ) )+ABS( A( 2, 2 ) ), SAFMIN )       BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 1, 2 ) )+ABS( B( 2, 2 ) ), SAFMIN )

      // Check for possible overflow.

      TEMP = ( SAFMIN*BNORM )*ABSW + ( SAFMIN*ANORM )*SCALES
      IF( TEMP.GE.ONE ) THEN

         // Scale down to avoid overflow

         INFO = 1
         TEMP = ONE / TEMP
         SCALES = SCALES*TEMP
         WRS = WRS*TEMP
         WIS = WIS*TEMP
         ABSW = ABS( WRS ) + ABS( WIS )
      END IF
      S1 = MAX( ULP*MAX( SCALES*ANORM, ABSW*BNORM ), SAFMIN*MAX( SCALES, ABSW ) )

      // Check for W and SCALE essentially zero.

      IF( S1.LT.SAFMIN ) THEN
         INFO = 2
         IF( SCALES.LT.SAFMIN .AND. ABSW.LT.SAFMIN ) THEN
            INFO = 3
            RESULT = ONE / ULP
            RETURN
         END IF

         // Scale up to avoid underflow

         TEMP = ONE / MAX( SCALES*ANORM+ABSW*BNORM, SAFMIN )
         SCALES = SCALES*TEMP
         WRS = WRS*TEMP
         WIS = WIS*TEMP
         ABSW = ABS( WRS ) + ABS( WIS )
         S1 = MAX( ULP*MAX( SCALES*ANORM, ABSW*BNORM ), SAFMIN*MAX( SCALES, ABSW ) )
         IF( S1.LT.SAFMIN ) THEN
            INFO = 3
            RESULT = ONE / ULP
            RETURN
         END IF
      END IF

      // Compute C = s A - w B

      CR11 = SCALES*A( 1, 1 ) - WRS*B( 1, 1 )
      CI11 = -WIS*B( 1, 1 )
      CR21 = SCALES*A( 2, 1 )
      CR12 = SCALES*A( 1, 2 ) - WRS*B( 1, 2 )
      CI12 = -WIS*B( 1, 2 )
      CR22 = SCALES*A( 2, 2 ) - WRS*B( 2, 2 )
      CI22 = -WIS*B( 2, 2 )

      // Compute the smallest singular value of s A - w B:

                  // |det( s A - w B )|
      // sigma_min = ------------------
                  // norm( s A - w B )

      CNORM = MAX( ABS( CR11 )+ABS( CI11 )+ABS( CR21 ), ABS( CR12 )+ABS( CI12 )+ABS( CR22 )+ABS( CI22 ), SAFMIN )
      CSCALE = ONE / SQRT( CNORM )
      DETR = ( CSCALE*CR11 )*( CSCALE*CR22 ) - ( CSCALE*CI11 )*( CSCALE*CI22 ) - ( CSCALE*CR12 )*( CSCALE*CR21 )       DETI = ( CSCALE*CR11 )*( CSCALE*CI22 ) + ( CSCALE*CI11 )*( CSCALE*CR22 ) - ( CSCALE*CI12 )*( CSCALE*CR21 )
      SIGMIN = ABS( DETR ) + ABS( DETI )
      RESULT = SIGMIN / S1
      RETURN

      // End of DGET53

      END
