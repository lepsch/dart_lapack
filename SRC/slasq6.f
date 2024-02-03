      void slasq6(I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2 ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I0, N0, PP;
      REAL               DMIN, DMIN1, DMIN2, DN, DNM1, DNM2;
      // ..
      // .. Array Arguments ..
      REAL               Z( * );
      // ..

// =====================================================================

      // .. Parameter ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J4, J4P2;
      REAL               D, EMIN, SAFMIN, TEMP;
      // ..
      // .. External Function ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      if( ( N0-I0-1 ) <= 0 ) return;

      SAFMIN = SLAMCH( 'Safe minimum' );
      J4 = 4*I0 + PP - 3;
      EMIN = Z( J4+4 );
      D = Z( J4 );
      DMIN = D;

      if ( PP == 0 ) {
         DO 10 J4 = 4*I0, 4*( N0-3 ), 4;
            Z( J4-2 ) = D + Z( J4-1 );
            if ( Z( J4-2 ) == ZERO ) {
               Z( J4 ) = ZERO;
               D = Z( J4+1 );
               DMIN = D;
               EMIN = ZERO;
            } else if ( SAFMIN*Z( J4+1 ) < Z( J4-2 ) && SAFMIN*Z( J4-2 ) < Z( J4+1 ) ) {
               TEMP = Z( J4+1 ) / Z( J4-2 );
               Z( J4 ) = Z( J4-1 )*TEMP;
               D = D*TEMP;
            } else {
               Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) );
               D = Z( J4+1 )*( D / Z( J4-2 ) );
            }
            DMIN = min( DMIN, D );
            EMIN = min( EMIN, Z( J4 ) );
         } // 10
      } else {
         DO 20 J4 = 4*I0, 4*( N0-3 ), 4;
            Z( J4-3 ) = D + Z( J4 );
            if ( Z( J4-3 ) == ZERO ) {
               Z( J4-1 ) = ZERO;
               D = Z( J4+2 );
               DMIN = D;
               EMIN = ZERO;
            } else if ( SAFMIN*Z( J4+2 ) < Z( J4-3 ) && SAFMIN*Z( J4-3 ) < Z( J4+2 ) ) {
               TEMP = Z( J4+2 ) / Z( J4-3 );
               Z( J4-1 ) = Z( J4 )*TEMP;
               D = D*TEMP;
            } else {
               Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) );
               D = Z( J4+2 )*( D / Z( J4-3 ) );
            }
            DMIN = min( DMIN, D );
            EMIN = min( EMIN, Z( J4-1 ) );
         } // 20
      }

      // Unroll last two steps.

      DNM2 = D;
      DMIN2 = DMIN;
      J4 = 4*( N0-2 ) - PP;
      J4P2 = J4 + 2*PP - 1;
      Z( J4-2 ) = DNM2 + Z( J4P2 );
      if ( Z( J4-2 ) == ZERO ) {
         Z( J4 ) = ZERO;
         DNM1 = Z( J4P2+2 );
         DMIN = DNM1;
         EMIN = ZERO;
      } else if ( SAFMIN*Z( J4P2+2 ) < Z( J4-2 ) && SAFMIN*Z( J4-2 ) < Z( J4P2+2 ) ) {
         TEMP = Z( J4P2+2 ) / Z( J4-2 );
         Z( J4 ) = Z( J4P2 )*TEMP;
         DNM1 = DNM2*TEMP;
      } else {
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) );
      }
      DMIN = min( DMIN, DNM1 );

      DMIN1 = DMIN;
      J4 = J4 + 4;
      J4P2 = J4 + 2*PP - 1;
      Z( J4-2 ) = DNM1 + Z( J4P2 );
      if ( Z( J4-2 ) == ZERO ) {
         Z( J4 ) = ZERO;
         DN = Z( J4P2+2 );
         DMIN = DN;
         EMIN = ZERO;
      } else if ( SAFMIN*Z( J4P2+2 ) < Z( J4-2 ) && SAFMIN*Z( J4-2 ) < Z( J4P2+2 ) ) {
         TEMP = Z( J4P2+2 ) / Z( J4-2 );
         Z( J4 ) = Z( J4P2 )*TEMP;
         DN = DNM1*TEMP;
      } else {
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) );
      }
      DMIN = min( DMIN, DN );

      Z( J4+2 ) = DN;
      Z( 4*N0-PP ) = EMIN;
      return;

      // End of SLASQ6

      }
