      void slasq4(final int I0, final int N0, final int Z, final int PP, final int N0IN, final int DMIN, final int DMIN1, final int DMIN2, final int DN, final int DN1, final int DN2, final int TAU, final int TTYPE, final int G,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                I0, N0, N0IN, PP, TTYPE;
      double               DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU;
      double               Z( * );
      // ..

      double               CNST1, CNST2, CNST3;
      const              CNST1 = 0.5630, CNST2 = 1.010, CNST3 = 1.050 ;
      double               QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD;
      const              QURTR = 0.250, THIRD = 0.3330, HALF = 0.50, ZERO = 0.0, ONE = 1.0, TWO = 2.0, HUNDRD = 100.0 ;
      int                I4, NN, NP;
      double               A2, B1, B2, GAM, GAP1, GAP2, S;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT

      // A negative DMIN forces the shift to take that absolute value
      // TTYPE records the type of shift.

      if ( DMIN <= ZERO ) {
         TAU = -DMIN;
         TTYPE = -1;
         return;
      }

      NN = 4*N0 + PP;
      if ( N0IN == N0 ) {

         // No eigenvalues deflated.

         if ( DMIN == DN || DMIN == DN1 ) {

            B1 = sqrt( Z( NN-3 ) )*sqrt( Z( NN-5 ) );
            B2 = sqrt( Z( NN-7 ) )*sqrt( Z( NN-9 ) );
            A2 = Z( NN-7 ) + Z( NN-5 );

            // Cases 2 and 3.

            if ( DMIN == DN && DMIN1 == DN1 ) {
               GAP2 = DMIN2 - A2 - DMIN2*QURTR;
               if ( GAP2 > ZERO && GAP2 > B2 ) {
                  GAP1 = A2 - DN - ( B2 / GAP2 )*B2;
               } else {
                  GAP1 = A2 - DN - ( B1+B2 );
               }
               if ( GAP1 > ZERO && GAP1 > B1 ) {
                  S = max( DN-( B1 / GAP1 )*B1, HALF*DMIN );
                  TTYPE = -2;
               } else {
                  S = ZERO;
                  if (DN > B1) S = DN - B1;
                  IF( A2 > ( B1+B2 ) ) S = min( S, A2-( B1+B2 ) );
                  S = max( S, THIRD*DMIN );
                  TTYPE = -3;
               }
            } else {

               // Case 4.

               TTYPE = -4;
               S = QURTR*DMIN;
               if ( DMIN == DN ) {
                  GAM = DN;
                  A2 = ZERO;
                  if( Z( NN-5 ) > Z( NN-7 ) ) return;
                  B2 = Z( NN-5 ) / Z( NN-7 );
                  NP = NN - 9;
               } else {
                  NP = NN - 2*PP;
                  GAM = DN1;
                  if( Z( NP-4 ) > Z( NP-2 ) ) return;
                  A2 = Z( NP-4 ) / Z( NP-2 );
                  if( Z( NN-9 ) > Z( NN-11 ) ) return;
                  B2 = Z( NN-9 ) / Z( NN-11 );
                  NP = NN - 13;
               }

               // Approximate contribution to norm squared from I < NN-1.

               A2 = A2 + B2;
               for (I4 = NP; I4 >= 4*I0 - 1 + PP; I4 -= 4) { // 10
                  if (B2 == ZERO) GO TO 20;
                  B1 = B2;
                  if( Z( I4 ) > Z( I4-2 ) ) return;
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) );
                  A2 = A2 + B2;
                  if( HUNDRD*max( B2, B1 ) < A2 || CNST1 < A2 ) GO TO 20;
               } // 10
               } // 20
               A2 = CNST3*A2;

               // Rayleigh quotient residual bound.

               if (A2 < CNST1) S = GAM*( ONE-sqrt( A2 ) ) / ( ONE+A2 );
            }
         } else if ( DMIN == DN2 ) {

            // Case 5.

            TTYPE = -5;
            S = QURTR*DMIN;

            // Compute contribution to norm squared from I > NN-2.

            NP = NN - 2*PP;
            B1 = Z( NP-2 );
            B2 = Z( NP-6 );
            GAM = DN2;
            if( Z( NP-8 ) > B2 || Z( NP-4 ) > B1 ) return;
            A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 );

            // Approximate contribution to norm squared from I < NN-2.

            if ( N0-I0 > 2 ) {
               B2 = Z( NN-13 ) / Z( NN-15 );
               A2 = A2 + B2;
               for (I4 = NN - 17; I4 >= 4*I0 - 1 + PP; I4 -= 4) { // 30
                  if (B2 == ZERO) GO TO 40;
                  B1 = B2;
                  if( Z( I4 ) > Z( I4-2 ) ) return;
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) );
                  A2 = A2 + B2;
                  if( HUNDRD*max( B2, B1 ) < A2 || CNST1 < A2 ) GO TO 40;
               } // 30
               } // 40
               A2 = CNST3*A2;
            }

            if (A2 < CNST1) S = GAM*( ONE-sqrt( A2 ) ) / ( ONE+A2 );
         } else {

            // Case 6, no information to guide us.

            if ( TTYPE == -6 ) {
               G = G + THIRD*( ONE-G );
            } else if ( TTYPE == -18 ) {
               G = QURTR*THIRD;
            } else {
               G = QURTR;
            }
            S = G*DMIN;
            TTYPE = -6;
         }

      } else if ( N0IN == ( N0+1 ) ) {

         // One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.

         if ( DMIN1 == DN1 && DMIN2 == DN2 ) {

            // Cases 7 and 8.

            TTYPE = -7;
            S = THIRD*DMIN1;
            if( Z( NN-5 ) > Z( NN-7 ) ) return;
            B1 = Z( NN-5 ) / Z( NN-7 );
            B2 = B1;
            if (B2 == ZERO) GO TO 60;
            for (I4 = 4*N0 - 9 + PP; I4 >= 4*I0 - 1 + PP; I4 -= 4) { // 50
               A2 = B1;
               if( Z( I4 ) > Z( I4-2 ) ) return;
               B1 = B1*( Z( I4 ) / Z( I4-2 ) );
               B2 = B2 + B1;
               if( HUNDRD*max( B1, A2 ) < B2 ) GO TO 60;
            } // 50
            } // 60
            B2 = sqrt( CNST3*B2 );
            A2 = DMIN1 / ( ONE+B2**2 );
            GAP2 = HALF*DMIN2 - A2;
            if ( GAP2 > ZERO && GAP2 > B2*A2 ) {
               S = max( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) );
            } else {
               S = max( S, A2*( ONE-CNST2*B2 ) );
               TTYPE = -8;
            }
         } else {

            // Case 9.

            S = QURTR*DMIN1;
            if (DMIN1 == DN1) S = HALF*DMIN1;
            TTYPE = -9;
         }

      } else if ( N0IN == ( N0+2 ) ) {

         // Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.

         // Cases 10 and 11.

         if ( DMIN2 == DN2 && TWO*Z( NN-5 ) < Z( NN-7 ) ) {
            TTYPE = -10;
            S = THIRD*DMIN2;
            if( Z( NN-5 ) > Z( NN-7 ) ) return;
            B1 = Z( NN-5 ) / Z( NN-7 );
            B2 = B1;
            if (B2 == ZERO) GO TO 80;
            for (I4 = 4*N0 - 9 + PP; I4 >= 4*I0 - 1 + PP; I4 -= 4) { // 70
               if( Z( I4 ) > Z( I4-2 ) ) return;
               B1 = B1*( Z( I4 ) / Z( I4-2 ) );
               B2 = B2 + B1;
               if (HUNDRD*B1 < B2) GO TO 80;
            } // 70
            } // 80
            B2 = sqrt( CNST3*B2 );
            A2 = DMIN2 / ( ONE+B2**2 );
            GAP2 = Z( NN-7 ) + Z( NN-9 ) - sqrt( Z( NN-11 ) )*sqrt( Z( NN-9 ) ) - A2;
            if ( GAP2 > ZERO && GAP2 > B2*A2 ) {
               S = max( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) );
            } else {
               S = max( S, A2*( ONE-CNST2*B2 ) );
            }
         } else {
            S = QURTR*DMIN2;
            TTYPE = -11;
         }
      } else if ( N0IN > ( N0+2 ) ) {

         // Case 12, more than two eigenvalues deflated. No information.

         S = ZERO;
         TTYPE = -12;
      }

      TAU = S;
      }
