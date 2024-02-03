      SUBROUTINE SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I0, N0, N0IN, PP, TTYPE;
      REAL               DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
      // ..
      // .. Array Arguments ..
      REAL               Z( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               CNST1, CNST2, CNST3
      const              CNST1 = 0.5630E0, CNST2 = 1.010E0, CNST3 = 1.050E0 ;
      REAL               QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD
      const              QURTR = 0.250E0, THIRD = 0.3330E0, HALF = 0.50E0, ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, HUNDRD = 100.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I4, NN, NP;
      REAL               A2, B1, B2, GAM, GAP1, GAP2, S
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // A negative DMIN forces the shift to take that absolute value
      // TTYPE records the type of shift.

      if ( DMIN.LE.ZERO ) {
         TAU = -DMIN
         TTYPE = -1
         RETURN
      }

      NN = 4*N0 + PP
      if ( N0IN.EQ.N0 ) {

         // No eigenvalues deflated.

         if ( DMIN.EQ.DN .OR. DMIN.EQ.DN1 ) {

            B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
            B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
            A2 = Z( NN-7 ) + Z( NN-5 )

            // Cases 2 and 3.

            if ( DMIN.EQ.DN .AND. DMIN1.EQ.DN1 ) {
               GAP2 = DMIN2 - A2 - DMIN2*QURTR
               if ( GAP2.GT.ZERO .AND. GAP2.GT.B2 ) {
                  GAP1 = A2 - DN - ( B2 / GAP2 )*B2
               } else {
                  GAP1 = A2 - DN - ( B1+B2 )
               }
               if ( GAP1.GT.ZERO .AND. GAP1.GT.B1 ) {
                  S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
                  TTYPE = -2
               } else {
                  S = ZERO
                  IF( DN.GT.B1 ) S = DN - B1                   IF( A2.GT.( B1+B2 ) ) S = MIN( S, A2-( B1+B2 ) )
                  S = MAX( S, THIRD*DMIN )
                  TTYPE = -3
               }
            } else {

               // Case 4.

               TTYPE = -4
               S = QURTR*DMIN
               if ( DMIN.EQ.DN ) {
                  GAM = DN
                  A2 = ZERO
                  IF( Z( NN-5 ) .GT. Z( NN-7 ) ) RETURN
                  B2 = Z( NN-5 ) / Z( NN-7 )
                  NP = NN - 9
               } else {
                  NP = NN - 2*PP
                  GAM = DN1
                  IF( Z( NP-4 ) .GT. Z( NP-2 ) ) RETURN
                  A2 = Z( NP-4 ) / Z( NP-2 )
                  IF( Z( NN-9 ) .GT. Z( NN-11 ) ) RETURN
                  B2 = Z( NN-9 ) / Z( NN-11 )
                  NP = NN - 13
               }

               // Approximate contribution to norm squared from I < NN-1.

               A2 = A2 + B2
               DO 10 I4 = NP, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO ) GO TO 20
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) ) RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 ) GO TO 20
               } // 10
               } // 20
               A2 = CNST3*A2

               // Rayleigh quotient residual bound.

               IF( A2.LT.CNST1 ) S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
            }
         } else if ( DMIN.EQ.DN2 ) {

            // Case 5.

            TTYPE = -5
            S = QURTR*DMIN

            // Compute contribution to norm squared from I > NN-2.

            NP = NN - 2*PP
            B1 = Z( NP-2 )
            B2 = Z( NP-6 )
            GAM = DN2
            IF( Z( NP-8 ).GT.B2 .OR. Z( NP-4 ).GT.B1 ) RETURN
            A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )

            // Approximate contribution to norm squared from I < NN-2.

            if ( N0-I0.GT.2 ) {
               B2 = Z( NN-13 ) / Z( NN-15 )
               A2 = A2 + B2
               DO 30 I4 = NN - 17, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO ) GO TO 40
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) ) RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 ) GO TO 40
               } // 30
               } // 40
               A2 = CNST3*A2
            }

            IF( A2.LT.CNST1 ) S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
         } else {

            // Case 6, no information to guide us.

            if ( TTYPE.EQ.-6 ) {
               G = G + THIRD*( ONE-G )
            } else if ( TTYPE.EQ.-18 ) {
               G = QURTR*THIRD
            } else {
               G = QURTR
            }
            S = G*DMIN
            TTYPE = -6
         }

      } else if ( N0IN.EQ.( N0+1 ) ) {

         // One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.

         if ( DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2 ) {

            // Cases 7 and 8.

            TTYPE = -7
            S = THIRD*DMIN1
            IF( Z( NN-5 ).GT.Z( NN-7 ) ) RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO ) GO TO 60
            DO 50 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               A2 = B1
               IF( Z( I4 ).GT.Z( I4-2 ) ) RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*MAX( B1, A2 ).LT.B2 ) GO TO 60
            } // 50
            } // 60
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN1 / ( ONE+B2**2 )
            GAP2 = HALF*DMIN2 - A2
            if ( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) {
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            } else {
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
               TTYPE = -8
            }
         } else {

            // Case 9.

            S = QURTR*DMIN1
            IF( DMIN1.EQ.DN1 ) S = HALF*DMIN1
            TTYPE = -9
         }

      } else if ( N0IN.EQ.( N0+2 ) ) {

         // Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.

         // Cases 10 and 11.

         if ( DMIN2.EQ.DN2 .AND. TWO*Z( NN-5 ).LT.Z( NN-7 ) ) {
            TTYPE = -10
            S = THIRD*DMIN2
            IF( Z( NN-5 ).GT.Z( NN-7 ) ) RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO ) GO TO 80
            DO 70 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               IF( Z( I4 ).GT.Z( I4-2 ) ) RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*B1.LT.B2 ) GO TO 80
            } // 70
            } // 80
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN2 / ( ONE+B2**2 )
            GAP2 = Z( NN-7 ) + Z( NN-9 ) - SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
            if ( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) {
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            } else {
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
            }
         } else {
            S = QURTR*DMIN2
            TTYPE = -11
         }
      } else if ( N0IN.GT.( N0+2 ) ) {

         // Case 12, more than two eigenvalues deflated. No information.

         S = ZERO
         TTYPE = -12
      }

      TAU = S
      RETURN

      // End of SLASQ4

      }
