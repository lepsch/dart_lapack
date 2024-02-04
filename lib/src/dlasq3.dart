      void dlasq3(I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, DN2, G, TAU ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               IEEE;
      int                I0, ITER, N0, NDIV, NFAIL, PP;
      double             DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, QMAX, SIGMA, TAU;
      // ..
      // .. Array Arguments ..
      double             Z( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             CBIAS;
      const              CBIAS = 1.50 ;
      double             ZERO, QURTR, HALF, ONE, TWO, HUNDRD;
      const              ZERO = 0.0, QURTR = 0.250, HALF = 0.5, ONE = 1.0, TWO = 2.0, HUNDRD = 100.0 ;
      // ..
      // .. Local Scalars ..
      int                IPN4, J4, N0IN, NN, TTYPE;
      double             EPS, S, T, TEMP, TOL, TOL2;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASQ4, DLASQ5, DLASQ6
      // ..
      // .. External Function ..
      double             DLAMCH;
      bool               DISNAN;
      // EXTERNAL DISNAN, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      N0IN = N0;
      EPS = DLAMCH( 'Precision' );
      TOL = EPS*HUNDRD;
      TOL2 = TOL**2;

      // Check for deflation.

      } // 10

      if (N0 < I0) return;
      IF( N0 == I0 ) GO TO 20;
      NN = 4*N0 + PP;
      if( N0 == ( I0+1 ) ) GO TO 40;

      // Check whether E(N0-1) is negligible, 1 eigenvalue.

      if( Z( NN-5 ) > TOL2*( SIGMA+Z( NN-3 ) ) && Z( NN-2*PP-4 ) > TOL2*Z( NN-7 ) ) GO TO 30;

      } // 20

      Z[4*N0-3] = Z( 4*N0+PP-3 ) + SIGMA;
      N0 = N0 - 1;
      GO TO 10;

      // Check  whether E(N0-2) is negligible, 2 eigenvalues.

      } // 30

      if( Z( NN-9 ) > TOL2*SIGMA && Z( NN-2*PP-8 ) > TOL2*Z( NN-11 ) ) GO TO 50;

      } // 40

      if ( Z( NN-3 ) > Z( NN-7 ) ) {
         S = Z( NN-3 );
         Z[NN-3] = Z( NN-7 );
         Z[NN-7] = S;
      }
      T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) );
      if ( Z( NN-5 ) > Z( NN-3 )*TOL2 && T != ZERO ) {
         S = Z( NN-3 )*( Z( NN-5 ) / T );
         if ( S <= T ) {
            S = Z( NN-3 )*( Z( NN-5 ) / ( T*( ONE+sqrt( ONE+S / T ) ) ) );
         } else {
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+sqrt( T )*sqrt( T+S ) ) );
         }
         T = Z( NN-7 ) + ( S+Z( NN-5 ) );
         Z[NN-3] = Z( NN-3 )*( Z( NN-7 ) / T );
         Z[NN-7] = T;
      }
      Z[4*N0-7] = Z( NN-7 ) + SIGMA;
      Z[4*N0-3] = Z( NN-3 ) + SIGMA;
      N0 = N0 - 2;
      GO TO 10;

      } // 50
      if (PP == 2) PP = 0;

      // Reverse the qd-array, if warranted.

      if ( DMIN <= ZERO || N0 < N0IN ) {
         if ( CBIAS*Z( 4*I0+PP-3 ) < Z( 4*N0+PP-3 ) ) {
            IPN4 = 4*( I0+N0 );
            for (J4 = 4*I0; J4 <= 2*( I0+N0-1 ); J4 += 4) { // 60
               TEMP = Z( J4-3 );
               Z[J4-3] = Z( IPN4-J4-3 );
               Z[IPN4-J4-3] = TEMP;
               TEMP = Z( J4-2 );
               Z[J4-2] = Z( IPN4-J4-2 );
               Z[IPN4-J4-2] = TEMP;
               TEMP = Z( J4-1 );
               Z[J4-1] = Z( IPN4-J4-5 );
               Z[IPN4-J4-5] = TEMP;
               TEMP = Z( J4 );
               Z[J4] = Z( IPN4-J4-4 );
               Z[IPN4-J4-4] = TEMP;
            } // 60
            if ( N0-I0 <= 4 ) {
               Z[4*N0+PP-1] = Z( 4*I0+PP-1 );
               Z[4*N0-PP] = Z( 4*I0-PP );
            }
            DMIN2 = min( DMIN2, Z( 4*N0+PP-1 ) );
            Z[4*N0+PP-1] = min( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ), Z( 4*I0+PP+3 ) )             Z( 4*N0-PP ) = min( Z( 4*N0-PP ), Z( 4*I0-PP ), Z( 4*I0-PP+4 ) );
            QMAX = max( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) );
            DMIN = -ZERO;
         }
      }

      // Choose a shift.

      dlasq4(I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G );

      // Call dqds until DMIN > 0.

      } // 70

      dlasq5(I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, DN1, DN2, IEEE, EPS );

      NDIV = NDIV + ( N0-I0+2 );
      ITER = ITER + 1;

      // Check status.

      if ( DMIN >= ZERO && DMIN1 >= ZERO ) {

         // Success.

         GO TO 90;

      } else if ( DMIN < ZERO && DMIN1 > ZERO && Z( 4*( N0-1 )-PP ) < TOL*( SIGMA+DN1 ) && ( DN ).abs() < TOL*SIGMA ) {

         // Convergence hidden by negative DN.

         Z[4*( N0-1 )-PP+2] = ZERO;
         DMIN = ZERO;
         GO TO 90;
      } else if ( DMIN < ZERO ) {

         // TAU too big. Select new TAU and try again.

         NFAIL = NFAIL + 1;
         if ( TTYPE < -22 ) {

            // Failed twice. Play it safe.

            TAU = ZERO;
         } else if ( DMIN1 > ZERO ) {

            // Late failure. Gives excellent shift.

            TAU = ( TAU+DMIN )*( ONE-TWO*EPS );
            TTYPE = TTYPE - 11;
         } else {

            // Early failure. Divide by 4.

            TAU = QURTR*TAU;
            TTYPE = TTYPE - 12;
         }
         GO TO 70;
      } else if ( DISNAN( DMIN ) ) {

         // NaN.

         if ( TAU == ZERO ) {
            GO TO 80;
         } else {
            TAU = ZERO;
            GO TO 70;
         }
      } else {

         // Possible underflow. Play it safe.

         GO TO 80;
      }

      // Risk of underflow.

      } // 80
      dlasq6(I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 );
      NDIV = NDIV + ( N0-I0+2 );
      ITER = ITER + 1;
      TAU = ZERO;

      } // 90
      if ( TAU < SIGMA ) {
         DESIG = DESIG + TAU;
         T = SIGMA + DESIG;
         DESIG = DESIG - ( T-SIGMA );
      } else {
         T = SIGMA + TAU;
         DESIG = SIGMA - ( T-TAU ) + DESIG;
      }
      SIGMA = T;

      return;
      }
