      SUBROUTINE SLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, DN2, G, TAU )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               IEEE;
      int                I0, ITER, N0, NDIV, NFAIL, PP;
      REAL               DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, QMAX, SIGMA, TAU
      // ..
      // .. Array Arguments ..
      REAL               Z( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               CBIAS
      const              CBIAS = 1.50E0 ;
      REAL               ZERO, QURTR, HALF, ONE, TWO, HUNDRD
      const              ZERO = 0.0E0, QURTR = 0.250E0, HALF = 0.5E0, ONE = 1.0E0, TWO = 2.0E0, HUNDRD = 100.0E0 ;
      // ..
      // .. Local Scalars ..
      int                IPN4, J4, N0IN, NN, TTYPE;
      REAL               EPS, S, T, TEMP, TOL, TOL2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASQ4, SLASQ5, SLASQ6
      // ..
      // .. External Function ..
      REAL               SLAMCH
      bool               SISNAN;
      // EXTERNAL SISNAN, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      N0IN = N0
      EPS = SLAMCH( 'Precision' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2

      // Check for deflation.

      } // 10

      if (N0.LT.I0) RETURN       IF( N0.EQ.I0 ) GO TO 20;
      NN = 4*N0 + PP
      IF( N0.EQ.( I0+1 ) ) GO TO 40

      // Check whether E(N0-1) is negligible, 1 eigenvalue.

      IF( Z( NN-5 ).GT.TOL2*( SIGMA+Z( NN-3 ) ) .AND. Z( NN-2*PP-4 ).GT.TOL2*Z( NN-7 ) ) GO TO 30

      } // 20

      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10

      // Check  whether E(N0-2) is negligible, 2 eigenvalues.

      } // 30

      IF( Z( NN-9 ).GT.TOL2*SIGMA .AND. Z( NN-2*PP-8 ).GT.TOL2*Z( NN-11 ) ) GO TO 50

      } // 40

      if ( Z( NN-3 ).GT.Z( NN-7 ) ) {
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      }
      T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
      if ( Z( NN-5 ).GT.Z( NN-3 )*TOL2.AND.T.NE.ZERO ) {
         S = Z( NN-3 )*( Z( NN-5 ) / T )
         if ( S.LE.T ) {
            S = Z( NN-3 )*( Z( NN-5 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
         } else {
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
         }
         T = Z( NN-7 ) + ( S+Z( NN-5 ) )
         Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
         Z( NN-7 ) = T
      }
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10

      } // 50
      if (PP.EQ.2) PP = 0;

      // Reverse the qd-array, if warranted.

      if ( DMIN.LE.ZERO .OR. N0.LT.N0IN ) {
         if ( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) {
            IPN4 = 4*( I0+N0 )
            DO 60 J4 = 4*I0, 2*( I0+N0-1 ), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
            } // 60
            if ( N0-I0.LE.4 ) {
               Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
               Z( 4*N0-PP ) = Z( 4*I0-PP )
            }
            DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
            Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ), Z( 4*I0+PP+3 ) )             Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ), Z( 4*I0-PP+4 ) )
            QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
            DMIN = -ZERO
         }
      }

      // Choose a shift.

      slasq4(I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G );

      // Call dqds until DMIN > 0.

      } // 70

      slasq5(I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, DN1, DN2, IEEE, EPS );

      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1

      // Check status.

      if ( DMIN.GE.ZERO .AND. DMIN1.GE.ZERO ) {

         // Success.

         GO TO 90

      } else if ( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND. Z( 4*( N0-1 )-PP ).LT.TOL*( SIGMA+DN1 ) .AND. ABS( DN ).LT.TOL*SIGMA ) {

         // Convergence hidden by negative DN.

         Z( 4*( N0-1 )-PP+2 ) = ZERO
         DMIN = ZERO
         GO TO 90
      } else if ( DMIN.LT.ZERO ) {

         // TAU too big. Select new TAU and try again.

         NFAIL = NFAIL + 1
         if ( TTYPE.LT.-22 ) {

            // Failed twice. Play it safe.

            TAU = ZERO
         } else if ( DMIN1.GT.ZERO ) {

            // Late failure. Gives excellent shift.

            TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
            TTYPE = TTYPE - 11
         } else {

            // Early failure. Divide by 4.

            TAU = QURTR*TAU
            TTYPE = TTYPE - 12
         }
         GO TO 70
      } else if ( SISNAN( DMIN ) ) {

         // NaN.

         if ( TAU.EQ.ZERO ) {
            GO TO 80
         } else {
            TAU = ZERO
            GO TO 70
         }
      } else {

         // Possible underflow. Play it safe.

         GO TO 80
      }

      // Risk of underflow.

      } // 80
      slasq6(I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 );
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
      TAU = ZERO

      } // 90
      if ( TAU.LT.SIGMA ) {
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - ( T-SIGMA )
      } else {
         T = SIGMA + TAU
         DESIG = SIGMA - ( T-TAU ) + DESIG
      }
      SIGMA = T

      RETURN

      // End of SLASQ3

      }
