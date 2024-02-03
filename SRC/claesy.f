      void claesy(A, B, C, RT1, RT2, EVSCAL, CS1, SN1 ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            A, B, C, CS1, EVSCAL, RT1, RT2, SN1;
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      REAL               ONE;
      const              ONE = 1.0 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      REAL               HALF;
      const              HALF = 0.5 ;
      REAL               THRESH;
      const              THRESH = 0.1 ;
      // ..
      // .. Local Scalars ..
      REAL               BABS, EVNORM, TABS, Z;
      COMPLEX            S, T, TMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..


      // Special case:  The matrix is actually diagonal.
      // To avoid divide by zero later, we treat this case separately.

      if ( ABS( B ) == ZERO ) {
         RT1 = A;
         RT2 = C;
         if ( ABS( RT1 ) < ABS( RT2 ) ) {
            TMP = RT1;
            RT1 = RT2;
            RT2 = TMP;
            CS1 = ZERO;
            SN1 = ONE;
         } else {
            CS1 = ONE;
            SN1 = ZERO;
         }
      } else {

         // Compute the eigenvalues and eigenvectors.
         // The characteristic equation is
            // lambda **2 - (A+C) lambda + (A*C - B*B)
         // and we solve it using the quadratic formula.

         S = ( A+C )*HALF;
         T = ( A-C )*HALF;

         // Take the square root carefully to avoid over/under flow.

         BABS = ABS( B );
         TABS = ABS( T );
         Z = max( BABS, TABS );
         if (Z > ZERO) T = Z*sqrt( ( T / Z )**2+( B / Z )**2 );

         // Compute the two eigenvalues.  RT1 and RT2 are exchanged
         // if necessary so that RT1 will have the greater magnitude.

         RT1 = S + T;
         RT2 = S - T;
         if ( ABS( RT1 ) < ABS( RT2 ) ) {
            TMP = RT1;
            RT1 = RT2;
            RT2 = TMP;
         }

         // Choose CS1 = 1 and SN1 to satisfy the first equation, then
         // scale the components of this eigenvector so that the matrix
         // of eigenvectors X satisfies  X * X**T = I .  (No scaling is
         // done if the norm of the eigenvalue matrix is less than THRESH.)

         SN1 = ( RT1-A ) / B;
         TABS = ABS( SN1 );
         if ( TABS > ONE ) {
            T = TABS*sqrt( ( ONE / TABS )**2+( SN1 / TABS )**2 );
         } else {
            T = sqrt( CONE+SN1*SN1 );
         }
         EVNORM = ABS( T );
         if ( EVNORM >= THRESH ) {
            EVSCAL = CONE / T;
            CS1 = EVSCAL;
            SN1 = SN1*EVSCAL;
         } else {
            EVSCAL = ZERO;
         }
      }
      return;
      }
