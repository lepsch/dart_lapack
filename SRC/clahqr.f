      void clahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, INFO ) {
      // IMPLICIT NONE

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * );
      // ..

// =========================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RZERO, RONE, HALF;
      const              RZERO = 0.0, RONE = 1.0, HALF = 0.5 ;
      REAL               DAT1;
      const              DAT1 = 3.0 / 4.0 ;
      int                KEXSH;
      const              KEXSH = 10 ;
      // ..
      // .. Local Scalars ..
      COMPLEX            CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U, V2, X, Y;
      REAL               AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX, SAFMIN, SMLNUM, SX, T2, TST, ULP;
      int                I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M, NH, NZ, KDEFL;
      // ..
      // .. Local Arrays ..
      COMPLEX            V( 2 );
      // ..
      // .. External Functions ..
      //- COMPLEX            CLADIV;
      //- REAL               SLAMCH;
      // EXTERNAL CLADIV, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLARFG, CSCAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, MIN, REAL, SQRT
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ( REAL( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      INFO = 0;

      // Quick return if possible

      if (N == 0) return;
      if ( ILO == IHI ) {
         W( ILO ) = H( ILO, ILO );
         return;
      }

      // ==== clear out the trash ====
      for (J = ILO; J <= IHI - 3; J++) { // 10
         H( J+2, J ) = ZERO;
         H( J+3, J ) = ZERO;
      } // 10
      if (ILO <= IHI-2) H( IHI, IHI-2 ) = ZERO;
      // ==== ensure that subdiagonal entries are real ====
      if ( WANTT ) {
         JLO = 1;
         JHI = N;
      } else {
         JLO = ILO;
         JHI = IHI;
      }
      for (I = ILO + 1; I <= IHI; I++) { // 20
         if ( AIMAG( H( I, I-1 ) ) != RZERO ) {
            // ==== The following redundant normalization
            // .    avoids problems with both gradual and
            // .    sudden underflow in ABS(H(I,I-1)) ====
            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) );
            SC = CONJG( SC ) / ( SC ).abs();
            H( I, I-1 ) = ( H( I, I-1 ) ).abs();
            cscal(JHI-I+1, SC, H( I, I ), LDH );
            cscal(min( JHI, I+1 )-JLO+1, CONJG( SC ), H( JLO, I ), 1 )             IF( WANTZ ) CALL CSCAL( IHIZ-ILOZ+1, CONJG( SC ), Z( ILOZ, I ), 1 );
         }
      } // 20

      NH = IHI - ILO + 1;
      NZ = IHIZ - ILOZ + 1;

      // Set machine-dependent constants for the stopping criterion.

      SAFMIN = SLAMCH( 'SAFE MINIMUM' );
      SAFMAX = RONE / SAFMIN;
      ULP = SLAMCH( 'PRECISION' );
      SMLNUM = SAFMIN*( REAL( NH ) / ULP );

      // I1 and I2 are the indices of the first row and last column of H
      // to which transformations must be applied. If eigenvalues only are
      // being computed, I1 and I2 are set inside the main loop.

      if ( WANTT ) {
         I1 = 1;
         I2 = N;
      }

      // ITMAX is the total number of QR iterations allowed.

      ITMAX = 30 * max( 10, NH );

      // KDEFL counts the number of iterations since a deflation

      KDEFL = 0;

      // The main loop begins here. I is the loop index and decreases from
      // IHI to ILO in steps of 1. Each iteration of the loop works
      // with the active submatrix in rows and columns L to I.
      // Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
      // H(L,L-1) is negligible so that the matrix splits.

      I = IHI;
      } // 30
      if (I < ILO) GO TO 150;

      // Perform QR iterations on rows and columns ILO to I until a
      // submatrix of order 1 splits off at the bottom because a
      // subdiagonal element has become negligible.

      L = ILO;
      for (ITS = 0; ITS <= ITMAX; ITS++) { // 130

         // Look for a single small subdiagonal element.

         DO 40 K = I, L + 1, -1;
            if( CABS1( H( K, K-1 ) ) <= SMLNUM ) GO TO 50;
            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) );
            if ( TST == ZERO ) {
               if (K-2 >= ILO) TST = TST + ABS( REAL( H( K-1, K-2 ) ) );
               IF( K+1 <= IHI ) TST = TST + ABS( REAL( H( K+1, K ) ) );
            }
            // ==== The following is a conservative small subdiagonal
            // .    deflation criterion due to Ahues & Tisseur (LAWN 122,
            // .    1997). It has better mathematical foundation and
            // .    improves accuracy in some examples.  ====
            if ( ABS( REAL( H( K, K-1 ) ) ) <= ULP*TST ) {
               AB = max( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) );
               BA = min( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) );
               AA = max( CABS1( H( K, K ) ), CABS1( H( K-1, K-1 )-H( K, K ) ) )                BB = min( CABS1( H( K, K ) ), CABS1( H( K-1, K-1 )-H( K, K ) ) );
               S = AA + AB;
               if( BA*( AB / S ) <= max( SMLNUM, ULP*( BB*( AA / S ) ) ) )GO TO 50;
            }
         } // 40
         } // 50
         L = K;
         if ( L > ILO ) {

            // H(L,L-1) is negligible

            H( L, L-1 ) = ZERO;
         }

         // Exit from loop if a submatrix of order 1 has split off.

         if (L >= I) GO TO 140;
         KDEFL = KDEFL + 1;

         // Now the active submatrix is in rows and columns L to I. If
         // eigenvalues only are being computed, only the active submatrix
         // need be transformed.

         if ( !WANTT ) {
            I1 = L;
            I2 = I;
         }

         if ( MOD(KDEFL,2*KEXSH) == 0 ) {

            // Exceptional shift.

            S = DAT1*ABS( REAL( H( I, I-1 ) ) );
            T = S + H( I, I );
         } else if ( MOD(KDEFL,KEXSH) == 0 ) {

            // Exceptional shift.

            S = DAT1*ABS( REAL( H( L+1, L ) ) );
            T = S + H( L, L );
         } else {

            // Wilkinson's shift.

            T = H( I, I );
            U = sqrt( H( I-1, I ) )*sqrt( H( I, I-1 ) );
            S = CABS1( U );
            if ( S != RZERO ) {
               X = HALF*( H( I-1, I-1 )-T );
               SX = CABS1( X );
               S = max( S, CABS1( X ) );
               Y = S*sqrt( ( X / S )**2+( U / S )**2 );
               if ( SX > RZERO ) {
                  if( REAL( X / SX )*REAL( Y )+AIMAG( X / SX )* AIMAG( Y ) < RZERO )Y = -Y;
               }
               T = T - U*CLADIV( U, ( X+Y ) );
            }
         }

         // Look for two consecutive small subdiagonal elements.

         DO 60 M = I - 1, L + 1, -1;

            // Determine the effect of starting the single-shift QR
            // iteration at row M, and see if this would make H(M,M-1)
            // negligible.

            H11 = H( M, M );
            H22 = H( M+1, M+1 );
            H11S = H11 - T;
            H21 = REAL( H( M+1, M ) );
            S = CABS1( H11S ) + ( H21 ).abs();
            H11S = H11S / S;
            H21 = H21 / S;
            V( 1 ) = H11S;
            V( 2 ) = H21;
            H10 = REAL( H( M, M-1 ) );
            if( ( H10 ).abs()*( H21 ).abs() <= ULP* ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) ) GO TO 70;
         } // 60
         H11 = H( L, L );
         H22 = H( L+1, L+1 );
         H11S = H11 - T;
         H21 = REAL( H( L+1, L ) );
         S = CABS1( H11S ) + ( H21 ).abs();
         H11S = H11S / S;
         H21 = H21 / S;
         V( 1 ) = H11S;
         V( 2 ) = H21;
         } // 70

         // Single-shift QR step

         for (K = M; K <= I - 1; K++) { // 120

            // The first iteration of this loop determines a reflection G
            // from the vector V and applies it from left and right to H,
            // thus creating a nonzero bulge below the subdiagonal.

            // Each subsequent iteration determines a reflection G to
            // restore the Hessenberg form in the (K-1)th column, and thus
            // chases the bulge one step toward the bottom of the active
            // submatrix.

            // V(2) is always real before the call to CLARFG, and hence
            // after the call T2 ( = T1*V(2) ) is also real.

            if (K > M) ccopy( 2, H( K, K-1 ), 1, V, 1 );
            clarfg(2, V( 1 ), V( 2 ), 1, T1 );
            if ( K > M ) {
               H( K, K-1 ) = V( 1 );
               H( K+1, K-1 ) = ZERO;
            }
            V2 = V( 2 );
            T2 = REAL( T1*V2 );

            // Apply G from the left to transform the rows of the matrix
            // in columns K to I2.

            for (J = K; J <= I2; J++) { // 80
               SUM = CONJG( T1 )*H( K, J ) + T2*H( K+1, J );
               H( K, J ) = H( K, J ) - SUM;
               H( K+1, J ) = H( K+1, J ) - SUM*V2;
            } // 80

            // Apply G from the right to transform the columns of the
            // matrix in rows I1 to min(K+2,I).

            for (J = I1; J <= min( K+2, I ); J++) { // 90
               SUM = T1*H( J, K ) + T2*H( J, K+1 );
               H( J, K ) = H( J, K ) - SUM;
               H( J, K+1 ) = H( J, K+1 ) - SUM*CONJG( V2 );
            } // 90

            if ( WANTZ ) {

               // Accumulate transformations in the matrix Z

               for (J = ILOZ; J <= IHIZ; J++) { // 100
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 );
                  Z( J, K ) = Z( J, K ) - SUM;
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*CONJG( V2 );
               } // 100
            }

            if ( K == M && M > L ) {

               // If the QR step was started at row M > L because two
               // consecutive small subdiagonals were found, then extra
               // scaling must be performed to ensure that H(M,M-1) remains
               // real.

               TEMP = ONE - T1;
               TEMP = TEMP / ( TEMP ).abs();
               H( M+1, M ) = H( M+1, M )*CONJG( TEMP );
               if (M+2 <= I) H( M+2, M+1 ) = H( M+2, M+1 )*TEMP;
               for (J = M; J <= I; J++) { // 110
                  if ( J != M+1 ) {
                     if (I2 > J) cscal( I2-J, TEMP, H( J, J+1 ), LDH );
                     cscal(J-I1, CONJG( TEMP ), H( I1, J ), 1 );
                     if ( WANTZ ) {
                        cscal(NZ, CONJG( TEMP ), Z( ILOZ, J ), 1 );
                     }
                  }
               } // 110
            }
         } // 120

         // Ensure that H(I,I-1) is real.

         TEMP = H( I, I-1 );
         if ( AIMAG( TEMP ) != RZERO ) {
            RTEMP = ( TEMP ).abs();
            H( I, I-1 ) = RTEMP;
            TEMP = TEMP / RTEMP;
            if (I2 > I) cscal( I2-I, CONJG( TEMP ), H( I, I+1 ), LDH );
            cscal(I-I1, TEMP, H( I1, I ), 1 );
            if ( WANTZ ) {
               cscal(NZ, TEMP, Z( ILOZ, I ), 1 );
            }
         }

      } // 130

      // Failure to converge in remaining number of iterations

      INFO = I;
      return;

      } // 140

      // H(I,I-1) is negligible: one eigenvalue has converged.

      W( I ) = H( I, I );
      // reset deflation counter
      KDEFL = 0;

      // return to start of the main loop with new value of I.

      I = L - 1;
      GO TO 30;

      } // 150
      return;
      }
