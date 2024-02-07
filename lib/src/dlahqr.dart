import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, INFO ) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N;
      bool               WANTT, WANTZ;
      double             H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * );
      // ..

// =========================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double             DAT1, DAT2;
      const              DAT1 = 3.0 / 4.0, DAT2 = -0.4375 ;
      int                KEXSH;
      const              KEXSH = 10 ;
      double             AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S, H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX, SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST, ULP, V2, V3;
      int                I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ, KDEFL;
      double             V( 3 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLANV2, DLARFG, DROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT

      INFO = 0;

      // Quick return if possible

      if (N == 0) return;
      if ( ILO == IHI ) {
         WR[ILO] = H( ILO, ILO );
         WI[ILO] = ZERO;
         return;
      }

      // ==== clear out the trash ====
      for (J = ILO; J <= IHI - 3; J++) { // 10
         H[J+2][J] = ZERO;
         H[J+3][J] = ZERO;
      } // 10
      if (ILO <= IHI-2) H( IHI, IHI-2 ) = ZERO;

      NH = IHI - ILO + 1;
      NZ = IHIZ - ILOZ + 1;

      // Set machine-dependent constants for the stopping criterion.

      SAFMIN = dlamch( 'SAFE MINIMUM' );
      SAFMAX = ONE / SAFMIN;
      ULP = dlamch( 'PRECISION' );
      SMLNUM = SAFMIN*( NH.toDouble() / ULP );

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
      // IHI to ILO in steps of 1 or 2. Each iteration of the loop works
      // with the active submatrix in rows and columns L to I.
      // Eigenvalues I+1 to IHI have already converged. Either L = ILO or
      // H(L,L-1) is negligible so that the matrix splits.

      I = IHI;
      } // 20
      L = ILO;
      if (I < ILO) GO TO 160;

      // Perform QR iterations on rows and columns ILO to I until a
      // submatrix of order 1 or 2 splits off at the bottom because a
      // subdiagonal element has become negligible.

      for (ITS = 0; ITS <= ITMAX; ITS++) { // 140

         // Look for a single small subdiagonal element.

         for (K = I; K >= L + 1; K--) { // 30
            if( ( H( K, K-1 ) ).abs() <= SMLNUM ) GO TO 40;
            TST = ( H( K-1, K-1 ) ).abs() + ( H( K, K ) ).abs();
            if ( TST == ZERO ) {
               if (K-2 >= ILO) TST = TST + ( H( K-1, K-2 ) ).abs();
               IF( K+1 <= IHI ) TST = TST + ( H( K+1, K ) ).abs();
            }
            // ==== The following is a conservative small subdiagonal
            // .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
            // .    1997). It has better mathematical foundation and
            // .    improves accuracy in some cases.  ====
            if ( ( H( K, K-1 ) ).abs() <= ULP*TST ) {
               AB = max( ( H( K, K-1 ) ).abs(), ( H( K-1, K ) ).abs() );
               BA = min( ( H( K, K-1 ) ).abs(), ( H( K-1, K ) ).abs() );
               AA = max( ( H( K, K ) ).abs(), ABS( H( K-1, K-1 )-H( K, K ) ) )                BB = min( ( H( K, K ) ).abs(), ABS( H( K-1, K-1 )-H( K, K ) ) );
               S = AA + AB;
               if( BA*( AB / S ) <= max( SMLNUM, ULP*( BB*( AA / S ) ) ) )GO TO 40;
            }
         } // 30
         } // 40
         L = K;
         if ( L > ILO ) {

            // H(L,L-1) is negligible

            H[L][L-1] = ZERO;
         }

         // Exit from loop if a submatrix of order 1 or 2 has split off.

         if (L >= I-1) GO TO 150;
         KDEFL = KDEFL + 1;

         // Now the active submatrix is in rows and columns L to I. If
         // eigenvalues only are being computed, only the active submatrix
         // need be transformed.

         if ( !WANTT ) {
            I1 = L;
            I2 = I;
         }

         if ( (KDEFL % 2*KEXSH) == 0 ) {

            // Exceptional shift.

            S = ( H( I, I-1 ) ).abs() + ( H( I-1, I-2 ) ).abs();
            H11 = DAT1*S + H( I, I );
            H12 = DAT2*S;
            H21 = S;
            H22 = H11;
         } else if ( (KDEFL % KEXSH) == 0 ) {

            // Exceptional shift.

            S = ( H( L+1, L ) ).abs() + ( H( L+2, L+1 ) ).abs();
            H11 = DAT1*S + H( L, L );
            H12 = DAT2*S;
            H21 = S;
            H22 = H11;
         } else {

            // Prepare to use Francis' double shift
            // (i.e. 2nd degree generalized Rayleigh quotient)

            H11 = H( I-1, I-1 );
            H21 = H( I, I-1 );
            H12 = H( I-1, I );
            H22 = H( I, I );
         }
         S = ( H11 ).abs() + ( H12 ).abs() + ( H21 ).abs() + ( H22 ).abs();
         if ( S == ZERO ) {
            RT1R = ZERO;
            RT1I = ZERO;
            RT2R = ZERO;
            RT2I = ZERO;
         } else {
            H11 = H11 / S;
            H21 = H21 / S;
            H12 = H12 / S;
            H22 = H22 / S;
            TR = ( H11+H22 ) / TWO;
            DET = ( H11-TR )*( H22-TR ) - H12*H21;
            RTDISC = sqrt( ( DET ).abs() );
            if ( DET >= ZERO ) {

               // ==== complex conjugate shifts ====

               RT1R = TR*S;
               RT2R = RT1R;
               RT1I = RTDISC*S;
               RT2I = -RT1I;
            } else {

               // ==== real shifts (use only one of them)  ====

               RT1R = TR + RTDISC;
               RT2R = TR - RTDISC;
               if ( ( RT1R-H22 ).abs() <= ( RT2R-H22 ).abs() ) {
                  RT1R = RT1R*S;
                  RT2R = RT1R;
               } else {
                  RT2R = RT2R*S;
                  RT1R = RT2R;
               }
               RT1I = ZERO;
               RT2I = ZERO;
            }
         }

         // Look for two consecutive small subdiagonal elements.

         for (M = I - 2; M >= L; M--) { // 50
            // Determine the effect of starting the double-shift QR
            // iteration at row M, and see if this would make H(M,M-1)
            // negligible.  (The following uses scaling to avoid
            // overflows and most underflows.)

            H21S = H( M+1, M );
            S = ABS( H( M, M )-RT2R ) + ( RT2I ).abs() + ( H21S ).abs();
            H21S = H( M+1, M ) / S;
            V[1] = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )* ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S );
            V[2] = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R );
            V[3] = H21S*H( M+2, M+1 );
            S = ( V( 1 ) ).abs() + ( V( 2 ) ).abs() + ( V( 3 ) ).abs();
            V[1] = V( 1 ) / S;
            V[2] = V( 2 ) / S;
            V[3] = V( 3 ) / S;
            if (M == L) GO TO 60;
            IF( ( H( M, M-1 ) ).abs()*( ( V( 2 ) ).abs()+( V( 3 ) ).abs() ) <= ULP*( V( 1 ) ).abs()*( ( H( M-1, M-1 ) ).abs()+( H( M, M ) ).abs()+( H( M+1, M+1 ) ).abs() ) )GO TO 60;
         } // 50
         } // 60

         // Double-shift QR step

         for (K = M; K <= I - 1; K++) { // 130

            // The first iteration of this loop determines a reflection G
            // from the vector V and applies it from left and right to H,
            // thus creating a nonzero bulge below the subdiagonal.

            // Each subsequent iteration determines a reflection G to
            // restore the Hessenberg form in the (K-1)th column, and thus
            // chases the bulge one step toward the bottom of the active
            // submatrix. NR is the order of G.

            NR = min( 3, I-K+1 );
            if (K > M) dcopy( NR, H( K, K-1 ), 1, V, 1 );
            dlarfg(NR, V( 1 ), V( 2 ), 1, T1 );
            if ( K > M ) {
               H[K][K-1] = V( 1 );
               H[K+1][K-1] = ZERO;
               if (K < I-1) H( K+2, K-1 ) = ZERO;
            } else if ( M > L ) {
                // ==== Use the following instead of
                // .    H( K, K-1 ) = -H( K, K-1 ) to
                // .    avoid a bug when v(2) and v(3)
                // .    underflow. ====
               H[K][K-1] = H( K, K-1 )*( ONE-T1 );
            }
            V2 = V( 2 );
            T2 = T1*V2;
            if ( NR == 3 ) {
               V3 = V( 3 );
               T3 = T1*V3;

               // Apply G from the left to transform the rows of the matrix
               // in columns K to I2.

               for (J = K; J <= I2; J++) { // 70
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J );
                  H[K][J] = H( K, J ) - SUM*T1;
                  H[K+1][J] = H( K+1, J ) - SUM*T2;
                  H[K+2][J] = H( K+2, J ) - SUM*T3;
               } // 70

               // Apply G from the right to transform the columns of the
               // matrix in rows I1 to min(K+3,I).

               for (J = I1; J <= min( K+3, I ); J++) { // 80
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 );
                  H[J][K] = H( J, K ) - SUM*T1;
                  H[J][K+1] = H( J, K+1 ) - SUM*T2;
                  H[J][K+2] = H( J, K+2 ) - SUM*T3;
               } // 80

               if ( WANTZ ) {

                  // Accumulate transformations in the matrix Z

                  for (J = ILOZ; J <= IHIZ; J++) { // 90
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 );
                     Z[J][K] = Z( J, K ) - SUM*T1;
                     Z[J][K+1] = Z( J, K+1 ) - SUM*T2;
                     Z[J][K+2] = Z( J, K+2 ) - SUM*T3;
                  } // 90
               }
            } else if ( NR == 2 ) {

               // Apply G from the left to transform the rows of the matrix
               // in columns K to I2.

               for (J = K; J <= I2; J++) { // 100
                  SUM = H( K, J ) + V2*H( K+1, J );
                  H[K][J] = H( K, J ) - SUM*T1;
                  H[K+1][J] = H( K+1, J ) - SUM*T2;
               } // 100

               // Apply G from the right to transform the columns of the
               // matrix in rows I1 to min(K+3,I).

               for (J = I1; J <= I; J++) { // 110
                  SUM = H( J, K ) + V2*H( J, K+1 );
                  H[J][K] = H( J, K ) - SUM*T1;
                  H[J][K+1] = H( J, K+1 ) - SUM*T2;
               } // 110

               if ( WANTZ ) {

                  // Accumulate transformations in the matrix Z

                  for (J = ILOZ; J <= IHIZ; J++) { // 120
                     SUM = Z( J, K ) + V2*Z( J, K+1 );
                     Z[J][K] = Z( J, K ) - SUM*T1;
                     Z[J][K+1] = Z( J, K+1 ) - SUM*T2;
                  } // 120
               }
            }
         } // 130

      } // 140

      // Failure to converge in remaining number of iterations

      INFO = I;
      return;

      } // 150

      if ( L == I ) {

         // H(I,I-1) is negligible: one eigenvalue has converged.

         WR[I] = H( I, I );
         WI[I] = ZERO;
      } else if ( L == I-1 ) {

         // H(I-1,I-2) is negligible: a pair of eigenvalues have converged.

         // Transform the 2-by-2 submatrix to standard Schur form,
         // and compute and store the eigenvalues.

         dlanv2(H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ), H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), CS, SN );

         if ( WANTT ) {

            // Apply the transformation to the rest of H.

            if (I2 > I) drot( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH, CS, SN );
            drot(I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN );
         }
         if ( WANTZ ) {

            // Apply the transformation to Z.

            drot(NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN );
         }
      }
      // reset deflation counter
      KDEFL = 0;

      // return to start of the main loop with new value of I.

      I = L - 1;
      GO TO 20;

      } // 160
      return;
      }
