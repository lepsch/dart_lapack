      SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
      // ..

*  =========================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO, RONE, HALF;
      const              RZERO = 0.0, RONE = 1.0, HALF = 0.5 ;
      double             DAT1;
      const              DAT1 = 3.0 / 4.0 ;
      int                KEXSH;
      const              KEXSH = 10 ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U, V2, X, Y;
      double             AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX, SAFMIN, SMLNUM, SX, T2, TST, ULP;
      int                I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M, NH, NZ, KDEFL;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         V( 2 )
      // ..
      // .. External Functions ..
      COMPLEX*16         ZLADIV
      double             DLAMCH;
      // EXTERNAL ZLADIV, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZLARFG, ZSCAL
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if (N == 0) RETURN;
      if ( ILO == IHI ) {
         W( ILO ) = H( ILO, ILO )
         RETURN
      }

      // ==== clear out the trash ====
      for (J = ILO; J <= IHI - 3; J++) { // 10
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
      } // 10
      if (ILO <= IHI-2) H( IHI, IHI-2 ) = ZERO;
      // ==== ensure that subdiagonal entries are real ====
      if ( WANTT ) {
         JLO = 1
         JHI = N
      } else {
         JLO = ILO
         JHI = IHI
      }
      for (I = ILO + 1; I <= IHI; I++) { // 20
         if ( DIMAG( H( I, I-1 ) ) != RZERO ) {
            // ==== The following redundant normalization
            // .    avoids problems with both gradual and
            // .    sudden underflow in ABS(H(I,I-1)) ====
            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
            SC = DCONJG( SC ) / ABS( SC )
            H( I, I-1 ) = ABS( H( I, I-1 ) )
            zscal(JHI-I+1, SC, H( I, I ), LDH );
            zscal(MIN( JHI, I+1 )-JLO+1, DCONJG( SC ), H( JLO, I ), 1 )             IF( WANTZ ) CALL ZSCAL( IHIZ-ILOZ+1, DCONJG( SC ), Z( ILOZ, I ), 1 );
         }
      } // 20

      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1

      // Set machine-dependent constants for the stopping criterion.

      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )

      // I1 and I2 are the indices of the first row and last column of H
      // to which transformations must be applied. If eigenvalues only are
      // being computed, I1 and I2 are set inside the main loop.

      if ( WANTT ) {
         I1 = 1
         I2 = N
      }

      // ITMAX is the total number of QR iterations allowed.

      ITMAX = 30 * MAX( 10, NH )

      // KDEFL counts the number of iterations since a deflation

      KDEFL = 0

      // The main loop begins here. I is the loop index and decreases from
      // IHI to ILO in steps of 1. Each iteration of the loop works
      // with the active submatrix in rows and columns L to I.
      // Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
      // H(L,L-1) is negligible so that the matrix splits.

      I = IHI
      } // 30
      if (I < ILO) GO TO 150;

      // Perform QR iterations on rows and columns ILO to I until a
      // submatrix of order 1 splits off at the bottom because a
      // subdiagonal element has become negligible.

      L = ILO
      for (ITS = 0; ITS <= ITMAX; ITS++) { // 130

         // Look for a single small subdiagonal element.

         DO 40 K = I, L + 1, -1
            IF( CABS1( H( K, K-1 ) ) <= SMLNUM ) GO TO 50
            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            if ( TST == ZERO ) {
               if (K-2 >= ILO) TST = TST + ABS( DBLE( H( K-1, K-2 ) ) )                IF( K+1 <= IHI ) TST = TST + ABS( DBLE( H( K+1, K ) ) );
            }
            // ==== The following is a conservative small subdiagonal
            // .    deflation criterion due to Ahues & Tisseur (LAWN 122,
            // .    1997). It has better mathematical foundation and
            // .    improves accuracy in some examples.  ====
            if ( ABS( DBLE( H( K, K-1 ) ) ) <= ULP*TST ) {
               AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               AA = MAX( CABS1( H( K, K ) ), CABS1( H( K-1, K-1 )-H( K, K ) ) )                BB = MIN( CABS1( H( K, K ) ), CABS1( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ) <= MAX( SMLNUM, ULP*( BB*( AA / S ) ) ) )GO TO 50
            }
         } // 40
         } // 50
         L = K
         if ( L > ILO ) {

            // H(L,L-1) is negligible

            H( L, L-1 ) = ZERO
         }

         // Exit from loop if a submatrix of order 1 has split off.

         if (L >= I) GO TO 140;
         KDEFL = KDEFL + 1

         // Now the active submatrix is in rows and columns L to I. If
         // eigenvalues only are being computed, only the active submatrix
         // need be transformed.

         if ( !WANTT ) {
            I1 = L
            I2 = I
         }

         if ( MOD(KDEFL,2*KEXSH) == 0 ) {

            // Exceptional shift.

            S = DAT1*ABS( DBLE( H( I, I-1 ) ) )
            T = S + H( I, I )
         } else if ( MOD(KDEFL,KEXSH) == 0 ) {

            // Exceptional shift.

            S = DAT1*ABS( DBLE( H( L+1, L ) ) )
            T = S + H( L, L )
         } else {

            // Wilkinson's shift.

            T = H( I, I )
            U = SQRT( H( I-1, I ) )*SQRT( H( I, I-1 ) )
            S = CABS1( U )
            if ( S != RZERO ) {
               X = HALF*( H( I-1, I-1 )-T )
               SX = CABS1( X )
               S = MAX( S, CABS1( X ) )
               Y = S*SQRT( ( X / S )**2+( U / S )**2 )
               if ( SX > RZERO ) {
                  IF( DBLE( X / SX )*DBLE( Y )+DIMAG( X / SX )* DIMAG( Y ) < RZERO )Y = -Y
               }
               T = T - U*ZLADIV( U, ( X+Y ) )
            }
         }

         // Look for two consecutive small subdiagonal elements.

         DO 60 M = I - 1, L + 1, -1

            // Determine the effect of starting the single-shift QR
            // iteration at row M, and see if this would make H(M,M-1)
            // negligible.

            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H11S = H11 - T
            H21 = DBLE( H( M+1, M ) )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            H10 = DBLE( H( M, M-1 ) )
            IF( ABS( H10 )*ABS( H21 ) <= ULP* ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) ) GO TO 70
         } // 60
         H11 = H( L, L )
         H22 = H( L+1, L+1 )
         H11S = H11 - T
         H21 = DBLE( H( L+1, L ) )
         S = CABS1( H11S ) + ABS( H21 )
         H11S = H11S / S
         H21 = H21 / S
         V( 1 ) = H11S
         V( 2 ) = H21
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

            // V(2) is always real before the call to ZLARFG, and hence
            // after the call T2 ( = T1*V(2) ) is also real.

            if (K > M) CALL ZCOPY( 2, H( K, K-1 ), 1, V, 1 );
            zlarfg(2, V( 1 ), V( 2 ), 1, T1 );
            if ( K > M ) {
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            }
            V2 = V( 2 )
            T2 = DBLE( T1*V2 )

            // Apply G from the left to transform the rows of the matrix
            // in columns K to I2.

            for (J = K; J <= I2; J++) { // 80
               SUM = DCONJG( T1 )*H( K, J ) + T2*H( K+1, J )
               H( K, J ) = H( K, J ) - SUM
               H( K+1, J ) = H( K+1, J ) - SUM*V2
            } // 80

            // Apply G from the right to transform the columns of the
            // matrix in rows I1 to min(K+2,I).

            DO 90 J = I1, MIN( K+2, I )
               SUM = T1*H( J, K ) + T2*H( J, K+1 )
               H( J, K ) = H( J, K ) - SUM
               H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
            } // 90

            if ( WANTZ ) {

               // Accumulate transformations in the matrix Z

               for (J = ILOZ; J <= IHIZ; J++) { // 100
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*DCONJG( V2 )
               } // 100
            }

            if ( K == M && M > L ) {

               // If the QR step was started at row M > L because two
               // consecutive small subdiagonals were found, then extra
               // scaling must be performed to ensure that H(M,M-1) remains
               // real.

               TEMP = ONE - T1
               TEMP = TEMP / ABS( TEMP )
               H( M+1, M ) = H( M+1, M )*DCONJG( TEMP )
               if (M+2 <= I) H( M+2, M+1 ) = H( M+2, M+1 )*TEMP;
               for (J = M; J <= I; J++) { // 110
                  if ( J != M+1 ) {
                     if (I2 > J) CALL ZSCAL( I2-J, TEMP, H( J, J+1 ), LDH );
                     zscal(J-I1, DCONJG( TEMP ), H( I1, J ), 1 );
                     if ( WANTZ ) {
                        zscal(NZ, DCONJG( TEMP ), Z( ILOZ, J ), 1 );
                     }
                  }
               } // 110
            }
         } // 120

         // Ensure that H(I,I-1) is real.

         TEMP = H( I, I-1 )
         if ( DIMAG( TEMP ) != RZERO ) {
            RTEMP = ABS( TEMP )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            if (I2 > I) CALL ZSCAL( I2-I, DCONJG( TEMP ), H( I, I+1 ), LDH );
            zscal(I-I1, TEMP, H( I1, I ), 1 );
            if ( WANTZ ) {
               zscal(NZ, TEMP, Z( ILOZ, I ), 1 );
            }
         }

      } // 130

      // Failure to converge in remaining number of iterations

      INFO = I
      RETURN

      } // 140

      // H(I,I-1) is negligible: one eigenvalue has converged.

      W( I ) = H( I, I )
      // reset deflation counter
      KDEFL = 0

      // return to start of the main loop with new value of I.

      I = L - 1
      GO TO 30

      } // 150
      RETURN

      // End of ZLAHQR

      }
