      void cupmtr(SIDE, UPLO, TRANS, M, N, AP, TAU, final Matrix<double> C, final int LDC, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS, UPLO;
      int                INFO, LDC, M, N;
      Complex            AP( * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               FORWRD, LEFT, NOTRAN, UPPER;
      int                I, I1, I2, I3, IC, II, JC, MI, NI, NQ;
      Complex            AII, TAUI;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX

      // Test the input arguments

      INFO = 0;
      LEFT = lsame( SIDE, 'L' );
      NOTRAN = lsame( TRANS, 'N' );
      UPPER = lsame( UPLO, 'U' );

      // NQ is the order of Q

      if ( LEFT ) {
         NQ = M;
      } else {
         NQ = N;
      }
      if ( !LEFT && !lsame( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( !NOTRAN && !lsame( TRANS, 'C' ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('CUPMTR', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      if ( UPPER ) {

         // Q was determined by a call to CHPTRD with UPLO = 'U'

         FORWRD = ( LEFT && NOTRAN ) || ( !LEFT && !NOTRAN );

         if ( FORWRD ) {
            I1 = 1;
            I2 = NQ - 1;
            I3 = 1;
            II = 2;
         } else {
            I1 = NQ - 1;
            I2 = 1;
            I3 = -1;
            II = NQ*( NQ+1 ) / 2 - 1;
         }

         if ( LEFT ) {
            NI = N;
         } else {
            MI = M;
         }

         for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) { // 10
            if ( LEFT ) {

               // H(i) or H(i)**H is applied to C(1:i,1:n)

               MI = I;
            } else {

               // H(i) or H(i)**H is applied to C(1:m,1:i)

               NI = I;
            }

            // Apply H(i) or H(i)**H

            if ( NOTRAN ) {
               TAUI = TAU( I );
            } else {
               TAUI = CONJG( TAU( I ) );
            }
            AII = AP( II );
            AP[II] = ONE;
            clarf(SIDE, MI, NI, AP( II-I+1 ), 1, TAUI, C, LDC, WORK );
            AP[II] = AII;

            if ( FORWRD ) {
               II = II + I + 2;
            } else {
               II = II - I - 1;
            }
         } // 10
      } else {

         // Q was determined by a call to CHPTRD with UPLO = 'L'.

         FORWRD = ( LEFT && !NOTRAN ) || ( !LEFT && NOTRAN );

         if ( FORWRD ) {
            I1 = 1;
            I2 = NQ - 1;
            I3 = 1;
            II = 2;
         } else {
            I1 = NQ - 1;
            I2 = 1;
            I3 = -1;
            II = NQ*( NQ+1 ) / 2 - 1;
         }

         if ( LEFT ) {
            NI = N;
            JC = 1;
         } else {
            MI = M;
            IC = 1;
         }

         for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) { // 20
            AII = AP( II );
            AP[II] = ONE;
            if ( LEFT ) {

               // H(i) or H(i)**H is applied to C(i+1:m,1:n)

               MI = M - I;
               IC = I + 1;
            } else {

               // H(i) or H(i)**H is applied to C(1:m,i+1:n)

               NI = N - I;
               JC = I + 1;
            }

            // Apply H(i) or H(i)**H

            if ( NOTRAN ) {
               TAUI = TAU( I );
            } else {
               TAUI = CONJG( TAU( I ) );
            }
            clarf(SIDE, MI, NI, AP( II ), 1, TAUI, C( IC, JC ), LDC, WORK );
            AP[II] = AII;

            if ( FORWRD ) {
               II = II + NQ - I + 1;
            } else {
               II = II - NQ + I - 2;
            }
         } // 20
      }
      }
