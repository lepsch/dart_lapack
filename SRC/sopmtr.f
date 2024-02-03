      void sopmtr(SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS, UPLO;
      int                INFO, LDC, M, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               FORWRD, LEFT, NOTRAN, UPPER;
      int                I, I1, I2, I3, IC, II, JC, MI, NI, NQ;
      REAL               AII;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      LEFT = LSAME( SIDE, 'L' );
      NOTRAN = LSAME( TRANS, 'N' );
      UPPER = LSAME( UPLO, 'U' );

      // NQ is the order of Q

      if ( LEFT ) {
         NQ = M;
      } else {
         NQ = N;
      }
      if ( !LEFT && !LSAME( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( !NOTRAN && !LSAME( TRANS, 'T' ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SOPMTR', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      if ( UPPER ) {

         // Q was determined by a call to SSPTRD with UPLO = 'U'

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

         DO 10 I = I1, I2, I3;
            if ( LEFT ) {

               // H(i) is applied to C(1:i,1:n)

               MI = I;
            } else {

               // H(i) is applied to C(1:m,1:i)

               NI = I;
            }

            // Apply H(i)

            AII = AP( II );
            AP( II ) = ONE;
            slarf(SIDE, MI, NI, AP( II-I+1 ), 1, TAU( I ), C, LDC, WORK );
            AP( II ) = AII;

            if ( FORWRD ) {
               II = II + I + 2;
            } else {
               II = II - I - 1;
            }
         } // 10
      } else {

         // Q was determined by a call to SSPTRD with UPLO = 'L'.

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

         DO 20 I = I1, I2, I3;
            AII = AP( II );
            AP( II ) = ONE;
            if ( LEFT ) {

               // H(i) is applied to C(i+1:m,1:n)

               MI = M - I;
               IC = I + 1;
            } else {

               // H(i) is applied to C(1:m,i+1:n)

               NI = N - I;
               JC = I + 1;
            }

            // Apply H(i)

            slarf(SIDE, MI, NI, AP( II ), 1, TAU( I ), C( IC, JC ), LDC, WORK );
            AP( II ) = AII;

            if ( FORWRD ) {
               II = II + NQ - I + 1;
            } else {
               II = II - NQ + I - 2;
            }
         } // 20
      }
      return;

      // End of SOPMTR

      }
