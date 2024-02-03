      void dladiv(A, B, C, D, P, Q ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             A, B, C, D, P, Q;
      // ..

// =====================================================================

      // .. Parameters ..
      double             BS;
      const              BS = 2.0 ;
      double             HALF;
      const              HALF = 0.5 ;
      double             TWO;
      const              TWO = 2.0 ;

      // .. Local Scalars ..
      double             AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLADIV1
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      AA = A;
      BB = B;
      CC = C;
      DD = D;
      AB = max( ABS(A), ABS(B) );
      CD = max( ABS(C), ABS(D) );
      S = 1.0;

      OV = DLAMCH( 'Overflow threshold' );
      UN = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Epsilon' );
      BE = BS / (EPS*EPS);

      if ( AB >= HALF*OV ) {
         AA = HALF * AA;
         BB = HALF * BB;
         S  = TWO * S;
      }
      if ( CD >= HALF*OV ) {
         CC = HALF * CC;
         DD = HALF * DD;
         S  = HALF * S;
      }
      if ( AB <= UN*BS/EPS ) {
         AA = AA * BE;
         BB = BB * BE;
         S  = S / BE;
      }
      if ( CD <= UN*BS/EPS ) {
         CC = CC * BE;
         DD = DD * BE;
         S  = S * BE;
      }
      if ( ABS( D ) <= ABS( C ) ) {
         dladiv1(AA, BB, CC, DD, P, Q);
      } else {
         dladiv1(BB, AA, DD, CC, P, Q);
         Q = -Q;
      }
      P = P * S;
      Q = Q * S;

      return;

      // End of DLADIV

      }

// > \ingroup ladiv


      void dladiv1(A, B, C, D, P, Q ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             A, B, C, D, P, Q;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;

      // .. Local Scalars ..
      double             R, T;
      // ..
      // .. External Functions ..
      double             DLADIV2;
      // EXTERNAL DLADIV2
      // ..
      // .. Executable Statements ..

      R = D / C;
      T = ONE / (C + D * R);
      P = DLADIV2(A, B, C, D, R, T);
      A = -A;
      Q = DLADIV2(B, A, C, D, R, T);

      return;

      // End of DLADIV1

      }

// > \ingroup ladiv

      double dladiv2(A, B, C, D, R, T ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             A, B, C, D, R, T;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;

      // .. Local Scalars ..
      double             BR;
      // ..
      // .. Executable Statements ..

      if ( R != ZERO ) {
         BR = B * R;
         if ( BR != ZERO ) {
            DLADIV2 = (A + BR) * T;
         } else {
            DLADIV2 = A * T + (B * T) * R;
         }
      } else {
         DLADIV2 = (A + D * (B / C)) * T;
      }

      return;

      // End of DLADIV2

      }
