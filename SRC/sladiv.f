      SUBROUTINE SLADIV( A, B, C, D, P, Q )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               BS
      const              BS = 2.0E0 ;
      REAL               HALF
      const              HALF = 0.5E0 ;
      REAL               TWO
      const              TWO = 2.0E0 ;

      // .. Local Scalars ..
      REAL               AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLADIV1
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0E0

      OV = SLAMCH( 'Overflow threshold' )
      UN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Epsilon' )
      BE = BS / (EPS*EPS)

      if ( AB >= HALF*OV ) {
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      }
      if ( CD >= HALF*OV ) {
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      }
      if ( AB <= UN*BS/EPS ) {
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      }
      if ( CD <= UN*BS/EPS ) {
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      }
      if ( ABS( D ) <= ABS( C ) ) {
         sladiv1(AA, BB, CC, DD, P, Q);
      } else {
         sladiv1(BB, AA, DD, CC, P, Q);
         Q = -Q
      }
      P = P * S
      Q = Q * S

      RETURN

      // End of SLADIV

      }

*> \ingroup ladiv


      SUBROUTINE SLADIV1( A, B, C, D, P, Q )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;

      // .. Local Scalars ..
      REAL               R, T
      // ..
      // .. External Functions ..
      REAL               SLADIV2
      // EXTERNAL SLADIV2
      // ..
      // .. Executable Statements ..

      R = D / C
      T = ONE / (C + D * R)
      P = SLADIV2(A, B, C, D, R, T)
      A = -A
      Q = SLADIV2(B, A, C, D, R, T)

      RETURN

      // End of SLADIV1

      }

*> \ingroup ladiv

      REAL FUNCTION SLADIV2( A, B, C, D, R, T )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               A, B, C, D, R, T
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;

      // .. Local Scalars ..
      REAL               BR
      // ..
      // .. Executable Statements ..

      if ( R != ZERO ) {
         BR = B * R
         if (BR != ZERO) {
            SLADIV2 = (A + BR) * T
         } else {
            SLADIV2 = A * T + (B * T) * R
         }
      } else {
         SLADIV2 = (A + D * (B / C)) * T
      }

      RETURN

      // End of SLADIV2

      }
