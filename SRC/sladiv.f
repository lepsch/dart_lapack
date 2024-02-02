      SUBROUTINE SLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               BS
      PARAMETER          ( BS = 2.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
*
*     .. Local Scalars ..
      REAL               AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLADIV1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
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

      IF( AB >= HALF*OV ) THEN
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      END IF
      IF( CD >= HALF*OV ) THEN
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      END IF
      IF( AB <= UN*BS/EPS ) THEN
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      END IF
      IF( CD <= UN*BS/EPS ) THEN
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      END IF
      IF( ABS( D ).LE.ABS( C ) ) THEN
         CALL SLADIV1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL SLADIV1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
*
      RETURN
*
*     End of SLADIV
*
      END

*> \ingroup ladiv


      SUBROUTINE SLADIV1( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*
*     .. Local Scalars ..
      REAL               R, T
*     ..
*     .. External Functions ..
      REAL               SLADIV2
      EXTERNAL           SLADIV2
*     ..
*     .. Executable Statements ..
*
      R = D / C
      T = ONE / (C + D * R)
      P = SLADIV2(A, B, C, D, R, T)
      A = -A
      Q = SLADIV2(B, A, C, D, R, T)
*
      RETURN
*
*     End of SLADIV1
*
      END

*> \ingroup ladiv

      REAL FUNCTION SLADIV2( A, B, C, D, R, T )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               A, B, C, D, R, T
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*
*     .. Local Scalars ..
      REAL               BR
*     ..
*     .. Executable Statements ..
*
      IF( R.NE.ZERO ) THEN
         BR = B * R
         if( BR.NE.ZERO ) THEN
            SLADIV2 = (A + BR) * T
         ELSE
            SLADIV2 = A * T + (B * T) * R
         END IF
      ELSE
         SLADIV2 = (A + D * (B / C)) * T
      END IF
*
      RETURN
*
*     End of SLADIV2
*
      END