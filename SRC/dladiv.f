      SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      double             A, B, C, D, P, Q;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             BS;
      PARAMETER          ( BS = 2.0D0 )
      double             HALF;
      PARAMETER          ( HALF = 0.5D0 )
      double             TWO;
      PARAMETER          ( TWO = 2.0D0 )
*
*     .. Local Scalars ..
      double             AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS;
*     ..
*     .. External Functions ..
      double             DLAMCH;
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV1
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
*     ..
*     .. Executable Statements ..
*
      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0D0

      OV = DLAMCH( 'Overflow threshold' )
      UN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Epsilon' )
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
         CALL DLADIV1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL DLADIV1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
*
      RETURN
*
*     End of DLADIV
*
      END

*> \ingroup ladiv


      SUBROUTINE DLADIV1( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      double             A, B, C, D, P, Q;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE;
      PARAMETER          ( ONE = 1.0D0 )
*
*     .. Local Scalars ..
      double             R, T;
*     ..
*     .. External Functions ..
      double             DLADIV2;
      EXTERNAL           DLADIV2
*     ..
*     .. Executable Statements ..
*
      R = D / C
      T = ONE / (C + D * R)
      P = DLADIV2(A, B, C, D, R, T)
      A = -A
      Q = DLADIV2(B, A, C, D, R, T)
*
      RETURN
*
*     End of DLADIV1
*
      END

*> \ingroup ladiv

      double           FUNCTION DLADIV2( A, B, C, D, R, T );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      double             A, B, C, D, R, T;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
*
*     .. Local Scalars ..
      double             BR;
*     ..
*     .. Executable Statements ..
*
      IF( R.NE.ZERO ) THEN
         BR = B * R
         IF( BR.NE.ZERO ) THEN
            DLADIV2 = (A + BR) * T
         ELSE
            DLADIV2 = A * T + (B * T) * R
         END IF
      ELSE
         DLADIV2 = (A + D * (B / C)) * T
      END IF
*
      RETURN
*
*     End of DLADIV2
*
      END
