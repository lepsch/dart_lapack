      SUBROUTINE CLAR2V( N, X, Y, Z, INCX, C, S, INCC )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INCC, INCX, N;
      // ..
      // .. Array Arguments ..
      REAL               C( * )
      COMPLEX            S( * ), X( * ), Y( * ), Z( * )
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      int                I, IC, IX;
      REAL               CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII, ZIR
      COMPLEX            SI, T2, T3, T4, ZI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, CONJG, REAL
      // ..
      // .. Executable Statements ..
*
      IX = 1
      IC = 1
      DO 10 I = 1, N
         XI = REAL( X( IX ) )
         YI = REAL( Y( IX ) )
         ZI = Z( IX )
         ZIR = REAL( ZI )
         ZII = AIMAG( ZI )
         CI = C( IC )
         SI = S( IC )
         SIR = REAL( SI )
         SII = AIMAG( SI )
         T1R = SIR*ZIR - SII*ZII
         T1I = SIR*ZII + SII*ZIR
         T2 = CI*ZI
         T3 = T2 - CONJG( SI )*XI
         T4 = CONJG( T2 ) + SI*YI
         T5 = CI*XI + T1R
         T6 = CI*YI - T1R
         X( IX ) = CI*T5 + ( SIR*REAL( T4 )+SII*AIMAG( T4 ) )
         Y( IX ) = CI*T6 - ( SIR*REAL( T3 )-SII*AIMAG( T3 ) )
         Z( IX ) = CI*T3 + CONJG( SI )*CMPLX( T6, T1I )
         IX = IX + INCX
         IC = IC + INCC
   10 CONTINUE
      RETURN
*
      // End of CLAR2V
*
      END
