      SUBROUTINE SLAMRG( N1, N2, A, STRD1, STRD2, INDEX )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N1, N2, STRD1, STRD2;
      // ..
      // .. Array Arguments ..
      int                INDEX( * );
      REAL               A( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, IND1, IND2, N1SV, N2SV;
      // ..
      // .. Executable Statements ..

      N1SV = N1
      N2SV = N2
      if ( STRD1.GT.0 ) {
         IND1 = 1
      } else {
         IND1 = N1
      }
      if ( STRD2.GT.0 ) {
         IND2 = 1 + N1
      } else {
         IND2 = N1 + N2
      }
      I = 1
      // while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      if ( N1SV.GT.0 .AND. N2SV.GT.0 ) {
         if ( A( IND1 ).LE.A( IND2 ) ) {
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + STRD1
            N1SV = N1SV - 1
         } else {
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + STRD2
            N2SV = N2SV - 1
         }
         GO TO 10
      }
      // end while
      if ( N1SV.EQ.0 ) {
         for (N1SV = 1; N1SV <= N2SV; N1SV++) { // 20
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + STRD2
   20    CONTINUE
      } else {
      // N2SV .EQ. 0
         for (N2SV = 1; N2SV <= N1SV; N2SV++) { // 30
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + STRD1
   30    CONTINUE
      }

      RETURN

      // End of SLAMRG

      }
