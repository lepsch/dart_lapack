      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, K1, K2, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * )
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, I1, I2, INC, IP, IX, IX0, J, K, N32;
      REAL               TEMP
      // ..
      // .. Executable Statements ..

      // Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
      // K1 through K2.

      if ( INCX.GT.0 ) {
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      } else if ( INCX.LT.0 ) {
         IX0 = K1 + ( K1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      } else {
         RETURN
      }

      N32 = ( N / 32 )*32
      if ( N32.NE.0 ) {
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               if ( IP.NE.I ) {
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               }
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      }
      if ( N32.NE.N ) {
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            if ( IP.NE.I ) {
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            }
            IX = IX + INCX
   50    CONTINUE
      }

      RETURN

      // End of SLASWP

      }
