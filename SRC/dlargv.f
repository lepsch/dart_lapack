      SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCC, INCX, INCY, N;
      // ..
      // .. Array Arguments ..
      double             C( * ), X( * ), Y( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IC, IX, IY;
      double             F, G, T, TT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      IX = 1
      IY = 1
      IC = 1
      DO 10 I = 1, N
         F = X( IX )
         G = Y( IY )
         if ( G.EQ.ZERO ) {
            C( IC ) = ONE
         } else if ( F.EQ.ZERO ) {
            C( IC ) = ZERO
            Y( IY ) = ONE
            X( IX ) = G
         } else if ( ABS( F ).GT.ABS( G ) ) {
            T = G / F
            TT = SQRT( ONE+T*T )
            C( IC ) = ONE / TT
            Y( IY ) = T*C( IC )
            X( IX ) = F*TT
         } else {
            T = F / G
            TT = SQRT( ONE+T*T )
            Y( IY ) = ONE / TT
            C( IC ) = T*Y( IY )
            X( IX ) = G*TT
         }
         IC = IC + INCC
         IY = IY + INCY
         IX = IX + INCX
   10 CONTINUE
      RETURN

      // End of DLARGV

      }
