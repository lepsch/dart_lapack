      SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC );

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
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IC, IX, IY;
      double             F, G, T, TT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      IX = 1;
      IY = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         F = X( IX );
         G = Y( IY );
         if ( G == ZERO ) {
            C( IC ) = ONE;
         } else if ( F == ZERO ) {
            C( IC ) = ZERO;
            Y( IY ) = ONE;
            X( IX ) = G;
         } else if ( ABS( F ) > ABS( G ) ) {
            T = G / F;
            TT = SQRT( ONE+T*T );
            C( IC ) = ONE / TT;
            Y( IY ) = T*C( IC );
            X( IX ) = F*TT;
         } else {
            T = F / G;
            TT = SQRT( ONE+T*T );
            Y( IY ) = ONE / TT;
            C( IC ) = T*Y( IY );
            X( IX ) = G*TT;
         }
         IC = IC + INCC;
         IY = IY + INCY;
         IX = IX + INCX;
      } // 10
      RETURN;

      // End of DLARGV

      }
