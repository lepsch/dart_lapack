      void clartv(final int N, final int X, final int INCX, final int Y, final int INCY, final int C, final int S, final int INCC,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCC, INCX, INCY, N;
      double               C( * );
      Complex            S( * ), X( * ), Y( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IC, IX, IY;
      Complex            XI, YI;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG

      IX = 1;
      IY = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         XI = X( IX );
         YI = Y( IY );
         X[IX] = C( IC )*XI + S( IC )*YI;
         Y[IY] = C( IC )*YI - CONJG( S( IC ) )*XI;
         IX = IX + INCX;
         IY = IY + INCY;
         IC = IC + INCC;
      } // 10
      }
