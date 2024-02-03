      void dlartv(N, X, INCX, Y, INCY, C, S, INCC ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCC, INCX, INCY, N;
      // ..
      // .. Array Arguments ..
      double             C( * ), S( * ), X( * ), Y( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IC, IX, IY;
      double             XI, YI;
      // ..
      // .. Executable Statements ..

      IX = 1;
      IY = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         XI = X( IX );
         YI = Y( IY );
         X( IX ) = C( IC )*XI + S( IC )*YI;
         Y( IY ) = C( IC )*YI - S( IC )*XI;
         IX = IX + INCX;
         IY = IY + INCY;
         IC = IC + INCC;
      } // 10
      return;
      }
