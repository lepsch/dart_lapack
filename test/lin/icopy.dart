      void icopy(N, SX, INCX, SY, final int INCY) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, INCY, N;
      int                SX( * ), SY( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IX, IY, M, MP1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD

      if (N <= 0) return;
      IF( INCX == 1 && INCY == 1 ) GO TO 20;

      // Code for unequal increments or equal increments not equal to 1

      IX = 1;
      IY = 1;
      if (INCX < 0) IX = ( -N+1 )*INCX + 1;
      IF( INCY < 0 ) IY = ( -N+1 )*INCY + 1;
      for (I = 1; I <= N; I++) { // 10
         SY[IY] = SX( IX );
         IX = IX + INCX;
         IY = IY + INCY;
      } // 10
      return;

      // Code for both increments equal to 1

      // Clean-up loop

      } // 20
      M = (N % 7);
      if (M == 0) GO TO 40;
      for (I = 1; I <= M; I++) { // 30
         SY[I] = SX( I );
      } // 30
      if (N < 7) return;
      } // 40
      MP1 = M + 1;
      for (I = MP1; I <= N; I += 7) { // 50
         SY[I] = SX( I );
         SY[I+1] = SX( I+1 );
         SY[I+2] = SX( I+2 );
         SY[I+3] = SX( I+3 );
         SY[I+4] = SX( I+4 );
         SY[I+5] = SX( I+5 );
         SY[I+6] = SX( I+6 );
      } // 50
      }
