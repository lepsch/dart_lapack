      void zlar2v(final int N, final int X, final int Y, final int Z, final int INCX, final int C, final int S, final int INCC) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCC, INCX, N;
      double             C( * );
      Complex         S( * ), X( * ), Y( * ), Z( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IC, IX;
      double             CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII, ZIR;
      Complex         SI, T2, T3, T4, ZI;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DCONJG, DIMAG

      IX = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         XI = (X( IX )).toDouble();
         YI = (Y( IX )).toDouble();
         ZI = Z( IX );
         ZIR = ZI.toDouble();
         ZII = DIMAG( ZI );
         CI = C( IC );
         SI = S( IC );
         SIR = SI.toDouble();
         SII = DIMAG( SI );
         T1R = SIR*ZIR - SII*ZII;
         T1I = SIR*ZII + SII*ZIR;
         T2 = CI*ZI;
         T3 = T2 - DCONJG( SI )*XI;
         T4 = DCONJG( T2 ) + SI*YI;
         T5 = CI*XI + T1R;
         T6 = CI*YI - T1R;
         X[IX] = CI*T5 + ( SIR*T4.toDouble()+SII*DIMAG( T4 ) );
         Y[IX] = CI*T6 - ( SIR*T3.toDouble()-SII*DIMAG( T3 ) );
         Z[IX] = CI*T3 + DCONJG( SI )*DCMPLX( T6, T1I );
         IX = IX + INCX;
         IC = IC + INCC;
      } // 10
      }
