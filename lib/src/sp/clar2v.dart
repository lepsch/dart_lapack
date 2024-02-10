      void clar2v(N, X, Y, Z, INCX, C, S, INCC ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCC, INCX, N;
      double               C( * );
      Complex            S( * ), X( * ), Y( * ), Z( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IC, IX;
      double               CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII, ZIR;
      Complex            SI, T2, T3, T4, ZI;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, CONJG, REAL

      IX = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         XI = double( X( IX ) );
         YI = double( Y( IX ) );
         ZI = Z( IX );
         ZIR = double( ZI );
         ZII = AIMAG( ZI );
         CI = C( IC );
         SI = S( IC );
         SIR = double( SI );
         SII = AIMAG( SI );
         T1R = SIR*ZIR - SII*ZII;
         T1I = SIR*ZII + SII*ZIR;
         T2 = CI*ZI;
         T3 = T2 - CONJG( SI )*XI;
         T4 = CONJG( T2 ) + SI*YI;
         T5 = CI*XI + T1R;
         T6 = CI*YI - T1R;
         X[IX] = CI*T5 + ( SIR*double( T4 )+SII*AIMAG( T4 ) );
         Y[IX] = CI*T6 - ( SIR*double( T3 )-SII*AIMAG( T3 ) );
         Z[IX] = CI*T3 + CONJG( SI )*CMPLX( T6, T1I );
         IX = IX + INCX;
         IC = IC + INCC;
      } // 10
      }
