import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlar2v(N, X, Y, Z, INCX, C, S, INCC ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCC, INCX, N;
      // ..
      // .. Array Arguments ..
      double             C( * ), S( * ), X( * ), Y( * ), Z( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IC, IX;
      double             CI, SI, T1, T2, T3, T4, T5, T6, XI, YI, ZI;
      // ..
      // .. Executable Statements ..

      IX = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         XI = X( IX );
         YI = Y( IX );
         ZI = Z( IX );
         CI = C( IC );
         SI = S( IC );
         T1 = SI*ZI;
         T2 = CI*ZI;
         T3 = T2 - SI*XI;
         T4 = T2 + SI*YI;
         T5 = CI*XI + T1;
         T6 = CI*YI - T1;
         X[IX] = CI*T5 + SI*T4;
         Y[IX] = CI*T6 - SI*T3;
         Z[IX] = CI*T4 - SI*T5;
         IX = IX + INCX;
         IC = IC + INCC;
      } // 10
      return;
      }
