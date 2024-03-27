import 'dart:io';
import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/drotg.dart';
import 'package:lapack/src/blas/drotm.dart';
import 'package:lapack/src/blas/drotmg.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsdot.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/epsilon.dart';
import 'package:lapack/src/intrinsics/huge.dart';
import 'package:lapack/src/intrinsics/random_number.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';

import '../test_driver.dart';
import 'common.dart';

void dblat1(final Nout nout, final TestDriver test) {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const SFAC = 9.765625e-4;

  nout.println(' Real BLAS Test Program Results');
  test.group('Double Precision Level 1 BLAS routines', () {
    for (final IC in 1.through(13)) {
      test(L[IC - 1].trim(), () {
        combla.ICASE = IC;
        _printHeader(nout);

        // .. Initialize  PASS,  INCX,  and INCY for a new case. ..
        // .. the value 9999 for INCX or INCY will appear in the ..
        // .. detailed  output, if any, for cases  that do not involve ..
        // .. these parameters ..

        combla.PASS = true;
        combla.INCX = 9999;
        combla.INCY = 9999;
        if (combla.ICASE == 3 || combla.ICASE == 11) {
          _check0(SFAC, nout);
        } else if (combla.ICASE == 7 ||
            combla.ICASE == 8 ||
            combla.ICASE == 9 ||
            combla.ICASE == 10) {
          _check1(SFAC, nout);
        } else if (combla.ICASE == 1 ||
            combla.ICASE == 2 ||
            combla.ICASE == 5 ||
            combla.ICASE == 6 ||
            combla.ICASE == 12 ||
            combla.ICASE == 13) {
          _check2(SFAC, nout);
        } else if (combla.ICASE == 4) {
          _check3(SFAC, nout);
        }

        if (test.expect(combla.PASS, true)) {
          nout.println('                                    ----- PASS -----');
        }
      });
    }
  });
}

const L = [
  ' DDOT ',
  'DAXPY ',
  'DROTG ',
  ' DROT ',
  'DCOPY ',
  'DSWAP ',
  'DNRM2 ',
  'DASUM ',
  'DSCAL ',
  'IDAMAX',
  'DROTMG',
  'DROTM ',
  'DSDOT ',
];

void _printHeader(final Nout nout) {
  nout.println(
      '\n Test of subprogram number${combla.ICASE.i3}${' ' * 12}${L[combla.ICASE - 1].a6}');
}

void _check0(final double SFAC, final Nout nout) {
  int I;
  final SA = Box(0.0), SB = Box(0.0), SC = Box(0.0), SS = Box(0.0);
  final DTEMP = Array<double>(9);
  final DA1 = Array.fromList([0.3, 0.4, -0.3, -0.4, -0.3, 0.0, 0.0, 1.0]);
  final DB1 = Array.fromList([0.4, 0.3, 0.4, 0.3, -0.4, 0.0, 1.0, 0.0]);
  final DC1 = Array.fromList([0.6, 0.8, -0.6, 0.8, 0.6, 1.0, 0.0, 1.0]);
  final DS1 = Array.fromList([0.8, 0.6, 0.8, -0.6, 0.8, 0.0, 1.0, 0.0]);
  final DATRUE = Array.fromList([0.5, 0.5, 0.5, -0.5, -0.5, 0.0, 1.0, 1.0]);
  final DBTRUE = Array.fromList([0.0, 0.6, 0.0, -0.6, 0.0, 0.0, 1.0, 0.0]);
  // INPUT FOR MODIFIED GIVENS
  final DAB = Matrix.fromData(
      [
        [.1, .3, 1.2, .2, .7, .2, .6, 4.2, 0.0],
        [0.0, 0.0, 0.0, 4.0, -1.0, 2.0, 4.0, 6e-10, 2e-2],
        [1e5, 10.0, 4e10, 2e-2, 1e-5, 10.0, 2e-10, 4e-2, 1e5],
        [10.0, 2e10, 4e-2, 1e-5, 10.0, 4.0, -2.0, 8.0, 4.0],
      ].flattened.toList(),
      (4, 9));
  // TRUE RESULTS FOR MODIFIED GIVENS
  final DTRUE = Matrix.fromData(
      [
        [0.0, 0.0, 1.3, 0.2, 0.0, 0.0, 0.0, 0.5, 0.0],
        [0.0, 0.0, 4.5, 4.2, 1.0, 0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 4.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 15e-3, 0.0, 10.0, -1.0, 0.0, -1e-4, 0.0, 1.0],
        [0.0, 0.0, 6144e-5, 10.0, -1.0, 4096.0, -1e6, 0.0, 1.0],
        [0.0, 0.0, 15.0, 10.0, -1.0, 5e-5, 0.0, 1.0, 0.0],
        [0.0, 0.0, 15.0, 10.0, -1.0, 5e5, -4096.0, 1.0, 4096e-6],
        [0.0, 0.0, 7.0, 4.0, 0.0, 0.0, -.5, -.25, 0.0],
      ].flattened.toList(),
      (9, 9));

  const D12 = 4096.0; // 4096 = 2 ** 12
  DTRUE[1][1] = 12.0 / 130.0;
  DTRUE[2][1] = 36.0 / 130.0;
  DTRUE[7][1] = -1.0 / 6.0;
  DTRUE[1][2] = 14.0 / 75.0;
  DTRUE[2][2] = 49.0 / 75.0;
  DTRUE[9][2] = 1.0 / 7.0;
  DTRUE[1][5] = 45e-11 * (D12 * D12);
  DTRUE[3][5] = 4e5 / (3.0 * D12);
  DTRUE[6][5] = 1.0 / D12;
  DTRUE[8][5] = 1e4 / (3.0 * D12);
  DTRUE[1][6] = 4e10 / (1.5 * D12 * D12);
  DTRUE[2][6] = 2e-2 / 1.5;
  DTRUE[8][6] = 5e-7 * D12;
  DTRUE[1][7] = 4.0 / 150.0;
  DTRUE[2][7] = (2e-10 / 1.5) * (D12 * D12);
  DTRUE[7][7] = -DTRUE[6][5];
  DTRUE[9][7] = 1e4 / D12;
  DTRUE[1][8] = DTRUE[1][7];
  DTRUE[2][8] = 2e10 / (1.5 * D12 * D12);
  DTRUE[1][9] = 32.0 / 7.0;
  DTRUE[2][9] = -16.0 / 7.0;

  // Compute true values which cannot be prestored
  // in decimal notation

  DBTRUE[1] = 1.0 / 0.6;
  DBTRUE[3] = -1.0 / 0.6;
  DBTRUE[5] = 1.0 / 0.6;

  for (var K = 1; K <= 8; K++) {
    // .. Set combla.N=K for identification in output if any ..
    combla.N = K;
    if (combla.ICASE == 3) {
      // .. DROTG ..
      if (K > 8) return;
      SA.value = DA1[K];
      SB.value = DB1[K];
      drotg(SA, SB, SC, SS);
      _stest1(SA.value, DATRUE[K], DATRUE(K), SFAC, nout);
      _stest1(SB.value, DBTRUE[K], DBTRUE(K), SFAC, nout);
      _stest1(SC.value, DC1[K], DC1(K), SFAC, nout);
      _stest1(SS.value, DS1[K], DS1(K), SFAC, nout);
    } else if (combla.ICASE == 11) {
      // .. DROTMG ..
      for (I = 1; I <= 4; I++) {
        DTEMP[I] = DAB[I][K];
        DTEMP[I + 4] = 0.0;
      }
      DTEMP[9] = 0.0;
      drotmg(DTEMP(1), DTEMP(2), DTEMP(3), DTEMP[4], DTEMP(5));
      _stest(
          9, DTEMP, DTRUE(1, K).asArray(), DTRUE(1, K).asArray(), SFAC, nout);
    } else {
      nout.println(' Shouldn\'t be here in CHECK0');
    }
  }
}

void _check1(final double SFAC, final Nout nout) {
  const THRESH = 10.0;
  int I, IX, LEN, NP1;
  final STEMP = Array<double>(1),
      STRUE = Array<double>(8),
      SX = Array<double>(8),
      SXR = Array<double>(15);
  final SA =
      Array.fromList([0.3, -1.0, 0.0, 1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]);
  final DV = Matrix3d.fromData(
      [
        [0.1, 2.0],
        [2.0, 2.0],
        [2.0, 2.0],
        [2.0, 2.0],
        [0.3, 3.0],
        [3.0, 3.0],
        [3.0, 3.0],
        [3.0, 3.0],
        [0.3, -0.4],
        [4.0, 4.0],
        [4.0, 4.0],
        [4.0, 4.0],
        [0.2, -0.6],
        [0.3, 5.0],
        [5.0, 5.0],
        [5.0, 5.0],
        [0.1, -0.3],
        [0.5, -0.1],
        [6.0, 6.0],
        [6.0, 6.0],
        [0.1, 8.0],
        [8.0, 8.0],
        [8.0, 8.0],
        [8.0, 8.0],
        [0.3, 9.0],
        [9.0, 9.0],
        [9.0, 9.0],
        [9.0, 9.0],
        [0.3, 2.0],
        [-0.4, 2.0],
        [2.0, 2.0],
        [2.0, 2.0],
        [0.2, 3.0],
        [-0.6, 5.0],
        [0.3, 2.0],
        [2.0, 2.0],
        [0.1, 4.0],
        [-0.3, 6.0],
        [-0.5, 7.0],
        [-0.1, 3.0],
      ].flattened.toList(),
      (8, 5, 2));
  final DVR = Array.fromList([8.0, -7.0, 9.0, 5.0, 9.0, 8.0, 7.0, 7.0]);
  final DTRUE1 = Array.fromList([0.0, 0.3, 0.5, 0.7, 0.6]);
  final DTRUE3 = Array.fromList([0.0, 0.3, 0.7, 1.1, 1.0]);
  final DTRUE5 = Matrix3d.fromData(
      [
        [0.10, 2.0],
        [2.0, 2.0],
        [2.0, 2.0],
        [2.0, 2.0],
        [-0.3, 3.0],
        [3.0, 3.0],
        [3.0, 3.0],
        [3.0, 3.0],
        [0.0, 0.0],
        [4.0, 4.0],
        [4.0, 4.0],
        [4.0, 4.0],
        [0.20, -0.60],
        [0.30, 5.0],
        [5.0, 5.0],
        [5.0, 5.0],
        [0.03, -0.09],
        [0.15, -0.03],
        [6.0, 6.0],
        [6.0, 6.0],
        [0.10, 8.0],
        [8.0, 8.0],
        [8.0, 8.0],
        [8.0, 8.0],
        [0.09, 9.0],
        [9.0, 9.0],
        [9.0, 9.0],
        [9.0, 9.0],
        [0.09, 2.0],
        [-0.12, 2.0],
        [2.0, 2.0],
        [2.0, 2.0],
        [0.06, 3.0],
        [-0.18, 5.0],
        [0.09, 2.0],
        [2.0, 2.0],
        [0.03, 4.0],
        [-0.09, 6.0],
        [-0.15, 7.0],
        [-0.03, 3.0],
      ].flattened.toList(),
      (8, 5, 2));
  final ITRUE2 = Array.fromList([0, 1, 2, 2, 3]);
  final ITRUEC = Array.fromList([0, 1, 1, 1, 1]);

  for (combla.INCX = 1; combla.INCX <= 2; combla.INCX++) {
    for (NP1 = 1; NP1 <= 5; NP1++) {
      combla.N = NP1 - 1;
      LEN = (2 * max(combla.N, 1));
      // .. Set vector arguments ..
      for (I = 1; I <= LEN; I++) {
        SX[I] = DV[I][NP1][combla.INCX];
      }

      if (combla.ICASE == 7) {
        // .. DNRM2 ..
        // Test scaling when some entries are tiny or huge
        _db1nrm2(combla.N, (combla.INCX - 2) * 2, THRESH, nout);
        _db1nrm2(combla.N, combla.INCX, THRESH, nout);
        // Test with hardcoded mid range entries
        STEMP[1] = DTRUE1[NP1];
        _stest1(dnrm2(combla.N, SX, combla.INCX), STEMP[1], STEMP, SFAC, nout);
      } else if (combla.ICASE == 8) {
        // .. DASUM ..
        STEMP[1] = DTRUE3[NP1];
        _stest1(dasum(combla.N, SX, combla.INCX), STEMP[1], STEMP, SFAC, nout);
      } else if (combla.ICASE == 9) {
        // .. DSCAL ..
        dscal(combla.N, SA[(combla.INCX - 1) * 5 + NP1], SX, combla.INCX);
        for (I = 1; I <= LEN; I++) {
          STRUE[I] = DTRUE5[I][NP1][combla.INCX];
        }
        _stest(LEN, SX, STRUE, STRUE, SFAC, nout);
      } else if (combla.ICASE == 10) {
        // .. idamax ..
        _itest1(idamax(combla.N, SX, combla.INCX), ITRUE2[NP1], nout);
        for (I = 1; I <= LEN; I++) {
          SX[I] = 42.0;
        }
        _itest1(idamax(combla.N, SX, combla.INCX), ITRUEC[NP1], nout);
      } else {
        throw UnimplementedError(' Shouldn\'t be here in CHECK1');
      }
    }

    if (combla.ICASE == 10) {
      combla.N = 8;
      IX = 1;
      for (I = 1; I <= combla.N; I++) {
        SXR[IX] = DVR[I];
        IX += combla.INCX;
      }
      _itest1(idamax(combla.N, SXR, combla.INCX), 3, nout);
    }
  }
}

void _check2(final double SFAC, final Nout nout) {
  int I, J, KI, KN, KNI, KPAR, KSIZE, LENX, LENY, LINCX, LINCY, MX, MY;
  final SSIZE = Array<double>(7),
      STX = Array<double>(7),
      STY = Array<double>(7),
      SX = Array<double>(7),
      SY = Array<double>(7),
      DTEMP = Array<double>(5),
      STY0 = Array<double>(1),
      SX0 = Array<double>(1),
      SY0 = Array<double>(1);

  const SA = 0.3;
  final INCXS = Array.fromList([1, 2, -2, -1]);
  final INCYS = Array.fromList([1, -2, 1, -2]);
  final LENS = Matrix.fromData([
    1, 1, 2, 4, 1, 1, 3, 7, //
  ], (
    4, 2 //
  ));
  final NS = Array.fromList([0, 1, 2, 4]);

  final DX1 = Array.fromList([0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4]);
  final DY1 = Array.fromList([0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8]);
  final DT7 = Matrix.fromData([
    0.0, 0.30, 0.21, 0.62, 0.0, 0.30, -0.07, 0.85, //
    0.0, 0.30, -0.79, -0.74, 0.0, 0.30, 0.33, 1.27, //
  ], (
    4, 4 //
  ));
  final DT8 = Matrix3d.fromData(
      [
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.68],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.68, -0.87],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.68, -0.87, 0.15],
        [0.94, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.68],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.35, -0.9],
        [0.48, 0.0, 0.0, 0.0],
        [0.0, 0.38, -0.9, 0.57],
        [0.7, -0.75, 0.2, 0.98],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.68],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.35, -0.72],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.38, -0.63, 0.15],
        [0.88, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.68],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.68, -0.9],
        [0.33, 0.0, 0.0, 0.0],
        [0.0, 0.68, -0.9, 0.33],
        [0.7, -0.75, 0.2, 1.04],
      ].flattened.toList(),
      (7, 4, 4));
  final DT10X = Matrix3d.fromData(
      [
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5, -0.9],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.5, -0.9, 0.3],
        [0.7, 0.0, 0.0, 0.0],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.3, 0.1],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.8, 0.1, -0.6],
        [0.8, 0.3, -0.3, 0.5],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -0.9, 0.1],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.7, 0.1, 0.3],
        [0.8, -0.9, -0.3, 0.5],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5, 0.3],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.5, 0.3, -0.6],
        [0.8, 0.0, 0.0, 0.0]
      ].flattened.toList(),
      (7, 4, 4));
  final DT10Y = Matrix3d.fromData(
      [
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.6],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.6, 0.1],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.6, 0.1, -0.5],
        [0.8, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.6],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -0.5, -0.9],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, -0.4, -0.9, 0.9],
        [0.7, -0.5, 0.2, 0.6],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.6],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -0.5, 0.6],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, -0.4, 0.9, -0.5],
        [0.6, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.6],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.6, -0.9],
        [0.1, 0.0, 0.0, 0.0],
        [0.0, 0.6, -0.9, 0.1],
        [0.7, -0.5, 0.2, 0.8],
      ].flattened.toList(),
      (7, 4, 4));
  final SSIZE1 = Array.fromList([0.0, 0.3, 1.6, 3.2]);
  final SSIZE2 = Matrix.fromData(
      [
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [1.17, 1.17],
        [1.17, 1.17],
        [1.17, 1.17],
        [1.17, 1.17],
        [1.17, 1.17],
        [1.17, 1.17],
        [1.17, 1.17],
      ].flattened.toList(),
      (14, 2));

  // FOR DROTM

  final DPAR = Matrix.fromData([
    -2.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, -3.0, -4.0, 5.0, //
    0.0, 0.0, 2.0, -3.0, 0.0, 1.0, 5.0, 2.0, 0.0, -4.0,
  ], (
    5, 4 //
  ));
  // TRUE X RESULTS F0R ROTATIONS DROTM
  final DT19X = Matrix3d.fromData(
      [
        [
          .6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, -.8, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, -.9, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .6, .1, 0.0, 0.0, 0.0, 0.0, 0.0, -.8
        ],
        [
          3.8, 0.0, 0.0, 0.0, 0.0, 0.0, -.9, 2.8, //
          0.0, 0.0, 0.0, 0.0, 0.0, 3.5, -.4, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, .6, .1, -.5, .8, //
          0.0, 0.0, 0.0, -.8, 3.8, -2.2, -1.2, 0.0
        ],
        [
          0.0, 0.0, -.9, 2.8, -1.4, -1.3, 0.0, 0.0, //
          0.0, 3.5, -.4, -2.2, 4.7, 0.0, 0.0, 0.0
        ],
        [
          .6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, -.8, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, -.9, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .6, .1, -.5, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          .1, -3.0, 0.0, 0.0, 0.0, 0.0, -.3, .1, //
          -2.0, 0.0, 0.0, 0.0, 0.0, 3.3, .1, -2.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, .6, .1, -.5, .8, //
          .9, -.3, -.4, -2.0, .1, 1.4, .8, .6
        ],
        [
          -.3, -2.8, -1.8, .1, 1.3, .8, 0.0, -.3, //
          -1.9, 3.8, .1, -3.1, .8, 4.8, -.3, -1.5
        ],
        [
          .6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, -.8, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, -.9, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .6, .1, -.5, 0.0, 0.0, 0.0, 0.0, 4.8
        ],
        [
          .1, -3.0, 0.0, 0.0, 0.0, 0.0, 3.3, .1, //
          -2.0, 0.0, 0.0, 0.0, 0.0, 2.1, .1, -2.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, .6, .1, -.5, .8, //
          .9, -.3, -.4, -1.6, .1, -2.2, .8, 5.4
        ],
        [
          -.3, -2.8, -1.5, .1, -1.4, .8, 3.6, -.3, //
          -1.9, 3.7, .1, -2.2, .8, 3.6, -.3, -1.5
        ],
        [
          .6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .6, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, -.8, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, -.9, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .6, .1, 0.0, 0.0, 0.0, 0.0, 0.0, -.8
        ],
        [
          -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -.9, -.8, //
          0.0, 0.0, 0.0, 0.0, 0.0, 3.5, .8, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, .6, .1, -.5, .8, //
          0.0, 0.0, 0.0, -.8, -1.0, 1.4, -1.6, 0.0
        ],
        [
          0.0, 0.0, -.9, -.8, 1.3, -1.6, 0.0, 0.0, //
          0.0, 3.5, .8, -3.1, 4.8, 0.0, 0.0, 0.0
        ],
      ].flattened.toList(),
      (7, 4, 16));
  // TRUE Y RESULTS FOR ROTATIONS DROTM
  final DT19Y = Matrix3d.fromData(
      [
        [
          .5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, .7, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, -2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .5, -.9, 0.0, 0.0, 0.0, 0.0, 0.0, .7
        ],
        [
          -4.8, 0.0, 0.0, 0.0, 0.0, 0.0, 1.7, -.7, //
          0.0, 0.0, 0.0, 0.0, 0.0, -2.6, 3.5, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, .5, -.9, .3, .7, //
          0.0, 0.0, 0.0, .7, -4.8, 3.0, 1.1, 0.0
        ],
        [
          0.0, 0.0, 1.7, -.7, -.7, 2.3, 0.0, 0.0, //
          0.0, -2.6, 3.5, -.7, -3.6, 0.0, 0.0, 0.0
        ],
        [
          .5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, .7, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, -2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .5, -.9, .3, 0.0, 0.0, 0.0, 0.0, 4.0
        ],
        [
          -.9, -.3, 0.0, 0.0, 0.0, 0.0, -.5, -.9, //
          1.5, 0.0, 0.0, 0.0, 0.0, -1.5, -.9, -1.8
        ],
        [
          0.0, 0.0, 0.0, 0.0, .5, -.9, .3, .7, //
          -.6, .2, .8, 3.7, -.9, -1.2, .7, -1.5
        ],
        [
          .2, 2.2, -.3, -.9, 2.1, .7, -1.6, .2, //
          2.0, -1.6, -.9, -2.1, .7, 2.9, .2, -3.8
        ],
        [
          .5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, .7, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, -2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .5, -.9, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0
        ],
        [
          -6.3, 0.0, 0.0, 0.0, 0.0, 0.0, -.5, .3, //
          0.0, 0.0, 0.0, 0.0, 0.0, -1.5, 3.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, .5, -.9, .3, .7, //
          0.0, 0.0, 0.0, 3.7, -7.2, 3.0, 1.7, 0.0
        ],
        [
          0.0, 0.0, -.3, .9, -.7, 1.9, 0.0, 0.0, //
          0.0, -1.6, 2.7, -.7, -3.4, 0.0, 0.0, 0.0
        ],
        [
          .5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, //
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0
        ],
        [
          0.0, 0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, //
          0.0, 0.0, 0.0, 0.0, .5, 0.0, 0.0, 0.0
        ],
        [
          0.0, 0.0, 0.0, .7, 0.0, 0.0, 0.0, 0.0, //
          0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 0.0, 0.0
        ],
        [
          0.0, -2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
          .5, -.9, .3, 0.0, 0.0, 0.0, 0.0, .7
        ],
        [
          -.9, 1.2, 0.0, 0.0, 0.0, 0.0, 1.7, -.9, //
          .5, 0.0, 0.0, 0.0, 0.0, -2.6, -.9, -1.3
        ],
        [
          0.0, 0.0, 0.0, 0.0, .5, -.9, .3, .7, //
          -.6, .2, .8, .7, -.9, 1.2, .7, -1.5
        ],
        [
          .2, 1.6, 1.7, -.9, .5, .7, -1.6, .2, //
          2.4, -2.6, -.9, -1.3, .7, 2.9, .2, -4.0
        ],
      ].flattened.toList(),
      (7, 4, 16));

  for (KI = 1; KI <= 4; KI++) {
    combla.INCX = INCXS[KI];
    combla.INCY = INCYS[KI];
    MX = combla.INCX.abs();
    MY = combla.INCY.abs();

    for (KN = 1; KN <= 4; KN++) {
      combla.N = NS[KN];
      KSIZE = min(2, KN);
      LENX = LENS[KN][MX];
      LENY = LENS[KN][MY];
      // .. Initialize all argument arrays ..
      for (I = 1; I <= 7; I++) {
        SX[I] = DX1[I];
        SY[I] = DY1[I];
      }

      if (combla.ICASE == 1) {
        // .. DDOT ..
        _stest1(ddot(combla.N, SX, combla.INCX, SY, combla.INCY), DT7[KN][KI],
            SSIZE1(KN), SFAC, nout);
      } else if (combla.ICASE == 2) {
        // .. DAXPY ..
        daxpy(combla.N, SA, SX, combla.INCX, SY, combla.INCY);
        for (J = 1; J <= LENY; J++) {
          STY[J] = DT8[J][KN][KI];
        }
        _stest(LENY, SY, STY, SSIZE2(1, KSIZE).asArray(), SFAC, nout);
      } else if (combla.ICASE == 5) {
        // .. DCOPY ..
        for (I = 1; I <= 7; I++) {
          STY[I] = DT10Y[I][KN][KI];
        }
        dcopy(combla.N, SX, combla.INCX, SY, combla.INCY);
        _stest(LENY, SY, STY, SSIZE2(1, 1).asArray(), 1.0, nout);
        if (KI == 1) {
          SX0[1] = 42.0;
          SY0[1] = 43.0;
          if (combla.N == 0) {
            STY0[1] = SY0[1];
          } else {
            STY0[1] = SX0[1];
          }
          LINCX = combla.INCX;
          combla.INCX = 0;
          LINCY = combla.INCY;
          combla.INCY = 0;
          dcopy(combla.N, SX0, combla.INCX, SY0, combla.INCY);
          _stest(1, SY0, STY0, SSIZE2(1, 1).asArray(), 1.0, nout);
          combla.INCX = LINCX;
          combla.INCY = LINCY;
        }
      } else if (combla.ICASE == 6) {
        // .. DSWAP ..
        dswap(combla.N, SX, combla.INCX, SY, combla.INCY);
        for (I = 1; I <= 7; I++) {
          STX[I] = DT10X[I][KN][KI];
          STY[I] = DT10Y[I][KN][KI];
        }
        _stest(LENX, SX, STX, SSIZE2(1, 1).asArray(), 1.0, nout);
        _stest(LENY, SY, STY, SSIZE2(1, 1).asArray(), 1.0, nout);
      } else if (combla.ICASE == 12) {
        // .. DROTM ..
        KNI = KN + 4 * (KI - 1);
        for (KPAR = 1; KPAR <= 4; KPAR++) {
          for (I = 1; I <= 7; I++) {
            SX[I] = DX1[I];
            SY[I] = DY1[I];
            STX[I] = DT19X[I][KPAR][KNI];
            STY[I] = DT19Y[I][KPAR][KNI];
          }

          for (I = 1; I <= 5; I++) {
            DTEMP[I] = DPAR[I][KPAR];
          }

          for (I = 1; I <= LENX; I++) {
            SSIZE[I] = STX[I];
          }
          // SEE REMARK ABOVE ABOUT DT11X(1,2,7)
          // AND DT11X(5,3,8).
          if ((KPAR == 2) && (KNI == 7)) SSIZE[1] = 2.4;
          if ((KPAR == 3) && (KNI == 8)) SSIZE[5] = 1.8;

          drotm(combla.N, SX, combla.INCX, SY, combla.INCY, DTEMP);
          _stest(LENX, SX, STX, SSIZE, SFAC, nout);
          _stest(LENY, SY, STY, STY, SFAC, nout);
        }
      } else if (combla.ICASE == 13) {
        // .. DSDOT ..
        _testdsdot(dsdot(combla.N, SX, combla.INCX, SY, combla.INCY),
            DT7[KN][KI], SSIZE1[KN], 0.3125E-1, nout);
      } else {
        throw UnsupportedError(' Shouldn' 't be here in CHECK2');
      }
    }
  }
}

void _check3(final double SFAC, final Nout nout) {
  int I, K, KI, KN, KSIZE, LENX, LENY, MX, MY;
  final MWPTX = Matrix<double>(11, 5), MWPTY = Matrix<double>(11, 5);
  final COPYX = Array<double>(5),
      COPYY = Array<double>(5),
      MWPC = Array<double>(11),
      MWPS = Array<double>(11),
      MWPSTX = Array<double>(5),
      MWPSTY = Array<double>(5),
      MWPX = Array<double>(5),
      MWPY = Array<double>(5),
      STX = Array<double>(7),
      STY = Array<double>(7),
      SX = Array<double>(7),
      SY = Array<double>(7);
  final MWPINX = Array<int>(11), MWPINY = Array<int>(11), MWPN = Array<int>(11);

  final INCXS = Array.fromList([1, 2, -2, -1]);
  final INCYS = Array.fromList([1, -2, 1, -2]);
  final LENS = Matrix.fromData([1, 1, 2, 4, 1, 1, 3, 7], (4, 2));
  final NS = Array.fromList([0, 1, 2, 4]);
  final DX1 = Array.fromList([0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4]);
  final DY1 = Array.fromList([0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8]);
  final (SC, SS) = (0.8, 0.6);
  final DT9X = Matrix3d.fromData(
      [
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.78],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.78, -0.46],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.78, -0.46, -0.22],
        [1.06, 0.0, 0.0, 0.0],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.78],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.66, 0.1],
        [-0.1, 0.0, 0.0, 0.0],
        [0.0, 0.96, 0.1, -0.76],
        [0.8, 0.90, -0.3, -0.02],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.78],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -0.06, 0.1],
        [-0.1, 0.0, 0.0, 0.0],
        [0.0, 0.90, 0.1, -0.22],
        [0.8, 0.18, -0.3, -0.02],
        [0.6, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.78],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.78, 0.26],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.78, 0.26, -0.76],
        [1.12, 0.0, 0.0, 0.0],
      ].flattened.toList(),
      (7, 4, 4));
  final DT9Y = Matrix3d.fromData(
      [
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.04],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.04, -0.78],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.04, -0.78, 0.54],
        [0.08, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.04],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.7, -0.9],
        [-0.12, 0.0, 0.0, 0.0],
        [0.0, 0.64, -0.9, -0.30],
        [0.7, -0.18, 0.2, 0.28],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.04],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.7, -1.08],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.64, -1.26, 0.54],
        [0.20, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.04],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.04, -0.9],
        [0.18, 0.0, 0.0, 0.0],
        [0.0, 0.04, -0.9, 0.18],
        [0.7, -0.18, 0.2, 0.16]
      ].flattened.toList(),
      (7, 4, 4));
  final SSIZE2 = Matrix.fromData([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
    1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, //
    1.17, 1.17, 1.17, 1.17, 1.17, 1.17, //
  ], (
    14, 2 //
  ));

  for (KI = 1; KI <= 4; KI++) {
    combla.INCX = INCXS[KI];
    combla.INCY = INCYS[KI];
    MX = combla.INCX.abs();
    MY = combla.INCY.abs();

    for (KN = 1; KN <= 4; KN++) {
      combla.N = NS[KN];
      KSIZE = min(2, KN);
      LENX = LENS[KN][MX];
      LENY = LENS[KN][MY];

      if (combla.ICASE == 4) {
        // .. DROT ..
        for (I = 1; I <= 7; I++) {
          SX[I] = DX1[I];
          SY[I] = DY1[I];
          STX[I] = DT9X[I][KN][KI];
          STY[I] = DT9Y[I][KN][KI];
        }
        drot(combla.N, SX, combla.INCX, SY, combla.INCY, SC, SS);
        _stest(LENX, SX, STX, SSIZE2(1, KSIZE).asArray(), SFAC, nout);
        _stest(LENY, SY, STY, SSIZE2(1, KSIZE).asArray(), SFAC, nout);
      } else {
        throw UnsupportedError(' Shouldn\'t be here in CHECK3');
      }
    }
  }

  MWPC[1] = 1;
  for (I = 2; I <= 11; I++) {
    MWPC[I] = 0;
  }
  MWPS[1] = 0;
  for (I = 2; I <= 6; I++) {
    MWPS[I] = 1;
  }
  for (I = 7; I <= 11; I++) {
    MWPS[I] = -1;
  }
  MWPINX[1] = 1;
  MWPINX[2] = 1;
  MWPINX[3] = 1;
  MWPINX[4] = -1;
  MWPINX[5] = 1;
  MWPINX[6] = -1;
  MWPINX[7] = 1;
  MWPINX[8] = 1;
  MWPINX[9] = -1;
  MWPINX[10] = 1;
  MWPINX[11] = -1;
  MWPINY[1] = 1;
  MWPINY[2] = 1;
  MWPINY[3] = -1;
  MWPINY[4] = -1;
  MWPINY[5] = 2;
  MWPINY[6] = 1;
  MWPINY[7] = 1;
  MWPINY[8] = -1;
  MWPINY[9] = -1;
  MWPINY[10] = 2;
  MWPINY[11] = 1;
  for (I = 1; I <= 11; I++) {
    MWPN[I] = 5;
  }
  MWPN[5] = 3;
  MWPN[10] = 3;
  for (I = 1; I <= 5; I++) {
    MWPX[I] = I.toDouble();
    MWPY[I] = I.toDouble();
    MWPTX[1][I] = I.toDouble();
    MWPTY[1][I] = I.toDouble();
    MWPTX[2][I] = I.toDouble();
    MWPTY[2][I] = -I.toDouble();
    MWPTX[3][I] = (6 - I).toDouble();
    MWPTY[3][I] = I - 6.toDouble();
    MWPTX[4][I] = I.toDouble();
    MWPTY[4][I] = -I.toDouble();
    MWPTX[6][I] = (6 - I).toDouble();
    MWPTY[6][I] = I - 6.toDouble();
    MWPTX[7][I] = -I.toDouble();
    MWPTY[7][I] = I.toDouble();
    MWPTX[8][I] = I - 6.toDouble();
    MWPTY[8][I] = (6 - I).toDouble();
    MWPTX[9][I] = -I.toDouble();
    MWPTY[9][I] = I.toDouble();
    MWPTX[11][I] = I - 6.toDouble();
    MWPTY[11][I] = (6 - I).toDouble();
  }
  MWPTX[5][1] = 1;
  MWPTX[5][2] = 3;
  MWPTX[5][3] = 5;
  MWPTX[5][4] = 4;
  MWPTX[5][5] = 5;
  MWPTY[5][1] = -1;
  MWPTY[5][2] = 2;
  MWPTY[5][3] = -2;
  MWPTY[5][4] = 4;
  MWPTY[5][5] = -3;
  MWPTX[10][1] = -1;
  MWPTX[10][2] = -3;
  MWPTX[10][3] = -5;
  MWPTX[10][4] = 4;
  MWPTX[10][5] = 5;
  MWPTY[10][1] = 1;
  MWPTY[10][2] = 2;
  MWPTY[10][3] = 2;
  MWPTY[10][4] = 4;
  MWPTY[10][5] = 3;
  for (I = 1; I <= 11; I++) {
    combla.INCX = MWPINX[I];
    combla.INCY = MWPINY[I];
    for (K = 1; K <= 5; K++) {
      COPYX[K] = MWPX[K];
      COPYY[K] = MWPY[K];
      MWPSTX[K] = MWPTX[I][K];
      MWPSTY[K] = MWPTY[I][K];
    }
    drot(MWPN[I], COPYX, combla.INCX, COPYY, combla.INCY, MWPC[I], MWPS[I]);
    _stest(5, COPYX, MWPSTX, MWPSTX, SFAC, nout);
    _stest(5, COPYY, MWPSTY, MWPSTY, SFAC, nout);
  }
}

void _stest(
  final int LEN,
  final Array<double> SCOMP_,
  final Array<double> STRUE_,
  final Array<double> SSIZE_,
  final double SFAC,
  final Nout nout,
) {
  // ********************************* STEST **************************

  // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
  // SEE if THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
  // NEGLIGIBLE.

  // C. L. LAWSON, JPL, 1974 DEC 10

  final SCOMP = SCOMP_.having();
  final STRUE = STRUE_.having();
  final SSIZE = SSIZE_.having();
  const ZERO = 0.0;

  for (var I = 1; I <= LEN; I++) {
    final SD = SCOMP[I] - STRUE[I];
    if ((SFAC * SD).abs() <= SSIZE[I].abs() * epsilon(ZERO)) continue;

    // HERE    SCOMP[I] IS NOT CLOSE TO STRUE[I].

    if (combla.PASS) {
      // PRINT FAIL MESSAGE AND HEADER.
      combla.PASS = false;
      nout.println('                                       FAIL');
      nout.println(
          '\n CASE  N INCX INCY  I                             COMP(I)                             TRUE(I)  DIFFERENCE     SIZE(I)\n ');
    }
    nout.println(
        ' ${combla.ICASE.i4}${combla.N.i3}${combla.INCX.i5}${combla.INCY.i5}${I.i3}${SCOMP[I].d36_8}${STRUE[I].d36_8}${SD.d12_4}${SSIZE[I].d12_4}');
  }
}

void _testdsdot(
  final double SCOMP,
  final double STRUE,
  final double SSIZE,
  final double SFAC,
  final Nout nout,
) {
  // ********************************* STEST **************************

  // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
  // SEE if THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
  // NEGLIGIBLE.

  // C. L. LAWSON, JPL, 1974 DEC 10

  const ZERO = 0.0;

  final SD = SCOMP - STRUE;
  if ((SFAC * SD).abs() <= SSIZE.abs() * epsilon(ZERO)) return;

  // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

  if (combla.PASS) {
    // PRINT FAIL MESSAGE AND HEADER.
    combla.PASS = false;
    nout.println('                                       FAIL');
    nout.println(
        '\n CASE  N INCX INCY                            COMP(I)                             TRUE(I)  DIFFERENCE     SIZE(I)\n ');
  }
  nout.println(
      ' ${combla.ICASE.i4}${combla.N.i3}${combla.INCX.i5}${combla.INCY.i3}${SCOMP.d36_8}${STRUE.d36_8}${SD.d12_4}${SSIZE.d12_4}');
}

void _stest1(
  final double SCOMP1,
  final double STRUE1,
  final Array<double> SSIZE_,
  final double SFAC,
  final Nout nout,
) {
  final SSIZE = SSIZE_.having();
  // ************************* STEST1 *****************************

  // THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
  // REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
  // ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.

  // C.L. LAWSON, JPL, 1978 DEC 6

  final SCOMP = Array<double>(1), STRUE = Array<double>(1);
  SCOMP[1] = SCOMP1;
  STRUE[1] = STRUE1;
  _stest(1, SCOMP, STRUE, SSIZE, SFAC, nout);
}

// double _sdiff(final double SA, final double SB) {
// // ********************************* SDIFF **************************
// // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
//   return SA - SB;
// }

void _itest1(final int ICOMP, final int ITRUE, final Nout nout) {
  // ********************************* ITEST1 *************************

  // THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
  // EQUALITY.
  // C. L. LAWSON, JPL, 1974 DEC 10

  if (ICOMP == ITRUE) return;

  // HERE ICOMP IS NOT EQUAL TO ITRUE.

  if (combla.PASS) {
    // PRINT FAIL MESSAGE AND HEADER.
    combla.PASS = false;
    nout.println('                                       FAIL');
    nout.println(
        '\n CASE  N INCX INCY                                COMP                                TRUE     DIFFERENCE\n ');
  }
  final ID = ICOMP - ITRUE;
  nout.println(
      ' ${combla.ICASE.i4}${combla.N.i3}${combla.INCX.i5}${combla.INCY.i5}${ICOMP.i36}${ITRUE.i36}${ID.i12}');
}

void _db1nrm2(
  final int N,
  final int INCX,
  final double THRESH,
  final Nout nout,
) {
  // Compare NRM2 with a reference computation using combinations
  // of the following values:

  // 0, very small, small, ulp, 1, 1/ulp, big, very big, infinity, NaN

  // one of these values is used to initialize x(1) and x(2:N) is
  // filled with random values from [-1,1] scaled by another of
  // these values.

  // This routine is adapted from the test suite provided by
  // Anderson E. (2017)
  // Algorithm 978: Safe Scaling in the Level 1 BLAS
  // ACM Trans Math Softw 44:1--28
  // https://doi.org/10.1145/3061665

  const NMAX = 20, NV = 10;
  const HALF = 0.5, ONE = 1.0, TWO = 2.0, ZERO = 0.0;
  const BIGNUM = 0.99792015476735990583e+292,
      SAFMAX = 0.44942328371557897693e+308,
      SAFMIN = 0.22250738585072013831e-307,
      SMLNUM = 0.10020841800044863890e-291,
      ULP = 0.22204460492503130808e-015;
  double ROGUE, SNRM, TRAT, V0 = 0, V1, WORKSSQ, Y1, Y2, YMAX, YMIN, YNRM, ZNRM;
  int I, IV, IW, IX;
  bool FIRST;
  final VALUES = Array<double>(NV),
      WORK = Array<double>(NMAX),
      X = Array<double>(NMAX),
      Z = Array<double>(NMAX);
  VALUES[1] = ZERO;
  VALUES[2] = TWO * SAFMIN;
  VALUES[3] = SMLNUM;
  VALUES[4] = ULP;
  VALUES[5] = ONE;
  VALUES[6] = ONE / ULP;
  VALUES[7] = BIGNUM;
  VALUES[8] = SAFMAX;
  VALUES[9] = _dxvals(V0, 2);
  VALUES[10] = _dxvals(V0, 3);
  ROGUE = -1234.5678;
  FIRST = true;

  // Check that the arrays are large enough

  if (N * INCX.abs() > NMAX) {
    nout.println(
        ' Not enough space to test DNRM2 : NMAX = ${NMAX.i6}, INCX = ${INCX.i6}\n   N = ${N.i6}, must be at least ${(N * INCX.abs()).i6}');
    return;
  }

  // Zero-sized inputs are tested in STEST1.
  if (N <= 0) return;

  // Generate (N-1) values in (-1,1).

  for (I = 2; I <= N; I++) {
    random_number(WORK.box(I));
    WORK[I] = ONE - TWO * WORK[I];
  }

  // Compute the sum of squares of the random values
  // by an unscaled algorithm.

  WORKSSQ = ZERO;
  for (I = 2; I <= N; I++) {
    WORKSSQ += WORK[I] * WORK[I];
  }

  // Construct the test vector with one known value
  // and the rest from the random work array multiplied
  // by a scaling factor.

  for (IV = 1; IV <= NV; IV++) {
    V0 = VALUES[IV];
    if (V0.abs() > ONE) {
      V0 *= HALF;
    }
    Z[1] = V0;
    for (IW = 1; IW <= NV; IW++) {
      V1 = VALUES[IW];
      if (V1.abs() > ONE) {
        V1 = (V1 * HALF) / sqrt(N);
      }
      for (I = 2; I <= N; I++) {
        Z[I] = V1 * WORK[I];
      }

      // Compute the expected value of the 2-norm

      Y1 = V0.abs();
      if (N > 1) {
        Y2 = V1.abs() * sqrt(WORKSSQ);
      } else {
        Y2 = ZERO;
      }
      YMIN = min(Y1, Y2);
      YMAX = max(Y1, Y2);

      // Expected value is NaN if either is NaN. The test
      // for YMIN == YMAX avoids further computation if both
      // are infinity.

      if ((Y1 != Y1) || (Y2 != Y2)) {
        // Add to propagate NaN
        YNRM = Y1 + Y2;
      } else if (YMAX == ZERO) {
        YNRM = ZERO;
      } else if (YMIN == YMAX) {
        YNRM = sqrt(TWO) * YMAX;
      } else {
        YNRM = YMAX * sqrt(ONE + pow((YMIN / YMAX), 2));
      }

      // Fill the input array to DNRM2 with steps of incx

      for (I = 1; I <= N; I++) {
        X[I] = ROGUE;
      }
      IX = 1;
      if (INCX < 0) IX = 1 - (N - 1) * INCX;
      for (I = 1; I <= N; I++) {
        X[IX] = Z[I];
        IX += INCX;
      }

      // Call DNRM2 to compute the 2-norm

      SNRM = dnrm2(N, X, INCX);

      // Compare SNRM and ZNRM.  Roundoff error grows like O(n)
      // in this implementation so we scale the test ratio accordingly.

      if (INCX == 0) {
        ZNRM = sqrt(N) * X[1].abs();
      } else {
        ZNRM = YNRM;
      }

      // The tests for NaN rely on the compiler not being overly
      // aggressive and removing the statements altogether.
      if ((SNRM != SNRM) || (ZNRM != ZNRM)) {
        if ((SNRM != SNRM) != (ZNRM != ZNRM)) {
          TRAT = ONE / ULP;
        } else {
          TRAT = ZERO;
        }
      } else if (SNRM == ZNRM) {
        TRAT = ZERO;
      } else if (ZNRM == ZERO) {
        TRAT = SNRM / ULP;
      } else {
        TRAT = ((SNRM - ZNRM).abs() / ZNRM) / (N * ULP);
      }
      if ((TRAT != TRAT) || (TRAT >= THRESH)) {
        if (FIRST) {
          FIRST = false;
          nout.println('                                       FAIL');
        }
        nout.println(
            ' DNRM2 : N=${N.i6}, INCX=${INCX.i4}, IV=${IV.i2}, IW=${IW.i2}, test=${TRAT.e15_8}');
      }
    }
  }
}

double _dxvals(final double XX, final int K) {
  final Y = huge(XX);
  final Z = Y * Y;
  if (K == 1) {
    return -Z;
  } else if (K == 2) {
    return Z;
  } else if (K == 3) {
    return Z / Z;
  }
  throw UnimplementedError();
}

void main() {
  final nout = Nout(stdout);
  dblat1(nout, lapackTestDriver);
  exit(lapackTestDriver.errors);
}
