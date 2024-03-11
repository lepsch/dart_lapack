import 'dart:io';
import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zdotu.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/intrinsics/epsilon.dart';
import 'package:lapack/src/intrinsics/huge.dart';
import 'package:lapack/src/intrinsics/random_number.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'common.dart';

void main() {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  final NOUT = Nout(stdout);
  const SFAC = 9.765625e-4;

  NOUT.println(' Complex BLAS Test Program Results\n ');
  for (var IC = 1; IC <= 10; IC++) {
    combla.ICASE = IC;
    _header(NOUT);

    // Initialize PASS, INCX, INCY, and MODE for a new case.
    // The value 9999 for INCX, INCY or MODE will appear in the
    // detailed  output, if any, for cases that do not involve
    // these parameters.

    combla.PASS = true;
    combla.INCX = 9999;
    combla.INCY = 9999;
    combla.MODE = 9999;
    if (combla.ICASE <= 5) {
      _check2(SFAC, NOUT);
    } else if (combla.ICASE >= 6) {
      _check1(SFAC, NOUT);
    }
    // -- Print
    if (combla.PASS) {
      NOUT.println('                                    ----- PASS -----');
    }
  }
}

void _header(final Nout NOUT) {
  const L = [
    'ZDOTC ',
    'ZDOTU ',
    'ZAXPY ',
    'ZCOPY ',
    'ZSWAP ',
    'DZNRM2',
    'DZASUM',
    'ZSCAL ',
    'ZDSCAL',
    'IZAMAX',
  ];

  NOUT.println(
      '\n Test of subprogram number${combla.ICASE.i3}${' ' * 12}${L[combla.ICASE - 1].a6}');
}

void _check1(final double SFAC, final Nout NOUT) {
  const THRESH = 10.0;
  int I, IX, LEN, NP1;
  var (SA, CA) = (0.3, Complex(0.4, -0.7));
  final CX = Array<Complex>(8),
      CXR = Array<Complex>(15),
      MWPCS = Array<Complex>(5),
      MWPCT = Array<Complex>(5);
  final CV = Matrix3d.fromData(
      [
        (0.1, 0.1), (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), //
        (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), //
        (0.3, -0.4), (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), //
        (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), //
        (0.1, -0.3), (0.5, -0.1), (5.0, 6.0), (5.0, 6.0), //
        (5.0, 6.0), (5.0, 6.0), (5.0, 6.0), (5.0, 6.0), //
        (0.1, 0.1), (-0.6, 0.1), (0.1, -0.3), (7.0, 8.0), //
        (7.0, 8.0), (7.0, 8.0), (7.0, 8.0), (7.0, 8.0), //
        (0.3, 0.1), (0.5, 0.0), (0.0, 0.5), (0.0, 0.2), //
        (2.0, 3.0), (2.0, 3.0), (2.0, 3.0), (2.0, 3.0), //
        (0.1, 0.1), (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), //
        (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), //
        (0.3, -0.4), (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), //
        (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), //
        (0.1, -0.3), (8.0, 9.0), (0.5, -0.1), (2.0, 5.0), //
        (2.0, 5.0), (2.0, 5.0), (2.0, 5.0), (2.0, 5.0), //
        (0.1, 0.1), (3.0, 6.0), (-0.6, 0.1), (4.0, 7.0), //
        (0.1, -0.3), (7.0, 2.0), (7.0, 2.0), (7.0, 2.0), //
        (0.3, 0.1), (5.0, 8.0), (0.5, 0.0), (6.0, 9.0), //
        (0.0, 0.5), (8.0, 3.0), (0.0, 0.2), (9.0, 4.0), //
      ].toComplexList(),
      (8, 5, 2));
  final CVR = Array.fromList([
    (8.0, 8.0), (-7.0, -7.0), (9.0, 9.0), (5.0, 5.0), //
    (9.0, 9.0), (8.0, 8.0), (7.0, 7.0), (7.0, 7.0)
  ].toComplexList());
  final STRUE2 = Array.fromList([0.0, 0.5, 0.6, 0.7, 0.8]);
  final STRUE4 = Array.fromList([0.0, 0.7, 1.0, 1.3, 1.6]);
  final CTRUE5 = Matrix3d.fromData(
      [
        (0.1, 0.1), (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), //
        (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), //
        (-0.16, -0.37), (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), //
        (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), //
        (-0.17, -0.19), (0.13, -0.39), (5.0, 6.0), (5.0, 6.0), //
        (5.0, 6.0), (5.0, 6.0), (5.0, 6.0), (5.0, 6.0), //
        (0.11, -0.03), (-0.17, 0.46), (-0.17, -0.19), (7.0, 8.0), //
        (7.0, 8.0), (7.0, 8.0), (7.0, 8.0), (7.0, 8.0), //
        (0.19, -0.17), (0.20, -0.35), (0.35, 0.20), (0.14, 0.08), //
        (2.0, 3.0), (2.0, 3.0), (2.0, 3.0), (2.0, 3.0),
        (0.1, 0.1), (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), //
        (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), //
        (-0.16, -0.37), (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), //
        (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), //
        (-0.17, -0.19), (8.0, 9.0), (0.13, -0.39), (2.0, 5.0), //
        (2.0, 5.0), (2.0, 5.0), (2.0, 5.0), (2.0, 5.0), //
        (0.11, -0.03), (3.0, 6.0), (-0.17, 0.46), (4.0, 7.0), //
        (-0.17, -0.19), (7.0, 2.0), (7.0, 2.0), (7.0, 2.0), //
        (0.19, -0.17), (5.0, 8.0), (0.20, -0.35), (6.0, 9.0), //
        (0.35, 0.20), (8.0, 3.0), (0.14, 0.08), (9.0, 4.0),
      ].toComplexList(),
      (8, 5, 2));
  final CTRUE6 = Matrix3d.fromData(
      [
        (0.1, 0.1), (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), //
        (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), (1.0, 2.0), //
        (0.09, -0.12), (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), //
        (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), (3.0, 4.0), //
        (0.03, -0.09), (0.15, -0.03), (5.0, 6.0), (5.0, 6.0), //
        (5.0, 6.0), (5.0, 6.0), (5.0, 6.0), (5.0, 6.0), //
        (0.03, 0.03), (-0.18, 0.03), (0.03, -0.09), (7.0, 8.0), //
        (7.0, 8.0), (7.0, 8.0), (7.0, 8.0), (7.0, 8.0), //
        (0.09, 0.03), (0.15, 0.00), (0.00, 0.15), (0.00, 0.06), //
        (2.0, 3.0), (2.0, 3.0), (2.0, 3.0), (2.0, 3.0),
        (0.1, 0.1), (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), //
        (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), (4.0, 5.0), //
        (0.09, -0.12), (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), //
        (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), (6.0, 7.0), //
        (0.03, -0.09), (8.0, 9.0), (0.15, -0.03), (2.0, 5.0), //
        (2.0, 5.0), (2.0, 5.0), (2.0, 5.0), (2.0, 5.0), //
        (0.03, 0.03), (3.0, 6.0), (-0.18, 0.03), (4.0, 7.0), //
        (0.03, -0.09), (7.0, 2.0), (7.0, 2.0), (7.0, 2.0), //
        (0.09, 0.03), (5.0, 8.0), (0.15, 0.00), (6.0, 9.0), //
        (0.00, 0.15), (8.0, 3.0), (0.00, 0.06), (9.0, 4.0),
      ].toComplexList(),
      (8, 5, 2));
  final ITRUE3 = Array.fromList([0, 1, 2, 2, 2]);
  final ITRUEC = Array.fromList([0, 1, 1, 1, 1]);

  for (combla.INCX = 1; combla.INCX <= 2; combla.INCX++) {
    for (NP1 = 1; NP1 <= 5; NP1++) {
      combla.N = NP1 - 1;
      LEN = 2 * max(combla.N, 1);
      // .. Set vector arguments ..
      for (I = 1; I <= LEN; I++) {
        CX[I] = CV[I][NP1][combla.INCX];
      }
      if (combla.ICASE == 6) {
        // .. dznrm2 ..
        // Test scaling when some entries are tiny or huge
        _zb1nrm2(combla.N, (combla.INCX - 2) * 2, THRESH, NOUT);
        _zb1nrm2(combla.N, combla.INCX, THRESH, NOUT);
        // Test with hardcoded mid range entries
        _stest1(dznrm2(combla.N, CX, combla.INCX), STRUE2[NP1], STRUE2(NP1),
            SFAC, NOUT);
      } else if (combla.ICASE == 7) {
        // .. dzasum ..
        _stest1(dzasum(combla.N, CX, combla.INCX), STRUE4[NP1], STRUE4(NP1),
            SFAC, NOUT);
      } else if (combla.ICASE == 8) {
        // .. ZSCAL ..
        zscal(combla.N, CA, CX, combla.INCX);
        _ctest(LEN, CX, CTRUE5(1, NP1, combla.INCX).asArray(),
            CTRUE5(1, NP1, combla.INCX).asArray(), SFAC, NOUT);
      } else if (combla.ICASE == 9) {
        // .. ZDSCAL ..
        zdscal(combla.N, SA, CX, combla.INCX);
        _ctest(LEN, CX, CTRUE6(1, NP1, combla.INCX).asArray(),
            CTRUE6(1, NP1, combla.INCX).asArray(), SFAC, NOUT);
      } else if (combla.ICASE == 10) {
        // .. izamax ..
        _itest1(izamax(combla.N, CX, combla.INCX), ITRUE3[NP1], NOUT);
        for (I = 1; I <= LEN; I++) {
          CX[I] = Complex(42.0, 43.0);
        }
        _itest1(izamax(combla.N, CX, combla.INCX), ITRUEC[NP1], NOUT);
      } else {
        NOUT.println(' Shouldn\'t be here in CHECK1');
        exit(1);
      }
    }
    if (combla.ICASE == 10) {
      combla.N = 8;
      IX = 1;
      for (I = 1; I <= combla.N; I++) {
        CXR[IX] = CVR[I];
        IX += combla.INCX;
      }
      _itest1(izamax(combla.N, CXR, combla.INCX), 3, NOUT);
    }
  }

  combla.INCX = 1;
  if (combla.ICASE == 8) {
    // ZSCAL
    // Add a test for alpha equal to zero.
    CA = Complex.zero;
    for (I = 1; I <= 5; I++) {
      MWPCT[I] = Complex.zero;
      MWPCS[I] = Complex(1.0, 1.0);
    }
    zscal(5, CA, CX, combla.INCX);
    _ctest(5, CX, MWPCT, MWPCS, SFAC, NOUT);
  } else if (combla.ICASE == 9) {
    // ZDSCAL
    // Add a test for alpha equal to zero.
    SA = 0.0;
    for (I = 1; I <= 5; I++) {
      MWPCT[I] = Complex.zero;
      MWPCS[I] = Complex(1.0, 1.0);
    }
    zdscal(5, SA, CX, combla.INCX);
    _ctest(5, CX, MWPCT, MWPCS, SFAC, NOUT);
    // Add a test for alpha equal to one.
    SA = 1.0;
    for (I = 1; I <= 5; I++) {
      MWPCT[I] = CX[I];
      MWPCS[I] = CX[I];
    }
    zdscal(5, SA, CX, combla.INCX);
    _ctest(5, CX, MWPCT, MWPCS, SFAC, NOUT);
    // Add a test for alpha equal to minus one.
    SA = -1.0;
    for (I = 1; I <= 5; I++) {
      MWPCT[I] = -CX[I];
      MWPCS[I] = -CX[I];
    }
    zdscal(5, SA, CX, combla.INCX);
    _ctest(5, CX, MWPCT, MWPCS, SFAC, NOUT);
  }
}

void _check2(final double SFAC, final Nout NOUT) {
  int I, KI, KN, KSIZE, LENX, LENY, LINCX, LINCY, MX, MY;
  final CDOT = Array<Complex>(1),
      CTY0 = Array<Complex>(1),
      CX = Array<Complex>(7),
      CX0 = Array<Complex>(1),
      CY = Array<Complex>(7),
      CY0 = Array<Complex>(1);
  const CA = Complex(0.4, -0.7);
  final INCXS = Array.fromList([1, 2, -2, -1]);
  final INCYS = Array.fromList([1, -2, 1, -2]);
  final LENS = Matrix.fromData([1, 1, 2, 4, 1, 1, 3, 7], (4, 2));
  final NS = Array.fromList([0, 1, 2, 4]);
  final CX1 = Array.fromList([
    (0.7, -0.8), (-0.4, -0.7), (-0.1, -0.9), (0.2, -0.8), //
    (-0.9, -0.4), (0.1, 0.4), (-0.6, 0.6)
  ].toComplexList());
  final CY1 = Array.fromList([
    (0.6, -0.6), (-0.9, 0.5), (0.7, -0.6), (0.1, -0.5), //
    (-0.1, -0.2), (-0.5, -0.3), (0.8, -0.7)
  ].toComplexList());

  final CT8 = Matrix3d.fromData(
      [
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.32, -1.41), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.32, -1.41), //
        (-1.55, 0.5), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.32, -1.41), (-1.55, 0.5), (0.03, -0.89),
        (-0.38, -0.96), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.32, -1.41), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (-0.07, -0.89), //
        (-0.9, 0.5), (0.42, -1.41), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.78, 0.06), (-0.9, 0.5), (0.06, -0.13), (0.1, -0.5), //
        (-0.77, -0.49), (-0.5, -0.3), (0.52, -1.51),
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.32, -1.41), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (-0.07, -0.89), //
        (-1.18, -0.31), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.78, 0.06), (-1.54, 0.97), (0.03, -0.89),
        (-0.18, -1.31), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.32, -1.41), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.32, -1.41), //
        (-0.9, 0.5), (0.05, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.32, -1.41), (-0.9, 0.5), (0.05, -0.6), (0.1, -0.5), //
        (-0.77, -0.49), (-0.5, -0.3), (0.32, -1.16),
      ].toComplexList(),
      (7, 4, 4));

  final CT7 = Matrix.fromData(
      [
        (0.0, 0.0), (-0.06, -0.90), (0.65, -0.47), (-0.34, -1.22), //
        (0.0, 0.0), (-0.06, -0.90), (-0.59, -1.46), (-1.04, -0.04), //
        (0.0, 0.0), (-0.06, -0.90), (-0.83, 0.59), (0.07, -0.37), //
        (0.0, 0.0), (-0.06, -0.90), (-0.76, -1.15), (-1.33, -1.82)
      ].toComplexList(),
      (4, 1));
  final CT6 = Matrix.fromData(
      [
        (0.0, 0.0), (0.90, 0.06), (0.91, -0.77), (1.80, -0.10), //
        (0.0, 0.0), (0.90, 0.06), (1.45, 0.74), (0.20, 0.90), //
        (0.0, 0.0), (0.90, 0.06), (-0.55, 0.23), (0.83, -0.39), //
        (0.0, 0.0), (0.90, 0.06), (1.04, 0.79), (1.95, 1.22)
      ].toComplexList(),
      (4, 1));

  final CT10X = Matrix3d.fromData(
      [
        (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.6, -0.6), //
        (-0.9, 0.5), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.6, -0.6), (-0.9, 0.5), (0.7, -0.6), (0.1, -0.5), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
        (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.7, -0.6), //
        (-0.4, -0.7), (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.8, -0.7), (-0.4, -0.7), (-0.1, -0.2), (0.2, -0.8), //
        (0.7, -0.6), (0.1, 0.4), (0.6, -0.6),
        (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (-0.9, 0.5), //
        (-0.4, -0.7), (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.1, -0.5), (-0.4, -0.7), (0.7, -0.6), (0.2, -0.8), //
        (-0.9, 0.5), (0.1, 0.4), (0.6, -0.6),
        (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.6, -0.6), //
        (0.7, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.6, -0.6), (0.7, -0.6), (-0.1, -0.2), (0.8, -0.7), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
      ].toComplexList(),
      (7, 4, 4));

  final CT10Y = Matrix3d.fromData(
      [
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.7, -0.8), //
        (-0.4, -0.7), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.7, -0.8), (-0.4, -0.7), (-0.1, -0.9), (0.2, -0.8), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (-0.1, -0.9), //
        (-0.9, 0.5), (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (-0.6, 0.6), (-0.9, 0.5), (-0.9, -0.4), (0.1, -0.5), //
        (-0.1, -0.9), (-0.5, -0.3), (0.7, -0.8),
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (-0.1, -0.9), //
        (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (-0.6, 0.6), (-0.9, -0.4), (-0.1, -0.9), (0.7, -0.8), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
        (0.6, -0.6), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.7, -0.8), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.7, -0.8), //
        (-0.9, 0.5), (-0.4, -0.7), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.7, -0.8), (-0.9, 0.5), (-0.4, -0.7), (0.1, -0.5), //
        (-0.1, -0.9), (-0.5, -0.3), (0.2, -0.8),
      ].toComplexList(),
      (7, 4, 4));

  final CSIZE1 = Array.fromList(
      [(0.0, 0.0), (0.9, 0.9), (1.63, 1.73), (2.90, 2.78)].toComplexList());
  final CSIZE3 = Array.fromList([
    (0.0, 0.0),
    (0.0, 0.0),
    (0.0, 0.0),
    (0.0, 0.0),
    (0.0, 0.0),
    (0.0, 0.0),
    (0.0, 0.0),
    (1.17, 1.17),
    (1.17, 1.17),
    (1.17, 1.17),
    (1.17, 1.17),
    (1.17, 1.17),
    (1.17, 1.17),
    (1.17, 1.17)
  ].toComplexList());
  final CSIZE2 = Matrix.fromData(
      [
        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), //
        (0.0, 0.0), (0.0, 0.0), (1.54, 1.54), (1.54, 1.54), (1.54, 1.54), //
        (1.54, 1.54), (1.54, 1.54), (1.54, 1.54), (1.54, 1.54)
      ].toComplexList(),
      (7, 1));

  for (KI = 1; KI <= 4; KI++) {
    combla.INCX = INCXS[KI];
    combla.INCY = INCYS[KI];
    MX = (combla.INCX).abs();
    MY = (combla.INCY).abs();

    for (KN = 1; KN <= 4; KN++) {
      combla.N = NS[KN];
      KSIZE = min(2, KN);
      LENX = LENS[KN][MX];
      LENY = LENS[KN][MY];
      // .. initialize all argument arrays ..
      for (I = 1; I <= 7; I++) {
        CX[I] = CX1[I];
        CY[I] = CY1[I];
      }
      if (combla.ICASE == 1) {
        // .. zdotc ..
        CDOT[1] = zdotc(combla.N, CX, combla.INCX, CY, combla.INCY);
        _ctest(1, CDOT, CT6(KN, KI).asArray(), CSIZE1(KN), SFAC, NOUT);
      } else if (combla.ICASE == 2) {
        // .. zdotu ..
        CDOT[1] = zdotu(combla.N, CX, combla.INCX, CY, combla.INCY);
        _ctest(1, CDOT, CT7(KN, KI).asArray(), CSIZE1(KN), SFAC, NOUT);
      } else if (combla.ICASE == 3) {
        // .. ZAXPY ..
        zaxpy(combla.N, CA, CX, combla.INCX, CY, combla.INCY);
        _ctest(LENY, CY, CT8(1, KN, KI).asArray(), CSIZE2(1, KSIZE).asArray(),
            SFAC, NOUT);
      } else if (combla.ICASE == 4) {
        // .. ZCOPY ..
        zcopy(combla.N, CX, combla.INCX, CY, combla.INCY);
        _ctest(LENY, CY, CT10Y(1, KN, KI).asArray(), CSIZE3, 1.0, NOUT);
        if (KI == 1) {
          CX0[1] = Complex(42.0, 43.0);
          CY0[1] = Complex(44.0, 45.0);
          if (combla.N == 0) {
            CTY0[1] = CY0[1];
          } else {
            CTY0[1] = CX0[1];
          }
          LINCX = combla.INCX;
          combla.INCX = 0;
          LINCY = combla.INCY;
          combla.INCY = 0;
          zcopy(combla.N, CX0, combla.INCX, CY0, combla.INCY);
          _ctest(1, CY0, CTY0, CSIZE3, 1.0, NOUT);
          combla.INCX = LINCX;
          combla.INCY = LINCY;
        }
      } else if (combla.ICASE == 5) {
        // .. ZSWAP ..
        zswap(combla.N, CX, combla.INCX, CY, combla.INCY);
        _ctest(LENX, CX, CT10X(1, KN, KI).asArray(), CSIZE3, 1.0, NOUT);
        _ctest(LENY, CY, CT10Y(1, KN, KI).asArray(), CSIZE3, 1.0, NOUT);
      } else {
        NOUT.println(' Shouldn\'t be here in CHECK2');
        exit(1);
      }
    }
  }
}

void _stest(
    final int LEN,
    final Array<double> SCOMP_,
    final Array<double> STRUE_,
    final Array<double> SSIZE_,
    final double SFAC,
    final Nout NOUT) {
  // ********************************* STEST **************************

  // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
  // SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
  // NEGLIGIBLE.

  // C. L. LAWSON, JPL, 1974 DEC 10
  final SCOMP = SCOMP_.having(length: LEN);
  final STRUE = STRUE_.having(length: LEN);
  final SSIZE = SSIZE_.having(length: LEN);
  const ZERO = 0.0;
  double SD;
  int I;

  for (I = 1; I <= LEN; I++) {
    SD = SCOMP[I] - STRUE[I];
    if ((SFAC * SD).abs() <= (SSIZE[I]).abs() * epsilon(ZERO)) continue;

    // HERE    SCOMP[I] IS NOT CLOSE TO STRUE[I].

    if (combla.PASS) {
      // PRINT FAIL MESSAGE AND HEADER.
      combla.PASS = false;
      NOUT.println('                                       FAIL');
      NOUT.println(
          '\n CASE  combla.N combla.INCX combla.INCY combla.MODE  I                             COMP[I]                             TRUE[I]  DIFFERENCE     SIZE[I]\n ');
    }
    NOUT.println(' ${combla.ICASE.i4}${combla.N.i3}${[
      combla.INCX,
      combla.INCY,
      combla.MODE
    ].i5()}${I.i3}${[SCOMP[I], STRUE[I]].d36_8}${[SD, SSIZE[I]].d12_4}');
  }
}

void _stest1(
  final double SCOMP1,
  final double STRUE1,
  final Array<double> SSIZE_,
  final double SFAC,
  Nout NOUT,
) {
  // ************************* STEST1 *****************************

  // THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
  // REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
  // ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.

  // C.L. LAWSON, JPL, 1978 DEC 6

  final SSIZE = SSIZE_.having();
  final SCOMP = Array<double>(1), STRUE = Array<double>(1);

  SCOMP[1] = SCOMP1;
  STRUE[1] = STRUE1;
  _stest(1, SCOMP, STRUE, SSIZE, SFAC, NOUT);
}

// double _sdiff(final double SA, final double SB) {
//   // ********************************* SDIFF **************************
//   // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

//   return SA - SB;
// }

void _ctest(
  final int LEN,
  final Array<Complex> CCOMP_,
  final Array<Complex> CTRUE_,
  final Array<Complex> CSIZE_,
  final double SFAC,
  final Nout NOUT,
) {
  // **************************** CTEST *****************************

  // C.L. LAWSON, JPL, 1978 DEC 6

  int I;
  final SCOMP = Array<double>(20),
      SSIZE = Array<double>(20),
      STRUE = Array<double>(20);
  final CCOMP = CCOMP_.having(length: LEN);
  final CTRUE = CTRUE_.having(length: LEN);
  final CSIZE = CSIZE_.having(length: LEN);

  for (I = 1; I <= LEN; I++) {
    SCOMP[2 * I - 1] = CCOMP[I].real;
    SCOMP[2 * I] = CCOMP[I].imaginary;
    STRUE[2 * I - 1] = CTRUE[I].real;
    STRUE[2 * I] = CTRUE[I].imaginary;
    SSIZE[2 * I - 1] = CSIZE[I].real;
    SSIZE[2 * I] = CSIZE[I].imaginary;
  }

  _stest(2 * LEN, SCOMP, STRUE, SSIZE, SFAC, NOUT);
}

void _itest1(final int ICOMP, final int ITRUE, Nout NOUT) {
  // ********************************* ITEST1 *************************

  // THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
  // EQUALITY.
  // C. L. LAWSON, JPL, 1974 DEC 10

  if (ICOMP == ITRUE) return;

  // HERE ICOMP IS NOT EQUAL TO ITRUE.

  if (combla.PASS) {
    // PRINT FAIL MESSAGE AND HEADER.
    combla.PASS = false;
    NOUT.println('                                       FAIL');
    NOUT.println(
        '\n CASE  N INCX INCY MODE                                COMP                                TRUE     DIFFERENCE\n ');
  }

  final ID = ICOMP - ITRUE;
  NOUT.println(' ${combla.ICASE.i4}${combla.N.i3}${[
    combla.INCX,
    combla.INCY,
    combla.MODE
  ].i5}${[ICOMP, ITRUE].i36()}${ID.i12}');
}

void _zb1nrm2(
  final int N,
  final int INCX,
  final double THRESH,
  final Nout NOUT,
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
  const HALF = 0.5, ONE = 1.0, TWO = 2.0, THREE = 3.0, ZERO = 0.0;
  const BIGNUM = 0.99792015476735990583e+292,
      SAFMAX = 0.44942328371557897693e+308,
      SAFMIN = 0.22250738585072013831e-307,
      SMLNUM = 0.10020841800044863890e-291,
      ULP = 0.22204460492503130808e-015;
  Complex ROGUE;
  double SNRM, TRAT, V0 = 0, V1, WORKSSQ, Y1, Y2, YMAX, YMIN, YNRM, ZNRM;
  int I, IV, IW, IX;
  bool FIRST;
  final X = Array<Complex>(NMAX), Z = Array<Complex>(NMAX);
  final VALUES = Array<double>(NV), WORK = Array<double>(NMAX);

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
  ROGUE = Complex(1234.5678, -1234.5678);
  FIRST = true;

  // Check that the arrays are large enough

  if (N * (INCX).abs() > NMAX) {
    NOUT.println(
        ' Not enough space to test DZNRM2: NMAX = ${NMAX.i6}, INCX = ${INCX.i6}\n   N = ${N.i6}, must be at least ${(N * INCX.abs()).i6}');
    return;
  }

  // Zero-sized inputs are tested in STEST1.
  if (N <= 0) return;

  // Generate 2*(N-1) values in (-1,1).

  final KS = 2 * (N - 1);
  for (I = 1; I <= KS; I++) {
    random_number(WORK(I));
    WORK[I] = ONE - TWO * WORK[I];
  }

  // Compute the sum of squares of the random values
  // by an unscaled algorithm.

  WORKSSQ = ZERO;
  for (I = 1; I <= KS; I++) {
    WORKSSQ += WORK[I] * WORK[I];
  }

  // Construct the test vector with one known value
  // and the rest from the random work array multiplied
  // by a scaling factor.

  for (IV = 1; IV <= NV; IV++) {
    V0 = VALUES[IV];
    if ((V0).abs() > ONE) {
      V0 = V0 * HALF * HALF;
    }
    Z[1] = Complex(V0, -THREE * V0);
    for (IW = 1; IW <= NV; IW++) {
      V1 = VALUES[IW];
      if ((V1).abs() > ONE) {
        V1 = (V1 * HALF) / sqrt((KS + 1));
      }
      for (I = 1; I <= N - 1; I++) {
        Z[I + 1] = Complex(V1 * WORK[2 * I - 1], V1 * WORK[2 * I]);
      }

      // Compute the expected value of the 2-norm

      Y1 = (V0).abs() * sqrt(10.0);
      if (N > 1) {
        Y2 = (V1).abs() * sqrt(WORKSSQ);
      } else {
        Y2 = ZERO;
      }
      YMIN = min(Y1, Y2);
      YMAX = max(Y1, Y2);

      // Expected value is NaN if either is NaN. The test
      // for YMIN == YMAX avoids further computation if both
      // are infinity.

      if (Y1.isNaN || Y2.isNaN) {
        // add to propagate NaN
        YNRM = _dxvals(V0, 3);
      } else if (YMIN == YMAX) {
        YNRM = sqrt(TWO) * YMAX;
      } else if (YMAX == ZERO) {
        YNRM = ZERO;
      } else {
        YNRM = YMAX * sqrt(ONE + pow(YMIN / YMAX, 2));
      }

      // Fill the input array to dznrm2 with steps of incx

      for (I = 1; I <= N; I++) {
        X[I] = ROGUE;
      }
      IX = 1;
      if (INCX < 0) IX = 1 - (N - 1) * INCX;
      for (I = 1; I <= N; I++) {
        X[IX] = Z[I];
        IX += INCX;
      }

      // Call dznrm2 to compute the 2-norm

      SNRM = dznrm2(N, X, INCX);

      // Compare SNRM and ZNRM.  Roundoff error grows like O(n)
      // in this implementation so we scale the test ratio accordingly.

      if (INCX == 0) {
        Y1 = X[1].real.abs();
        Y2 = X[1].imaginary.abs();
        YMIN = min(Y1, Y2);
        YMAX = max(Y1, Y2);
        if (Y1.isNaN || Y2.isNaN) {
          // add to propagate NaN
          ZNRM = _dxvals(V0, 3);
        } else if (YMIN == YMAX) {
          ZNRM = sqrt(TWO) * YMAX;
        } else if (YMAX == ZERO) {
          ZNRM = ZERO;
        } else {
          ZNRM = YMAX * sqrt(ONE + pow(YMIN / YMAX, 2));
        }
        ZNRM = sqrt(N.toDouble()) * ZNRM;
      } else {
        ZNRM = YNRM;
      }

      // The tests for NaN rely on the compiler not being overly
      // aggressive and removing the statements altogether.
      if (SNRM.isNaN || ZNRM.isNaN) {
        if (SNRM.isNaN != ZNRM.isNaN) {
          TRAT = ONE / ULP;
        } else {
          TRAT = ZERO;
        }
      } else if (ZNRM == ZERO) {
        TRAT = SNRM / ULP;
      } else {
        TRAT = ((SNRM - ZNRM).abs() / ZNRM) / (TWO * N * ULP);
      }
      if (TRAT >= THRESH) {
        if (FIRST) {
          FIRST = false;
          NOUT.println('                                       FAIL');
        }
        NOUT.println(
            ' DZNRM2: N=${N.i6}, INCX=${INCX.i4}, IV=${IV.i2}, IW=${IW.i2}, test=${TRAT.e15_8}');
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
