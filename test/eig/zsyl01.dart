// -- LAPACK test routine --
import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/ztrsyl.dart';
import 'package:lapack/src/ztrsyl3.dart';

import '../matgen/zlatmr.dart';

void zsyl01(
  final double THRESH,
  final Array<int> NFAIL_,
  final Array<double> RMAX_,
  final Array<int> NINFO_,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NFAIL = NFAIL_.having(length: 3);
  final NINFO = NINFO_.having(length: 2);
  final RMAX = RMAX_.having(length: 2);
  const ZERO = 0.0, ONE = 1.0;
  const MAXM = 185, MAXN = 192, LDSWORK = 36;
  String TRANA = '', TRANB = '';
  int I, ISGN, ITRANA, ITRANB, J, KLA, KUA, KLB, KUB, M, N;
  double ANRM, BNRM, BIGNUM, EPS, RES, RES1, SMLNUM, TNRM, XNRM;
  Complex RMUL;
  final DUML = Array<Complex>(MAXM),
      DUMR = Array<Complex>(MAXN),
      D = Array<Complex>(max(MAXM, MAXN));
  final DUM = Array<double>(MAXN), VM = Array<double>(2);
  final ISEED = Array<int>(4), IWORK = Array<int>(MAXM + MAXN + 2);
  final A = Matrix<Complex>(MAXM, MAXM),
      B = Matrix<Complex>(MAXN, MAXN),
      C = Matrix<Complex>(MAXM, MAXN),
      CC = Matrix<Complex>(MAXM, MAXN),
      X = Matrix<Complex>(MAXM, MAXN);
  final SWORK = Matrix<double>(LDSWORK, 103);
  final INFO = Box(0), IINFO = Box(0);
  final SCALE = Box(0.0), SCALE3 = Box(0.0);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // Expect INFO.value = 0
  VM[1] = ONE;
  // Expect INFO.value = 1
  VM[2] = 0.05;

  // Begin test loop

  NINFO[1] = 0;
  NINFO[2] = 0;
  NFAIL[1] = 0;
  NFAIL[2] = 0;
  NFAIL[3] = 0;
  RMAX[1] = ZERO;
  RMAX[2] = ZERO;
  KNT.value = 0;
  ISEED[1] = 1;
  ISEED[2] = 1;
  ISEED[3] = 1;
  ISEED[4] = 1;
  SCALE.value = ONE;
  SCALE3.value = ONE;
  for (J = 1; J <= 2; J++) {
    for (ISGN = -1; 2 < 0 ? ISGN >= 1 : ISGN <= 1; ISGN += 2) {
      // Reset seed (overwritten by LATMR)
      ISEED[1] = 1;
      ISEED[2] = 1;
      ISEED[3] = 1;
      ISEED[4] = 1;
      for (M = 32; M <= MAXM; M += 51) {
        KLA = 0;
        KUA = M - 1;
        zlatmr(
            M,
            M,
            'S',
            ISEED,
            'N',
            D,
            6,
            ONE,
            Complex.one,
            'T',
            'N',
            DUML,
            1,
            ONE,
            DUMR,
            1,
            ONE,
            'N',
            IWORK,
            KLA,
            KUA,
            ZERO,
            ONE,
            'NO',
            A,
            MAXM,
            IWORK,
            IINFO);
        for (I = 1; I <= M; I++) {
          A[I][I] = A[I][I] * VM[J].toComplex();
        }
        ANRM = zlange('M', M, M, A, MAXM, DUM);
        for (N = 51; N <= MAXN; N += 47) {
          KLB = 0;
          KUB = N - 1;
          zlatmr(
              N,
              N,
              'S',
              ISEED,
              'N',
              D,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              DUML,
              1,
              ONE,
              DUMR,
              1,
              ONE,
              'N',
              IWORK,
              KLB,
              KUB,
              ZERO,
              ONE,
              'NO',
              B,
              MAXN,
              IWORK,
              IINFO);
          for (I = 1; I <= N; I++) {
            B[I][I] = B[I][I] * VM[J].toComplex();
          }
          BNRM = zlange('M', N, N, B, MAXN, DUM);
          TNRM = max(ANRM, BNRM);
          zlatmr(
              M,
              N,
              'S',
              ISEED,
              'N',
              D,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              DUML,
              1,
              ONE,
              DUMR,
              1,
              ONE,
              'N',
              IWORK,
              M,
              N,
              ZERO,
              ONE,
              'NO',
              C,
              MAXM,
              IWORK,
              IINFO);
          for (ITRANA = 1; ITRANA <= 2; ITRANA++) {
            if (ITRANA == 1) TRANA = 'N';
            if (ITRANA == 2) TRANA = 'C';
            for (ITRANB = 1; ITRANB <= 2; ITRANB++) {
              if (ITRANB == 1) TRANB = 'N';
              if (ITRANB == 2) TRANB = 'C';
              KNT.value++;

              zlacpy('All', M, N, C, MAXM, X, MAXM);
              zlacpy('All', M, N, C, MAXM, CC, MAXM);
              ztrsyl(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM, SCALE,
                  IINFO);
              if (IINFO.value != 0) NINFO[1]++;
              XNRM = zlange('M', M, N, X, MAXM, DUM);
              RMUL = Complex.one;
              if (XNRM > ONE && TNRM > ONE) {
                if (XNRM > BIGNUM / TNRM) {
                  RMUL = Complex.one / max(XNRM, TNRM).toComplex();
                }
              }
              zgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM,
                  -SCALE.value.toComplex() * RMUL, CC, MAXM);
              zgemm('N', TRANB, M, N, N, ISGN.toComplex() * RMUL, X, MAXM, B,
                  MAXN, Complex.one, CC, MAXM);
              RES1 = zlange('M', M, N, CC, MAXM, DUM);
              RES = RES1 /
                  max(SMLNUM,
                      max(SMLNUM * XNRM, (((RMUL).abs() * TNRM) * EPS) * XNRM));
              if (RES > THRESH) NFAIL[1]++;
              if (RES > RMAX[1]) RMAX[1] = RES;

              zlacpy('All', M, N, C, MAXM, X, MAXM);
              zlacpy('All', M, N, C, MAXM, CC, MAXM);
              ztrsyl3(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM,
                  SCALE3, SWORK, Box(LDSWORK), INFO);
              if (INFO.value != 0) NINFO[2]++;
              XNRM = zlange('M', M, N, X, MAXM, DUM);
              RMUL = Complex.one;
              if (XNRM > ONE && TNRM > ONE) {
                if (XNRM > BIGNUM / TNRM) {
                  RMUL = Complex.one / max(XNRM, TNRM).toComplex();
                }
              }
              zgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM,
                  -SCALE3.value.toComplex() * RMUL, CC, MAXM);
              zgemm('N', TRANB, M, N, N, ISGN.toComplex() * RMUL, X, MAXM, B,
                  MAXN, Complex.one, CC, MAXM);
              RES1 = zlange('M', M, N, CC, MAXM, DUM);
              RES = RES1 /
                  max(SMLNUM,
                      max(SMLNUM * XNRM, (((RMUL).abs() * TNRM) * EPS) * XNRM));
              // Verify that TRSYL3 only flushes if TRSYL flushes (but
              // there may be cases where TRSYL3 avoid flushing).
              if (SCALE3.value == ZERO && SCALE.value > ZERO ||
                  IINFO.value != INFO.value) {
                NFAIL[3]++;
              }
              if (RES > THRESH || disnan(RES)) NFAIL[2]++;
              if (RES > RMAX[2]) RMAX[2] = RES;
            }
          }
        }
      }
    }
  }
}
