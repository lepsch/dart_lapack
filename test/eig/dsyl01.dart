import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dtrsyl.dart';
import 'package:lapack/src/dtrsyl3.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import '../matgen/dlatmr.dart';

void dsyl01(
  final double THRESH,
  final Array<int> NFAIL,
  final Array<double> RMAX,
  final Array<int> NINFO,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const MAXM = 245, MAXN = 192, LDSWORK = 36;
  String TRANA = '', TRANB = '';
  int I, ISGN, ITRANA, ITRANB, J, KLA, KUA, KLB, KUB, LIWORK, M, N;
  double ANRM, BNRM, BIGNUM, EPS, RES, RES1, RMUL, SMLNUM, TNRM, XNRM;
  final DUML = Array<double>(MAXM),
      DUMR = Array<double>(MAXN),
      D = Array<double>(max(MAXM, MAXN)),
      DUM = Array<double>(MAXN),
      VM = Array<double>(2);
  final ISEED = Array<int>(4), IWORK = Array<int>(MAXM + MAXN + 2);
  final A = Matrix<double>(MAXM, MAXM),
      B = Matrix<double>(MAXN, MAXN),
      C = Matrix<double>(MAXM, MAXN),
      CC = Matrix<double>(MAXM, MAXN),
      X = Matrix<double>(MAXM, MAXN),
      SWORK = Matrix<double>(LDSWORK, 126);
  final INFO = Box(0), IINFO = Box(0);
  final SCALE = Box(0.0), SCALE3 = Box(0.0);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  VM[1] = ONE;
  VM[2] = 0.000001;

  // Begin test loop

  NINFO[1] = 0;
  NINFO[2] = 0;
  NFAIL[1] = 0;
  NFAIL[2] = 0;
  NFAIL[3] = 0;
  RMAX[1] = ZERO;
  RMAX[2] = ZERO;
  KNT.value = 0;
  for (I = 1; I <= 4; I++) {
    ISEED[I] = 1;
  }
  SCALE.value = ONE;
  SCALE3.value = ONE;
  LIWORK = MAXM + MAXN + 2;
  for (J = 1; J <= 2; J++) {
    for (ISGN = -1; 2 < 0 ? ISGN >= 1 : ISGN <= 1; ISGN += 2) {
      // Reset seed (overwritten by LATMR)
      for (I = 1; I <= 4; I++) {
        ISEED[I] = 1;
      }
      for (M = 32; M <= MAXM; M += 71) {
        KLA = 0;
        KUA = M - 1;
        dlatmr(
            M,
            M,
            'S',
            ISEED,
            'N',
            D,
            6,
            ONE,
            ONE,
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
          A[I][I] = A[I][I] * VM[J];
        }
        ANRM = dlange('M', M, M, A, MAXM, DUM);
        for (N = 51; N <= MAXN; N += 47) {
          KLB = 0;
          KUB = N - 1;
          dlatmr(
              N,
              N,
              'S',
              ISEED,
              'N',
              D,
              6,
              ONE,
              ONE,
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
          BNRM = dlange('M', N, N, B, MAXN, DUM);
          TNRM = max(ANRM, BNRM);
          dlatmr(
              M,
              N,
              'S',
              ISEED,
              'N',
              D,
              6,
              ONE,
              ONE,
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
            if (ITRANA == 1) {
              TRANA = 'N';
            }
            if (ITRANA == 2) {
              TRANA = 'T';
            }
            for (ITRANB = 1; ITRANB <= 2; ITRANB++) {
              if (ITRANB == 1) {
                TRANB = 'N';
              }
              if (ITRANB == 2) {
                TRANB = 'T';
              }
              KNT.value = KNT.value + 1;

              dlacpy('All', M, N, C, MAXM, X, MAXM);
              dlacpy('All', M, N, C, MAXM, CC, MAXM);
              dtrsyl(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM, SCALE,
                  IINFO);
              if (IINFO.value != 0) NINFO[1] = NINFO[1] + 1;
              XNRM = dlange('M', M, N, X, MAXM, DUM);
              RMUL = ONE;
              if (XNRM > ONE && TNRM > ONE) {
                if (XNRM > BIGNUM / TNRM) {
                  RMUL = ONE / max(XNRM, TNRM);
                }
              }
              dgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM,
                  -SCALE.value * RMUL, CC, MAXM);
              dgemm('N', TRANB, M, N, N, ISGN.toDouble() * RMUL, X, MAXM, B,
                  MAXN, ONE, CC, MAXM);
              RES1 = dlange('M', M, N, CC, MAXM, DUM);
              RES = RES1 /
                  max(SMLNUM, max(SMLNUM * XNRM, ((RMUL * TNRM) * EPS) * XNRM));
              if (RES > THRESH) NFAIL[1] = NFAIL[1] + 1;
              if (RES > RMAX[1]) RMAX[1] = RES;

              dlacpy('All', M, N, C, MAXM, X, MAXM);
              dlacpy('All', M, N, C, MAXM, CC, MAXM);
              final LDSWORKINOUT = Box(LDSWORK);
              dtrsyl3(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM,
                  SCALE3, IWORK, LIWORK, SWORK, LDSWORKINOUT, INFO);
              if (INFO.value != 0) NINFO[2] = NINFO[2] + 1;
              XNRM = dlange('M', M, N, X, MAXM, DUM);
              RMUL = ONE;
              if (XNRM > ONE && TNRM > ONE) {
                if (XNRM > BIGNUM / TNRM) {
                  RMUL = ONE / max(XNRM, TNRM);
                }
              }
              dgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM,
                  -SCALE3.value * RMUL, CC, MAXM);
              dgemm('N', TRANB, M, N, N, ISGN.toDouble() * RMUL, X, MAXM, B,
                  MAXN, ONE, CC, MAXM);
              RES1 = dlange('M', M, N, CC, MAXM, DUM);
              RES = RES1 /
                  max(SMLNUM, max(SMLNUM * XNRM, ((RMUL * TNRM) * EPS) * XNRM));
              // Verify that TRSYL3 only flushes if TRSYL flushes (but
              // there may be cases where TRSYL3 avoid flushing).
              if (SCALE3.value == ZERO && SCALE.value > ZERO ||
                  IINFO.value != INFO.value) {
                NFAIL[3] = NFAIL[3] + 1;
              }
              if (RES > THRESH || disnan(RES)) NFAIL[2] = NFAIL[2] + 1;
              if (RES > RMAX[2]) RMAX[2] = RES;
            }
          }
        }
      }
    }
  }
}
