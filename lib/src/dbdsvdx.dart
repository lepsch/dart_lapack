import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dstevx.dart';
import 'package:lapack/src/f2c/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dbdsvdx(
  final String UPLO,
  final String JOBZ,
  final String RANGE,
  final int N,
  final Array<double> D,
  final Array<double> E,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final Box<int> NS,
  final Array<double> S,
  final Matrix<double> Z,
  final int LDZ,
  final Array<double> WORK,
  final Array<int> IWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0, HNDRD = 100.0, MEIGTH = -0.1250;
  const FUDGE = 2.0;
  String RNGVX = '';
  bool ALLSV, INDSV, LOWER, SPLIT, SVEQ0, VALSV, WANTZ;
  int I,
      ICOLZ,
      IDBEG,
      IDEND,
      IDTGK,
      IDPTR,
      IEPTR,
      IETGK,
      IIFAIL,
      IIWORK,
      ILTGK,
      IROWU,
      IROWV,
      IROWZ,
      ISBEG,
      ISPLT,
      ITEMP,
      IUTGK,
      J,
      K,
      NTGK = 0,
      NRU,
      NRV;
  double ABSTOL,
      EPS,
      EMIN,
      MU,
      NRMU,
      NRMV,
      ORTOL,
      SMAX,
      SMIN,
      SQRT2,
      THRESH,
      TOL,
      ULP,
      VLTGK,
      VUTGK,
      ZJTJI;
  final NSL = Box(0);

  // Test the input parameters.

  ALLSV = lsame(RANGE, 'A');
  VALSV = lsame(RANGE, 'V');
  INDSV = lsame(RANGE, 'I');
  WANTZ = lsame(JOBZ, 'V');
  LOWER = lsame(UPLO, 'L');

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !LOWER) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(ALLSV || VALSV || INDSV)) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (N > 0) {
    if (VALSV) {
      if (VL < ZERO) {
        INFO.value = -7;
      } else if (VU <= VL) {
        INFO.value = -8;
      }
    } else if (INDSV) {
      if (IL < 1 || IL > max(1, N)) {
        INFO.value = -9;
      } else if (IU < min(N, IL) || IU > N) {
        INFO.value = -10;
      }
    }
  }
  if (INFO.value == 0) {
    if (LDZ < 1 || (WANTZ && LDZ < N * 2)) INFO.value = -14;
  }

  if (INFO.value != 0) {
    xerbla('DBDSVDX', -INFO.value);
    return;
  }

  // Quick return if possible (N <= 1)

  NS.value = 0;
  if (N == 0) return;

  if (N == 1) {
    if (ALLSV || INDSV) {
      NS.value = 1;
      S[1] = (D[1]).abs();
    } else {
      if (VL < (D[1]).abs() && VU >= (D[1]).abs()) {
        NS.value = 1;
        S[1] = (D[1]).abs();
      }
    }
    if (WANTZ) {
      Z[1][1] = sign(ONE, D[1]).toDouble();
      Z[2][1] = ONE;
    }
    return;
  }

  ABSTOL = 2 * dlamch('Safe Minimum');
  ULP = dlamch('Precision');
  EPS = dlamch('Epsilon');
  SQRT2 = sqrt(2.0);
  ORTOL = sqrt(ULP);

  // Criterion for splitting is taken from DBDSQR when singular
  // values are computed to relative accuracy TOL. (See J. Demmel and
  // W. Kahan, Accurate singular values of bidiagonal matrices, SIAM
  // J. Sci. and Stat. Comput., 11:873–912, 1990.)

  TOL = max(TEN, min(HNDRD, pow(EPS, MEIGTH))) * EPS;

  // Compute approximate maximum, minimum singular values.

  I = idamax(N, D, 1);
  SMAX = (D[I]).abs();
  I = idamax(N - 1, E, 1);
  SMAX = max(SMAX, (E[I]).abs());

  // Compute threshold for neglecting D's and E's.

  SMIN = (D[1]).abs();
  if (SMIN != ZERO) {
    MU = SMIN;
    for (I = 2; I <= N; I++) {
      MU = (D[I]).abs() * (MU / (MU + (E[I - 1]).abs()));
      SMIN = min(SMIN, MU);
      if (SMIN == ZERO) break;
    }
  }
  SMIN = SMIN / sqrt(N.toDouble());
  THRESH = TOL * SMIN;

  // Check for zeros in D and E (splits), i.e. submatrices.

  for (I = 1; I <= N - 1; I++) {
    if ((D[I]).abs() <= THRESH) D[I] = ZERO;
    if ((E[I]).abs() <= THRESH) E[I] = ZERO;
  }
  if ((D[N]).abs() <= THRESH) D[N] = ZERO;

  // Pointers for arrays used by DSTEVX.

  IDTGK = 1;
  IETGK = IDTGK + N * 2;
  ITEMP = IETGK + N * 2;
  IIFAIL = 1;
  IIWORK = IIFAIL + N * 2;

  // Set RNGVX, which corresponds to RANGE for DSTEVX in TGK mode.
  // VL,VU or IL,IU are redefined to conform to implementation a)
  // described in the leading comments.

  ILTGK = 0;
  IUTGK = 0;
  VLTGK = ZERO;
  VUTGK = ZERO;

  if (ALLSV) {
    // All singular values will be found. We aim at -s (see
    // leading comments) with RNGVX = 'I'. IL and IU are set
    // later (as ILTGK and IUTGK) according to the dimension
    // of the active submatrix.

    RNGVX = 'I';
    if (WANTZ) dlaset('F', N * 2, N + 1, ZERO, ZERO, Z, LDZ);
  } else if (VALSV) {
    // Find singular values in a half-open interval. We aim
    // at -s (see leading comments) and we swap VL and VU
    // (as VUTGK and VLTGK), changing their signs.

    RNGVX = 'V';
    VLTGK = -VU;
    VUTGK = -VL;
    for (var i = IDTGK; i <= IDTGK + 2 * N - 1; i++) {
      WORK[i] = ZERO;
    }
    dcopy(N, D, 1, WORK(IETGK), 2);
    dcopy(N - 1, E, 1, WORK(IETGK + 1), 2);
    dstevx(
      'N',
      'V',
      N * 2,
      WORK(IDTGK),
      WORK(IETGK),
      VLTGK,
      VUTGK,
      ILTGK,
      ILTGK,
      ABSTOL,
      NS,
      S,
      Z,
      LDZ,
      WORK(ITEMP),
      IWORK(IIWORK),
      IWORK(IIFAIL),
      INFO,
    );
    if (NS.value == 0) {
      return;
    } else {
      if (WANTZ) dlaset('F', N * 2, NS.value, ZERO, ZERO, Z, LDZ);
    }
  } else if (INDSV) {
    // Find the IL-th through the IU-th singular values. We aim
    // at -s (see leading comments) and indices are mapped into
    // values, therefore mimicking DSTEBZ, where

    // GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
    // GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN

    ILTGK = IL;
    IUTGK = IU;
    RNGVX = 'V';
    for (var i = IDTGK; i <= IDTGK + 2 * N - 1; i++) {
      WORK[i] = ZERO;
    }
    dcopy(N, D, 1, WORK(IETGK), 2);
    dcopy(N - 1, E, 1, WORK(IETGK + 1), 2);
    dstevx(
      'N',
      'I',
      N * 2,
      WORK(IDTGK),
      WORK(IETGK),
      VLTGK,
      VLTGK,
      ILTGK,
      ILTGK,
      ABSTOL,
      NS,
      S,
      Z,
      LDZ,
      WORK(ITEMP),
      IWORK(IIWORK),
      IWORK(IIFAIL),
      INFO,
    );
    VLTGK = S[1] - FUDGE * SMAX * ULP * N;
    for (var i = IDTGK; i <= IDTGK + 2 * N - 1; i++) {
      WORK[i] = ZERO;
    }
    dcopy(N, D, 1, WORK(IETGK), 2);
    dcopy(N - 1, E, 1, WORK(IETGK + 1), 2);
    dstevx(
      'N',
      'I',
      N * 2,
      WORK(IDTGK),
      WORK(IETGK),
      VUTGK,
      VUTGK,
      IUTGK,
      IUTGK,
      ABSTOL,
      NS,
      S,
      Z,
      LDZ,
      WORK(ITEMP),
      IWORK(IIWORK),
      IWORK(IIFAIL),
      INFO,
    );
    VUTGK = S[1] + FUDGE * SMAX * ULP * N;
    VUTGK = min(VUTGK, ZERO);

    // If VLTGK=VUTGK, DSTEVX returns an error message,
    // so if needed we change VUTGK slightly.

    if (VLTGK == VUTGK) VLTGK = VLTGK - TOL;

    if (WANTZ) dlaset('F', N * 2, IU - IL + 1, ZERO, ZERO, Z, LDZ);
  }

  // Initialize variables and pointers for S, Z, and WORK.

  // NRU, NRV: number of rows in U and V for the active submatrix
  // IDBEG, ISBEG: offsets for the entries of D and S
  // IROWZ, ICOLZ: offsets for the rows and columns of Z
  // IROWU, IROWV: offsets for the rows of U and V

  NS.value = 0;
  NRU = 0;
  NRV = 0;
  IDBEG = 1;
  ISBEG = 1;
  IROWZ = 1;
  ICOLZ = 1;
  IROWU = 2;
  IROWV = 1;
  SPLIT = false;
  SVEQ0 = false;

  // Form the tridiagonal TGK matrix.
  for (var i = 1; i <= N; i++) {
    S[i] = ZERO;
  }
  WORK[IETGK + 2 * N - 1] = ZERO;
  for (var i = IDTGK; i <= IDTGK + 2 * N - 1; i++) {
    WORK[i] = ZERO;
  }
  dcopy(N, D, 1, WORK(IETGK), 2);
  dcopy(N - 1, E, 1, WORK(IETGK + 1), 2);

  // Check for splits in two levels, outer level
  // in E and inner level in D.

  for (IEPTR = 2; 2 < 0 ? IEPTR >= N * 2 : IEPTR <= N * 2; IEPTR += 2) {
    if (WORK[IETGK + IEPTR - 1] == ZERO) {
      // Split in E (this piece of B is square) or bottom
      // of the (input bidiagonal) matrix.

      ISPLT = IDBEG;
      IDEND = IEPTR - 1;
      for (IDPTR = IDBEG; IDPTR <= IDEND; IDPTR += 2) {
        if (WORK[IETGK + IDPTR - 1] == ZERO) {
          // Split in D (rectangular submatrix). Set the number
          // of rows in U and V (NRU and NRV) accordingly.

          if (IDPTR == IDBEG) {
            // D=0 at the top.

            SVEQ0 = true;
            if (IDBEG == IDEND) {
              NRU = 1;
              NRV = 1;
            }
          } else if (IDPTR == IDEND) {
            // D=0 at the bottom.

            SVEQ0 = true;
            NRU = (IDEND - ISPLT) ~/ 2 + 1;
            NRV = NRU;
            if (ISPLT != IDBEG) {
              NRU = NRU + 1;
            }
          } else {
            if (ISPLT == IDBEG) {
              // Split: top rectangular submatrix.

              NRU = (IDPTR - IDBEG) ~/ 2;
              NRV = NRU + 1;
            } else {
              // Split: middle square submatrix.

              NRU = (IDPTR - ISPLT) ~/ 2 + 1;
              NRV = NRU;
            }
          }
        } else if (IDPTR == IDEND) {
          // Last entry of D in the active submatrix.

          if (ISPLT == IDBEG) {
            // No split (trivial case).

            NRU = (IDEND - IDBEG) ~/ 2 + 1;
            NRV = NRU;
          } else {
            // Split: bottom rectangular submatrix.

            NRV = (IDEND - ISPLT) ~/ 2 + 1;
            NRU = NRV + 1;
          }
        }

        NTGK = NRU + NRV;

        if (NTGK > 0) {
          // Compute eigenvalues/vectors of the active
          // submatrix according to RANGE:
          // if RANGE='A' (ALLSV) then RNGVX = 'I'
          // if RANGE='V' (VALSV) then RNGVX = 'V'
          // if RANGE='I' (INDSV) then RNGVX = 'V'

          ILTGK = 1;
          IUTGK = NTGK ~/ 2;
          if (ALLSV || VUTGK == ZERO) {
            if (SVEQ0 || SMIN < EPS || (NTGK % 2) > 0) {
              // Special case: eigenvalue equal to zero or very
              // small, additional eigenvector is needed.
              IUTGK = IUTGK + 1;
            }
          }

          // Workspace needed by DSTEVX:
          // WORK[ ITEMP: ]: 2*5*NTGK
          // IWORK[ 1: ]: 2*6*NTGK

          dstevx(
            JOBZ,
            RNGVX,
            NTGK,
            WORK(IDTGK + ISPLT - 1),
            WORK(IETGK + ISPLT - 1),
            VLTGK,
            VUTGK,
            ILTGK,
            IUTGK,
            ABSTOL,
            NSL,
            S(ISBEG),
            Z(IROWZ, ICOLZ),
            LDZ,
            WORK(ITEMP),
            IWORK(IIWORK),
            IWORK(IIFAIL),
            INFO,
          );
          if (INFO.value != 0) {
            // Exit with the error code from DSTEVX.
            return;
          }
          EMIN = S.maxval(ISBEG, ISBEG + NSL.value - 1).abs();

          if (NSL.value > 0 && WANTZ) {
            // Normalize u=Z[[2,4,...]][:] and v=Z[[1,3,...]][:],
            // changing the sign of v as discussed in the leading
            // comments. The norms of u and v may be (slightly)
            // different from 1/sqrt(2) if the corresponding
            // eigenvalues are very small or too close. We check
            // those norms and, if needed, reorthogonalize the
            // vectors.

            if (NSL.value > 1 &&
                VUTGK == ZERO &&
                (NTGK % 2) == 0 &&
                EMIN == 0 &&
                !SPLIT) {
              // D=0 at the top or bottom of the active submatrix:
              // one eigenvalue is equal to zero; concatenate the
              // eigenvectors corresponding to the two smallest
              // eigenvalues.

              for (var i = IROWZ; i <= IROWZ + NTGK - 1; i++) {
                Z[i][ICOLZ + NSL.value - 2] =
                    Z[i][ICOLZ + NSL.value - 2] + Z[i][ICOLZ + NSL.value - 1];
                Z[i][ICOLZ + NSL.value - 1] = ZERO;
              }

              // if( IUTGK*2 > NTGK ) THEN
              // Eigenvalue equal to zero or very small.
              // NSL.value = NSL.value - 1
              // END if
            }

            for (I = 0; I <= min(NSL.value - 1, NRU - 1); I++) {
              NRMU = dnrm2(NRU, Z(IROWU, ICOLZ + I).asArray(), 2);
              if (NRMU == ZERO) {
                INFO.value = N * 2 + 1;
                return;
              }
              dscal(NRU, ONE / NRMU, Z(IROWU, ICOLZ + I).asArray(), 2);
              if (NRMU != ONE && (NRMU - ORTOL).abs() * SQRT2 > ONE) {
                for (J = 0; J <= I - 1; J++) {
                  ZJTJI = -ddot(
                    NRU,
                    Z(IROWU, ICOLZ + J).asArray(),
                    2,
                    Z(IROWU, ICOLZ + I).asArray(),
                    2,
                  );
                  daxpy(
                    NRU,
                    ZJTJI,
                    Z(IROWU, ICOLZ + J).asArray(),
                    2,
                    Z(IROWU, ICOLZ + I).asArray(),
                    2,
                  );
                }
                NRMU = dnrm2(NRU, Z(IROWU, ICOLZ + I).asArray(), 2);
                dscal(NRU, ONE / NRMU, Z(IROWU, ICOLZ + I).asArray(), 2);
              }
            }
            for (I = 0; I <= min(NSL.value - 1, NRV - 1); I++) {
              NRMV = dnrm2(NRV, Z(IROWV, ICOLZ + I).asArray(), 2);
              if (NRMV == ZERO) {
                INFO.value = N * 2 + 1;
                return;
              }
              dscal(NRV, -ONE / NRMV, Z(IROWV, ICOLZ + I).asArray(), 2);
              if (NRMV != ONE && (NRMV - ORTOL).abs() * SQRT2 > ONE) {
                for (J = 0; J <= I - 1; J++) {
                  ZJTJI = -ddot(
                    NRV,
                    Z(IROWV, ICOLZ + J).asArray(),
                    2,
                    Z(IROWV, ICOLZ + I).asArray(),
                    2,
                  );
                  daxpy(
                    NRU,
                    ZJTJI,
                    Z(IROWV, ICOLZ + J).asArray(),
                    2,
                    Z(IROWV, ICOLZ + I).asArray(),
                    2,
                  );
                }
                NRMV = dnrm2(NRV, Z(IROWV, ICOLZ + I).asArray(), 2);
                dscal(NRV, ONE / NRMV, Z(IROWV, ICOLZ + I).asArray(), 2);
              }
            }
            if (VUTGK == ZERO && IDPTR < IDEND && (NTGK % 2) > 0) {
              // D=0 in the middle of the active submatrix (one
              // eigenvalue is equal to zero): save the corresponding
              // eigenvector for later use (when bottom of the
              // active submatrix is reached).

              SPLIT = true;

              for (var i = IROWZ; i <= IROWZ + NTGK - 1; i++) {
                Z[i][N + 1] = Z[i][NS.value + NSL.value];
                Z[i][NS.value + NSL.value] = ZERO;
              }
            }
          } // !** WANTZ **!;

          NSL.value = min(NSL.value, NRU);
          SVEQ0 = false;

          // Absolute values of the eigenvalues of TGK.

          for (I = 0; I <= NSL.value - 1; I++) {
            S[ISBEG + I] = (S[ISBEG + I]).abs();
          }

          // Update pointers for TGK, S and Z.

          ISBEG = ISBEG + NSL.value;
          IROWZ = IROWZ + NTGK;
          ICOLZ = ICOLZ + NSL.value;
          IROWU = IROWZ;
          IROWV = IROWZ + 1;
          ISPLT = IDPTR + 1;
          NS.value = NS.value + NSL.value;
          NRU = 0;
          NRV = 0;
        } // !** NTGK > 0 **!;
        if (IROWZ < N * 2 && WANTZ) {
          for (var i = 1; i <= IROWZ - 1; i++) {
            Z[i][ICOLZ] = ZERO;
          }
        }
      } // !** IDPTR loop **!;
      if (SPLIT && WANTZ) {
        // Bring back eigenvector corresponding
        // to eigenvalue equal to zero.

        for (var i = IDBEG; i <= IDEND - NTGK + 1; i++) {
          Z[i][ISBEG - 1] = Z[i][ISBEG - 1] + Z[i][N + 1];
          Z[i][N + 1] = ZERO;
        }
      }
      IROWV = IROWV - 1;
      IROWU = IROWU + 1;
      IDBEG = IEPTR + 1;
      SVEQ0 = false;
      SPLIT = false;
    } // !** Check for split in E **!;
  } // !** IEPTR loop **!;

  // Sort the singular values into decreasing order (insertion sort on
  // singular values, but only one transposition per singular vector)

  for (I = 1; I <= NS.value - 1; I++) {
    K = 1;
    SMIN = S[1];
    for (J = 2; J <= NS.value + 1 - I; J++) {
      if (S[J] <= SMIN) {
        K = J;
        SMIN = S[J];
      }
    }
    if (K != NS.value + 1 - I) {
      S[K] = S[NS.value + 1 - I];
      S[NS.value + 1 - I] = SMIN;
      if (WANTZ) {
        dswap(N * 2, Z(1, K).asArray(), 1, Z(1, NS.value + 1 - I).asArray(), 1);
      }
    }
  }

  // If RANGE=I, check for singular values/vectors to be discarded.

  if (INDSV) {
    K = IU - IL + 1;
    if (K < NS.value) {
      for (var i = K + 1; i <= NS.value; i++) {
        S[i] = ZERO;
      }

      if (WANTZ) {
        for (var i = 1; i <= N * 2; i++) {
          for (var j = K + 1; K <= NS.value; j++) {
            Z[i][j] = ZERO;
          }
        }
      }
      NS.value = K;
    }
  }

  // Reorder Z: U = Z[ 1:N][1:NS.value ], V = Z[ N+1:N*2][1:NS.value ].
  // If B is a lower diagonal, swap U and V.

  if (WANTZ) {
    for (I = 1; I <= NS.value; I++) {
      dcopy(N * 2, Z(1, I).asArray(), 1, WORK, 1);
      if (LOWER) {
        dcopy(N, WORK(2), 2, Z(N + 1, I).asArray(), 1);
        dcopy(N, WORK(1), 2, Z(1, I).asArray(), 1);
      } else {
        dcopy(N, WORK(2), 2, Z(1, I).asArray(), 1);
        dcopy(N, WORK(1), 2, Z(N + 1, I).asArray(), 1);
      }
    }
  }
}
