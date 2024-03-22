import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/nint.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/qr/ll/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgelqf.dart';
import 'package:lapack/src/zgeqp3.dart';
import 'package:lapack/src/zgesvj.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlapmr.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zlassq.dart';
import 'package:lapack/src/zlaswp.dart';
import 'package:lapack/src/zpocon.dart';
import 'package:lapack/src/zungqr.dart';
import 'package:lapack/src/zunmlq.dart';
import 'package:lapack/src/zunmqr.dart';

void zgejsv(
  final String JOBA,
  final String JOBU,
  final String JOBV,
  final String JOBR,
  final String JOBT,
  final String JOBP,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> SVA_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> CWORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final SVA = SVA_.having(length: N);
  final CWORK = CWORK_.having();
  final RWORK = RWORK_.having(length: LRWORK);
  final IWORK = IWORK_.having();

  const ZERO = 0.0, ONE = 1.0;
  Complex CTEMP;
  double AATMAX,
      AATMIN,
      BIG,
      BIG1,
      COND_OK,
      CONDR1,
      CONDR2,
      ENTRA,
      ENTRAT,
      EPSLN,
      MAXPRJ,
      SCALEM,
      SCONDA,
      SFMIN,
      SMALL,
      USCAL1,
      USCAL2;
  int N1 = 0, NR, NUMRANK, p, q, WARNING;
  bool ALMORT,
      DEFR,
      ERREST,
      GOSCAL,
      JRACC,
      KILL,
      LQUERY,
      LSVEC,
      L2ABER,
      L2KILL,
      L2PERT,
      L2RANK,
      L2TRAN,
      NOSCAL,
      ROWPIV,
      RSVEC,
      TRANSP;

  int OPTWRK = 0, MINWRK = 0, MINRWRK = 0, MINIWRK = 0;
  int LWCON,
      LWLQF,
      LWQP3,
      LWQRF,
      LWUNMLQ,
      LWUNMQR,
      LWUNMQRM,
      LWSVDJ,
      LWSVDJV,
      LRWQP3,
      LRWCON,
      LRWSVDJ,
      IWOFF = 0;
  int LWRK_ZGELQF = 0,
      LWRK_ZGEQP3 = 0,
      LWRK_ZGEQP3N,
      LWRK_ZGEQRF = 0,
      LWRK_ZGESVJ = 0,
      LWRK_ZGESVJV,
      LWRK_ZGESVJU,
      LWRK_ZUNMLQ,
      LWRK_ZUNMQR,
      LWRK_ZUNMQRM;
  final CDUMMY = Array<Complex>(1);
  final RDUMMY = Array<double>(1);
  final IERR = Box(0);
  final AAPP = Box(0.0), AAQQ = Box(0.0), XSC = Box(0.0), TEMP1 = Box(0.0);

  // Test the input arguments

  LSVEC = lsame(JOBU, 'U') || lsame(JOBU, 'F');
  JRACC = lsame(JOBV, 'J');
  RSVEC = lsame(JOBV, 'V') || JRACC;
  ROWPIV = lsame(JOBA, 'F') || lsame(JOBA, 'G');
  L2RANK = lsame(JOBA, 'R');
  L2ABER = lsame(JOBA, 'A');
  ERREST = lsame(JOBA, 'E') || lsame(JOBA, 'G');
  L2TRAN = lsame(JOBT, 'T') && (M == N);
  L2KILL = lsame(JOBR, 'R');
  DEFR = lsame(JOBR, 'N');
  L2PERT = lsame(JOBP, 'P');

  LQUERY = (LWORK == -1) || (LRWORK == -1);

  if (!(ROWPIV || L2RANK || L2ABER || ERREST || lsame(JOBA, 'C'))) {
    INFO.value = -1;
  } else if (!(LSVEC ||
      lsame(JOBU, 'N') ||
      (lsame(JOBU, 'W') && RSVEC && L2TRAN))) {
    INFO.value = -2;
  } else if (!(RSVEC ||
      lsame(JOBV, 'N') ||
      (lsame(JOBV, 'W') && LSVEC && L2TRAN))) {
    INFO.value = -3;
  } else if (!(L2KILL || DEFR)) {
    INFO.value = -4;
  } else if (!(lsame(JOBT, 'T') || lsame(JOBT, 'N'))) {
    INFO.value = -5;
  } else if (!(L2PERT || lsame(JOBP, 'N'))) {
    INFO.value = -6;
  } else if (M < 0) {
    INFO.value = -7;
  } else if ((N < 0) || (N > M)) {
    INFO.value = -8;
  } else if (LDA < M) {
    INFO.value = -10;
  } else if (LSVEC && (LDU < M)) {
    INFO.value = -13;
  } else if (RSVEC && (LDV < N)) {
    INFO.value = -15;
  } else {
    INFO.value = 0;
  }

  if (INFO.value == 0) {
    // .. compute the minimal and the optimal workspace lengths
    // [[The expressions for computing the minimal and the optimal
    // values of LCWORK, LRWORK are written with a lot of redundancy and
    // can be simplified. However, this verbose form is useful for
    // maintenance and modifications of the code.]]

    // .. minimal workspace length for ZGEQP3 of an M x N matrix,
    //  ZGEQRF of an N x N matrix, ZGELQF of an N x N matrix,
    //  ZUNMLQ for computing N x N matrix, ZUNMQR for computing N x N
    //  matrix, ZUNMQR for computing M x N matrix, respectively.
    LWQP3 = N + 1;
    LWQRF = max(1, N);
    LWLQF = max(1, N);
    LWUNMLQ = max(1, N);
    LWUNMQR = max(1, N);
    LWUNMQRM = max(1, M);
    // .. minimal workspace length for ZPOCON of an N x N matrix
    LWCON = 2 * N;
    // .. minimal workspace length for ZGESVJ of an N x N matrix,
    //  without and with explicit accumulation of Jacobi rotations
    LWSVDJ = max(2 * N, 1);
    LWSVDJV = max(2 * N, 1);
    // .. minimal REAL workspace length for ZGEQP3, ZPOCON, ZGESVJ
    LRWQP3 = 2 * N;
    LRWCON = N;
    LRWSVDJ = N;
    if (LQUERY) {
      zgeqp3(M, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1, RDUMMY, IERR);
      LWRK_ZGEQP3 = CDUMMY[1].toInt();
      zgeqrf(N, N, A, LDA, CDUMMY, CDUMMY, -1, IERR);
      LWRK_ZGEQRF = CDUMMY[1].toInt();
      zgelqf(N, N, A, LDA, CDUMMY, CDUMMY, -1, IERR);
      LWRK_ZGELQF = CDUMMY[1].toInt();
    }
    MINWRK = 2;
    OPTWRK = 2;
    MINIWRK = N;
    if (!(LSVEC || RSVEC)) {
      // .. minimal and optimal sizes of the complex workspace if
      // only the singular values are requested
      if (ERREST) {
        MINWRK = max(max(N + LWQP3, pow(N, 2).toInt().toInt() + LWCON),
            max(N + LWQRF, LWSVDJ));
      } else {
        MINWRK = max(N + LWQP3, max(N + LWQRF, LWSVDJ));
      }
      if (LQUERY) {
        zgesvj('L', 'N', 'N', N, N, A, LDA, SVA, N, V, LDV, CDUMMY, -1, RDUMMY,
            -1, IERR);
        LWRK_ZGESVJ = CDUMMY[1].toInt();
        if (ERREST) {
          OPTWRK = max(max(N + LWRK_ZGEQP3, pow(N, 2).toInt() + LWCON),
              max(N + LWRK_ZGEQRF, LWRK_ZGESVJ));
        } else {
          OPTWRK = max(N + LWRK_ZGEQP3, max(N + LWRK_ZGEQRF, LWRK_ZGESVJ));
        }
      }
      if (L2TRAN || ROWPIV) {
        if (ERREST) {
          MINRWRK = max(7, max(max(2 * M, LRWQP3), max(LRWCON, LRWSVDJ)));
        } else {
          MINRWRK = max(max(7, 2 * M), max(LRWQP3, LRWSVDJ));
        }
      } else {
        if (ERREST) {
          MINRWRK = max(max(7, LRWQP3), max(LRWCON, LRWSVDJ));
        } else {
          MINRWRK = max(7, max(LRWQP3, LRWSVDJ));
        }
      }
      if (ROWPIV || L2TRAN) MINIWRK += M;
    } else if (RSVEC && !LSVEC) {
      // .. minimal and optimal sizes of the complex workspace if the
      // singular values and the right singular vectors are requested
      if (ERREST) {
        MINWRK = max(max(N + LWQP3, max(LWCON, LWSVDJ)),
            max(max(N + LWLQF, 2 * N + LWQRF), max(N + LWSVDJ, N + LWUNMLQ)));
      } else {
        MINWRK = max(max(N + LWQP3, LWSVDJ),
            max(max(N + LWLQF, 2 * N + LWQRF), max(N + LWSVDJ, N + LWUNMLQ)));
      }
      if (LQUERY) {
        zgesvj('L', 'U', 'N', N, N, U, LDU, SVA, N, A, LDA, CDUMMY, -1, RDUMMY,
            -1, IERR);
        LWRK_ZGESVJ = CDUMMY[1].toInt();
        zunmlq('L', 'C', N, N, N, A, LDA, CDUMMY, V, LDV, CDUMMY, -1, IERR);
        LWRK_ZUNMLQ = CDUMMY[1].toInt();
        if (ERREST) {
          OPTWRK = max(
              max(N + LWRK_ZGEQP3, max(LWCON, LWRK_ZGESVJ)),
              max(max(N + LWRK_ZGELQF, 2 * N + LWRK_ZGEQRF),
                  max(N + LWRK_ZGESVJ, N + LWRK_ZUNMLQ)));
        } else {
          OPTWRK = max(
              max(N + LWRK_ZGEQP3, LWRK_ZGESVJ),
              max(max(N + LWRK_ZGELQF, 2 * N + LWRK_ZGEQRF),
                  max(N + LWRK_ZGESVJ, N + LWRK_ZUNMLQ)));
        }
      }
      if (L2TRAN || ROWPIV) {
        if (ERREST) {
          MINRWRK = max(7, max(max(2 * M, LRWQP3), max(LRWSVDJ, LRWCON)));
        } else {
          MINRWRK = max(max(7, 2 * M), max(LRWQP3, LRWSVDJ));
        }
      } else {
        if (ERREST) {
          MINRWRK = max(max(7, LRWQP3), max(LRWSVDJ, LRWCON));
        } else {
          MINRWRK = max(7, max(LRWQP3, LRWSVDJ));
        }
      }
      if (ROWPIV || L2TRAN) MINIWRK += M;
    } else if (LSVEC && !RSVEC) {
      // .. minimal and optimal sizes of the complex workspace if the
      // singular values and the left singular vectors are requested
      if (ERREST) {
        MINWRK =
            N + max(LWQP3, max(max(LWCON, N + LWQRF), max(LWSVDJ, LWUNMQRM)));
      } else {
        MINWRK = N + max(max(LWQP3, N + LWQRF), max(LWSVDJ, LWUNMQRM));
      }
      if (LQUERY) {
        zgesvj('L', 'U', 'N', N, N, U, LDU, SVA, N, A, LDA, CDUMMY, -1, RDUMMY,
            -1, IERR);
        LWRK_ZGESVJ = CDUMMY[1].toInt();
        zunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR);
        LWRK_ZUNMQRM = CDUMMY[1].toInt();
        if (ERREST) {
          OPTWRK = N +
              max(
                  LWRK_ZGEQP3,
                  max(max(LWCON, N + LWRK_ZGEQRF),
                      max(LWRK_ZGESVJ, LWRK_ZUNMQRM)));
        } else {
          OPTWRK = N +
              max(max(LWRK_ZGEQP3, N + LWRK_ZGEQRF),
                  max(LWRK_ZGESVJ, LWRK_ZUNMQRM));
        }
      }
      if (L2TRAN || ROWPIV) {
        if (ERREST) {
          MINRWRK = max(7, max(max(2 * M, LRWQP3), max(LRWSVDJ, LRWCON)));
        } else {
          MINRWRK = max(max(7, 2 * M), max(LRWQP3, LRWSVDJ));
        }
      } else {
        if (ERREST) {
          MINRWRK = max(max(7, LRWQP3), max(LRWSVDJ, LRWCON));
        } else {
          MINRWRK = max(7, max(LRWQP3, LRWSVDJ));
        }
      }
      if (ROWPIV || L2TRAN) MINIWRK += M;
    } else {
      // .. minimal and optimal sizes of the complex workspace if the
      // full SVD is requested
      if (!JRACC) {
        if (ERREST) {
          MINWRK = max(
              max(
                  N + LWQP3,
                  max(max(N + LWCON, 2 * N + pow(N, 2).toInt() + LWCON),
                      max(2 * N + LWQRF, 2 * N + LWQP3))),
              max(
                  max(
                      max(
                          2 * N + pow(N, 2).toInt() + N + LWLQF,
                          2 * N +
                              pow(N, 2).toInt() +
                              N +
                              pow(N, 2).toInt() +
                              LWCON),
                      max(2 * N + pow(N, 2).toInt() + N + LWSVDJ,
                          2 * N + pow(N, 2).toInt() + N + LWSVDJV)),
                  max(
                      max(2 * N + pow(N, 2).toInt() + N + LWUNMQR,
                          2 * N + pow(N, 2).toInt() + N + LWUNMLQ),
                      max(N + pow(N, 2).toInt() + LWSVDJ, N + LWUNMQRM))));
        } else {
          MINWRK = max(
              max(max(N + LWQP3, 2 * N + pow(N, 2).toInt() + LWCON),
                  max(2 * N + LWQRF, 2 * N + LWQP3)),
              max(
                  max(
                      max(
                          2 * N + pow(N, 2).toInt() + N + LWLQF,
                          2 * N +
                              pow(N, 2).toInt() +
                              N +
                              pow(N, 2).toInt() +
                              LWCON),
                      max(2 * N + pow(N, 2).toInt() + N + LWSVDJ,
                          2 * N + pow(N, 2).toInt() + N + LWSVDJV)),
                  max(
                      max(2 * N + pow(N, 2).toInt() + N + LWUNMQR,
                          2 * N + pow(N, 2).toInt() + N + LWUNMLQ),
                      max(N + pow(N, 2).toInt() + LWSVDJ, N + LWUNMQRM))));
        }
        MINIWRK += N;
        if (ROWPIV || L2TRAN) MINIWRK += M;
      } else {
        if (ERREST) {
          MINWRK = max(
              max(N + LWQP3, N + LWCON),
              max(max(2 * N + LWQRF, 2 * N + pow(N, 2).toInt() + LWSVDJV),
                  max(2 * N + pow(N, 2).toInt() + N + LWUNMQR, N + LWUNMQRM)));
        } else {
          MINWRK = max(
              N + LWQP3,
              max(max(2 * N + LWQRF, 2 * N + pow(N, 2).toInt() + LWSVDJV),
                  max(2 * N + pow(N, 2).toInt() + N + LWUNMQR, N + LWUNMQRM)));
        }
        if (ROWPIV || L2TRAN) MINIWRK += M;
      }
      if (LQUERY) {
        zunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR);
        LWRK_ZUNMQRM = CDUMMY[1].toInt();
        zunmqr('L', 'N', N, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR);
        LWRK_ZUNMQR = CDUMMY[1].toInt();
        if (!JRACC) {
          zgeqp3(N, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1, RDUMMY, IERR);
          LWRK_ZGEQP3N = CDUMMY[1].toInt();
          zgesvj('L', 'U', 'N', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1,
              RDUMMY, -1, IERR);
          LWRK_ZGESVJ = CDUMMY[1].toInt();
          zgesvj('U', 'U', 'N', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1,
              RDUMMY, -1, IERR);
          LWRK_ZGESVJU = CDUMMY[1].toInt();
          zgesvj('L', 'U', 'V', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1,
              RDUMMY, -1, IERR);
          LWRK_ZGESVJV = CDUMMY[1].toInt();
          zunmlq('L', 'C', N, N, N, A, LDA, CDUMMY, V, LDV, CDUMMY, -1, IERR);
          LWRK_ZUNMLQ = CDUMMY[1].toInt();
          if (ERREST) {
            OPTWRK = max(
                max(
                    N + LWRK_ZGEQP3,
                    max(max(N + LWCON, 2 * N + pow(N, 2).toInt() + LWCON),
                        max(2 * N + LWRK_ZGEQRF, 2 * N + LWRK_ZGEQP3N))),
                max(
                    max(
                        max(
                            2 * N + pow(N, 2).toInt() + N + LWRK_ZGELQF,
                            2 * N +
                                pow(N, 2).toInt() +
                                N +
                                pow(N, 2).toInt() +
                                LWCON),
                        max(2 * N + pow(N, 2).toInt() + N + LWRK_ZGESVJ,
                            2 * N + pow(N, 2).toInt() + N + LWRK_ZGESVJV)),
                    max(
                        max(2 * N + pow(N, 2).toInt() + N + LWRK_ZUNMQR,
                            2 * N + pow(N, 2).toInt() + N + LWRK_ZUNMLQ),
                        max(N + pow(N, 2).toInt() + LWRK_ZGESVJU,
                            N + LWRK_ZUNMQRM))));
          } else {
            OPTWRK = max(
                max(max(N + LWRK_ZGEQP3, 2 * N + pow(N, 2).toInt() + LWCON),
                    max(2 * N + LWRK_ZGEQRF, 2 * N + LWRK_ZGEQP3N)),
                max(
                    max(
                        max(
                            2 * N + pow(N, 2).toInt() + N + LWRK_ZGELQF,
                            2 * N +
                                pow(N, 2).toInt() +
                                N +
                                pow(N, 2).toInt() +
                                LWCON),
                        max(2 * N + pow(N, 2).toInt() + N + LWRK_ZGESVJ,
                            2 * N + pow(N, 2).toInt() + N + LWRK_ZGESVJV)),
                    max(
                        max(2 * N + pow(N, 2).toInt() + N + LWRK_ZUNMQR,
                            2 * N + pow(N, 2).toInt() + N + LWRK_ZUNMLQ),
                        max(N + pow(N, 2).toInt() + LWRK_ZGESVJU,
                            N + LWRK_ZUNMQRM))));
          }
        } else {
          zgesvj('L', 'U', 'V', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1,
              RDUMMY, -1, IERR);
          LWRK_ZGESVJV = CDUMMY[1].toInt();
          zunmqr('L', 'N', N, N, N, CDUMMY.asMatrix(N), N, CDUMMY, V, LDV,
              CDUMMY, -1, IERR);
          LWRK_ZUNMQR = CDUMMY[1].toInt();
          zunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR);
          LWRK_ZUNMQRM = CDUMMY[1].toInt();
          if (ERREST) {
            OPTWRK = max(
                max(N + LWRK_ZGEQP3, max(N + LWCON, 2 * N + LWRK_ZGEQRF)),
                max(
                    max(2 * N + pow(N, 2).toInt(),
                        2 * N + pow(N, 2).toInt() + LWRK_ZGESVJV),
                    max(2 * N + pow(N, 2).toInt() + N + LWRK_ZUNMQR,
                        N + LWRK_ZUNMQRM)));
          } else {
            OPTWRK = max(
                max(N + LWRK_ZGEQP3, 2 * N + LWRK_ZGEQRF),
                max(
                    max(2 * N + pow(N, 2).toInt(),
                        2 * N + pow(N, 2).toInt() + LWRK_ZGESVJV),
                    max(2 * N + pow(N, 2).toInt() + N + LWRK_ZUNMQR,
                        N + LWRK_ZUNMQRM)));
          }
        }
      }
      if (L2TRAN || ROWPIV) {
        MINRWRK = max(7, max(max(2 * M, LRWQP3), max(LRWSVDJ, LRWCON)));
      } else {
        MINRWRK = max(max(7, LRWQP3), max(LRWSVDJ, LRWCON));
      }
    }
    MINWRK = max(2, MINWRK);
    OPTWRK = max(MINWRK, OPTWRK);
    if (LWORK < MINWRK && !LQUERY) INFO.value = -17;
    if (LRWORK < MINRWRK && !LQUERY) INFO.value = -19;
  }

  if (INFO.value != 0) {
    xerbla('ZGEJSV', -INFO.value);
    return;
  } else if (LQUERY) {
    CWORK[1] = OPTWRK.toComplex();
    CWORK[2] = MINWRK.toComplex();
    RWORK[1] = MINRWRK.toDouble();
    IWORK[1] = max(4, MINIWRK);
    return;
  }

  // Quick return for void matrix (Y3K safe)
  if ((M == 0) || (N == 0)) {
    for (var i = 1; i <= 4; i++) {
      IWORK[i] = 0;
    }
    for (var i = 1; i <= 7; i++) {
      RWORK[i] = 0;
    }
    return;
  }

  // Determine whether the matrix U should be M x N or M x M

  if (LSVEC) {
    N1 = N;
    if (lsame(JOBU, 'F')) N1 = M;
  }

  // Set numerical parameters

  // !    NOTE: Make sure dlamch() does not fail on the target architecture.

  EPSLN = dlamch('Epsilon');
  SFMIN = dlamch('SafeMinimum');
  SMALL = SFMIN / EPSLN;
  BIG = dlamch('O');
  // BIG   = ONE / SFMIN

  // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N

  // (!)  If necessary, scale SVA() to protect the largest norm from
  // overflow. It is possible that this scaling pushes the smallest
  // column norm left from the underflow threshold (extreme case).

  SCALEM = ONE / sqrt(M.toDouble() * N.toDouble());
  NOSCAL = true;
  GOSCAL = true;
  for (p = 1; p <= N; p++) {
    AAPP.value = ZERO;
    AAQQ.value = ONE;
    zlassq(M, A(1, p).asArray(), 1, AAPP, AAQQ);
    if (AAPP.value > BIG) {
      INFO.value = -9;
      xerbla('ZGEJSV', -INFO.value);
      return;
    }
    AAQQ.value = sqrt(AAQQ.value);
    if ((AAPP.value < (BIG / AAQQ.value)) && NOSCAL) {
      SVA[p] = AAPP.value * AAQQ.value;
    } else {
      NOSCAL = false;
      SVA[p] = AAPP.value * (AAQQ.value * SCALEM);
      if (GOSCAL) {
        GOSCAL = false;
        dscal(p - 1, SCALEM, SVA, 1);
      }
    }
  }

  if (NOSCAL) SCALEM = ONE;

  AAPP.value = ZERO;
  AAQQ.value = BIG;
  for (p = 1; p <= N; p++) {
    AAPP.value = max(AAPP.value, SVA[p]);
    if (SVA[p] != ZERO) AAQQ.value = min(AAQQ.value, SVA[p]);
  }

  // Quick return for zero M x N matrix
// #:)
  if (AAPP.value == ZERO) {
    if (LSVEC) zlaset('G', M, N1, Complex.zero, Complex.one, U, LDU);
    if (RSVEC) zlaset('G', N, N, Complex.zero, Complex.one, V, LDV);
    RWORK[1] = ONE;
    RWORK[2] = ONE;
    if (ERREST) RWORK[3] = ONE;
    if (LSVEC && RSVEC) {
      RWORK[4] = ONE;
      RWORK[5] = ONE;
    }
    if (L2TRAN) {
      RWORK[6] = ZERO;
      RWORK[7] = ZERO;
    }
    IWORK[1] = 0;
    IWORK[2] = 0;
    IWORK[3] = 0;
    IWORK[4] = -1;
    return;
  }

  // Issue warning if denormalized column norms detected. Override the
  // high relative accuracy request. Issue licence to kill nonzero columns
  // (set them to zero) whose norm is less than sigma_max / BIG (roughly).
// #:(
  WARNING = 0;
  if (AAQQ.value <= SFMIN) {
    L2RANK = true;
    L2KILL = true;
    WARNING = 1;
  }

  // Quick return for one-column matrix
// #:)
  if (N == 1) {
    if (LSVEC) {
      zlascl('G', 0, 0, SVA[1], SCALEM, M, 1, A(1, 1), LDA, IERR);
      zlacpy('A', M, 1, A, LDA, U, LDU);
      // computing all M left singular vectors of the M x 1 matrix
      if (N1 != N) {
        zgeqrf(M, N, U, LDU, CWORK, CWORK(N + 1), LWORK - N, IERR);
        zungqr(M, N1, 1, U, LDU, CWORK, CWORK(N + 1), LWORK - N, IERR);
        zcopy(M, A(1, 1).asArray(), 1, U(1, 1).asArray(), 1);
      }
    }
    if (RSVEC) {
      V[1][1] = Complex.one;
    }
    if (SVA[1] < (BIG * SCALEM)) {
      SVA[1] /= SCALEM;
      SCALEM = ONE;
    }
    RWORK[1] = ONE / SCALEM;
    RWORK[2] = ONE;
    if (SVA[1] != ZERO) {
      IWORK[1] = 1;
      if ((SVA[1] / SCALEM) >= SFMIN) {
        IWORK[2] = 1;
      } else {
        IWORK[2] = 0;
      }
    } else {
      IWORK[1] = 0;
      IWORK[2] = 0;
    }
    IWORK[3] = 0;
    IWORK[4] = -1;
    if (ERREST) RWORK[3] = ONE;
    if (LSVEC && RSVEC) {
      RWORK[4] = ONE;
      RWORK[5] = ONE;
    }
    if (L2TRAN) {
      RWORK[6] = ZERO;
      RWORK[7] = ZERO;
    }
    return;
  }

  TRANSP = false;

  AATMAX = -ONE;
  AATMIN = BIG;
  if (ROWPIV || L2TRAN) {
    // Compute the row norms, needed to determine row pivoting sequence
    // (in the case of heavily row weighted A, row pivoting is strongly
    // advised) and to collect information needed to compare the
    // structures of A * A^* and A^* * A (in the case L2TRAN == true ).

    if (L2TRAN) {
      for (p = 1; p <= M; p++) {
        XSC.value = ZERO;
        TEMP1.value = ONE;
        zlassq(N, A(p, 1).asArray(), LDA, XSC, TEMP1);
        // ZLASSQ gets both the ell_2 and the ell_infinity norm
        // in one pass through the vector
        RWORK[M + p] = XSC.value * SCALEM;
        RWORK[p] = XSC.value * (SCALEM * sqrt(TEMP1.value));
        AATMAX = max(AATMAX, RWORK[p]);
        if (RWORK[p] != ZERO) AATMIN = min(AATMIN, RWORK[p]);
      }
    } else {
      for (p = 1; p <= M; p++) {
        RWORK[M + p] = SCALEM * A[p][izamax(N, A(p, 1).asArray(), LDA)].abs();
        AATMAX = max(AATMAX, RWORK[M + p]);
        AATMIN = min(AATMIN, RWORK[M + p]);
      }
    }
  }

  // For square matrix A try to determine whether A^*  would be better
  // input for the preconditioned Jacobi SVD, with faster convergence.
  // The decision is based on an O(N) function of the vector of column
  // and row norms of A, based on the Shannon entropy. This should give
  // the right choice in most cases when the difference actually matters.
  // It may fail and pick the slower converging side.

  ENTRA = ZERO;
  ENTRAT = ZERO;
  if (L2TRAN) {
    XSC.value = ZERO;
    TEMP1.value = ONE;
    dlassq(N, SVA, 1, XSC, TEMP1);
    TEMP1.value = ONE / TEMP1.value;

    ENTRA = ZERO;
    for (p = 1; p <= N; p++) {
      BIG1 = (pow(SVA[p] / XSC.value, 2)) * TEMP1.value;
      if (BIG1 != ZERO) ENTRA += BIG1 * log(BIG1);
    }
    ENTRA = -ENTRA / log(N);

    // Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex.
    // It is derived from the diagonal of  A^* * A.  Do the same with the
    // diagonal of A * A^*, compute the entropy of the corresponding
    // probability distribution. Note that A * A^* and A^* * A have the
    // same trace.

    ENTRAT = ZERO;
    for (p = 1; p <= M; p++) {
      BIG1 = (pow(RWORK[p] / XSC.value, 2)) * TEMP1.value;
      if (BIG1 != ZERO) ENTRAT += BIG1 * log(BIG1);
    }
    ENTRAT = -ENTRAT / log(M);

    // Analyze the entropies and decide A or A^*. Smaller entropy
    // usually means better input for the algorithm.

    TRANSP = (ENTRAT < ENTRA);

    // If A^* is better than A, take the adjoint of A. This is allowed
    // only for square matrices, M=N.
    if (TRANSP) {
      // In an optimal implementation, this trivial transpose
      // should be replaced with faster transpose.
      for (p = 1; p <= N - 1; p++) {
        A[p][p] = A[p][p].conjugate();
        for (q = p + 1; q <= N; q++) {
          CTEMP = A[q][p].conjugate();
          A[q][p] = A[p][q].conjugate();
          A[p][q] = CTEMP;
        }
      }
      A[N][N] = A[N][N].conjugate();
      for (p = 1; p <= N; p++) {
        RWORK[M + p] = SVA[p];
        SVA[p] = RWORK[p];
        // previously computed row 2-norms are now column 2-norms
        // of the transposed matrix
      }
      TEMP1.value = AAPP.value;
      AAPP.value = AATMAX;
      AATMAX = TEMP1.value;
      TEMP1.value = AAQQ.value;
      AAQQ.value = AATMIN;
      AATMIN = TEMP1.value;
      KILL = LSVEC;
      LSVEC = RSVEC;
      RSVEC = KILL;
      if (LSVEC) N1 = N;

      ROWPIV = true;
    }
  }
  // END IF L2TRAN

  // Scale the matrix so that its maximal singular value remains less
  // than sqrt(BIG) -- the matrix is scaled so that its maximal column
  // has Euclidean norm equal to sqrt(BIG/N). The only reason to keep
  // sqrt(BIG) instead of BIG is the fact that ZGEJSV uses LAPACK and
  // BLAS routines that, in some implementations, are not capable of
  // working in the full interval [SFMIN,BIG] and that they may provoke
  // overflows in the intermediate results. If the singular values spread
  // from SFMIN to BIG, then ZGESVJ will compute them. So, in that case,
  // one should use ZGESVJ instead of ZGEJSV.
  // >> change in the April 2016 update: allow bigger range, i.e. the
  // largest column is allowed up to BIG/N and ZGESVJ will do the rest.
  BIG1 = sqrt(BIG);
  TEMP1.value = sqrt(BIG / N.toDouble());
  // TEMP1.value  = BIG/N.toDouble()

  dlascl('G', 0, 0, AAPP.value, TEMP1.value, N, 1, SVA.asMatrix(N), N, IERR);
  if (AAQQ.value > (AAPP.value * SFMIN)) {
    AAQQ.value = (AAQQ.value / AAPP.value) * TEMP1.value;
  } else {
    AAQQ.value = (AAQQ.value * TEMP1.value) / AAPP.value;
  }
  TEMP1.value *= SCALEM;
  zlascl('G', 0, 0, AAPP.value, TEMP1.value, M, N, A, LDA, IERR);

  // To undo scaling at the end of this procedure, multiply the
  // computed singular values with USCAL2 / USCAL1.

  USCAL1 = TEMP1.value;
  USCAL2 = AAPP.value;

  if (L2KILL) {
    // L2KILL enforces computation of nonzero singular values in
    // the restricted range of condition number of the initial A,
    // sigma_max(A) / sigma_min(A) approx. sqrt(BIG)/sqrt(SFMIN).
    XSC.value = sqrt(SFMIN);
  } else {
    XSC.value = SMALL;

    // Now, if the condition number of A is too big,
    // sigma_max(A) / sigma_min(A) > sqrt(BIG/N) * EPSLN / SFMIN,
    // as a precaution measure, the full SVD is computed using ZGESVJ
    // with accumulated Jacobi rotations. This provides numerically
    // more robust computation, at the cost of slightly increased run
    // time. Depending on the concrete implementation of BLAS and LAPACK
    // (i.e. how they behave in presence of extreme ill-conditioning) the
    // implementor may decide to remove this switch.
    if ((AAQQ.value < sqrt(SFMIN)) && LSVEC && RSVEC) {
      JRACC = true;
    }
  }
  if (AAQQ.value < XSC.value) {
    for (p = 1; p <= N; p++) {
      if (SVA[p] < XSC.value) {
        zlaset('A', M, 1, Complex.zero, Complex.zero, A(1, p), LDA);
        SVA[p] = ZERO;
      }
    }
  }

  // Preconditioning using QR factorization with pivoting

  if (ROWPIV) {
    // Optional row permutation (Bjoerck row pivoting):
    // A result by Cox and Higham shows that the Bjoerck's
    // row pivoting combined with standard column pivoting
    // has similar effect as Powell-Reid complete pivoting.
    // The ell-infinity norms of A are made nonincreasing.
    if ((LSVEC && RSVEC) && !JRACC) {
      IWOFF = 2 * N;
    } else {
      IWOFF = N;
    }
    for (p = 1; p <= M - 1; p++) {
      q = idamax(M - p + 1, RWORK(M + p), 1) + p - 1;
      IWORK[IWOFF + p] = q;
      if (p != q) {
        TEMP1.value = RWORK[M + p];
        RWORK[M + p] = RWORK[M + q];
        RWORK[M + q] = TEMP1.value;
      }
    }
    zlaswp(N, A, LDA, 1, M - 1, IWORK(IWOFF + 1), 1);
  }

  // End of the preparation phase (scaling, optional sorting and
  // transposing, optional flushing of small columns).

  // Preconditioning

  // If the full SVD is needed, the right singular vectors are computed
  // from a matrix equation, and for that we need theoretical analysis
  // of the Businger-Golub pivoting. So we use ZGEQP3 as the first RR QRF.
  // In all other cases the first RR QRF can be chosen by other criteria
  // (eg speed by replacing global with restricted window pivoting, such
  // as in xGEQPX from TOMS # 782). Good results will be obtained using
  // xGEQPX with properly (!) chosen numerical parameters.
  // Any improvement of ZGEQP3 improves overall performance of ZGEJSV.

  // A * P1 = Q1 * [ R1^* 0]^*:
  for (p = 1; p <= N; p++) {
    // .. all columns are free columns
    IWORK[p] = 0;
  }
  zgeqp3(M, N, A, LDA, IWORK, CWORK, CWORK(N + 1), LWORK - N, RWORK, IERR);

  // The upper triangular matrix R1 from the first QRF is inspected for
  // rank deficiency and possibilities for deflation, or possible
  // ill-conditioning. Depending on the user specified flag L2RANK,
  // the procedure explores possibilities to reduce the numerical
  // rank by inspecting the computed upper triangular factor. If
  // L2RANK or L2ABER are up, then ZGEJSV will compute the SVD of
  // A + dA, where ||dA|| <= f(M,N)*EPSLN.

  NR = 1;
  if (L2ABER) {
    // Standard absolute error bound suffices. All sigma_i with
    // sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
    // aggressive enforcement of lower numerical rank by introducing a
    // backward error of the order of N*EPSLN*||A||.
    TEMP1.value = sqrt(N.toDouble()) * EPSLN;
    for (p = 2; p <= N; p++) {
      if (A[p][p].abs() >= (TEMP1.value * A[1][1].abs())) {
        NR++;
      } else {
        break;
      }
    }
  } else if (L2RANK) {
    // .. similarly as above, only slightly more gentle (less aggressive).
    // Sudden drop on the diagonal of R1 is used as the criterion for
    // close-to-rank-deficient.
    TEMP1.value = sqrt(SFMIN);
    for (p = 2; p <= N; p++) {
      if ((A[p][p].abs() < (EPSLN * A[p - 1][p - 1].abs())) ||
          (A[p][p].abs() < SMALL) ||
          (L2KILL && (A[p][p].abs() < TEMP1.value))) break;
      NR++;
    }
  } else {
    // The goal is high relative accuracy. However, if the matrix
    // has high scaled condition number the relative accuracy is in
    // general not feasible. Later on, a condition number estimator
    // will be deployed to estimate the scaled condition number.
    // Here we just remove the underflowed part of the triangular
    // factor. This prevents the situation in which the code is
    // working hard to get the accuracy not warranted by the data.
    TEMP1.value = sqrt(SFMIN);
    for (p = 2; p <= N; p++) {
      if ((A[p][p].abs() < SMALL) ||
          (L2KILL && (A[p][p].abs() < TEMP1.value))) {
        break;
      }
      NR++;
    }
  }

  ALMORT = false;
  if (NR == N) {
    MAXPRJ = ONE;
    for (p = 2; p <= N; p++) {
      TEMP1.value = A[p][p].abs() / SVA[IWORK[p]];
      MAXPRJ = min(MAXPRJ, TEMP1.value);
    }
    if (pow(MAXPRJ, 2) >= ONE - N.toDouble() * EPSLN) ALMORT = true;
  }

  SCONDA = -ONE;
  CONDR1 = -ONE;
  CONDR2 = -ONE;

  if (ERREST) {
    if (N == NR) {
      if (RSVEC) {
        // .. V is available as workspace
        zlacpy('U', N, N, A, LDA, V, LDV);
        for (p = 1; p <= N; p++) {
          TEMP1.value = SVA[IWORK[p]];
          zdscal(p, ONE / TEMP1.value, V(1, p).asArray(), 1);
        }
        if (LSVEC) {
          zpocon('U', N, V, LDV, ONE, TEMP1, CWORK(N + 1), RWORK, IERR);
        } else {
          zpocon('U', N, V, LDV, ONE, TEMP1, CWORK, RWORK, IERR);
        }
      } else if (LSVEC) {
        // .. U is available as workspace
        zlacpy('U', N, N, A, LDA, U, LDU);
        for (p = 1; p <= N; p++) {
          TEMP1.value = SVA[IWORK[p]];
          zdscal(p, ONE / TEMP1.value, U(1, p).asArray(), 1);
        }
        zpocon('U', N, U, LDU, ONE, TEMP1, CWORK(N + 1), RWORK, IERR);
      } else {
        zlacpy('U', N, N, A, LDA, CWORK.asMatrix(N), N);
        // []            CALL ZLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
        // Change: here index shifted by N to the left, CWORK(1:N)
        // not needed for SIGMA only computation
        for (p = 1; p <= N; p++) {
          TEMP1.value = SVA[IWORK[p]];
          // []               CALL ZDSCAL( p, ONE/TEMP1.value, CWORK(N+(p-1)*N+1), 1 )
          zdscal(p, ONE / TEMP1.value, CWORK((p - 1) * N + 1), 1);
        }
        // .. the columns of R are scaled to have unit Euclidean lengths.
        // []               CALL ZPOCON( 'U', N, CWORK(N+1), N, ONE, TEMP1.value,
        // []     $              CWORK(N+N*N+1), RWORK, IERR.value )
        zpocon('U', N, CWORK.asMatrix(N), N, ONE, TEMP1, CWORK(N * N + 1),
            RWORK, IERR);
      }
      if (TEMP1.value != ZERO) {
        SCONDA = ONE / sqrt(TEMP1.value);
      } else {
        SCONDA = -ONE;
      }
      // SCONDA is an estimate of sqrt(||(R^* * R)^(-1)||_1).
      // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
    } else {
      SCONDA = -ONE;
    }
  }

  L2PERT = L2PERT && ((A[1][1] / A[NR][NR]).abs() > sqrt(BIG1));
  // If there is no violent scaling, artificial perturbation is not needed.

  // Phase 3:

  if (!(RSVEC || LSVEC)) {
    // Singular Values only

    // .. transpose A(1:NR,1:N)
    for (p = 1; p <= min(N - 1, NR); p++) {
      zcopy(N - p, A(p, p + 1).asArray(), LDA, A(p + 1, p).asArray(), 1);
      zlacgv(N - p + 1, A(p, p).asArray(), 1);
    }
    if (NR == N) A[N][N] = A[N][N].conjugate();

    // The following two DO-loops introduce small relative perturbation
    // into the strict upper triangle of the lower triangular matrix.
    // Small entries below the main diagonal are also changed.
    // This modification is useful if the computing environment does not
    // provide/allow FLUSH TO ZERO underflow, for it prevents many
    // annoying denormalized numbers in case of strongly scaled matrices.
    // The perturbation is structured so that it does not introduce any
    // new perturbation of the singular values, and it does not destroy
    // the job done by the preconditioner.
    // The licence for this perturbation is in the variable L2PERT, which
    // should be false if FLUSH TO ZERO underflow is active.

    if (!ALMORT) {
      if (L2PERT) {
        // XSC.value = sqrt(SMALL)
        XSC.value = EPSLN / N.toDouble();
        for (q = 1; q <= NR; q++) {
          CTEMP = Complex(XSC.value * A[q][q].abs(), ZERO);
          for (p = 1; p <= N; p++) {
            if (((p > q) && (A[p][q].abs() <= TEMP1.value)) || (p < q)) {
              // $                     A[p][q] = TEMP1.value * ( A[p][q] / ABS(A[p][q]) )
              A[p][q] = CTEMP;
            }
          }
        }
      } else {
        zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, A(1, 2), LDA);
      }

      // .. second preconditioning using the QR factorization

      zgeqrf(N, NR, A, LDA, CWORK, CWORK(N + 1), LWORK - N, IERR);

      // .. and transpose upper to lower triangular
      for (p = 1; p <= NR - 1; p++) {
        zcopy(NR - p, A(p, p + 1).asArray(), LDA, A(p + 1, p).asArray(), 1);
        zlacgv(NR - p + 1, A(p, p).asArray(), 1);
      }
    }

    // Row-cyclic Jacobi SVD algorithm with column pivoting

    // .. again some perturbation (a "background noise") is added
    // to drown denormals
    if (L2PERT) {
      // XSC.value = sqrt(SMALL)
      XSC.value = EPSLN / N.toDouble();
      for (q = 1; q <= NR; q++) {
        CTEMP = Complex(XSC.value * A[q][q].abs(), ZERO);
        for (p = 1; p <= NR; p++) {
          if (((p > q) && (A[p][q].abs() <= TEMP1.value)) || (p < q)) {
            // $                   A[p][q] = TEMP1.value * ( A[p][q] / ABS(A[p][q]) )
            A[p][q] = CTEMP;
          }
        }
      }
    } else {
      zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, A(1, 2), LDA);
    }

    // .. and one-sided Jacobi rotations are started on a lower
    // triangular matrix (plus perturbation which is ignored in
    // the part which destroys triangular form (confusing?!))

    zgesvj('L', 'N', 'N', NR, NR, A, LDA, SVA, N, V, LDV, CWORK, LWORK, RWORK,
        LRWORK, INFO);

    SCALEM = RWORK[1];
    NUMRANK = nint(RWORK[2]);
  } else if ((RSVEC && !LSVEC && !JRACC) || (JRACC && !LSVEC && (NR != N))) {
    // -> Singular Values and Right Singular Vectors <-

    if (ALMORT) {
      // .. in this case NR equals N
      for (p = 1; p <= NR; p++) {
        zcopy(N - p + 1, A(p, p).asArray(), LDA, V(p, p).asArray(), 1);
        zlacgv(N - p + 1, V(p, p).asArray(), 1);
      }
      zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, V(1, 2), LDV);

      zgesvj('L', 'U', 'N', N, NR, V, LDV, SVA, NR, A, LDA, CWORK, LWORK, RWORK,
          LRWORK, INFO);
      SCALEM = RWORK[1];
      NUMRANK = nint(RWORK[2]);
    } else {
      // .. two more QR factorizations ( one QRF is not enough, two require
      // accumulated product of Jacobi rotations, three are perfect )

      zlaset('L', NR - 1, NR - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
      zgelqf(NR, N, A, LDA, CWORK, CWORK(N + 1), LWORK - N, IERR);
      zlacpy('L', NR, NR, A, LDA, V, LDV);
      zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, V(1, 2), LDV);
      zgeqrf(
          NR, NR, V, LDV, CWORK(N + 1), CWORK(2 * N + 1), LWORK - 2 * N, IERR);
      for (p = 1; p <= NR; p++) {
        zcopy(NR - p + 1, V(p, p).asArray(), LDV, V(p, p).asArray(), 1);
        zlacgv(NR - p + 1, V(p, p).asArray(), 1);
      }
      zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, V(1, 2), LDV);

      zgesvj('L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU, CWORK(N + 1),
          LWORK - N, RWORK, LRWORK, INFO);
      SCALEM = RWORK[1];
      NUMRANK = nint(RWORK[2]);
      if (NR < N) {
        zlaset('A', N - NR, NR, Complex.zero, Complex.zero, V(NR + 1, 1), LDV);
        zlaset('A', NR, N - NR, Complex.zero, Complex.zero, V(1, NR + 1), LDV);
        zlaset('A', N - NR, N - NR, Complex.zero, Complex.one,
            V(NR + 1, NR + 1), LDV);
      }

      zunmlq('L', 'C', N, N, NR, A, LDA, CWORK, V, LDV, CWORK(N + 1), LWORK - N,
          IERR);
    }
    // .. permute the rows of V
    // DO 8991 p = 1, N
    //    CALL ZCOPY( N, V(p,1), LDV, A(IWORK[p],1), LDA )
    // 8991    CONTINUE
    // CALL ZLACPY( 'All', N, N, A, LDA, V, LDV )
    zlapmr(false, N, N, V, LDV, IWORK);

    if (TRANSP) {
      zlacpy('A', N, N, V, LDV, U, LDU);
    }
  } else if (JRACC && !LSVEC && (NR == N)) {
    zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);

    zgesvj('U', 'N', 'V', N, N, A, LDA, SVA, N, V, LDV, CWORK, LWORK, RWORK,
        LRWORK, INFO);
    SCALEM = RWORK[1];
    NUMRANK = nint(RWORK[2]);
    zlapmr(false, N, N, V, LDV, IWORK);
  } else if (LSVEC && !RSVEC) {
    // .. Singular Values and Left Singular Vectors                 ..

    // .. second preconditioning step to avoid need to accumulate
    // Jacobi rotations in the Jacobi iterations.
    for (p = 1; p <= NR; p++) {
      zcopy(N - p + 1, A(p, p).asArray(), LDA, U(p, p).asArray(), 1);
      zlacgv(N - p + 1, U(p, p).asArray(), 1);
    }
    zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, U(1, 2), LDU);

    zgeqrf(N, NR, U, LDU, CWORK(N + 1), CWORK(2 * N + 1), LWORK - 2 * N, IERR);

    for (p = 1; p <= NR - 1; p++) {
      zcopy(NR - p, U(p, p + 1).asArray(), LDU, U(p + 1, p).asArray(), 1);
      zlacgv(N - p + 1, U(p, p).asArray(), 1);
    }
    zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, U(1, 2), LDU);

    zgesvj('L', 'U', 'N', NR, NR, U, LDU, SVA, NR, A, LDA, CWORK(N + 1),
        LWORK - N, RWORK, LRWORK, INFO);
    SCALEM = RWORK[1];
    NUMRANK = nint(RWORK[2]);

    if (NR < M) {
      zlaset('A', M - NR, NR, Complex.zero, Complex.zero, U(NR + 1, 1), LDU);
      if (NR < N1) {
        zlaset('A', NR, N1 - NR, Complex.zero, Complex.zero, U(1, NR + 1), LDU);
        zlaset('A', M - NR, N1 - NR, Complex.zero, Complex.one,
            U(NR + 1, NR + 1), LDU);
      }
    }

    zunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N + 1), LWORK - N,
        IERR);

    if (ROWPIV) zlaswp(N1, U, LDU, 1, M - 1, IWORK(IWOFF + 1), -1);

    for (p = 1; p <= N1; p++) {
      XSC.value = ONE / dznrm2(M, U(1, p).asArray(), 1);
      zdscal(M, XSC.value, U(1, p).asArray(), 1);
    }

    if (TRANSP) {
      zlacpy('A', N, N, U, LDU, V, LDV);
    }
  } else {
    // .. Full SVD ..

    if (!JRACC) {
      if (!ALMORT) {
        // Second Preconditioning Step (QRF [with pivoting])
        // Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
        // equivalent to an LQF CALL. Since in many libraries the QRF
        // seems to be better optimized than the LQF, we do explicit
        // transpose and use the QRF. This is subject to changes in an
        // optimized implementation of ZGEJSV.

        for (p = 1; p <= NR; p++) {
          zcopy(N - p + 1, A(p, p).asArray(), LDA, V(p, p).asArray(), 1);
          zlacgv(N - p + 1, V(p, p).asArray(), 1);
        }

        // .. the following two loops perturb small entries to avoid
        // denormals in the second QR factorization, where they are
        // as good as zeros. This is done to avoid painfully slow
        // computation with denormals. The relative size of the perturbation
        // is a parameter that can be changed by the implementer.
        // This perturbation device will be obsolete on machines with
        // properly implemented arithmetic.
        // To switch it off, set L2PERT= false To remove it from  the
        // code, remove the action under L2PERT= true , leave the ELSE part.
        // The following two loops should be blocked and fused with the
        // transposed copy above.

        if (L2PERT) {
          XSC.value = sqrt(SMALL);
          for (q = 1; q <= NR; q++) {
            CTEMP = Complex(XSC.value * V[q][q].abs(), ZERO);
            for (p = 1; p <= N; p++) {
              if ((p > q) && (V[p][q].abs() <= TEMP1.value) || (p < q)) {
                // $                   V[p][q] = TEMP1.value * ( V[p][q] / ABS(V[p][q]) )
                V[p][q] = CTEMP;
              }
              if (p < q) V[p][q] = -V[p][q];
            }
          }
        } else {
          zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, V(1, 2), LDV);
        }

        // Estimate the row scaled condition number of R1
        // (If R1 is rectangular, N > NR, then the condition number
        // of the leading NR x NR submatrix is estimated.)

        zlacpy('L', NR, NR, V, LDV, CWORK(2 * N + 1).asMatrix(NR), NR);
        for (p = 1; p <= NR; p++) {
          TEMP1.value = dznrm2(NR - p + 1, CWORK(2 * N + (p - 1) * NR + p), 1);
          zdscal(NR - p + 1, ONE / TEMP1.value, CWORK(2 * N + (p - 1) * NR + p),
              1);
        }
        zpocon('L', NR, CWORK(2 * N + 1).asMatrix(NR), NR, ONE, TEMP1,
            CWORK(2 * N + NR * NR + 1), RWORK, IERR);
        CONDR1 = ONE / sqrt(TEMP1.value);
        // .. here need a second opinion on the condition number
        // .. then assume worst case scenario
        // R1 is OK for inverse <=> CONDR1 < N.toDouble()
        // more conservative    <=> CONDR1 < sqrt(N.toDouble())

        COND_OK = sqrt(sqrt(NR.toDouble()));
        // [TP]       COND_OK is a tuning parameter.

        if (CONDR1 < COND_OK) {
          // .. the second QRF without pivoting. Note: in an optimized
          // implementation, this QRF should be implemented as the QRF
          // of a lower triangular matrix.
          // R1^* = Q2 * R2
          zgeqrf(N, NR, V, LDV, CWORK(N + 1), CWORK(2 * N + 1), LWORK - 2 * N,
              IERR);

          if (L2PERT) {
            XSC.value = sqrt(SMALL) / EPSLN;
            for (p = 2; p <= NR; p++) {
              for (q = 1; q <= p - 1; q++) {
                CTEMP = Complex(
                    XSC.value * min(V[p][p].abs(), V[q][q].abs()), ZERO);
                if (V[q][p].abs() <= TEMP1.value) {
                  // $                     V[q][p] = TEMP1.value * ( V[q][p] / ABS(V[q][p]) )
                  V[q][p] = CTEMP;
                }
              }
            }
          }

          if (NR != N) {
            zlacpy('A', N, NR, V, LDV, CWORK(2 * N + 1).asMatrix(N), N);
          }
          // .. save ...

          // .. this transposed copy should be better than naive
          for (p = 1; p <= NR - 1; p++) {
            zcopy(NR - p, V(p, p + 1).asArray(), LDV, V(p + 1, p).asArray(), 1);
            zlacgv(NR - p + 1, V(p, p).asArray(), 1);
          }
          V[NR][NR] = V[NR][NR].conjugate();

          CONDR2 = CONDR1;
        } else {
          // .. ill-conditioned case: second QRF with pivoting
          // Note that windowed pivoting would be equally good
          // numerically, and more run-time efficient. So, in
          // an optimal implementation, the next call to ZGEQP3
          // should be replaced with eg. CALL ZGEQPX (ACM TOMS #782)
          // with properly (carefully) chosen parameters.

          // R1^* * P2 = Q2 * R2
          for (p = 1; p <= NR; p++) {
            IWORK[N + p] = 0;
          }
          zgeqp3(N, NR, V, LDV, IWORK(N + 1), CWORK(N + 1), CWORK(2 * N + 1),
              LWORK - 2 * N, RWORK, IERR);
// *               CALL ZGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1),
// *     $              LWORK-2*N, IERR.value )
          if (L2PERT) {
            XSC.value = sqrt(SMALL);
            for (p = 2; p <= NR; p++) {
              for (q = 1; q <= p - 1; q++) {
                CTEMP = Complex(
                    XSC.value * min(V[p][p].abs(), V[q][q].abs()), ZERO);
                if (V[q][p].abs() <= TEMP1.value) {
                  // V[q][p] = TEMP1.value * ( V[q][p] / ABS(V[q][p]) )
                  V[q][p] = CTEMP;
                }
              }
            }
          }

          zlacpy('A', N, NR, V, LDV, CWORK(2 * N + 1).asMatrix(N), N);

          if (L2PERT) {
            XSC.value = sqrt(SMALL);
            for (p = 2; p <= NR; p++) {
              for (q = 1; q <= p - 1; q++) {
                CTEMP = Complex(
                    XSC.value * min(V[p][p].abs(), V[q][q].abs()), ZERO);
                // V[p][q] = - TEMP1.value*( V[q][p] / ABS(V[q][p]) )
                V[p][q] = -CTEMP;
              }
            }
          } else {
            zlaset(
                'L', NR - 1, NR - 1, Complex.zero, Complex.zero, V(2, 1), LDV);
          }
          // Now, compute R2 = L3 * Q3, the LQ factorization.
          zgelqf(
              NR,
              NR,
              V,
              LDV,
              CWORK(2 * N + N * NR + 1),
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);
          // .. and estimate the condition number
          zlacpy('L', NR, NR, V, LDV,
              CWORK(2 * N + N * NR + NR + 1).asMatrix(NR), NR);
          for (p = 1; p <= NR; p++) {
            TEMP1.value = dznrm2(p, CWORK(2 * N + N * NR + NR + p), NR);
            zdscal(p, ONE / TEMP1.value, CWORK(2 * N + N * NR + NR + p), NR);
          }
          zpocon('L', NR, CWORK(2 * N + N * NR + NR + 1).asMatrix(NR), NR, ONE,
              TEMP1, CWORK(2 * N + N * NR + NR + NR * NR + 1), RWORK, IERR);
          CONDR2 = ONE / sqrt(TEMP1.value);

          if (CONDR2 >= COND_OK) {
            // .. save the Householder vectors used for Q3
            // (this overwrites the copy of R2, as it will not be
            // needed in this branch, but it does not overwrite the
            // Huseholder vectors of Q2.).
            zlacpy('U', NR, NR, V, LDV, CWORK(2 * N + 1).asMatrix(N), N);
            // .. and the rest of the information on Q3 is in
            // WORK(2*N+N*NR+1:2*N+N*NR+N)
          }
        }

        if (L2PERT) {
          XSC.value = sqrt(SMALL);
          for (q = 2; q <= NR; q++) {
            CTEMP = XSC.value.toComplex() * V[q][q];
            for (p = 1; p <= q - 1; p++) {
              // V[p][q] = - TEMP1.value*( V[p][q] / ABS(V[p][q]) )
              V[p][q] = -CTEMP;
            }
          }
        } else {
          zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, V(1, 2), LDV);
        }

        // Second preconditioning finished; continue with Jacobi SVD
        // The input matrix is lower triangular.

        // Recover the right singular vectors as solution of a well
        // conditioned triangular matrix equation.

        if (CONDR1 < COND_OK) {
          zgesvj(
              'L',
              'U',
              'N',
              NR,
              NR,
              V,
              LDV,
              SVA,
              NR,
              U,
              LDU,
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              RWORK,
              LRWORK,
              INFO);
          SCALEM = RWORK[1];
          NUMRANK = nint(RWORK[2]);
          for (p = 1; p <= NR; p++) {
            zcopy(NR, V(1, p).asArray(), 1, U(1, p).asArray(), 1);
            zdscal(NR, SVA[p], V(1, p).asArray(), 1);
          }

          // .. pick the right matrix equation and solve it

          if (NR == N) {
            // .. best case, R1 is inverted. The solution of this matrix
            // equation is Q2*V2 = the product of the Jacobi rotations
            // used in ZGESVJ, premultiplied with the orthogonal matrix
            // from the second QR factorization.
            ztrsm('L', 'U', 'N', 'N', NR, NR, Complex.one, A, LDA, V, LDV);
          } else {
            // .. R1 is well conditioned, but non-square. Adjoint of R2
            // is inverted to get the product of the Jacobi rotations
            // used in ZGESVJ. The Q-factor from the second QR
            // factorization is then built in explicitly.
            ztrsm('L', 'U', 'C', 'N', NR, NR, Complex.one,
                CWORK(2 * N + 1).asMatrix(N), N, V, LDV);
            if (NR < N) {
              zlaset('A', N - NR, NR, Complex.zero, Complex.zero, V(NR + 1, 1),
                  LDV);
              zlaset('A', NR, N - NR, Complex.zero, Complex.zero, V(1, NR + 1),
                  LDV);
              zlaset('A', N - NR, N - NR, Complex.zero, Complex.one,
                  V(NR + 1, NR + 1), LDV);
            }
            zunmqr(
                'L',
                'N',
                N,
                N,
                NR,
                CWORK(2 * N + 1).asMatrix(N),
                N,
                CWORK(N + 1),
                V,
                LDV,
                CWORK(2 * N + N * NR + NR + 1),
                LWORK - 2 * N - N * NR - NR,
                IERR);
          }
        } else if (CONDR2 < COND_OK) {
          // The matrix R2 is inverted. The solution of the matrix equation
          // is Q3^* * V3 = the product of the Jacobi rotations (applied to
          // the lower triangular L3 from the LQ factorization of
          // R2=L3*Q3), pre-multiplied with the transposed Q3.
          zgesvj(
              'L',
              'U',
              'N',
              NR,
              NR,
              V,
              LDV,
              SVA,
              NR,
              U,
              LDU,
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              RWORK,
              LRWORK,
              INFO);
          SCALEM = RWORK[1];
          NUMRANK = nint(RWORK[2]);
          for (p = 1; p <= NR; p++) {
            zcopy(NR, V(1, p).asArray(), 1, U(1, p).asArray(), 1);
            zdscal(NR, SVA[p], U(1, p).asArray(), 1);
          }
          ztrsm('L', 'U', 'N', 'N', NR, NR, Complex.one,
              CWORK(2 * N + 1).asMatrix(N), N, U, LDU);
          // .. apply the permutation from the second QR factorization
          for (q = 1; q <= NR; q++) {
            for (p = 1; p <= NR; p++) {
              CWORK[2 * N + N * NR + NR + IWORK[N + p]] = U[p][q];
            }
            for (p = 1; p <= NR; p++) {
              U[p][q] = CWORK[2 * N + N * NR + NR + p];
            }
          }
          if (NR < N) {
            zlaset(
                'A', N - NR, NR, Complex.zero, Complex.zero, V(NR + 1, 1), LDV);
            zlaset(
                'A', NR, N - NR, Complex.zero, Complex.zero, V(1, NR + 1), LDV);
            zlaset('A', N - NR, N - NR, Complex.zero, Complex.one,
                V(NR + 1, NR + 1), LDV);
          }
          zunmqr(
              'L',
              'N',
              N,
              N,
              NR,
              CWORK(2 * N + 1).asMatrix(N),
              N,
              CWORK(N + 1),
              V,
              LDV,
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);
        } else {
          // Last line of defense.
// #:(          This is a rather pathological case: no scaled condition
          // improvement after two pivoted QR factorizations. Other
          // possibility is that the rank revealing QR factorization
          // or the condition estimator has failed, or the COND_OK
          // is set very close to ONE (which is unnecessary). Normally,
          // this branch should never be executed, but in rare cases of
          // failure of the RRQR or condition estimator, the last line of
          // defense ensures that ZGEJSV completes the task.
          // Compute the full SVD of L3 using ZGESVJ with explicit
          // accumulation of Jacobi rotations.
          zgesvj(
              'L',
              'U',
              'V',
              NR,
              NR,
              V,
              LDV,
              SVA,
              NR,
              U,
              LDU,
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              RWORK,
              LRWORK,
              INFO);
          SCALEM = RWORK[1];
          NUMRANK = nint(RWORK[2]);
          if (NR < N) {
            zlaset(
                'A', N - NR, NR, Complex.zero, Complex.zero, V(NR + 1, 1), LDV);
            zlaset(
                'A', NR, N - NR, Complex.zero, Complex.zero, V(1, NR + 1), LDV);
            zlaset('A', N - NR, N - NR, Complex.zero, Complex.one,
                V(NR + 1, NR + 1), LDV);
          }
          zunmqr(
              'L',
              'N',
              N,
              N,
              NR,
              CWORK(2 * N + 1).asMatrix(N),
              N,
              CWORK(N + 1),
              V,
              LDV,
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);

          zunmlq(
              'L',
              'C',
              NR,
              NR,
              NR,
              CWORK(2 * N + 1).asMatrix(N),
              N,
              CWORK(2 * N + N * NR + 1),
              U,
              LDU,
              CWORK(2 * N + N * NR + NR + 1),
              LWORK - 2 * N - N * NR - NR,
              IERR);
          for (q = 1; q <= NR; q++) {
            for (p = 1; p <= NR; p++) {
              CWORK[2 * N + N * NR + NR + IWORK[N + p]] = U[p][q];
            }
            for (p = 1; p <= NR; p++) {
              U[p][q] = CWORK[2 * N + N * NR + NR + p];
            }
          }
        }

        // Permute the rows of V using the (column) permutation from the
        // first QRF. Also, scale the columns to make them unit in
        // Euclidean norm. This applies to all cases.

        TEMP1.value = sqrt(N.toDouble()) * EPSLN;
        for (q = 1; q <= N; q++) {
          for (p = 1; p <= N; p++) {
            CWORK[2 * N + N * NR + NR + IWORK[p]] = V[p][q];
          }
          for (p = 1; p <= N; p++) {
            V[p][q] = CWORK[2 * N + N * NR + NR + p];
          }
          XSC.value = ONE / dznrm2(N, V(1, q).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            zdscal(N, XSC.value, V(1, q).asArray(), 1);
          }
        }
        // At this moment, V contains the right singular vectors of A.
        // Next, assemble the left singular vector matrix U (M x N).
        if (NR < M) {
          zlaset(
              'A', M - NR, NR, Complex.zero, Complex.zero, U(NR + 1, 1), LDU);
          if (NR < N1) {
            zlaset('A', NR, N1 - NR, Complex.zero, Complex.zero, U(1, NR + 1),
                LDU);
            zlaset('A', M - NR, N1 - NR, Complex.zero, Complex.one,
                U(NR + 1, NR + 1), LDU);
          }
        }

        // The Q matrix from the first QRF is built into the left singular
        // matrix U. This applies to all cases.

        zunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N + 1),
            LWORK - N, IERR);

        // The columns of U are normalized. The cost is O(M*N) flops.
        TEMP1.value = sqrt(M.toDouble()) * EPSLN;
        for (p = 1; p <= NR; p++) {
          XSC.value = ONE / dznrm2(M, U(1, p).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            zdscal(M, XSC.value, U(1, p).asArray(), 1);
          }
        }

        // If the initial QRF is computed with row pivoting, the left
        // singular vectors must be adjusted.

        if (ROWPIV) zlaswp(N1, U, LDU, 1, M - 1, IWORK(IWOFF + 1), -1);
      } else {
        // .. the initial matrix A has almost orthogonal columns and
        // the second QRF is not needed

        zlacpy('U', N, N, A, LDA, CWORK(N + 1).asMatrix(N), N);
        if (L2PERT) {
          XSC.value = sqrt(SMALL);
          for (p = 2; p <= N; p++) {
            CTEMP = XSC.value.toComplex() * CWORK[N + (p - 1) * N + p];
            for (q = 1; q <= p - 1; q++) {
              // CWORK(N+(q-1)*N+p)=-TEMP1.value * ( CWORK(N+(p-1)*N+q) /
              // $                                        ABS(CWORK(N+(p-1)*N+q)) )
              CWORK[N + (q - 1) * N + p] = -CTEMP;
            }
          }
        } else {
          zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
              CWORK(N + 2).asMatrix(N), N);
        }

        zgesvj('U', 'U', 'N', N, N, CWORK(N + 1).asMatrix(N), N, SVA, N, U, LDU,
            CWORK(N + N * N + 1), LWORK - N - N * N, RWORK, LRWORK, INFO);

        SCALEM = RWORK[1];
        NUMRANK = nint(RWORK[2]);
        for (p = 1; p <= N; p++) {
          zcopy(N, CWORK(N + (p - 1) * N + 1), 1, U(1, p).asArray(), 1);
          zdscal(N, SVA[p], CWORK(N + (p - 1) * N + 1), 1);
        }

        ztrsm('L', 'U', 'N', 'N', N, N, Complex.one, A, LDA,
            CWORK(N + 1).asMatrix(N), N);
        for (p = 1; p <= N; p++) {
          zcopy(N, CWORK(N + p), N, V(IWORK[p], 1).asArray(), LDV);
        }
        TEMP1.value = sqrt(N.toDouble()) * EPSLN;
        for (p = 1; p <= N; p++) {
          XSC.value = ONE / dznrm2(N, V(1, p).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            zdscal(N, XSC.value, V(1, p).asArray(), 1);
          }
        }

        // Assemble the left singular vector matrix U (M x N).

        if (N < M) {
          zlaset('A', M - N, N, Complex.zero, Complex.zero, U(N + 1, 1), LDU);
          if (N < N1) {
            zlaset(
                'A', N, N1 - N, Complex.zero, Complex.zero, U(1, N + 1), LDU);
            zlaset('A', M - N, N1 - N, Complex.zero, Complex.one,
                U(N + 1, N + 1), LDU);
          }
        }
        zunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N + 1),
            LWORK - N, IERR);
        TEMP1.value = sqrt(M.toDouble()) * EPSLN;
        for (p = 1; p <= N1; p++) {
          XSC.value = ONE / dznrm2(M, U(1, p).asArray(), 1);
          if ((XSC.value < (ONE - TEMP1.value)) ||
              (XSC.value > (ONE + TEMP1.value))) {
            zdscal(M, XSC.value, U(1, p).asArray(), 1);
          }
        }

        if (ROWPIV) zlaswp(N1, U, LDU, 1, M - 1, IWORK(IWOFF + 1), -1);
      }

      // end of the  >> almost orthogonal case <<  in the full SVD
    } else {
      // This branch deploys a preconditioned Jacobi SVD with explicitly
      // accumulated rotations. It is included as optional, mainly for
      // experimental purposes. It does perform well, and can also be used.
      // In this implementation, this branch will be automatically activated
      // if the  condition number sigma_max(A) / sigma_min(A) is predicted
      // to be greater than the overflow threshold. This is because the
      // a posteriori computation of the singular vectors assumes robust
      // implementation of BLAS and some LAPACK procedures, capable of working
      // in presence of extreme values, e.g. when the singular values spread from
      // the underflow to the overflow threshold.

      for (p = 1; p <= NR; p++) {
        zcopy(N - p + 1, A(p, p).asArray(), LDA, V(p, p).asArray(), 1);
        zlacgv(N - p + 1, V(p, p).asArray(), 1);
      }

      if (L2PERT) {
        XSC.value = sqrt(SMALL / EPSLN);
        for (q = 1; q <= NR; q++) {
          CTEMP = Complex(XSC.value * V[q][q].abs(), ZERO);
          for (p = 1; p <= N; p++) {
            if ((p > q) && (V[p][q].abs() <= TEMP1.value) || (p < q)) {
              // $                V[p][q] = TEMP1.value * ( V[p][q] / ABS(V[p][q]) )
              V[p][q] = CTEMP;
            }
            if (p < q) V[p][q] = -V[p][q];
          }
        }
      } else {
        zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, V(1, 2), LDV);
      }
      zgeqrf(
          N, NR, V, LDV, CWORK(N + 1), CWORK(2 * N + 1), LWORK - 2 * N, IERR);
      zlacpy('L', N, NR, V, LDV, CWORK(2 * N + 1).asMatrix(N), N);

      for (p = 1; p <= NR; p++) {
        zcopy(NR - p + 1, V(p, p).asArray(), LDV, U(p, p).asArray(), 1);
        zlacgv(NR - p + 1, U(p, p).asArray(), 1);
      }

      if (L2PERT) {
        XSC.value = sqrt(SMALL / EPSLN);
        for (q = 2; q <= NR; q++) {
          for (p = 1; p <= q - 1; p++) {
            CTEMP =
                Complex(XSC.value * min(U[p][p].abs(), U[q][q].abs()), ZERO);
            // U[p][q] = - TEMP1.value * ( U[q][p] / ABS(U[q][p]) )
            U[p][q] = -CTEMP;
          }
        }
      } else {
        zlaset('U', NR - 1, NR - 1, Complex.zero, Complex.zero, U(1, 2), LDU);
      }
      zgesvj(
          'L',
          'U',
          'V',
          NR,
          NR,
          U,
          LDU,
          SVA,
          N,
          V,
          LDV,
          CWORK(2 * N + N * NR + 1),
          LWORK - 2 * N - N * NR,
          RWORK,
          LRWORK,
          INFO);
      SCALEM = RWORK[1];
      NUMRANK = nint(RWORK[2]);

      if (NR < N) {
        zlaset('A', N - NR, NR, Complex.zero, Complex.zero, V(NR + 1, 1), LDV);
        zlaset('A', NR, N - NR, Complex.zero, Complex.zero, V(1, NR + 1), LDV);
        zlaset('A', N - NR, N - NR, Complex.zero, Complex.one,
            V(NR + 1, NR + 1), LDV);
      }
      zunmqr(
          'L',
          'N',
          N,
          N,
          NR,
          CWORK(2 * N + 1).asMatrix(N),
          N,
          CWORK(N + 1),
          V,
          LDV,
          CWORK(2 * N + N * NR + NR + 1),
          LWORK - 2 * N - N * NR - NR,
          IERR);

      // Permute the rows of V using the (column) permutation from the
      // first QRF. Also, scale the columns to make them unit in
      // Euclidean norm. This applies to all cases.

      TEMP1.value = sqrt(N.toDouble()) * EPSLN;
      for (q = 1; q <= N; q++) {
        for (p = 1; p <= N; p++) {
          CWORK[2 * N + N * NR + NR + IWORK[p]] = V[p][q];
        }
        for (p = 1; p <= N; p++) {
          V[p][q] = CWORK[2 * N + N * NR + NR + p];
        }
        XSC.value = ONE / dznrm2(N, V(1, q).asArray(), 1);
        if ((XSC.value < (ONE - TEMP1.value)) ||
            (XSC.value > (ONE + TEMP1.value))) {
          zdscal(N, XSC.value, V(1, q).asArray(), 1);
        }
      }

      // At this moment, V contains the right singular vectors of A.
      // Next, assemble the left singular vector matrix U (M x N).

      if (NR < M) {
        zlaset('A', M - NR, NR, Complex.zero, Complex.zero, U(NR + 1, 1), LDU);
        if (NR < N1) {
          zlaset(
              'A', NR, N1 - NR, Complex.zero, Complex.zero, U(1, NR + 1), LDU);
          zlaset('A', M - NR, N1 - NR, Complex.zero, Complex.one,
              U(NR + 1, NR + 1), LDU);
        }
      }

      zunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N + 1), LWORK - N,
          IERR);

      if (ROWPIV) zlaswp(N1, U, LDU, 1, M - 1, IWORK(IWOFF + 1), -1);
    }
    if (TRANSP) {
      // .. swap U and V because the procedure worked on A^*
      for (p = 1; p <= N; p++) {
        zswap(N, U(1, p).asArray(), 1, V(1, p).asArray(), 1);
      }
    }
  }
  // end of the full SVD

  // Undo scaling, if necessary (and possible)

  if (USCAL2 <= (BIG / SVA[1]) * USCAL1) {
    dlascl('G', 0, 0, USCAL1, USCAL2, NR, 1, SVA.asMatrix(N), N, IERR);
    USCAL1 = ONE;
    USCAL2 = ONE;
  }

  if (NR < N) {
    for (p = NR + 1; p <= N; p++) {
      SVA[p] = ZERO;
    }
  }

  RWORK[1] = USCAL2 * SCALEM;
  RWORK[2] = USCAL1;
  if (ERREST) RWORK[3] = SCONDA;
  if (LSVEC && RSVEC) {
    RWORK[4] = CONDR1;
    RWORK[5] = CONDR2;
  }
  if (L2TRAN) {
    RWORK[6] = ENTRA;
    RWORK[7] = ENTRAT;
  }

  IWORK[1] = NR;
  IWORK[2] = NUMRANK;
  IWORK[3] = WARNING;
  if (TRANSP) {
    IWORK[4] = 1;
  } else {
    IWORK[4] = -1;
  }
}
