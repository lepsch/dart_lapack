import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhb2st_kernels.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

void zhetrd_hb2st(
  final String STAGE1,
  final String VECT,
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> D_,
  final Array<double> E_,
  final Array<Complex> HOUS_,
  final int LHOUS,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// #if defined(_OPENMP)
  // use omp_lib;
// #endif

  final AB = AB_.having(ld: LDAB);
  final WORK = WORK_.having();
  final HOUS = HOUS_.having();
  final D = D_.having();
  final E = E_.having();
  const RZERO = 0.0;
  bool LQUERY, WANTQ, UPPER, AFTERS1;
  int I,
      M,
      K,
      IB,
      SWEEPID,
      MYID,
      SHIFT,
      STT,
      ST,
      ED,
      STIND,
      EDIND,
      BLKLASTIND,
      COLPT,
      THED,
      STEPERCOL,
      GRSIZ,
      THGRSIZ,
      THGRNB,
      THGRID,
      // NBTILES,
      TTYPE,
      // TID,
      // NTHREADS,
      ABDPOS,
      ABOFDPOS,
      DPOS,
      OFDPOS,
      AWPOS,
      INDA,
      INDW,
      APOS,
      SIZEA,
      LDA,
      INDV,
      INDTAU,
      // SIZEV,
      SIZETAU,
      LDV,
      LHMIN,
      LWMIN;
  double ABSTMP;
  Complex TMP;

  // Determine the minimal workspace size required.
  // Test the input parameters

  INFO.value = 0;
  AFTERS1 = lsame(STAGE1, 'Y');
  WANTQ = lsame(VECT, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1) || (LHOUS == -1);

  // Determine the block size, the workspace size and the hous size.

  IB = ilaenv2stage(2, 'ZHETRD_HB2ST', VECT, N, KD, -1, -1);
  if (N == 0 || KD <= 1) {
    LHMIN = 1;
    LWMIN = 1;
  } else {
    LHMIN = ilaenv2stage(3, 'ZHETRD_HB2ST', VECT, N, KD, IB, -1);
    LWMIN = ilaenv2stage(4, 'ZHETRD_HB2ST', VECT, N, KD, IB, -1);
  }

  if (!AFTERS1 && !lsame(STAGE1, 'N')) {
    INFO.value = -1;
  } else if (!lsame(VECT, 'N')) {
    INFO.value = -2;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (KD < 0) {
    INFO.value = -5;
  } else if (LDAB < (KD + 1)) {
    INFO.value = -7;
  } else if (LHOUS < LHMIN && !LQUERY) {
    INFO.value = -11;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    HOUS[1] = LHMIN.toComplex();
    WORK[1] = LWMIN.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZHETRD_HB2ST', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    HOUS[1] = Complex.one;
    WORK[1] = Complex.one;
    return;
  }

  // Determine pointer position

  LDV = KD + IB;
  SIZETAU = 2 * N;
  // SIZEV = 2 * N;
  INDTAU = 1;
  INDV = INDTAU + SIZETAU;
  LDA = 2 * KD + 1;
  SIZEA = LDA * N;
  INDA = 1;
  INDW = INDA + SIZEA;
  // NTHREADS = 1;
  // TID = 0;

  if (UPPER) {
    APOS = INDA + KD;
    AWPOS = INDA;
    DPOS = APOS + KD;
    OFDPOS = DPOS - 1;
    ABDPOS = KD + 1;
    ABOFDPOS = KD;
  } else {
    APOS = INDA;
    AWPOS = INDA + KD + 1;
    DPOS = APOS;
    OFDPOS = DPOS + 1;
    ABDPOS = 1;
    ABOFDPOS = 2;
  }

  // Case KD=0:
  // The matrix is diagonal. We just copy it (convert to "real" for
  // complex because D is double and the imaginary part should be 0)
  // and store it in D. A sequential code here is better or
  // in a parallel environment it might need two cores for D and E

  if (KD == 0) {
    for (I = 1; I <= N; I++) {
      D[I] = (AB[ABDPOS][I]).toDouble();
    }
    for (I = 1; I <= N - 1; I++) {
      E[I] = RZERO;
    }

    HOUS[1] = Complex.one;
    WORK[1] = Complex.one;
    return;
  }

  // Case KD=1:
  // The matrix is already Tridiagonal. We have to make diagonal
  // and offdiagonal elements real, and store them in D and E.
  // For that, for real precision just copy the diag and offdiag
  // to D and E while for the COMPLEX case the bulge chasing is
  // performed to convert the hermetian tridiagonal to symmetric
  // tridiagonal. A simpler conversion formula might be used, but then
  // updating the Q matrix will be required and based if Q is generated
  // or not this might complicate the story.

  if (KD == 1) {
    for (I = 1; I <= N; I++) {
      D[I] = (AB[ABDPOS][I]).toDouble();
    }

    // make off-diagonal elements real and copy them to E

    if (UPPER) {
      for (I = 1; I <= N - 1; I++) {
        TMP = AB[ABOFDPOS][I + 1];
        ABSTMP = TMP.abs();
        AB[ABOFDPOS][I + 1] = ABSTMP.toComplex();
        E[I] = ABSTMP;
        if (ABSTMP != RZERO) {
          TMP /= ABSTMP.toComplex();
        } else {
          TMP = Complex.one;
        }
        if (I < N - 1) AB[ABOFDPOS][I + 2] *= TMP;
        // IF( WANTZ ) THEN
        //    CALL ZSCAL( N, DCONJG( TMP ), Q( 1, I+1 ), 1 )
        // END IF
      }
    } else {
      for (I = 1; I <= N - 1; I++) {
        TMP = AB[ABOFDPOS][I];
        ABSTMP = TMP.abs();
        AB[ABOFDPOS][I] = ABSTMP.toComplex();
        E[I] = ABSTMP;
        if (ABSTMP != RZERO) {
          TMP /= ABSTMP.toComplex();
        } else {
          TMP = Complex.one;
        }
        if (I < N - 1) AB[ABOFDPOS][I + 1] *= TMP;
        // IF( WANTQ ) THEN
        //    CALL ZSCAL( N, TMP, Q( 1, I+1 ), 1 )
        // END IF
      }
    }

    HOUS[1] = Complex.one;
    WORK[1] = Complex.one;
    return;
  }

  // Main code start here.
  // Reduce the hermitian band of A to a tridiagonal matrix.

  THGRSIZ = N;
  GRSIZ = 1;
  SHIFT = 3;
  // NBTILES = (N / KD).ceil();
  STEPERCOL = (SHIFT / GRSIZ).ceil();
  THGRNB = (N - 1 / THGRSIZ).ceil();

  zlacpy('A', KD + 1, N, AB, LDAB, WORK(APOS).asMatrix(LDA), LDA);
  zlaset(
      'A', KD, N, Complex.zero, Complex.zero, WORK(AWPOS).asMatrix(LDA), LDA);

  // openMP parallelisation start here

  // #if defined(_OPENMP)
  // $OMP PARALLEL PRIVATE( TID, THGRID, BLKLASTIND )
  // $OMP$         PRIVATE( THED, I, M, K, ST, ED, STT, SWEEPID )
  // $OMP$         PRIVATE( MYID, TTYPE, COLPT, STIND, EDIND )
  // $OMP$         SHARED ( UPLO, WANTQ, INDV, INDTAU, HOUS, WORK)
  // $OMP$         SHARED ( N, KD, IB, NBTILES, LDA, LDV, INDA )
  // $OMP$         SHARED ( STEPERCOL, THGRNB, THGRSIZ, GRSIZ, SHIFT )
  // $OMP MASTER
  // #endif

  // main bulge chasing loop

  for (THGRID = 1; THGRID <= THGRNB; THGRID++) {
    STT = (THGRID - 1) * THGRSIZ + 1;
    THED = min((STT + THGRSIZ - 1), (N - 1));
    for (I = STT; I <= N - 1; I++) {
      ED = min(I, THED);
      if (STT > ED) break;
      for (M = 1; M <= STEPERCOL; M++) {
        ST = STT;
        for (SWEEPID = ST; SWEEPID <= ED; SWEEPID++) {
          for (K = 1; K <= GRSIZ; K++) {
            MYID = (I - SWEEPID) * (STEPERCOL * GRSIZ) + (M - 1) * GRSIZ + K;
            if (MYID == 1) {
              TTYPE = 1;
            } else {
              TTYPE = (MYID % 2) + 2;
            }

            if (TTYPE == 2) {
              COLPT = (MYID ~/ 2) * KD + SWEEPID;
              STIND = COLPT - KD + 1;
              EDIND = min(COLPT, N);
              BLKLASTIND = COLPT;
            } else {
              COLPT = ((MYID + 1) ~/ 2) * KD + SWEEPID;
              STIND = COLPT - KD + 1;
              EDIND = min(COLPT, N);
              if ((STIND >= EDIND - 1) && (EDIND == N)) {
                BLKLASTIND = N;
              } else {
                BLKLASTIND = 0;
              }
            }

            // Call the kernel

// // #if defined(_OPENMP) &&  _OPENMP >= 201307

//             if (TTYPE != 1) {
// // $OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
// // $OMP$     DEPEND(in:WORK(MYID-1))
// // $OMP$     DEPEND(out:WORK(MYID))
//               TID = OMP_GET_THREAD_NUM();
//               zhb2st_kernels(
//                   UPLO,
//                   WANTQ,
//                   TTYPE,
//                   STIND,
//                   EDIND,
//                   SWEEPID,
//                   N,
//                   KD,
//                   IB,
//                   WORK(INDA),
//                   LDA,
//                   HOUS(INDV),
//                   HOUS(INDTAU),
//                   LDV,
//                   WORK(INDW + TID * KD));
// // $OMP END TASK
//             } else {
// // $OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
// // $OMP$     DEPEND(out:WORK(MYID))
//               TID = OMP_GET_THREAD_NUM();
//               zhb2st_kernels(
//                   UPLO,
//                   WANTQ,
//                   TTYPE,
//                   STIND,
//                   EDIND,
//                   SWEEPID,
//                   N,
//                   KD,
//                   IB,
//                   WORK(INDA),
//                   LDA,
//                   HOUS(INDV),
//                   HOUS(INDTAU),
//                   LDV,
//                   WORK(INDW + TID * KD));
// // $OMP END TASK
//             }
// // #else
            zhb2st_kernels(
                UPLO,
                WANTQ,
                TTYPE,
                STIND,
                EDIND,
                SWEEPID,
                N,
                KD,
                IB,
                WORK(INDA).asMatrix(LDA),
                LDA,
                HOUS(INDV),
                HOUS(INDTAU),
                LDV,
                WORK(INDW));
// #endif
            if (BLKLASTIND >= (N - 1)) {
              STT++;
              break;
            }
          }
        }
      }
    }
  }

// #if defined(_OPENMP)
// $OMP END MASTER
// $OMP END PARALLEL
// #endif

  // Copy the diagonal from A to D. Note that D is REAL thus only
  // the Real part is needed, the imaginary part should be zero.

  for (I = 1; I <= N; I++) {
    D[I] = WORK[DPOS + (I - 1) * LDA].real;
  }

  // Copy the off diagonal from A to E. Note that E is REAL thus only
  // the Real part is needed, the imaginary part should be zero.

  if (UPPER) {
    for (I = 1; I <= N - 1; I++) {
      E[I] = WORK[OFDPOS + I * LDA].real;
    }
  } else {
    for (I = 1; I <= N - 1; I++) {
      E[I] = WORK[OFDPOS + (I - 1) * LDA].real;
    }
  }

  WORK[1] = LWMIN.toComplex();
}
