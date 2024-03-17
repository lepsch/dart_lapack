import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zladiv.dart';
import 'package:lapack/src/zlatrs.dart';

void zlaein(
  final bool RIGHTV,
  final bool NOINIT,
  final int N,
  final Matrix<Complex> H_,
  final int LDH,
  final Complex W,
  final Array<Complex> V_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> RWORK_,
  final double EPS3,
  final double SMLNUM,
  final Box<int> INFO,
) {
  final H = H_.having(ld: LDH);
  final B = B_.having(ld: LDB);
  final V = V_.having();
  final RWORK = RWORK_.having();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, TENTH = 1.0e-1;
  String NORMIN, TRANS;
  int I, ITS, J;
  double GROWTO, NRMSML, ROOTN, RTEMP, VNORM;
  Complex EI, EJ, TEMP, X;
  final SCALE = Box(0.0);
  final IERR = Box(0);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  INFO.value = 0;

  // GROWTO is the threshold used in the acceptance test for an
  // eigenvector.

  ROOTN = sqrt(N.toDouble());
  GROWTO = TENTH / ROOTN;
  NRMSML = max(ONE, EPS3 * ROOTN) * SMLNUM;

  // Form B = H - W*I (except that the subdiagonal elements are not
  // stored).

  for (J = 1; J <= N; J++) {
    // 20
    for (I = 1; I <= J - 1; I++) {
      // 10
      B[I][J] = H[I][J];
    } // 10
    B[J][J] = H[J][J] - W;
  } // 20

  if (NOINIT) {
    // Initialize V.

    for (I = 1; I <= N; I++) {
      // 30
      V[I] = EPS3.toComplex();
    } // 30
  } else {
    // Scale supplied initial vector.

    VNORM = dznrm2(N, V, 1);
    zdscal(N, (EPS3 * ROOTN) / max(VNORM, NRMSML), V, 1);
  }

  if (RIGHTV) {
    // LU decomposition with partial pivoting of B, replacing zero
    // pivots by EPS3.

    for (I = 1; I <= N - 1; I++) {
      // 60
      EI = H[I + 1][I];
      if (CABS1(B[I][I]) < CABS1(EI)) {
        // Interchange rows and eliminate.

        X = zladiv(B[I][I], EI);
        B[I][I] = EI;
        for (J = I + 1; J <= N; J++) {
          // 40
          TEMP = B[I + 1][J];
          B[I + 1][J] = B[I][J] - X * TEMP;
          B[I][J] = TEMP;
        } // 40
      } else {
        // Eliminate without interchange.

        if (B[I][I] == Complex.zero) B[I][I] = EPS3.toComplex();
        X = zladiv(EI, B[I][I]);
        if (X != Complex.zero) {
          for (J = I + 1; J <= N; J++) {
            // 50
            B[I + 1][J] -= X * B[I][J];
          } // 50
        }
      }
    } // 60
    if (B[N][N] == Complex.zero) B[N][N] = EPS3.toComplex();

    TRANS = 'N';
  } else {
    // UL decomposition with partial pivoting of B, replacing zero
    // pivots by EPS3.

    for (J = N; J >= 2; J--) {
      // 90
      EJ = H[J][J - 1];
      if (CABS1(B[J][J]) < CABS1(EJ)) {
        // Interchange columns and eliminate.

        X = zladiv(B[J][J], EJ);
        B[J][J] = EJ;
        for (I = 1; I <= J - 1; I++) {
          // 70
          TEMP = B[I][J - 1];
          B[I][J - 1] = B[I][J] - X * TEMP;
          B[I][J] = TEMP;
        } // 70
      } else {
        // Eliminate without interchange.

        if (B[J][J] == Complex.zero) B[J][J] = EPS3.toComplex();
        X = zladiv(EJ, B[J][J]);
        if (X != Complex.zero) {
          for (I = 1; I <= J - 1; I++) {
            // 80
            B[I][J - 1] -= X * B[I][J];
          } // 80
        }
      }
    } // 90
    if (B[1][1] == Complex.zero) B[1][1] = EPS3.toComplex();

    TRANS = 'C';
  }

  var solved = false;
  NORMIN = 'N';
  for (ITS = 1; ITS <= N; ITS++) {
    // 110

    // Solve U*x = scale*v for a right eigenvector
    //   or U**H *x = scale*v for a left eigenvector,
    // overwriting x on v.

    zlatrs('Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, V, SCALE, RWORK, IERR);
    NORMIN = 'Y';

    // Test for sufficient growth in the norm of v.

    VNORM = dzasum(N, V, 1);
    if (VNORM >= GROWTO * SCALE.value) {
      solved = true;
      break;
    }

    // Choose new orthogonal starting vector and try again.

    RTEMP = EPS3 / (ROOTN + ONE);
    V[1] = EPS3.toComplex();
    for (I = 2; I <= N; I++) {
      // 100
      V[I] = RTEMP.toComplex();
    } // 100
    V[N - ITS + 1] = V[N - ITS + 1] - (EPS3 * ROOTN).toComplex();
  } // 110

  if (!solved) {
    // Failure to find eigenvector in N iterations.
    INFO.value = 1;
  }

  // Normalize eigenvector.

  I = izamax(N, V, 1);
  zdscal(N, ONE / CABS1(V[I]), V, 1);
}
