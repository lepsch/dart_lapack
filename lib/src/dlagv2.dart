import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlag2.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlasv2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlagv2(
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final int LDB,
  final Array<double> ALPHAR,
  final Array<double> ALPHAI,
  final Array<double> BETA,
  final Box<double> CSL,
  final Box<double> SNL,
  final Box<double> CSR,
  final Box<double> SNR,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  double ANORM, ASCALE, BNORM, BSCALE, H1, H2, H3, QQ, RR, SAFMIN = 0, ULP;
  final R = Box(0.0),
      T = Box(0.0),
      SCALE1 = Box(0.0),
      SCALE2 = Box(0.0),
      WI = Box(0.0),
      WR1 = Box(0.0),
      WR2 = Box(0.0);

  SAFMIN = dlamch('S');
  ULP = dlamch('P');

  // Scale A

  ANORM = max(max(A[1][1].abs() + A[2][1].abs(), A[1][2].abs() + A[2][2].abs()),
      SAFMIN);
  ASCALE = ONE / ANORM;
  A[1][1] = ASCALE * A[1][1];
  A[1][2] = ASCALE * A[1][2];
  A[2][1] = ASCALE * A[2][1];
  A[2][2] = ASCALE * A[2][2];

  // Scale B

  BNORM = max(max(B[1][1].abs(), B[1][2].abs() + B[2][2].abs()), SAFMIN);
  BSCALE = ONE / BNORM;
  B[1][1] = BSCALE * B[1][1];
  B[1][2] = BSCALE * B[1][2];
  B[2][2] = BSCALE * B[2][2];

  // Check if A can be deflated

  if (A[2][1].abs() <= ULP) {
    CSL.value = ONE;
    SNL.value = ZERO;
    CSR.value = ONE;
    SNR.value = ZERO;
    A[2][1] = ZERO;
    B[2][1] = ZERO;
    WI.value = ZERO;

    // Check if B is singular
  } else if (B[1][1].abs() <= ULP) {
    dlartg(A[1][1], A[2][1], CSL, SNL, R);
    CSR.value = ONE;
    SNR.value = ZERO;
    drot(2, A(1, 1).asArray(), LDA, A(2, 1).asArray(), LDA, CSL.value,
        SNL.value);
    drot(2, B(1, 1).asArray(), LDB, B(2, 1).asArray(), LDB, CSL.value,
        SNL.value);
    A[2][1] = ZERO;
    B[1][1] = ZERO;
    B[2][1] = ZERO;
    WI.value = ZERO;
  } else if (B[2][2].abs() <= ULP) {
    dlartg(A[2][2], A[2][1], CSR, SNR, T);
    SNR.value = -SNR.value;
    drot(2, A(1, 1).asArray(), 1, A(1, 2).asArray(), 1, CSR.value, SNR.value);
    drot(2, B(1, 1).asArray(), 1, B(1, 2).asArray(), 1, CSR.value, SNR.value);
    CSL.value = ONE;
    SNL.value = ZERO;
    A[2][1] = ZERO;
    B[2][1] = ZERO;
    B[2][2] = ZERO;
    WI.value = ZERO;
  } else {
    // B is nonsingular, first compute the eigenvalues of (A,B)

    dlag2(A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2, WI);

    if (WI.value == ZERO) {
      // two real eigenvalues, compute s*A-w*B

      H1 = SCALE1.value * A[1][1] - WR1.value * B[1][1];
      H2 = SCALE1.value * A[1][2] - WR1.value * B[1][2];
      H3 = SCALE1.value * A[2][2] - WR1.value * B[2][2];

      RR = dlapy2(H1, H2);
      QQ = dlapy2(SCALE1.value * A[2][1], H3);

      if (RR > QQ) {
        // find right rotation matrix to zero 1,1 element of
        // (sA - wB)

        dlartg(H2, H1, CSR, SNR, T);
      } else {
        // find right rotation matrix to zero 2,1 element of
        // (sA - wB)

        dlartg(H3, SCALE1.value * A[2][1], CSR, SNR, T);
      }

      SNR.value = -SNR.value;
      drot(2, A(1, 1).asArray(), 1, A(1, 2).asArray(), 1, CSR.value, SNR.value);
      drot(2, B(1, 1).asArray(), 1, B(1, 2).asArray(), 1, CSR.value, SNR.value);

      // compute inf norms of A and B

      H1 = max(A[1][1].abs() + A[1][2].abs(), A[2][1].abs() + A[2][2].abs());

      H2 = max(B[1][1].abs() + B[1][2].abs(), B[2][1].abs() + B[2][2].abs());

      if ((SCALE1.value * H1) >= (WR1.value).abs() * H2) {
        // find left rotation matrix Q to zero out B(2,1)

        dlartg(B[1][1], B[2][1], CSL, SNL, R);
      } else {
        // find left rotation matrix Q to zero out A(2,1)

        dlartg(A[1][1], A[2][1], CSL, SNL, R);
      }

      drot(2, A(1, 1).asArray(), LDA, A(2, 1).asArray(), LDA, CSL.value,
          SNL.value);
      drot(2, B(1, 1).asArray(), LDB, B(2, 1).asArray(), LDB, CSL.value,
          SNL.value);

      A[2][1] = ZERO;
      B[2][1] = ZERO;
    } else {
      // a pair of complex conjugate eigenvalues
      // first compute the SVD of the matrix B

      dlasv2(B[1][1], B[1][2], B[2][2], R, T, SNR, CSR, SNL, CSL);

      // Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and
      // Z is right rotation matrix computed from DLASV2

      drot(2, A(1, 1).asArray(), LDA, A(2, 1).asArray(), LDA, CSL.value,
          SNL.value);
      drot(2, B(1, 1).asArray(), LDB, B(2, 1).asArray(), LDB, CSL.value,
          SNL.value);
      drot(2, A(1, 1).asArray(), 1, A(1, 2).asArray(), 1, CSR.value, SNR.value);
      drot(2, B(1, 1).asArray(), 1, B(1, 2).asArray(), 1, CSR.value, SNR.value);

      B[2][1] = ZERO;
      B[1][2] = ZERO;
    }
  }

  // Unscaling

  A[1][1] = ANORM * A[1][1];
  A[2][1] = ANORM * A[2][1];
  A[1][2] = ANORM * A[1][2];
  A[2][2] = ANORM * A[2][2];
  B[1][1] = BNORM * B[1][1];
  B[2][1] = BNORM * B[2][1];
  B[1][2] = BNORM * B[1][2];
  B[2][2] = BNORM * B[2][2];

  if (WI.value == ZERO) {
    ALPHAR[1] = A[1][1];
    ALPHAR[2] = A[2][2];
    ALPHAI[1] = ZERO;
    ALPHAI[2] = ZERO;
    BETA[1] = B[1][1];
    BETA[2] = B[2][2];
  } else {
    ALPHAR[1] = ANORM * WR1.value / SCALE1.value / BNORM;
    ALPHAI[1] = ANORM * WI.value / SCALE1.value / BNORM;
    ALPHAR[2] = ALPHAR[1];
    ALPHAI[2] = -ALPHAI[1];
    BETA[1] = ONE;
    BETA[2] = ONE;
  }
}
