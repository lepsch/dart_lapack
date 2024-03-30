import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlassq.dart';
import 'package:lapack/src/zrot.dart';

void ztgex2(
  final bool WANTQ,
  final bool WANTZ,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final int J1,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWENTY = 2.0e+1;
  const LDST = 2;
  const WANDS = true;
  bool STRONG, WEAK;
  int I, M;
  double EPS, SA, SB, SMLNUM, THRESHA, THRESHB;
  Complex F, G;
  final S = Matrix<Complex>(LDST, LDST),
      T = Matrix<Complex>(LDST, LDST),
      WORK = Array<Complex>(8);
  final SCALE = Box(0.0), SUM = Box(0.0), CQ = Box(0.0), CZ = Box(0.0);
  final SQ = Box(Complex.zero),
      CDUM = Box(Complex.zero),
      SZ = Box(Complex.zero);

  INFO.value = 0;

  // Quick return if possible

  if (N <= 1) return;

  M = LDST;
  WEAK = false;
  STRONG = false;

  // Make a local copy of selected block in (A, B)

  zlacpy('Full', M, M, A(J1, J1), LDA, S, LDST);
  zlacpy('Full', M, M, B(J1, J1), LDB, T, LDST);

  // Compute the threshold for testing the acceptance of swapping.

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  SCALE.value = ZERO;
  SUM.value = ONE;
  zlacpy('Full', M, M, S, LDST, WORK.asMatrix(M), M);
  zlacpy('Full', M, M, T, LDST, WORK(M * M + 1).asMatrix(M), M);
  zlassq(M * M, WORK, 1, SCALE, SUM);
  SA = SCALE.value * sqrt(SUM.value);
  SCALE.value = ZERO;
  SUM.value = ONE;
  zlassq(M * M, WORK(M * M + 1), 1, SCALE, SUM);
  SB = SCALE.value * sqrt(SUM.value);

  // THRES has been changed from
  //    THRESH = max( TEN*EPS*SA, SMLNUM )
  // to
  //    THRESH = max( TWENTY*EPS*SA, SMLNUM )
  // on 04/01/10.
  // "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
  // Jim Demmel and Guillaume Revy. See forum post 1783.

  THRESHA = max(TWENTY * EPS * SA, SMLNUM);
  THRESHB = max(TWENTY * EPS * SB, SMLNUM);

  // Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
  // using Givens rotations and perform the swap tentatively.

  F = S[2][2] * T[1][1] - T[2][2] * S[1][1];
  G = S[2][2] * T[1][2] - T[2][2] * S[1][2];
  SA = S[2][2].abs() * T[1][1].abs();
  SB = S[1][1].abs() * T[2][2].abs();
  zlartg(G, F, CZ, SZ, CDUM);
  SZ.value = -SZ.value;
  zrot(2, S(1, 1).asArray(), 1, S(1, 2).asArray(), 1, CZ.value,
      SZ.value.conjugate());
  zrot(2, T(1, 1).asArray(), 1, T(1, 2).asArray(), 1, CZ.value,
      SZ.value.conjugate());
  if (SA >= SB) {
    zlartg(S[1][1], S[2][1], CQ, SQ, CDUM);
  } else {
    zlartg(T[1][1], T[2][1], CQ, SQ, CDUM);
  }
  zrot(2, S(1, 1).asArray(), LDST, S(2, 1).asArray(), LDST, CQ.value, SQ.value);
  zrot(2, T(1, 1).asArray(), LDST, T(2, 1).asArray(), LDST, CQ.value, SQ.value);

  // Weak stability test: |S21| <= O(EPS F-norm((A)))
  //                      and  |T21| <= O(EPS F-norm((B)))

  WEAK = S[2][1].abs() <= THRESHA && T[2][1].abs() <= THRESHB;
  if (!WEAK) {
    // Exit with INFO = 1 if swap was rejected.
    INFO.value = 1;
    return;
  }

  if (WANDS) {
    // Strong stability test:
    //    F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
    //    and
    //    F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))

    zlacpy('Full', M, M, S, LDST, WORK.asMatrix(M), M);
    zlacpy('Full', M, M, T, LDST, WORK(M * M + 1).asMatrix(M), M);
    zrot(2, WORK, 1, WORK(3), 1, CZ.value, -SZ.value.conjugate());
    zrot(2, WORK(5), 1, WORK(7), 1, CZ.value, -SZ.value.conjugate());
    zrot(2, WORK, 2, WORK(2), 2, CQ.value, -SQ.value);
    zrot(2, WORK(5), 2, WORK(6), 2, CQ.value, -SQ.value);
    for (I = 1; I <= 2; I++) {
      WORK[I] -= A[J1 + I - 1][J1];
      WORK[I + 2] -= A[J1 + I - 1][J1 + 1];
      WORK[I + 4] -= B[J1 + I - 1][J1];
      WORK[I + 6] -= B[J1 + I - 1][J1 + 1];
    }
    SCALE.value = ZERO;
    SUM.value = ONE;
    zlassq(M * M, WORK, 1, SCALE, SUM);
    SA = SCALE.value * sqrt(SUM.value);
    SCALE.value = ZERO;
    SUM.value = ONE;
    zlassq(M * M, WORK(M * M + 1), 1, SCALE, SUM);
    SB = SCALE.value * sqrt(SUM.value);
    STRONG = SA <= THRESHA && SB <= THRESHB;
    if (!STRONG) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }
  }

  // If the swap is accepted ("weakly" and "strongly"), apply the
  // equivalence transformations to the original matrix pair (A,B)

  zrot(J1 + 1, A(1, J1).asArray(), 1, A(1, J1 + 1).asArray(), 1, CZ.value,
      SZ.value.conjugate());
  zrot(J1 + 1, B(1, J1).asArray(), 1, B(1, J1 + 1).asArray(), 1, CZ.value,
      SZ.value.conjugate());
  zrot(N - J1 + 1, A(J1, J1).asArray(), LDA, A(J1 + 1, J1).asArray(), LDA,
      CQ.value, SQ.value);
  zrot(N - J1 + 1, B(J1, J1).asArray(), LDB, B(J1 + 1, J1).asArray(), LDB,
      CQ.value, SQ.value);

  // Set  N1 by N2 (2,1) blocks to 0

  A[J1 + 1][J1] = Complex.zero;
  B[J1 + 1][J1] = Complex.zero;

  // Accumulate transformations into Q and Z if requested.

  if (WANTZ) {
    zrot(N, Z(1, J1).asArray(), 1, Z(1, J1 + 1).asArray(), 1, CZ.value,
        SZ.value.conjugate());
  }
  if (WANTQ) {
    zrot(N, Q(1, J1).asArray(), 1, Q(1, J1 + 1).asArray(), 1, CQ.value,
        SQ.value.conjugate());
  }

  // Exit with INFO = 0 if swap was successfully performed.
}
