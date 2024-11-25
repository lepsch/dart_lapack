// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/nio.dart';

void alahd(final Nout IOUNIT, final String PATH) {
  String SUBNAM;
  String SYM;

  final C1 = PATH.substring(0, 1);
  final C3 = PATH.substring(2, 3);
  final P2 = PATH.substring(1, 3);
  final SORD = lsame(C1, 'S') || lsame(C1, 'D');
  final CORZ = lsame(C1, 'C') || lsame(C1, 'Z');
  if (!(SORD || CORZ)) return;

  if (lsamen(2, P2, 'GE')) {
    // GE: General dense

    IOUNIT.println('\n ${PATH.a3}:  General dense matrices');
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9979();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9962(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9959(4);
    IOUNIT.print9958(5);
    IOUNIT.print9957(6);
    IOUNIT.print9956(7);
    IOUNIT.print9955(8);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'GB')) {
    // GB: General band

    IOUNIT.print9998(PATH);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9978();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9962(1);
    IOUNIT.print9960(2);
    IOUNIT.print9959(3);
    IOUNIT.print9958(4);
    IOUNIT.print9957(5);
    IOUNIT.print9956(6);
    IOUNIT.print9955(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'GT')) {
    // GT: General tridiagonal

    IOUNIT.print9997(PATH);
    IOUNIT.print9977();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9962(1);
    IOUNIT.print9960(2);
    IOUNIT.print9959(3);
    IOUNIT.print9958(4);
    IOUNIT.print9957(5);
    IOUNIT.print9956(6);
    IOUNIT.print9955(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PO') || lsamen(2, P2, 'PP')) {
    // PO: Positive definite full
    // PP: Positive definite packed

    if (SORD) {
      SYM = 'Symmetric';
    } else {
      SYM = 'Hermitian';
    }
    if (lsame(C3, 'O')) {
      IOUNIT.print9996(PATH, SYM);
    } else {
      IOUNIT.print9995(PATH, SYM);
    }
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9975(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9954(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9959(4);
    IOUNIT.print9958(5);
    IOUNIT.print9957(6);
    IOUNIT.print9956(7);
    IOUNIT.print9955(8);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PS')) {
    // PS: Positive semi-definite full

    if (SORD) {
      SYM = 'Symmetric';
    } else {
      SYM = 'Hermitian';
    }
    final EIGCNM = (lsame(C1, 'S') || lsame(C1, 'C')) ? '1E04' : '1D12';

    IOUNIT.print9995(PATH, SYM);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print8973(EIGCNM);
    IOUNIT.println(' Difference:');
    IOUNIT.print8972(C1);
    IOUNIT.println(' Test ratio:');
    IOUNIT.print8950();
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PB')) {
    // PB: Positive definite band

    if (SORD) {
      IOUNIT.print9994(PATH, 'Symmetric');
    } else {
      IOUNIT.print9994(PATH, 'Hermitian');
    }
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9973(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9954(1);
    IOUNIT.print9960(2);
    IOUNIT.print9959(3);
    IOUNIT.print9958(4);
    IOUNIT.print9957(5);
    IOUNIT.print9956(6);
    IOUNIT.print9955(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PT')) {
    // PT: Positive definite tridiagonal

    if (SORD) {
      IOUNIT.print9993(PATH, 'Symmetric');
    } else {
      IOUNIT.print9993(PATH, 'Hermitian');
    }
    IOUNIT.print9976();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9952(1);
    IOUNIT.print9960(2);
    IOUNIT.print9959(3);
    IOUNIT.print9958(4);
    IOUNIT.print9957(5);
    IOUNIT.print9956(6);
    IOUNIT.print9955(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'SY')) {
    // SY: Symmetric indefinite full,
    //     with partial (Bunch-Kaufman) pivoting algorithm

    if (lsame(C3, 'Y')) {
      IOUNIT.print9992(PATH, 'Symmetric');
    } else {
      IOUNIT.print9991(PATH, 'Symmetric');
    }
    IOUNIT.println(' Matrix types:');
    if (SORD) {
      IOUNIT.print9972();
    } else {
      IOUNIT.print9971();
    }
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9960(4);
    IOUNIT.print9959(5);
    IOUNIT.print9958(6);
    IOUNIT.print9956(7);
    IOUNIT.print9957(8);
    IOUNIT.print9955(9);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'SR') || lsamen(2, P2, 'SK')) {
    // SR: Symmetric indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm

    // SK: Symmetric indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm,
    //     ( new storage format for factors:
    //       L and diagonal of D is stored in A,
    //       subdiagonal of D is stored in E )

    IOUNIT.print9892(PATH, 'Symmetric');

    IOUNIT.println(' Matrix types:');
    if (SORD) {
      IOUNIT.print9972();
    } else {
      IOUNIT.print9971();
    }

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9927(3);
    IOUNIT.print9928();
    IOUNIT.print9926(4);
    IOUNIT.print9928();
    IOUNIT.print9960(5);
    IOUNIT.print9959(6);
    IOUNIT.print9955(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'SP')) {
    // SP: Symmetric indefinite packed,
    //     with partial (Bunch-Kaufman) pivoting algorithm

    if (lsame(C3, 'Y')) {
      IOUNIT.print9992(PATH, 'Symmetric');
    } else {
      IOUNIT.print9991(PATH, 'Symmetric');
    }
    IOUNIT.println(' Matrix types:');
    if (SORD) {
      IOUNIT.print9972();
    } else {
      IOUNIT.print9971();
    }
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9959(4);
    IOUNIT.print9958(5);
    IOUNIT.print9956(6);
    IOUNIT.print9957(7);
    IOUNIT.print9955(8);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HA')) {
    // HA: Hermitian,
    //     with Assen Algorithm

    IOUNIT.print9992(PATH, 'Hermitian');

    IOUNIT.println(' Matrix types:');
    IOUNIT.print9972();

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9960(4);
    IOUNIT.print9959(5);
    IOUNIT.print9958(6);
    IOUNIT.print9956(7);
    IOUNIT.print9957(8);
    IOUNIT.print9955(9);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HE')) {
    // HE: Hermitian indefinite full,
    //     with partial (Bunch-Kaufman) pivoting algorithm

    IOUNIT.print9992(PATH, 'Hermitian');

    IOUNIT.println(' Matrix types:');
    IOUNIT.print9972();

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9960(4);
    IOUNIT.print9959(5);
    IOUNIT.print9958(6);
    IOUNIT.print9956(7);
    IOUNIT.print9957(8);
    IOUNIT.print9955(9);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HR') || lsamen(2, P2, 'HR')) {
    // HR: Hermitian indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm

    // HK: Hermitian indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm,
    //     ( new storage format for factors:
    //       L and diagonal of D is stored in A,
    //       subdiagonal of D is stored in E )

    IOUNIT.print9892(PATH, 'Hermitian');

    IOUNIT.println(' Matrix types:');
    IOUNIT.print9972();

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9927(3);
    IOUNIT.print9928();
    IOUNIT.print9926(4);
    IOUNIT.print9928();
    IOUNIT.print9960(5);
    IOUNIT.print9959(6);
    IOUNIT.print9955(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HP')) {
    // HP: Hermitian indefinite packed,
    //     with partial (Bunch-Kaufman) pivoting algorithm

    if (lsame(C3, 'E')) {
      IOUNIT.print9992(PATH, 'Hermitian');
    } else {
      IOUNIT.print9991(PATH, 'Hermitian');
    }
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9972();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9953(1);
    IOUNIT.print9961(2);
    IOUNIT.print9960(3);
    IOUNIT.print9959(4);
    IOUNIT.print9958(5);
    IOUNIT.print9956(6);
    IOUNIT.print9957(7);
    IOUNIT.print9955(8);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'TR') || lsamen(2, P2, 'TP')) {
    // TR: Triangular full
    // TP: Triangular packed

    if (lsame(C3, 'R')) {
      IOUNIT.print9990(PATH);
      SUBNAM = '${PATH[0]}LATRS';
    } else {
      IOUNIT.print9989(PATH);
      SUBNAM = '${PATH[0]}LATPS';
    }
    IOUNIT.print9966(PATH);
    IOUNIT.print9965(SUBNAM.trim());
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9961(1);
    IOUNIT.print9960(2);
    IOUNIT.print9959(3);
    IOUNIT.print9958(4);
    IOUNIT.print9957(5);
    IOUNIT.print9956(6);
    IOUNIT.print9955(7);
    IOUNIT.print9951(SUBNAM.trim(), 8);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'TB')) {
    // TB: Triangular band

    IOUNIT.print9988(PATH);
    SUBNAM = '${PATH[0]}LATBS';
    IOUNIT.print9964(PATH);
    IOUNIT.print9963(SUBNAM.trim());
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9960(1);
    IOUNIT.print9959(2);
    IOUNIT.print9958(3);
    IOUNIT.print9957(4);
    IOUNIT.print9956(5);
    IOUNIT.print9955(6);
    IOUNIT.print9951(SUBNAM.trim(), 7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'QR')) {
    // QR decomposition of rectangular matrices

    IOUNIT.print9987(PATH, 'QR');
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9970();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9950(1);
    IOUNIT.print6950(8);
    IOUNIT.print9946(2);
    IOUNIT.print9944(3, 'M');
    IOUNIT.print9943(4, 'M');
    IOUNIT.print9942(5, 'M');
    IOUNIT.print9941(6, 'M');
    IOUNIT.print9960(7);
    IOUNIT.print6660(9);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'LQ')) {
    // LQ decomposition of rectangular matrices

    IOUNIT.print9987(PATH, 'LQ');
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9970();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9949(1);
    IOUNIT.print9945(2);
    IOUNIT.print9944(3, 'N');
    IOUNIT.print9943(4, 'N');
    IOUNIT.print9942(5, 'N');
    IOUNIT.print9941(6, 'N');
    IOUNIT.print9960(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'QL')) {
    // QL decomposition of rectangular matrices

    IOUNIT.print9987(PATH, 'QL');
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9970();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9948(1);
    IOUNIT.print9946(2);
    IOUNIT.print9944(3, 'M');
    IOUNIT.print9943(4, 'M');
    IOUNIT.print9942(5, 'M');
    IOUNIT.print9941(6, 'M');
    IOUNIT.print9960(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'RQ')) {
    // RQ decomposition of rectangular matrices

    IOUNIT.print9987(PATH, 'RQ');
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9970();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9947(1);
    IOUNIT.print9945(2);
    IOUNIT.print9944(3, 'N');
    IOUNIT.print9943(4, 'N');
    IOUNIT.print9942(5, 'N');
    IOUNIT.print9941(6, 'N');
    IOUNIT.print9960(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'QP')) {
    // QR decomposition with column pivoting

    IOUNIT.print8006(PATH);
    IOUNIT.print9969();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9940(1);
    IOUNIT.print9939(2);
    IOUNIT.print9938(3);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'QK')) {
    // truncated QR decomposition with column pivoting

    IOUNIT.print8006(PATH);
    IOUNIT.print9871();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8060(1);
    IOUNIT.print8061(2);
    IOUNIT.print8062(3);
    IOUNIT.print8063(4);
    IOUNIT.print8064(5);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'TZ')) {
    // TZ:  Trapezoidal

    IOUNIT.print9985(PATH);
    IOUNIT.print9968();
    IOUNIT.print9929(C1);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9940(1);
    IOUNIT.print9937(2);
    IOUNIT.print9938(3);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'LS')) {
    // LS:  Least Squares driver routines for
    //      LS, LST, TSLS, LSD, LSS, LSX and LSY.

    IOUNIT.print9984(PATH);
    IOUNIT.print9967();
    IOUNIT.print9921(C1);
    IOUNIT.print9935(1);
    IOUNIT.print9931(2);
    IOUNIT.print9919();
    IOUNIT.print9933(7);
    IOUNIT.print9935(8);
    IOUNIT.print9934(9);
    IOUNIT.print9932(10);
    IOUNIT.print9920();
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'LU')) {
    // LU factorization variants

    IOUNIT.print9983(PATH);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9979();
    IOUNIT.println(' Test ratio:');
    IOUNIT.print9962(1);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'CH')) {
    // Cholesky factorization variants

    IOUNIT.print9982(PATH);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9974();
    IOUNIT.println(' Test ratio:');
    IOUNIT.print9954(1);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'QS')) {
    // QR factorization variants

    IOUNIT.print9981(PATH);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9970();
    IOUNIT.println(' Test ratios:');
  } else if (lsamen(2, P2, 'QT')) {
    // QRT (general matrices)

    IOUNIT.print8000(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8011(1);
    IOUNIT.print8012(2);
    IOUNIT.print8013(3);
    IOUNIT.print8014(4);
    IOUNIT.print8015(5);
    IOUNIT.print8016(6);
  } else if (lsamen(2, P2, 'QX')) {
    // QRT (triangular-pentagonal)

    IOUNIT.print8001(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8017(1);
    IOUNIT.print8018(2);
    IOUNIT.print8019(3);
    IOUNIT.print8020(4);
    IOUNIT.print8021(5);
    IOUNIT.print8022(6);
  } else if (lsamen(2, P2, 'TQ')) {
    // QRT (triangular-pentagonal)

    IOUNIT.print8002(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8023(1);
    IOUNIT.print8024(2);
    IOUNIT.print8025(3);
    IOUNIT.print8026(4);
    IOUNIT.print8027(5);
    IOUNIT.print8028(6);
  } else if (lsamen(2, P2, 'XQ')) {
    // QRT (triangular-pentagonal)

    IOUNIT.print8003(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8029(1);
    IOUNIT.print8030(2);
    IOUNIT.print8031(3);
    IOUNIT.print8032(4);
    IOUNIT.print8033(5);
    IOUNIT.print8034(6);
  } else if (lsamen(2, P2, 'TS')) {
    // TS:  QR routines for tall-skinny and short-wide matrices

    IOUNIT.print8004(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8035(1);
    IOUNIT.print8036(2);
    IOUNIT.print8037(3);
    IOUNIT.print8038(4);
    IOUNIT.print8039(5);
    IOUNIT.print8040(6);
  } else if (lsamen(2, P2, 'HH')) {
    // HH:  Householder reconstruction for tall-skinny matrices

    IOUNIT.print8005(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print8050(1);
    IOUNIT.print8051(2);
    IOUNIT.print8052(3);
    IOUNIT.print8053(4);
    IOUNIT.print8054(5);
    IOUNIT.print8055(6);
  } else {
    // void Print error message if no header is available.

    IOUNIT.print9980(PATH);
  }

  // First line of header
}

extension on Nout {
  // GE matrix types

  void print9979() {
    println(
        '    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero');
  }

  void print9998(String PATH) {
    println('\n ${PATH.a3}:  General band matrices');
  }

  void print9997(String PATH) {
    println('\n ${PATH.a3}:  General tridiagonal');
  }

  void print9996(String PATH, String SYM) {
    println('\n ${PATH.a3}:  ${SYM.a9} positive definite matrices');
  }

  void print9995(String PATH, String SYM) {
    println('\n ${PATH.a3}:  ${SYM.a9} positive definite packed matrices');
  }

  void print9994(String PATH, String SYM) {
    println('\n ${PATH.a3}:  ${SYM.a9} positive definite band matrices');
  }

  void print9993(String PATH, String SYM) {
    println('\n ${PATH.a3}:  ${SYM.a9} positive definite tridiagonal');
  }

  void print9992(String PATH, String SYM) {
    println(
        '\n ${PATH.a3}:  ${SYM.a9} indefinite matrices, partial (Bunch-Kaufman) pivoting');
  }

  void print9991(String PATH, String SYM) {
    println(
        '\n ${PATH.a3}:  ${SYM.a9} indefinite packed matrices, partial (Bunch-Kaufman) pivoting');
  }

  void print9892(String PATH, String SYM) {
    println(
        '\n ${PATH.a3}:  ${SYM.a9} indefinite matrices, "rook" (bounded Bunch-Kaufman) pivoting');
  }

  // void print9891(String PATH, String SYM) {
  //   println(
  //       '\n ${PATH.a3}:  ${SYM.a9} indefinite packed matrices, "rook" (bounded Bunch-Kaufman) pivoting');
  // }

  void print9990(String PATH) {
    println('\n ${PATH.a3}:  Triangular matrices');
  }

  void print9989(String PATH) {
    println('\n ${PATH.a3}:  Triangular packed matrices');
  }

  void print9988(String PATH) {
    println('\n ${PATH.a3}:  Triangular band matrices');
  }

  void print9987(String PATH, String TYPE) {
    println('\n ${PATH.a3}:  ${TYPE.a2} factorization of general matrices');
  }

  // void print9986(String PATH) {
  //   println('\n ${PATH.a3}:  QR factorization with column pivoting');
  // }

  void print9985(String PATH) {
    println('\n ${PATH.a3}:  RQ factorization of trapezoidal matrix');
  }

  void print9984(String PATH) {
    println('\n ${PATH.a3}:  Least squares driver routines');
  }

  void print9983(String PATH) {
    println('\n ${PATH.a3}:  LU factorization variants');
  }

  void print9982(String PATH) {
    println('\n ${PATH.a3}:  Cholesky factorization variants');
  }

  void print9981(String PATH) {
    println('\n ${PATH.a3}:  QR factorization variants');
  }

  void print9980(String PATH) {
    println('\n ${PATH.a3}:  No header available');
  }

  void print8000(String PATH) {
    println('\n ${PATH.a3}:  QRT factorization for general matrices');
  }

  void print8001(String PATH) {
    println(
        '\n ${PATH.a3}:  QRT factorization for triangular-pentagonal matrices');
  }

  void print8002(String PATH) {
    println('\n ${PATH.a3}:  LQT factorization for general matrices');
  }

  void print8003(String PATH) {
    println(
        '\n ${PATH.a3}:  LQT factorization for triangular-pentagonal matrices');
  }

  void print8004(String PATH) {
    println(
        '\n ${PATH.a3}:  TS factorization for tall-skinny or short-wide matrices');
  }

  void print8005(String PATH) {
    println(
        '\n ${PATH.a3}:  Householder reconstruction from TSQR factorization output \n for tall-skinny matrices.');
  }

  void print8006(String PATH) {
    println('\n ${PATH.a3}:  truncated QR factorization with column pivoting');
  }

  // GB matrix types

  void print9978() {
    println(
        '    1. Random, CNDNUM = 2${' ' * 14}5. Random, CNDNUM = sqrt(0.1/EPS)\n    2. First column zero${' ' * 15}6. Random, CNDNUM = .01/EPS\n    3. Last column zero${' ' * 16}7. Scaled near underflow\n    4. Last n/2 columns zero${' ' * 11}8. Scaled near overflow');
  }

  // GT matrix types

  void print9977() {
    println(
        ' Matrix types (1-6 have specified condition numbers):\n    1. Diagonal${' ' * 24}7. Random, unspecified CNDNUM\n    2. Random, CNDNUM = 2${' ' * 14}8. First column zero\n    3. Random, CNDNUM = sqrt(0.1/EPS)  9. Last column zero\n    4. Random, CNDNUM = 0.1/EPS${' ' * 7}10. Last n/2 columns zero\n    5. Scaled near underflow${' ' * 10}11. Scaled near underflow\n    6. Scaled near overflow${' ' * 11}12. Scaled near overflow');
  }

  // PT matrix types

  void print9976() {
    println(
        ' Matrix types (1-6 have specified condition numbers):\n    1. Diagonal${' ' * 24}7. Random, unspecified CNDNUM\n    2. Random, CNDNUM = 2${' ' * 14}8. First row and column zero\n    3. Random, CNDNUM = sqrt(0.1/EPS)  9. Last row and column zero\n    4. Random, CNDNUM = 0.1/EPS${' ' * 7}10. Middle row and column zero\n    5. Scaled near underflow${' ' * 10}11. Scaled near underflow\n    6. Scaled near overflow${' ' * 11}12. Scaled near overflow');
  }

  // PO, PP matrix types

  void print9975(String PATH) {
    println(
        '    1. Diagonal${' ' * 24}6. Random, CNDNUM = sqrt(0.1/EPS)\n    2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = 0.1/EPS\n   *3. First row and column zero${' ' * 7}8. Scaled near underflow\n   *4. Last row and column zero${' ' * 8}9. Scaled near overflow\n   *5. Middle row and column zero\n   (* - tests error exits from ${PATH.a3}TRF, no test ratios are computed)');
  }

  // CH matrix types

  void print9974() {
    println(
        '    1. Diagonal${' ' * 24}6. Random, CNDNUM = sqrt(0.1/EPS)\n    2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = 0.1/EPS\n   *3. First row and column zero${' ' * 7}8. Scaled near underflow\n   *4. Last row and column zero${' ' * 8}9. Scaled near overflow\n   *5. Middle row and column zero\n   (* - tests error exits, no test ratios are computed)');
  }

  // PS matrix types

  void print8973(String EIGCNM) {
    println(
        '    1. Diagonal\n    2. Random, CNDNUM = 2${' ' * 14}\n   *3. Nonzero eigenvalues of: D(1:RANK-1)=1 and D(RANK) = 1.0/${EIGCNM.a4}\n   *4. Nonzero eigenvalues of: D(1)=1 and  D(2:RANK) = 1.0/${EIGCNM.a4}\n   *5. Nonzero eigenvalues of: D(I) = ${EIGCNM.a4}**(-(I-1)/(RANK-1))  I=1:RANK\n    6. Random, CNDNUM = sqrt(0.1/EPS)\n    7. Random, CNDNUM = 0.1/EPS\n    8. Scaled near underflow\n    9. Scaled near overflow\n   (* - Semi-definite tests )');
  }

  void print8972(String prefix) {
    println('   RANK minus computed rank, returned by ${prefix}PSTRF');
  }

  // PB matrix types

  void print9973(String PATH) {
    println(
        '    1. Random, CNDNUM = 2${' ' * 14}5. Random, CNDNUM = sqrt(0.1/EPS)\n   *2. First row and column zero${' ' * 7}6. Random, CNDNUM = 0.1/EPS\n   *3. Last row and column zero${' ' * 8}7. Scaled near underflow\n   *4. Middle row and column zero${' ' * 6}8. Scaled near overflow\n   (* - tests error exits from ${PATH.a3}TRF, no test ratios are computed)');
  }

  // SSY, SSR, SSP, CHE, CHR, CHP matrix types

  void print9972() {
    println(
        '    1. Diagonal${' ' * 24}6. Last n/2 rows and columns zero\n    2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = sqrt(0.1/EPS)\n    3. First row and column zero${' ' * 7}8. Random, CNDNUM = 0.1/EPS\n    4. Last row and column zero${' ' * 8}9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow');
  }

  // CSY, CSR, CSP matrix types

  void print9971() {
    println(
        '    1. Diagonal${' ' * 24}7. Random, CNDNUM = sqrt(0.1/EPS)\n    2. Random, CNDNUM = 2${' ' * 14}8. Random, CNDNUM = 0.1/EPS\n    3. First row and column zero${' ' * 7}9. Scaled near underflow\n    4. Last row and column zero${' ' * 7}10. Scaled near overflow\n    5. Middle row and column zero     11. Block diagonal matrix\n    6. Last n/2 rows and columns zero');
  }

  // QR matrix types

  void print9970() {
    println(
        '    1. Diagonal${' ' * 24}5. Random, CNDNUM = sqrt(0.1/EPS)\n    2. Upper triangular${' ' * 16}6. Random, CNDNUM = 0.1/EPS\n    3. Lower triangular${' ' * 16}7. Scaled near underflow\n    4. Random, CNDNUM = 2${' ' * 14}8. Scaled near overflow');
  }

  // QP matrix types

  void print9969() {
    println(
        ' Matrix types (2-6 have condition 1/EPS):\n    1. Zero matrix${' ' * 21}4. First n/2 columns fixed\n    2. One small eigenvalue${' ' * 12}5. Last n/2 columns fixed\n    3. Geometric distribution${' ' * 10}6. Every second column fixed');
  }

  // QK matrix types

  void print9871() {
    println(
        '     1. Zero matrix\n     2. Random, Diagonal, CNDNUM = 2\n     3. Random, Upper triangular, CNDNUM = 2\n     4. Random, Lower triangular, CNDNUM = 2\n     5. Random, First column is zero, CNDNUM = 2\n     6. Random, Last MINMN column is zero, CNDNUM = 2\n     7. Random, Last N column is zero, CNDNUM = 2\n     8. Random, Middle column in MINMN is zero, CNDNUM = 2\n     9. Random, First half of MINMN columns are zero, CNDNUM = 2\n    10. Random, Last columns are zero starting from MINMN/2+1, CNDNUM = 2\n    11. Random, Half MINMN columns in the middle are zero starting from MINMN/2-(MINMN/2)/2+1, CNDNUM = 2\n    12. Random, Odd columns are ZERO, CNDNUM = 2\n    13. Random, Even columns are ZERO, CNDNUM = 2\n    14. Random, CNDNUM = 2\n    15. Random, CNDNUM = sqrt(0.1/EPS)\n    16. Random, CNDNUM = 0.1/EPS\n    17. Random, CNDNUM = 0.1/EPS, one small singular value S(N)=1/CNDNUM\n    18. Random, CNDNUM = 2, scaled near underflow, NORM = SMALL = SAFMIN\n    19. Random, CNDNUM = 2, scaled near overflow, NORM = LARGE = 1.0/( 0.25 * ( SAFMIN / EPS ) )');
  }

  // TZ matrix types

  void print9968() {
    println(
        ' Matrix types (2-3 have condition 1/EPS):\n    1. Zero matrix\n    2. One small eigenvalue\n    3. Geometric distribution');
  }

  // LS matrix types

  void print9967() {
    println(
        ' Matrix types (1-3: full rank, 4-6: rank deficient):\n    1 and 4. Normal scaling\n    2 and 5. Scaled near overflow\n    3 and 6. Scaled near underflow');
  }

  // TR, TP matrix types

  void print9966(String PATH) {
    println(
        ' Matrix types for ${PATH.a3} routines:\n    1. Diagonal${' ' * 24}6. Scaled near overflow\n    2. Random, CNDNUM = 2${' ' * 14}7. Identity\n    3. Random, CNDNUM = sqrt(0.1/EPS)  8. Unit triangular, CNDNUM = 2\n    4. Random, CNDNUM = 0.1/EPS${' ' * 8}9. Unit, CNDNUM = sqrt(0.1/EPS)\n    5. Scaled near underflow${' ' * 10}10. Unit, CNDNUM = 0.1/EPS');
  }

  void print9965(String SUBNAM) {
    println(
        ' Special types for testing $SUBNAM:\n   11. Matrix elements are O(1), large right hand side\n   12. First diagonal causes overflow, offdiagonal column norms < 1\n   13. First diagonal causes overflow, offdiagonal column norms > 1\n   14. Growth factor underflows, solution does not overflow\n   15. Small diagonal causes gradual overflow\n   16. One zero diagonal element\n   17. Large offdiagonals cause overflow when adding a column\n   18. Unit triangular with large right hand side');
  }

  // TB matrix types

  void print9964(String PATH) {
    println(
        ' Matrix types for ${PATH.a3} routines:\n    1. Random, CNDNUM = 2${' ' * 14}6. Identity\n    2. Random, CNDNUM = sqrt(0.1/EPS)  7. Unit triangular, CNDNUM = 2\n    3. Random, CNDNUM = 0.1/EPS${' ' * 8}8. Unit, CNDNUM = sqrt(0.1/EPS)\n    4. Scaled near underflow${' ' * 11}9. Unit, CNDNUM = 0.1/EPS\n    5. Scaled near overflow');
  }

  void print9963(String SUBNAM) {
    println(
        ' Special types for testing $SUBNAM:\n   10. Matrix elements are O(1), large right hand side\n   11. First diagonal causes overflow, offdiagonal column norms < 1\n   12. First diagonal causes overflow, offdiagonal column norms > 1\n   13. Growth factor underflows, solution does not overflow\n   14. Small diagonal causes gradual overflow\n   15. One zero diagonal element\n   16. Large offdiagonals cause overflow when adding a column\n   17. Unit triangular with large right hand side');
  }

  // Test ratios

  void print9962(int i) {
    println('   ${i.i2}: norm( L * U - A )  / ( N * norm(A) * EPS )');
  }

  void print9961(int i) {
    println(
        '   ${i.i2}: norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS )');
  }

  void print9960(int i) {
    println('   ${i.i2}: norm( B - A * X )  / ( norm(A) * norm(X) * EPS )');
  }

  void print6660(int i) {
    println('   ${i.i2}: diagonal is not non-negative');
  }

  void print9959(int i) {
    println('   ${i.i2}: norm( X - XACT )   / ( norm(XACT) * CNDNUM * EPS )');
  }

  void print9958(int i) {
    println(
        '   ${i.i2}: norm( X - XACT )   / ( norm(XACT) * CNDNUM * EPS ), refined');
  }

  void print9957(int i) {
    println('   ${i.i2}: norm( X - XACT )   / ( norm(XACT) * (error bound) )');
  }

  void print9956(int i) {
    println('   ${i.i2}: (backward error)   / EPS');
  }

  void print9955(int i) {
    println('   ${i.i2}: RCOND * CNDNUM - 1.0');
  }

  void print9954(int i) {
    println(
        '   ${i.i2}: norm( U\' * U - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L * L\' - A ) / ( N * norm(A) * EPS )');
  }

  void print8950() {
    println(
        '   norm( P * U\' * U * P\' - A ) / ( N * norm(A) * EPS ), or\n   norm( P * L * L\' * P\' - A ) / ( N * norm(A) * EPS )');
  }

  void print9953(int i) {
    println(
        '   ${i.i2}: norm( U*D*U\' - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L*D*L\' - A ) / ( N * norm(A) * EPS )');
  }

  void print9952(int i) {
    println(
        '   ${i.i2}: norm( U\'*D*U - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L*D*L\' - A ) / ( N * norm(A) * EPS )');
  }

  void print9951(String SUBNAM, int i) {
    println(
        ' Test ratio for $SUBNAM:\n   ${i.i2}: norm( s*b - A*x )  / ( norm(A) * norm(x) * EPS )');
  }

  void print9950(int i) {
    println('   ${i.i2}: norm( R - Q\' * A ) / ( M * norm(A) * EPS )');
  }

  void print6950(int i) {
    println('   ${i.i2}: norm( R - Q\' * A ) / ( M * norm(A) * EPS ) [RFPG]');
  }

  void print9949(int i) {
    println('   ${i.i2}: norm( L - A * Q\' ) / ( N * norm(A) * EPS )');
  }

  void print9948(int i) {
    println('   ${i.i2}: norm( L - Q\' * A ) / ( M * norm(A) * EPS )');
  }

  void print9947(int i) {
    println('   ${i.i2}: norm( R - A * Q\' ) / ( N * norm(A) * EPS )');
  }

  void print9946(int i) {
    println('   ${i.i2}: norm( I - Q\'*Q )   / ( M * EPS )');
  }

  void print9945(int i) {
    println('   ${i.i2}: norm( I - Q*Q\' )   / ( N * EPS )');
  }

  void print9944(int i, String name) {
    println('   ${i.i2}: norm( Q*C - Q*C )  / ( ${name.a1} * norm(C) * EPS )');
  }

  void print9943(int i, String name) {
    println('   ${i.i2}: norm( C*Q - C*Q )  / ( ${name.a1} * norm(C) * EPS )');
  }

  void print9942(int i, String name) {
    println(
        '   ${i.i2}: norm( Q\'*C - Q\'*C )/ ( ${name.a1} * norm(C) * EPS )');
  }

  void print9941(int i, String name) {
    println(
        '   ${i.i2}: norm( C*Q\' - C*Q\' )/ ( ${name.a1} * norm(C) * EPS )');
  }

  void print9940(int i) {
    println('   ${i.i2}: norm(svd(A) - svd(R)) / ( M * norm(svd(R)) * EPS )');
  }

  void print9939(int i) {
    println('   ${i.i2}: norm( A*P - Q*R ) / ( M * norm(A) * EPS )');
  }

  void print9938(int i) {
    println('   ${i.i2}: norm( I - Q\'*Q ) / ( M * EPS )');
  }

  void print9937(int i) {
    println('   ${i.i2}: norm( A - R*Q )       / ( M * norm(A) * EPS )');
  }

  void print9935(int i) {
    println(
        '   ${i.i2}: norm( B - A * X )   / ( max(M,N) * norm(A) * norm(X) * EPS )');
  }

  void print9934(int i) {
    println(
        '   ${i.i2}: norm( (A*X-B)\' *A ) / ( max(M,N,NRHS) * norm(A) * norm(B) * EPS )');
  }

  void print9933(int i) {
    println(
        '   ${i.i2}: norm(svd(A)-svd(R)) / ( min(M,N) * norm(svd(R)) * EPS )');
  }

  void print9932(int i) {
    println('   ${i.i2}: Check if X is in the row space of A or A\'');
  }

  void print9931(int i) {
    println(
        '   ${i.i2}: norm( (A*X-B)\' *A ) / ( max(M,N,NRHS) * norm(A) * norm(B) * EPS )\n${' ' * 7}if TRANS=\'N\' and M >= N or TRANS=\'T\' and M < N, otherwise\n${' ' * 7}check if X is in the row space of A or A\' (overdetermined case)');
  }

  void print9929(String prefix) {
    println(' Test ratios (1-3: ${prefix.a1}TZRZF):');
  }

  void print9919() {
    println('    3-4: same as 1-2    5-6: same as 1-2');
  }

  void print9920() {
    println('    11-14: same as 7-10    15-18: same as 7-10');
  }

  void print9921(String prefix) {
    println(
        ' Test ratios:\n    (1-2: ${prefix.a1}GELS, 3-4: ${prefix.a1}GELST, 5-6: ${prefix.a1}GETSLS, 7-10: ${prefix.a1}GELSY, 11-14: ${prefix.a1}GETSS, 15-18: ${prefix.a1}GELSD)');
  }

  void print9928() {
    println('${' ' * 7}where ALPHA = ( 1 + sqrt( 17 ) ) / 8');
  }

  void print9927(int i) {
    println(
        '   ${i.i2}: ( Largest element in L ).abs()\n${' ' * 12} - ( 1 / ( 1 - ALPHA ) ) + THRESH');
  }

  void print9926(int i) {
    println(
        '   ${i.i2}: Largest 2-Norm of 2-by-2 pivots\n${' ' * 12} - ( ( 1 + ALPHA ) / ( 1 - ALPHA ) ) + THRESH');
  }

  void print8011(int i) {
    println('   ${i.i2}: norm( R - Q\'*A ) / ( M * norm(A) * EPS )');
  }

  void print8012(int i) {
    println('   ${i.i2}: norm( I - Q\'*Q ) / ( M * EPS )');
  }

  void print8013(int i) {
    println('   ${i.i2}: norm( Q*C - Q*C ) / ( M * norm(C) * EPS )');
  }

  void print8014(int i) {
    println('   ${i.i2}: norm( Q\'*C - Q\'*C ) / ( M * norm(C) * EPS )');
  }

  void print8015(int i) {
    println('   ${i.i2}: norm( C*Q - C*Q ) / ( M * norm(C) * EPS )');
  }

  void print8016(int i) {
    println('   ${i.i2}: norm( C*Q\' - C*Q\' ) / ( M * norm(C) * EPS )');
  }

  void print8017(int i) {
    println('   ${i.i2}: norm( R - Q\'*A ) / ( (M+N) * norm(A) * EPS )');
  }

  void print8018(int i) {
    println('   ${i.i2}: norm( I - Q\'*Q ) / ( (M+N) * EPS )');
  }

  void print8019(int i) {
    println('   ${i.i2}: norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8020(int i) {
    println('   ${i.i2}: norm( Q\'*C - Q\'*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8021(int i) {
    println('   ${i.i2}: norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8022(int i) {
    println('   ${i.i2}: norm( C*Q\' - C*Q\' ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8023(int i) {
    println('   ${i.i2}: norm( L - A*Q\' ) / ( (M+N) * norm(A) * EPS )');
  }

  void print8024(int i) {
    println('   ${i.i2}: norm( I - Q*Q\' ) / ( (M+N) * EPS )');
  }

  void print8025(int i) {
    println('   ${i.i2}: norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8026(int i) {
    println('   ${i.i2}: norm( Q\'*C - Q\'*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8027(int i) {
    println('   ${i.i2}: norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8028(int i) {
    println('   ${i.i2}: norm( C*Q\' - C*Q\' ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8029(int i) {
    println('   ${i.i2}: norm( L - A*Q\' ) / ( (M+N) * norm(A) * EPS )');
  }

  void print8030(int i) {
    println('   ${i.i2}: norm( I - Q*Q\' ) / ( (M+N) * EPS )');
  }

  void print8031(int i) {
    println('   ${i.i2}: norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8032(int i) {
    println('   ${i.i2}: norm( Q\'*C - Q\'*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8033(int i) {
    println('   ${i.i2}: norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8034(int i) {
    println('   ${i.i2}: norm( C*Q\' - C*Q\' ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8035(int i) {
    println('   ${i.i2}: norm( R - Q\'*A ) / ( (M+N) * norm(A) * EPS )');
  }

  void print8036(int i) {
    println('   ${i.i2}: norm( I - Q\'*Q ) / ( (M+N) * EPS )');
  }

  void print8037(int i) {
    println('   ${i.i2}: norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8038(int i) {
    println('   ${i.i2}: norm( Q\'*C - Q\'*C ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8039(int i) {
    println('   ${i.i2}: norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8040(int i) {
    println('   ${i.i2}: norm( C*Q\' - C*Q\' ) / ( (M+N) * norm(C) * EPS )');
  }

  void print8050(int i) {
    println('   ${i.i2}: norm( R - Q\'*A ) / ( M * norm(A) * EPS )');
  }

  void print8051(int i) {
    println('   ${i.i2}: norm( I - Q\'*Q ) / ( M * EPS )');
  }

  void print8052(int i) {
    println('   ${i.i2}: norm( Q*C - Q*C ) / ( M * norm(C) * EPS )');
  }

  void print8053(int i) {
    println('   ${i.i2}: norm( Q\'*C - Q\'*C ) / ( M * norm(C) * EPS )');
  }

  void print8054(int i) {
    println('   ${i.i2}: norm( C*Q - C*Q ) / ( M * norm(C) * EPS )');
  }

  void print8055(int i) {
    println('   ${i.i2}: norm( C*Q\' - C*Q\' ) / ( M * norm(C) * EPS )');
  }

  void print8060(int i) {
    println(
        '   ${i.i2}: 2-norm(svd(A) - svd(R)) / ( max(M,N) * 2-norm(svd(R)) * EPS )');
  }

  void print8061(int i) {
    println('   ${i.i2}: 1-norm( A*P - Q*R ) / ( max(M,N) * 1-norm(A) * EPS )');
  }

  void print8062(int i) {
    println('   ${i.i2}: 1-norm( I - Q\'*Q ) / ( M * EPS )');
  }

  void print8063(int i) {
    println(
        '   ${i.i2}: Returns 1.0e+100, if abs(R(K+1,K+1)) > abs(R(K,K)), where K=1:KFACT-1');
  }

  void print8064(int i) {
    println('   ${i.i2}: 1-norm(Q**T * B - Q**T * B ) / ( M * EPS )');
  }
}
