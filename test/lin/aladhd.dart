// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/nio.dart';

void aladhd(
  final Nout IOUNIT,
  final String PATH,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  String SYM;

  final C1 = PATH.substring(0, 1);
  final C3 = PATH.substring(2, 3);
  final P2 = PATH.substring(1, 3);
  final SORD = lsame(C1, 'S') || lsame(C1, 'D');
  final CORZ = lsame(C1, 'C') || lsame(C1, 'Z');
  if (!(SORD || CORZ)) return;

  if (lsamen(2, P2, 'GE')) {
    // GE: General dense

    IOUNIT.print9999(PATH);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9989();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9981(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9978(4);
    IOUNIT.print9977(5);
    IOUNIT.print9976(6);
    IOUNIT.print9972(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'GB')) {
    // GB: General band

    IOUNIT.print9998(PATH);
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9988();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9981(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9978(4);
    IOUNIT.print9977(5);
    IOUNIT.print9976(6);
    IOUNIT.print9972(7);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'GT')) {
    // GT: General tridiagonal

    IOUNIT.print9997(PATH);
    IOUNIT.print9987();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9981(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9978(4);
    IOUNIT.print9977(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PO') ||
      lsamen(2, P2, 'PP') ||
      lsamen(2, P2, 'PS')) {
    // PO: Positive definite full
    // PS: Positive definite full
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
    IOUNIT.print9985(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9975(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9978(4);
    IOUNIT.print9977(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PB')) {
    // PB: Positive definite band

    if (SORD) {
      IOUNIT.print9994(PATH, 'Symmetric');
    } else {
      IOUNIT.print9994(PATH, 'Hermitian');
    }
    IOUNIT.println(' Matrix types:');
    IOUNIT.print9984(PATH);
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9975(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9978(4);
    IOUNIT.print9977(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'PT')) {
    // PT: Positive definite tridiagonal

    if (SORD) {
      IOUNIT.print9993(PATH, 'Symmetric');
    } else {
      IOUNIT.print9993(PATH, 'Hermitian');
    }
    IOUNIT.print9986();
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9973(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9978(4);
    IOUNIT.print9977(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'SY') || lsamen(2, P2, 'SP')) {
    // SY: Symmetric indefinite full
    //     with partial (Bunch-Kaufman) pivoting algorithm
    // SP: Symmetric indefinite packed
    //     with partial (Bunch-Kaufman) pivoting algorithm

    if (lsame(C3, 'Y')) {
      IOUNIT.print9992(PATH, 'Symmetric');
    } else {
      IOUNIT.print9991(PATH, 'Symmetric');
    }
    IOUNIT.println(' Matrix types:');
    if (SORD) {
      IOUNIT.print9983();
    } else {
      IOUNIT.print9982();
    }
    IOUNIT.println(' Test ratios:');
    IOUNIT.print9974(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9977(4);
    IOUNIT.print9978(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'SR') || lsamen(2, P2, 'SK')) {
    // SR: Symmetric indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm

    // SK: Symmetric indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm,
    //     ( new storage format for factors:
    //       L and diagonal of D is stored in A,
    //       subdiagonal of D is stored in E )

    IOUNIT.print9992(PATH, 'Symmetric');

    IOUNIT.println(' Matrix types:');
    if (SORD) {
      IOUNIT.print9983();
    } else {
      IOUNIT.print9982();
    }

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9974(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HA')) {
    // HA: Hermitian
    //     Aasen algorithm
    IOUNIT.print9971(PATH, 'Hermitian');

    IOUNIT.println(' Matrix types:');
    IOUNIT.print9983();

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9974(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9977(4);
    IOUNIT.print9978(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HE') || lsamen(2, P2, 'HP')) {
    // HE: Hermitian indefinite full
    //     with partial (Bunch-Kaufman) pivoting algorithm
    // HP: Hermitian indefinite packed
    //     with partial (Bunch-Kaufman) pivoting algorithm

    if (lsame(C3, 'E')) {
      IOUNIT.print9992(PATH, 'Hermitian');
    } else {
      IOUNIT.print9991(PATH, 'Hermitian');
    }

    IOUNIT.println(' Matrix types:');
    IOUNIT.print9983();

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9974(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.print9977(4);
    IOUNIT.print9978(5);
    IOUNIT.print9976(6);
    IOUNIT.println(' Messages:');
  } else if (lsamen(2, P2, 'HR') || lsamen(2, P2, 'HK')) {
    // HR: Hermitian indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm

    // HK: Hermitian indefinite full,
    //     with rook (bounded Bunch-Kaufman) pivoting algorithm,
    //     ( new storage format for factors:
    //       L and diagonal of D is stored in A,
    //       subdiagonal of D is stored in E )

    IOUNIT.print9992(PATH, 'Hermitian');

    IOUNIT.println(' Matrix types:');
    IOUNIT.print9983();

    IOUNIT.println(' Test ratios:');
    IOUNIT.print9974(1);
    IOUNIT.print9980(2);
    IOUNIT.print9979(3);
    IOUNIT.println(' Messages:');
  } else {
    // Print error message if no header is available.

    IOUNIT.print9990(PATH);
  }
}

extension on Nout {
  // First line of header

  void print9999(String PATH) {
    println('\n ${PATH.a3} drivers:  General dense matrices');
  }

  void print9998(String PATH) {
    println('\n ${PATH.a3} drivers:  General band matrices');
  }

  void print9997(String PATH) {
    println('\n ${PATH.a3} drivers:  General tridiagonal');
  }

  void print9996(String PATH, String SYM) {
    println('\n ${PATH.a3} drivers:  ${SYM.a9} positive definite matrices');
  }

  void print9995(String PATH, String SYM) {
    println(
        '\n ${PATH.a3} drivers:  ${SYM.a9} positive definite packed matrices');
  }

  void print9994(String PATH, String SYM) {
    println(
        '\n ${PATH.a3} drivers:  ${SYM.a9} positive definite band matrices');
  }

  void print9993(String PATH, String SYM) {
    println('\n ${PATH.a3} drivers:  ${SYM.a9} positive definite tridiagonal');
  }

  void print9971(String PATH, String SYM) {
    println(
        '\n ${PATH.a3} drivers:  ${SYM.a9} indefinite matrices, "Aasen" Algorithm');
  }

  void print9992(String PATH, String SYM) {
    println(
        '\n ${PATH.a3} drivers:  ${SYM.a9} indefinite matrices, "rook" (bounded Bunch-Kaufman) pivoting');
  }

  void print9991(String PATH, String SYM) {
    println(
        '\n ${PATH.a3} drivers:  ${SYM.a9} indefinite packed matrices, partial (Bunch-Kaufman) pivoting');
  }

//  void print9891(String PATH, String SYM){println('\n ${PATH.a3} drivers:  ${SYM.a9} indefinite packed matrices, "rook" (bounded Bunch-Kaufman) pivoting' );}
  void print9990(String PATH) {
    println('\n ${PATH.a3}:  No header available');
  }

  // GE matrix types

  void print9989() {
    println(
        '    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero');
  }

  // GB matrix types

  void print9988() {
    println(
        '    1. Random, CNDNUM = 2${' ' * 14}5. Random, CNDNUM = sqrt(0.1/EPS)\n    2. First column zero${' ' * 15}6. Random, CNDNUM = 0.1/EPS\n    3. Last column zero${' ' * 16}7. Scaled near underflow\n    4. Last n/2 columns zero${' ' * 11}8. Scaled near overflow');
  }

  // GT matrix types

  void print9987() {
    println(
        ' Matrix types (1-6 have specified condition numbers):\n    1. Diagonal${' ' * 24}7. Random, unspecified CNDNUM\n    2. Random, CNDNUM = 2${' ' * 14}8. First column zero\n    3. Random, CNDNUM = sqrt(0.1/EPS)  9. Last column zero\n    4. Random, CNDNUM = 0.1/EPS${' ' * 7}10. Last n/2 columns zero\n    5. Scaled near underflow${' ' * 10}11. Scaled near underflow\n    6. Scaled near overflow${' ' * 11}12. Scaled near overflow');
  }

  // PT matrix types

  void print9986() {
    println(
        ' Matrix types (1-6 have specified condition numbers):\n    1. Diagonal${' ' * 24}7. Random, unspecified CNDNUM\n    2. Random, CNDNUM = 2${' ' * 14}8. First row and column zero\n    3. Random, CNDNUM = sqrt(0.1/EPS)  9. Last row and column zero\n    4. Random, CNDNUM = 0.1/EPS${' ' * 7}10. Middle row and column zero\n    5. Scaled near underflow${' ' * 10}11. Scaled near underflow\n    6. Scaled near overflow${' ' * 11}12. Scaled near overflow');
  }

  // PO, PP matrix types

  void print9985(String PATH) {
    println(
        '    1. Diagonal${' ' * 24}6. Random, CNDNUM = sqrt(0.1/EPS)\n    2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = 0.1/EPS\n   *3. First row and column zero${' ' * 7}8. Scaled near underflow\n   *4. Last row and column zero${' ' * 8}9. Scaled near overflow\n   *5. Middle row and column zero\n   (* - tests error exits from ${PATH.a3}TRF, no test ratios are computed)');
  }

  // PB matrix types

  void print9984(String PATH) {
    println(
        '    1. Random, CNDNUM = 2${' ' * 14}5. Random, CNDNUM = sqrt(0.1/EPS)\n   *2. First row and column zero${' ' * 7}6. Random, CNDNUM = 0.1/EPS\n   *3. Last row and column zero${' ' * 8}7. Scaled near underflow\n   *4. Middle row and column zero${' ' * 6}8. Scaled near overflow\n   (* - tests error exits from ${PATH.a3}TRF, no test ratios are computed)');
  }

  // SSY, SSP, CHE, CHP matrix types

  void print9983() {
    println(
        '    1. Diagonal${' ' * 24}6. Last n/2 rows and columns zero\n    2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = sqrt(0.1/EPS)\n    3. First row and column zero${' ' * 7}8. Random, CNDNUM = 0.1/EPS\n    4. Last row and column zero${' ' * 8}9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow');
  }

  // CSY, CSP matrix types

  void print9982() {
    println(
        '    1. Diagonal${' ' * 24}7. Random, CNDNUM = sqrt(0.1/EPS)\n    2. Random, CNDNUM = 2${' ' * 14}8. Random, CNDNUM = 0.1/EPS\n    3. First row and column zero${' ' * 7}9. Scaled near underflow\n    4. Last row and column zero${' ' * 7}10. Scaled near overflow\n    5. Middle row and column zero     11. Block diagonal matrix\n    6. Last n/2 rows and columns zero');
  }

  // Test ratios

  void print9981(int i) {
    println('   ${i.i2}: norm( L * U - A )  / ( N * norm(A) * EPS )');
  }

  void print9980(int i) {
    println('   ${i.i2}: norm( B - A * X )  / ( norm(A) * norm(X) * EPS )');
  }

  void print9979(int i) {
    println('   ${i.i2}: norm( X - XACT )   / ( norm(XACT) * CNDNUM * EPS )');
  }

  void print9978(int i) {
    println('   ${i.i2}: norm( X - XACT )   / ( norm(XACT) * (error bound) )');
  }

  void print9977(int i) {
    println('   ${i.i2}: (backward error)   / EPS');
  }

  void print9976(int i) {
    println('   ${i.i2}: RCOND * CNDNUM - 1.0');
  }

  void print9975(int i) {
    println(
        '   ${i.i2}: norm( U\' * U - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L * L\' - A ) / ( N * norm(A) * EPS )');
  }

  void print9974(int i) {
    println(
        '   ${i.i2}: norm( U*D*U\' - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L*D*L\' - A ) / ( N * norm(A) * EPS )');
  }

  void print9973(int i) {
    println(
        '   ${i.i2}: norm( U\'*D*U - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L*D*L\' - A ) / ( N * norm(A) * EPS )');
  }

  void print9972(int i) {
    println(
        '   ${i.i2}: abs( WORK(1) - RPVGRW ) / ( max( WORK(1), RPVGRW ) * EPS )');
  }
}
