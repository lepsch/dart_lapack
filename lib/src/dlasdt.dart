// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlasdt(
  final int N,
  final Box<int> LVL,
  final Box<int> ND,
  final Array<int> INODE_,
  final Array<int> NDIML_,
  final Array<int> NDIMR_,
  final int MSUB,
) {
  final INODE = INODE_.having();
  final NDIML = NDIML_.having();
  final NDIMR = NDIMR_.having();
  const TWO = 2.0;
  int I, IL, IR, LLST, MAXN, NCRNT, NLVL;
  double TEMP;

  // Find the number of levels on the tree.

  MAXN = max(1, N);
  TEMP = log(MAXN / (MSUB + 1)) / log(TWO);
  LVL.value = TEMP.toInt() + 1;

  I = N ~/ 2;
  INODE[1] = I + 1;
  NDIML[1] = I;
  NDIMR[1] = N - I - 1;
  IL = 0;
  IR = 1;
  LLST = 1;
  for (NLVL = 1; NLVL <= LVL.value - 1; NLVL++) {
    // Constructing the tree at (NLVL+1)-st level. The number of
    // nodes created on this level is LLST * 2.

    for (I = 0; I <= LLST - 1; I++) {
      IL += 2;
      IR += 2;
      NCRNT = LLST + I;
      NDIML[IL] = NDIML[NCRNT] ~/ 2;
      NDIMR[IL] = NDIML[NCRNT] - NDIML[IL] - 1;
      INODE[IL] = INODE[NCRNT] - NDIMR[IL] - 1;
      NDIML[IR] = NDIMR[NCRNT] ~/ 2;
      NDIMR[IR] = NDIMR[NCRNT] - NDIML[IR] - 1;
      INODE[IR] = INODE[NCRNT] + NDIML[IR] + 1;
    }
    LLST *= 2;
  }
  ND.value = LLST * 2 - 1;
}
