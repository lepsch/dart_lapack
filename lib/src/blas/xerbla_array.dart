import 'dart:math';

import 'package:lapack/src/blas/xerbla.dart';

void xerbla_array(
  final List<int> SRNAME_ARRAY,
  final int SRNAME_LEN,
  final int INFO,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I;
  String SRNAME;

  SRNAME = ' ';
  for (I = 1; I <= min(SRNAME_LEN, SRNAME.length); I++) {
    SRNAME += String.fromCharCode(SRNAME_ARRAY[I]);
  }

  xerbla(SRNAME, INFO);
}
