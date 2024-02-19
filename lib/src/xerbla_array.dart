import 'dart:math';

import 'package:lapack/src/xerbla.dart';

void xerbla_array(
  final String SRNAME_ARRAY,
  final int SRNAME_LEN,
  final int INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I;
  String SRNAME;

  SRNAME = ' ';
  for (I = 1; I <= min(SRNAME_LEN, SRNAME.length); I++) {
    SRNAME += SRNAME_ARRAY[I - 1];
  }

  xerbla(SRNAME, INFO);
}
