import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'dlahd2.dart';

void dlafts(
  final String TYPE,
  final int M,
  final int N,
  final int IMAT,
  final int NTESTS,
  final Array<double> RESULT,
  final Array<int> ISEED,
  final double THRESH,
  final Nout IOUNIT,
  final Box<int> IE,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int K;

  if (M == N) {
    // Output for square matrices:

    for (K = 1; K <= NTESTS; K++) {
      if (RESULT[K] >= THRESH) {
        // If this is the first test to fail, call DLAHD2
        // to print a header to the data file.

        if (IE.value == 0) dlahd2(IOUNIT, TYPE);
        IE.value += 1;
        if (RESULT[K] < 10000.0) {
          IOUNIT.println(
              ' Matrix order=${N.i5}, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')} result ${K.i3} is${RESULT[K].f8_2}');
        } else {
          IOUNIT.println(
              ' Matrix order=${N.i5}, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')} result ${K.i3} is${(RESULT[K] * 10).d10_3}');
        }
      }
    }
  } else {
    // Output for rectangular matrices

    for (K = 1; K <= NTESTS; K++) {
      if (RESULT[K] >= THRESH) {
        // If this is the first test to fail, call DLAHD2
        // to print a header to the data file.

        if (IE.value == 0) dlahd2(IOUNIT, TYPE);
        IE.value = IE.value + 1;
        if (RESULT[K] < 10000.0) {
          IOUNIT.println(
              ' ${M.i5} x${N.i5} matrix, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')}: result ${K.i3} is${RESULT[K].f8_2}');
        } else {
          IOUNIT.println(
              ' ${M.i5} x${N.i5} matrix, type=${IMAT.i2}, seed=${ISEED.i4(4, ',')}: result ${K.i3} is${(RESULT[K] * 10).d10_3}');
        }
      }
    }
  }
}
