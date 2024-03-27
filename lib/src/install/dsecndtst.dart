import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dsecnd.dart';
import 'package:lapack/src/matrix.dart';

void main() {
  const NMAX = 1000, ITS = 50000;
  int I, J;
  double ALPHA, AVG, T1, T2, TNOSEC, TOTAL;
  final X = Array<double>(NMAX), Y = Array<double>(NMAX);

  // .. Figure TOTAL flops ..
  TOTAL = NMAX * ITS * 2.0;

  // Initialize X and Y

  for (I = 1; I <= NMAX; I++) {
    X[I] = 1 / I;
    Y[I] = (NMAX - I) / NMAX;
  }
  ALPHA = 0.315;

  // Time TOTAL SAXPY operations

  T1 = dsecnd();
  for (J = 1; J <= ITS; J++) {
    for (I = 1; I <= NMAX; I++) {
      Y[I] += ALPHA * X[I];
    }
    ALPHA = -ALPHA;
  }
  T2 = dsecnd();
  TNOSEC = T2 - T1;
  print(' Time for ${TOTAL.g10_3} DAXPY ops = ${TNOSEC.g10_3} seconds');
  if (TNOSEC > 0.0) {
    print(
        ' DAXPY performance rate        = ${((TOTAL / 1.0e6) / TNOSEC).g10_3} mflops ');
  } else {
    print(
        ' *** Warning:  Time for operations was less or equal than zero => timing in TESTING might be dubious');
  }

  // Time TOTAL DAXPY operations with dsecnd in the outer loop

  T1 = dsecnd();
  for (J = 1; J <= ITS; J++) {
    for (I = 1; I <= NMAX; I++) {
      Y[I] += ALPHA * X[I];
    }
    ALPHA = -ALPHA;
    T2 = dsecnd();
  }

  // Compute the time used in milliseconds used by an average call
  // to dsecnd.

  print(' Including dsecnd, time        = ${(T2 - T1).g10_3} seconds');
  AVG = ((T2 - T1) - TNOSEC) * 1000.0e+00 / ITS;
  if (AVG > 0.0) {
    print(' Average time for dsecnd       = ${AVG.g10_3} milliseconds');
  }

  // Compute the equivalent number of floating point operations used
  // by an average call to dsecnd.

  if ((AVG > 0.0) && (TNOSEC > 0.0)) {
    print(
        ' Equivalent floating point ops = ${((AVG / 1000) * TOTAL / TNOSEC).g10_3} ops');
  }
}
