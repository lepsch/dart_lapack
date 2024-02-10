double second() {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

  // .. Local Scalars ..

  double T;

// .. Intrinsic Functions ..

  // INTRINSIC CPU_TIME

// .. Executable Statements .. *

  cpu_time(T);
  SECOND = T;
}
