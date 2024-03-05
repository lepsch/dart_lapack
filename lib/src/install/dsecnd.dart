final _timer = Stopwatch()..start();

double dsecnd() {
  // -- LAPACK auxiliary routine --
  // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  // =====================================================================

  return _timer.elapsed.inMilliseconds / 1000;
}
