void Function(String SRNAME, int INFO) xerbla = _xerbla;

void _xerbla(final String SRNAME, final int INFO) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  throw ArgumentError(
      ' ** On entry to $SRNAME parameter number $INFO had an illegal value');
}
