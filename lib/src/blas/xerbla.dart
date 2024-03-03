void Function(String SRNAME, int INFO) xerbla = _xerbla;

void _xerbla(final String SRNAME, final int INFO) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  throw ArgumentError(
      ' ** On entry to ${SRNAME.trim()} parameter number $INFO had an illegal value');
}
