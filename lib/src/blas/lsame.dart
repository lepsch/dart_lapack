bool lsame(final String CA, final String CB) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  return CA.isNotEmpty &&
      CB.isNotEmpty &&
      CA.toUpperCase()[0] == CB.toUpperCase()[0];
}
