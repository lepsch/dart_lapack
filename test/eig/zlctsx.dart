bool zlctsx(ALPHA, BETA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex ALPHA, BETA;
  // ..

// =====================================================================

  // .. Parameters ..
  // double                         ZERO;
  // PARAMETER          ( ZERO = 0.0 )
  // Complex            CZERO
  // PARAMETER          ( CZERO = ( 0.0, 0.0 ) )
  // ..
  // .. Scalars in Common ..
  bool FS;
  int I, M, MPLUSN, N;
  // ..
  // .. Common blocks ..
  // COMMON / MN / M, N, MPLUSN, I, FS
  // ..
  // .. Save statement ..
  SAVE;
  // ..
  // .. Executable Statements ..

  if (FS) {
    I = I + 1;
    if (I <= M) {
      ZLCTSX = false;
    } else {
      ZLCTSX = true;
    }
    if (I == MPLUSN) {
      FS = false;
      I = 0;
    }
  } else {
    I = I + 1;
    if (I <= N) {
      ZLCTSX = true;
    } else {
      ZLCTSX = false;
    }
    if (I == MPLUSN) {
      FS = true;
      I = 0;
    }
  }

  // IF( BETA == CZERO ) THEN
  // ZLCTSX = ( ALPHA.toDouble() > ZERO )
  // ELSE
  // ZLCTSX = ( (ALPHA/BETA).toDouble() > ZERO )
  // END IF

  return;
}