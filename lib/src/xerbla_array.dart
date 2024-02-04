      void xerbla_array(SRNAME_ARRAY, SRNAME_LEN, INFO) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     SRNAME_LEN, INFO;
      // ..
      // .. Array Arguments ..
      String   (1) SRNAME_ARRAY(SRNAME_LEN);
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      int     I;
      // ..
      // .. Local Arrays ..
      String       SRNAME;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, LEN
      // ..
      // .. External Functions ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..
      SRNAME = ' ';
      DO I = 1, min( SRNAME_LEN, SRNAME.length );
         SRNAME[I:I] = SRNAME_ARRAY( I );
      }

      xerbla(SRNAME, INFO );
      }
