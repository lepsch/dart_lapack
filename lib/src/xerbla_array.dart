      void xerbla_array(SRNAME_ARRAY, SRNAME_LEN, final Box<int> INFO) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     SRNAME_LEN, INFO;
      String   (1) SRNAME_ARRAY(SRNAME_LEN);
      // ..

// =====================================================================

      int     I;
      String       SRNAME;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, LEN
      // ..
      // .. External Functions ..
      // EXTERNAL XERBLA
      SRNAME = ' ';
      DO I = 1, min( SRNAME_LEN, SRNAME.length );
         SRNAME[I:I] = SRNAME_ARRAY( I );
      }

      xerbla(SRNAME, INFO );
      }
