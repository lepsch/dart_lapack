      void xerbla(SRNAME, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>         SRNAME;
      int                INFO;
      // ..

// =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. Executable Statements ..

      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO;

      STOP;

 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', 'an illegal value' );

      // End of XERBLA

      }
