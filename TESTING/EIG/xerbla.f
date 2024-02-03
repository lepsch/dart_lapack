      SUBROUTINE XERBLA( SRNAME, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>         SRNAME;
      int                INFO;
      // ..

*  =====================================================================

      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      LERR = true;
      if ( INFO.NE.INFOT ) {
         if ( INFOT.NE.0 ) {
            WRITE( NOUT, FMT = 9999 ) SRNAMT( 1:LEN_TRIM( SRNAMT ) ), INFO, INFOT
         } else {
            WRITE( NOUT, FMT = 9997 ) SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
         }
         OK = false;
      }
      if ( SRNAME.NE.SRNAMT ) {
         WRITE( NOUT, FMT = 9998 ) SRNAME( 1:LEN_TRIM( SRNAME ) ), SRNAMT( 1:LEN_TRIM( SRNAMT ) )
         OK = false;
      }
      RETURN

 9999 FORMAT( ' *** XERBLA was called from ', A, ' with INFO = ', I6, ' instead of ', I2, ' ***' )
 9998 FORMAT( ' *** XERBLA was called with SRNAME = ', A, ' instead of ', A6, ' ***' )
 9997 FORMAT( ' *** On entry to ', A, ' parameter number ', I6, ' had an illegal value ***' )

      // End of XERBLA

      }
