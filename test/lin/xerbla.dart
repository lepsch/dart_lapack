import 'common.dart';
      void xerbla(SRNAME, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      List<String>         SRNAME;
      int                INFO;
      // ..

// =====================================================================

      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      infoc.LERR = true;
      if ( INFO != infoc.INFOT ) {
         if ( infoc.INFOT != 0 ) {
            WRITE( infoc.NOUT, FMT = 9999 )srnamc.SRNAMT( 1:LEN_TRIM(srnamc.SRNAMT ) ), INFO, infoc.INFOT;
         } else {
            WRITE( infoc.NOUT, FMT = 9997 ) SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO;
         }
         infoc.OK = false;
      }
      if ( SRNAME !=srnamc.SRNAMT ) {
         WRITE( infoc.NOUT, FMT = 9998 ) SRNAME( 1:LEN_TRIM( SRNAME ) ),srnamc.SRNAMT( 1:LEN_TRIM(srnamc.SRNAMT ) );
         infoc.OK = false;
      }
      return;

 9999 FORMAT( ' *** XERBLA was called from ', A, ' with INFO = ', I6, ' instead of ', I2, ' ***' );
 9998 FORMAT( ' *** XERBLA was called with SRNAME = ', A, ' instead of ', A9, ' ***' );
 9997 FORMAT( ' *** On entry to ', A, ' parameter number ', I6, ' had an illegal value ***' );
      }
