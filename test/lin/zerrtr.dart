      void zerrtr(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                infoc.NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                INFO;
      double             RCOND, SCALE, SCALES(0);
      // ..
      // .. Local Arrays ..
      double             R1( NMAX ), R2( NMAX ), RW( NMAX );
      Complex         A( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZLATBS, ZLATPS, ZLATRS, ZLATRS3, ZTBCON, ZTBRFS, ZTBTRS, ZTPCON, ZTPRFS, ZTPTRI, ZTPTRS, ZTRCON, ZTRRFS, ZTRTI2, ZTRTRI, ZTRTRS
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A[1, 1] = 1.0;
      A[1, 2] = 2.0;
      A[2, 2] = 3.0;
      A[2, 1] = 4.0;
      infoc.OK = true;

      // Test error exits for the general triangular routines.

      if ( lsamen( 2, C2, 'TR' ) ) {

         // ZTRTRI

        srnamc.SRNAMT = 'ZTRTRI';
         infoc.INFOT = 1;
         ztrtri('/', 'N', 0, A, 1, INFO );
         chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztrtri('U', '/', 0, A, 1, INFO );
         chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztrtri('U', 'N', -1, A, 1, INFO );
         chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztrtri('U', 'N', 2, A, 1, INFO );
         chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTRTI2

        srnamc.SRNAMT = 'ZTRTI2';
         infoc.INFOT = 1;
         ztrti2('/', 'N', 0, A, 1, INFO );
         chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztrti2('U', '/', 0, A, 1, INFO );
         chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztrti2('U', 'N', -1, A, 1, INFO );
         chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztrti2('U', 'N', 2, A, 1, INFO );
         chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );


         // ZTRTRS

        srnamc.SRNAMT = 'ZTRTRS';
         infoc.INFOT = 1;
         ztrtrs('/', 'N', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztrtrs('U', '/', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztrtrs('U', 'N', '/', 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztrtrs('U', 'N', 'N', -1, 0, A, 1, X, 1, INFO );
         chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztrtrs('U', 'N', 'N', 0, -1, A, 1, X, 1, INFO );
         chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;

         // ZTRRFS

        srnamc.SRNAMT = 'ZTRRFS';
         infoc.INFOT = 1;
         ztrrfs('/', 'N', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztrrfs('U', '/', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztrrfs('U', 'N', '/', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztrrfs('U', 'N', 'N', -1, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztrrfs('U', 'N', 'N', 0, -1, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         ztrrfs('U', 'N', 'N', 2, 1, A, 1, B, 2, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         ztrrfs('U', 'N', 'N', 2, 1, A, 2, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         ztrrfs('U', 'N', 'N', 2, 1, A, 2, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTRCON

        srnamc.SRNAMT = 'ZTRCON';
         infoc.INFOT = 1;
         ztrcon('/', 'U', 'N', 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztrcon('1', '/', 'N', 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztrcon('1', 'U', '/', 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztrcon('1', 'U', 'N', -1, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         ztrcon('1', 'U', 'N', 2, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZLATRS

        srnamc.SRNAMT = 'ZLATRS';
         infoc.INFOT = 1;
         zlatrs('/', 'N', 'N', 'N', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zlatrs('U', '/', 'N', 'N', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zlatrs('U', 'N', '/', 'N', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zlatrs('U', 'N', 'N', '/', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zlatrs('U', 'N', 'N', 'N', -1, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zlatrs('U', 'N', 'N', 'N', 2, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZLATRS3

        srnamc.SRNAMT = 'ZLATRS3';
         infoc.INFOT = 1;
         zlatrs3('/', 'N', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zlatrs3('U', '/', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zlatrs3('U', 'N', '/', 'N', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zlatrs3('U', 'N', 'N', '/', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zlatrs3('U', 'N', 'N', 'N', -1, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zlatrs3('U', 'N', 'N', 'N', 0, -1, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zlatrs3('U', 'N', 'N', 'N', 2, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zlatrs3('U', 'N', 'N', 'N', 2, 0, A, 2, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         zlatrs3('U', 'N', 'N', 'N', 1, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 0, INFO );
         chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Test error exits for the packed triangular routines.

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // ZTPTRI

        srnamc.SRNAMT = 'ZTPTRI';
         infoc.INFOT = 1;
         ztptri('/', 'N', 0, A, INFO );
         chkxer('ZTPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztptri('U', '/', 0, A, INFO );
         chkxer('ZTPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztptri('U', 'N', -1, A, INFO );
         chkxer('ZTPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTPTRS

        srnamc.SRNAMT = 'ZTPTRS';
         infoc.INFOT = 1;
         ztptrs('/', 'N', 'N', 0, 0, A, X, 1, INFO );
         chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztptrs('U', '/', 'N', 0, 0, A, X, 1, INFO );
         chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztptrs('U', 'N', '/', 0, 0, A, X, 1, INFO );
         chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztptrs('U', 'N', 'N', -1, 0, A, X, 1, INFO );
         chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztptrs('U', 'N', 'N', 0, -1, A, X, 1, INFO );
         chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         ztptrs('U', 'N', 'N', 2, 1, A, X, 1, INFO );
         chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTPRFS

        srnamc.SRNAMT = 'ZTPRFS';
         infoc.INFOT = 1;
         ztprfs('/', 'N', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztprfs('U', '/', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztprfs('U', 'N', '/', 0, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztprfs('U', 'N', 'N', -1, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztprfs('U', 'N', 'N', 0, -1, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         ztprfs('U', 'N', 'N', 2, 1, A, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         ztprfs('U', 'N', 'N', 2, 1, A, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTPCON

        srnamc.SRNAMT = 'ZTPCON';
         infoc.INFOT = 1;
         ztpcon('/', 'U', 'N', 0, A, RCOND, W, RW, INFO );
         chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztpcon('1', '/', 'N', 0, A, RCOND, W, RW, INFO );
         chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztpcon('1', 'U', '/', 0, A, RCOND, W, RW, INFO );
         chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztpcon('1', 'U', 'N', -1, A, RCOND, W, RW, INFO );
         chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZLATPS

        srnamc.SRNAMT = 'ZLATPS';
         infoc.INFOT = 1;
         zlatps('/', 'N', 'N', 'N', 0, A, X, SCALE, RW, INFO );
         chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zlatps('U', '/', 'N', 'N', 0, A, X, SCALE, RW, INFO );
         chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zlatps('U', 'N', '/', 'N', 0, A, X, SCALE, RW, INFO );
         chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zlatps('U', 'N', 'N', '/', 0, A, X, SCALE, RW, INFO );
         chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zlatps('U', 'N', 'N', 'N', -1, A, X, SCALE, RW, INFO );
         chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Test error exits for the banded triangular routines.

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // ZTBTRS

        srnamc.SRNAMT = 'ZTBTRS';
         infoc.INFOT = 1;
         ztbtrs('/', 'N', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztbtrs('U', '/', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztbtrs('U', 'N', '/', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztbtrs('U', 'N', 'N', -1, 0, 0, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztbtrs('U', 'N', 'N', 0, -1, 0, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         ztbtrs('U', 'N', 'N', 0, 0, -1, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         ztbtrs('U', 'N', 'N', 2, 1, 1, A, 1, X, 2, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         ztbtrs('U', 'N', 'N', 2, 0, 1, A, 1, X, 1, INFO );
         chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTBRFS

        srnamc.SRNAMT = 'ZTBRFS';
         infoc.INFOT = 1;
         ztbrfs('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztbrfs('U', '/', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztbrfs('U', 'N', '/', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztbrfs('U', 'N', 'N', -1, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztbrfs('U', 'N', 'N', 0, -1, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         ztbrfs('U', 'N', 'N', 0, 0, -1, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         ztbrfs('U', 'N', 'N', 2, 1, 1, A, 1, B, 2, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         ztbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         ztbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZTBCON

        srnamc.SRNAMT = 'ZTBCON';
         infoc.INFOT = 1;
         ztbcon('/', 'U', 'N', 0, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztbcon('1', '/', 'N', 0, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         ztbcon('1', 'U', '/', 0, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztbcon('1', 'U', 'N', -1, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         ztbcon('1', 'U', 'N', 0, -1, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         ztbcon('1', 'U', 'N', 2, 1, A, 1, RCOND, W, RW, INFO );
         chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZLATBS

        srnamc.SRNAMT = 'ZLATBS';
         infoc.INFOT = 1;
         zlatbs('/', 'N', 'N', 'N', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zlatbs('U', '/', 'N', 'N', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zlatbs('U', 'N', '/', 'N', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zlatbs('U', 'N', 'N', '/', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zlatbs('U', 'N', 'N', 'N', -1, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zlatbs('U', 'N', 'N', 'N', 1, -1, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zlatbs('U', 'N', 'N', 'N', 2, 1, A, 1, X, SCALE, RW, INFO );
         chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
