import 'common.dart';
      void derrtr(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
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
      int                IW( NMAX );
      double             A( NMAX, NMAX ), B( NMAX ), R1( NMAX ), R2( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DLATBS, DLATPS, DLATRS, DLATRS3, DTBCON, DTBRFS, DTBTRS, DTPCON, DTPRFS, DTPTRI, DTPTRS, DTRCON, DTRRFS, DTRTI2, DTRTRI, DTRTRS
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String            srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A[1, 1] = 1.0;
      A[1, 2] = 2.0;
      A[2, 2] = 3.0;
      A[2, 1] = 4.0;
      infoc.OK = true;

      if ( lsamen( 2, C2, 'TR' ) ) {

         // Test error exits for the general triangular routines.

         // DTRTRI

        srnamc.SRNAMT = 'DTRTRI';
         infoc.INFOT = 1;
         dtrtri('/', 'N', 0, A, 1, INFO );
         chkxer('DTRTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrtri('U', '/', 0, A, 1, INFO );
         chkxer('DTRTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtrtri('U', 'N', -1, A, 1, INFO );
         chkxer('DTRTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtrtri('U', 'N', 2, A, 1, INFO );
         chkxer('DTRTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTRTI2

        srnamc.SRNAMT = 'DTRTI2';
         infoc.INFOT = 1;
         dtrti2('/', 'N', 0, A, 1, INFO );
         chkxer('DTRTI2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrti2('U', '/', 0, A, 1, INFO );
         chkxer('DTRTI2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtrti2('U', 'N', -1, A, 1, INFO );
         chkxer('DTRTI2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtrti2('U', 'N', 2, A, 1, INFO );
         chkxer('DTRTI2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTRTRS

        srnamc.SRNAMT = 'DTRTRS';
         infoc.INFOT = 1;
         dtrtrs('/', 'N', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrtrs('U', '/', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtrtrs('U', 'N', '/', 0, 0, A, 1, X, 1, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtrtrs('U', 'N', 'N', -1, 0, A, 1, X, 1, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtrtrs('U', 'N', 'N', 0, -1, A, 1, X, 1, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dtrtrs('U', 'N', 'N', 2, 1, A, 1, X, 2, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dtrtrs('U', 'N', 'N', 2, 1, A, 2, X, 1, INFO );
         chkxer('DTRTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTRRFS

        srnamc.SRNAMT = 'DTRRFS';
         infoc.INFOT = 1;
         dtrrfs('/', 'N', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrrfs('U', '/', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtrrfs('U', 'N', '/', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtrrfs('U', 'N', 'N', -1, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtrrfs('U', 'N', 'N', 0, -1, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dtrrfs('U', 'N', 'N', 2, 1, A, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dtrrfs('U', 'N', 'N', 2, 1, A, 2, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dtrrfs('U', 'N', 'N', 2, 1, A, 2, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTRRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTRCON

        srnamc.SRNAMT = 'DTRCON';
         infoc.INFOT = 1;
         dtrcon('/', 'U', 'N', 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTRCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrcon('1', '/', 'N', 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTRCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtrcon('1', 'U', '/', 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTRCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtrcon('1', 'U', 'N', -1, A, 1, RCOND, W, IW, INFO );
         chkxer('DTRCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dtrcon('1', 'U', 'N', 2, A, 1, RCOND, W, IW, INFO );
         chkxer('DTRCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DLATRS

        srnamc.SRNAMT = 'DLATRS';
         infoc.INFOT = 1;
         dlatrs('/', 'N', 'N', 'N', 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dlatrs('U', '/', 'N', 'N', 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dlatrs('U', 'N', '/', 'N', 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dlatrs('U', 'N', 'N', '/', 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dlatrs('U', 'N', 'N', 'N', -1, A, 1, X, SCALE, W, INFO );
         chkxer('DLATRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dlatrs('U', 'N', 'N', 'N', 2, A, 1, X, SCALE, W, INFO );
         chkxer('DLATRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DLATRS3

        srnamc.SRNAMT = 'DLATRS3';
         infoc.INFOT = 1;
         dlatrs3('/', 'N', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dlatrs3('U', '/', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dlatrs3('U', 'N', '/', 'N', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dlatrs3('U', 'N', 'N', '/', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dlatrs3('U', 'N', 'N', 'N', -1, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dlatrs3('U', 'N', 'N', 'N', 0, -1, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dlatrs3('U', 'N', 'N', 'N', 2, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dlatrs3('U', 'N', 'N', 'N', 2, 0, A, 2, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         dlatrs3('U', 'N', 'N', 'N', 1, 0, A, 1, X, 1, SCALES, W, W( 2 ), 0, INFO );
         chkxer('DLATRS3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // Test error exits for the packed triangular routines.

         // DTPTRI

        srnamc.SRNAMT = 'DTPTRI';
         infoc.INFOT = 1;
         dtptri('/', 'N', 0, A, INFO );
         chkxer('DTPTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtptri('U', '/', 0, A, INFO );
         chkxer('DTPTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtptri('U', 'N', -1, A, INFO );
         chkxer('DTPTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTPTRS

        srnamc.SRNAMT = 'DTPTRS';
         infoc.INFOT = 1;
         dtptrs('/', 'N', 'N', 0, 0, A, X, 1, INFO );
         chkxer('DTPTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtptrs('U', '/', 'N', 0, 0, A, X, 1, INFO );
         chkxer('DTPTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtptrs('U', 'N', '/', 0, 0, A, X, 1, INFO );
         chkxer('DTPTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtptrs('U', 'N', 'N', -1, 0, A, X, 1, INFO );
         chkxer('DTPTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtptrs('U', 'N', 'N', 0, -1, A, X, 1, INFO );
         chkxer('DTPTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dtptrs('U', 'N', 'N', 2, 1, A, X, 1, INFO );
         chkxer('DTPTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTPRFS

        srnamc.SRNAMT = 'DTPRFS';
         infoc.INFOT = 1;
         dtprfs('/', 'N', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtprfs('U', '/', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtprfs('U', 'N', '/', 0, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtprfs('U', 'N', 'N', -1, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtprfs('U', 'N', 'N', 0, -1, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dtprfs('U', 'N', 'N', 2, 1, A, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dtprfs('U', 'N', 'N', 2, 1, A, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTPRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTPCON

        srnamc.SRNAMT = 'DTPCON';
         infoc.INFOT = 1;
         dtpcon('/', 'U', 'N', 0, A, RCOND, W, IW, INFO );
         chkxer('DTPCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtpcon('1', '/', 'N', 0, A, RCOND, W, IW, INFO );
         chkxer('DTPCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtpcon('1', 'U', '/', 0, A, RCOND, W, IW, INFO );
         chkxer('DTPCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtpcon('1', 'U', 'N', -1, A, RCOND, W, IW, INFO );
         chkxer('DTPCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DLATPS

        srnamc.SRNAMT = 'DLATPS';
         infoc.INFOT = 1;
         dlatps('/', 'N', 'N', 'N', 0, A, X, SCALE, W, INFO );
         chkxer('DLATPS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dlatps('U', '/', 'N', 'N', 0, A, X, SCALE, W, INFO );
         chkxer('DLATPS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dlatps('U', 'N', '/', 'N', 0, A, X, SCALE, W, INFO );
         chkxer('DLATPS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dlatps('U', 'N', 'N', '/', 0, A, X, SCALE, W, INFO );
         chkxer('DLATPS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dlatps('U', 'N', 'N', 'N', -1, A, X, SCALE, W, INFO );
         chkxer('DLATPS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // Test error exits for the banded triangular routines.

         // DTBTRS

        srnamc.SRNAMT = 'DTBTRS';
         infoc.INFOT = 1;
         dtbtrs('/', 'N', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtbtrs('U', '/', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtbtrs('U', 'N', '/', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtbtrs('U', 'N', 'N', -1, 0, 0, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtbtrs('U', 'N', 'N', 0, -1, 0, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dtbtrs('U', 'N', 'N', 0, 0, -1, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dtbtrs('U', 'N', 'N', 2, 1, 1, A, 1, X, 2, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dtbtrs('U', 'N', 'N', 2, 0, 1, A, 1, X, 1, INFO );
         chkxer('DTBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTBRFS

        srnamc.SRNAMT = 'DTBRFS';
         infoc.INFOT = 1;
         dtbrfs('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtbrfs('U', '/', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtbrfs('U', 'N', '/', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtbrfs('U', 'N', 'N', -1, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtbrfs('U', 'N', 'N', 0, -1, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dtbrfs('U', 'N', 'N', 0, 0, -1, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dtbrfs('U', 'N', 'N', 2, 1, 1, A, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dtbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         dtbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DTBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DTBCON

        srnamc.SRNAMT = 'DTBCON';
         infoc.INFOT = 1;
         dtbcon('/', 'U', 'N', 0, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtbcon('1', '/', 'N', 0, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dtbcon('1', 'U', '/', 0, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtbcon('1', 'U', 'N', -1, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('DTBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dtbcon('1', 'U', 'N', 0, -1, A, 1, RCOND, W, IW, INFO );
         chkxer('DTBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dtbcon('1', 'U', 'N', 2, 1, A, 1, RCOND, W, IW, INFO );
         chkxer('DTBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DLATBS

        srnamc.SRNAMT = 'DLATBS';
         infoc.INFOT = 1;
         dlatbs('/', 'N', 'N', 'N', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dlatbs('U', '/', 'N', 'N', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dlatbs('U', 'N', '/', 'N', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dlatbs('U', 'N', 'N', '/', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dlatbs('U', 'N', 'N', 'N', -1, 0, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dlatbs('U', 'N', 'N', 'N', 1, -1, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dlatbs('U', 'N', 'N', 'N', 2, 1, A, 1, X, SCALE, W, INFO );
         chkxer('DLATBS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      return;
      }
