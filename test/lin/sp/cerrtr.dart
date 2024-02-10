      void cerrtr(final int PATH, final int NUNIT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      String             C2;
      int                INFO;
      double               RCOND, SCALE, SCALES(0);
      double               R1( NMAX ), R2( NMAX ), RW( NMAX );
      Complex            A( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CLATBS, CLATPS, CLATRS, CLATRS3, CTBCON, CTBRFS, CTBTRS, CTPCON, CTPRFS, CTPTRI, CTPTRS, CTRCON, CTRRFS, CTRTI2, CTRTRI, CTRTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A[1][1] = 1.;
      A[1][2] = 2.;
      A[2][2] = 3.;
      A[2][1] = 4.;
      OK = true;

      // Test error exits for the general triangular routines.

      if ( lsamen( 2, C2, 'TR' ) ) {

         // CTRTRI

        srnamc.SRNAMT = 'CTRTRI';
         INFOT = 1;
         ctrtri('/', 'N', 0, A, 1, INFO );
         chkxer('CTRTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctrtri('U', '/', 0, A, 1, INFO );
         chkxer('CTRTRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctrtri('U', 'N', -1, A, 1, INFO );
         chkxer('CTRTRI', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctrtri('U', 'N', 2, A, 1, INFO );
         chkxer('CTRTRI', INFOT, NOUT, LERR, OK );

         // CTRTI2

        srnamc.SRNAMT = 'CTRTI2';
         INFOT = 1;
         ctrti2('/', 'N', 0, A, 1, INFO );
         chkxer('CTRTI2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctrti2('U', '/', 0, A, 1, INFO );
         chkxer('CTRTI2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctrti2('U', 'N', -1, A, 1, INFO );
         chkxer('CTRTI2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctrti2('U', 'N', 2, A, 1, INFO );
         chkxer('CTRTI2', INFOT, NOUT, LERR, OK );


         // CTRTRS

        srnamc.SRNAMT = 'CTRTRS';
         INFOT = 1;
         ctrtrs('/', 'N', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('CTRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctrtrs('U', '/', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('CTRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctrtrs('U', 'N', '/', 0, 0, A, 1, X, 1, INFO );
         chkxer('CTRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctrtrs('U', 'N', 'N', -1, 0, A, 1, X, 1, INFO );
         chkxer('CTRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctrtrs('U', 'N', 'N', 0, -1, A, 1, X, 1, INFO );
         chkxer('CTRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;

         // CTRRFS

        srnamc.SRNAMT = 'CTRRFS';
         INFOT = 1;
         ctrrfs('/', 'N', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctrrfs('U', '/', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctrrfs('U', 'N', '/', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctrrfs('U', 'N', 'N', -1, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctrrfs('U', 'N', 'N', 0, -1, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ctrrfs('U', 'N', 'N', 2, 1, A, 1, B, 2, X, 2, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ctrrfs('U', 'N', 'N', 2, 1, A, 2, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ctrrfs('U', 'N', 'N', 2, 1, A, 2, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTRRFS', INFOT, NOUT, LERR, OK );

         // CTRCON

        srnamc.SRNAMT = 'CTRCON';
         INFOT = 1;
         ctrcon('/', 'U', 'N', 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTRCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctrcon('1', '/', 'N', 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTRCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctrcon('1', 'U', '/', 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTRCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctrcon('1', 'U', 'N', -1, A, 1, RCOND, W, RW, INFO );
         chkxer('CTRCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctrcon('1', 'U', 'N', 2, A, 1, RCOND, W, RW, INFO );
         chkxer('CTRCON', INFOT, NOUT, LERR, OK );

         // CLATRS

        srnamc.SRNAMT = 'CLATRS';
         INFOT = 1;
         clatrs('/', 'N', 'N', 'N', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         clatrs('U', '/', 'N', 'N', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         clatrs('U', 'N', '/', 'N', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         clatrs('U', 'N', 'N', '/', 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         clatrs('U', 'N', 'N', 'N', -1, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         clatrs('U', 'N', 'N', 'N', 2, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATRS', INFOT, NOUT, LERR, OK );

         // CLATRS3

        srnamc.SRNAMT = 'CLATRS3';
         INFOT = 1;
         clatrs3('/', 'N', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         clatrs3('U', '/', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         clatrs3('U', 'N', '/', 'N', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         clatrs3('U', 'N', 'N', '/', 0, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         clatrs3('U', 'N', 'N', 'N', -1, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         clatrs3('U', 'N', 'N', 'N', 0, -1, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         clatrs3('U', 'N', 'N', 'N', 2, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         clatrs3('U', 'N', 'N', 'N', 2, 0, A, 2, X, 1, SCALES, RW, RW( 2 ), 1, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         clatrs3('U', 'N', 'N', 'N', 1, 0, A, 1, X, 1, SCALES, RW, RW( 2 ), 0, INFO );
         chkxer('CLATRS3', INFOT, NOUT, LERR, OK );

      // Test error exits for the packed triangular routines.

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // CTPTRI

        srnamc.SRNAMT = 'CTPTRI';
         INFOT = 1;
         ctptri('/', 'N', 0, A, INFO );
         chkxer('CTPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctptri('U', '/', 0, A, INFO );
         chkxer('CTPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctptri('U', 'N', -1, A, INFO );
         chkxer('CTPTRI', INFOT, NOUT, LERR, OK );

         // CTPTRS

        srnamc.SRNAMT = 'CTPTRS';
         INFOT = 1;
         ctptrs('/', 'N', 'N', 0, 0, A, X, 1, INFO );
         chkxer('CTPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctptrs('U', '/', 'N', 0, 0, A, X, 1, INFO );
         chkxer('CTPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctptrs('U', 'N', '/', 0, 0, A, X, 1, INFO );
         chkxer('CTPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctptrs('U', 'N', 'N', -1, 0, A, X, 1, INFO );
         chkxer('CTPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctptrs('U', 'N', 'N', 0, -1, A, X, 1, INFO );
         chkxer('CTPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctptrs('U', 'N', 'N', 2, 1, A, X, 1, INFO );
         chkxer('CTPTRS', INFOT, NOUT, LERR, OK );

         // CTPRFS

        srnamc.SRNAMT = 'CTPRFS';
         INFOT = 1;
         ctprfs('/', 'N', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctprfs('U', '/', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctprfs('U', 'N', '/', 0, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctprfs('U', 'N', 'N', -1, 0, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctprfs('U', 'N', 'N', 0, -1, A, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctprfs('U', 'N', 'N', 2, 1, A, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctprfs('U', 'N', 'N', 2, 1, A, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTPRFS', INFOT, NOUT, LERR, OK );

         // CTPCON

        srnamc.SRNAMT = 'CTPCON';
         INFOT = 1;
         ctpcon('/', 'U', 'N', 0, A, RCOND, W, RW, INFO );
         chkxer('CTPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctpcon('1', '/', 'N', 0, A, RCOND, W, RW, INFO );
         chkxer('CTPCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctpcon('1', 'U', '/', 0, A, RCOND, W, RW, INFO );
         chkxer('CTPCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctpcon('1', 'U', 'N', -1, A, RCOND, W, RW, INFO );
         chkxer('CTPCON', INFOT, NOUT, LERR, OK );

         // CLATPS

        srnamc.SRNAMT = 'CLATPS';
         INFOT = 1;
         clatps('/', 'N', 'N', 'N', 0, A, X, SCALE, RW, INFO );
         chkxer('CLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         clatps('U', '/', 'N', 'N', 0, A, X, SCALE, RW, INFO );
         chkxer('CLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         clatps('U', 'N', '/', 'N', 0, A, X, SCALE, RW, INFO );
         chkxer('CLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         clatps('U', 'N', 'N', '/', 0, A, X, SCALE, RW, INFO );
         chkxer('CLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         clatps('U', 'N', 'N', 'N', -1, A, X, SCALE, RW, INFO );
         chkxer('CLATPS', INFOT, NOUT, LERR, OK );

      // Test error exits for the banded triangular routines.

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // CTBTRS

        srnamc.SRNAMT = 'CTBTRS';
         INFOT = 1;
         ctbtrs('/', 'N', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctbtrs('U', '/', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctbtrs('U', 'N', '/', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctbtrs('U', 'N', 'N', -1, 0, 0, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctbtrs('U', 'N', 'N', 0, -1, 0, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctbtrs('U', 'N', 'N', 0, 0, -1, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctbtrs('U', 'N', 'N', 2, 1, 1, A, 1, X, 2, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctbtrs('U', 'N', 'N', 2, 0, 1, A, 1, X, 1, INFO );
         chkxer('CTBTRS', INFOT, NOUT, LERR, OK );

         // CTBRFS

        srnamc.SRNAMT = 'CTBRFS';
         INFOT = 1;
         ctbrfs('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctbrfs('U', '/', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctbrfs('U', 'N', '/', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctbrfs('U', 'N', 'N', -1, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctbrfs('U', 'N', 'N', 0, -1, 0, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctbrfs('U', 'N', 'N', 0, 0, -1, A, 1, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctbrfs('U', 'N', 'N', 2, 1, 1, A, 1, B, 2, X, 2, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ctbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('CTBRFS', INFOT, NOUT, LERR, OK );

         // CTBCON

        srnamc.SRNAMT = 'CTBCON';
         INFOT = 1;
         ctbcon('/', 'U', 'N', 0, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctbcon('1', '/', 'N', 0, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctbcon('1', 'U', '/', 0, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTBCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctbcon('1', 'U', 'N', -1, 0, A, 1, RCOND, W, RW, INFO );
         chkxer('CTBCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctbcon('1', 'U', 'N', 0, -1, A, 1, RCOND, W, RW, INFO );
         chkxer('CTBCON', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ctbcon('1', 'U', 'N', 2, 1, A, 1, RCOND, W, RW, INFO );
         chkxer('CTBCON', INFOT, NOUT, LERR, OK );

         // CLATBS

        srnamc.SRNAMT = 'CLATBS';
         INFOT = 1;
         clatbs('/', 'N', 'N', 'N', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         clatbs('U', '/', 'N', 'N', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         clatbs('U', 'N', '/', 'N', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         clatbs('U', 'N', 'N', '/', 0, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         clatbs('U', 'N', 'N', 'N', -1, 0, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         clatbs('U', 'N', 'N', 'N', 1, -1, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         clatbs('U', 'N', 'N', 'N', 2, 1, A, 1, X, SCALE, RW, INFO );
         chkxer('CLATBS', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
