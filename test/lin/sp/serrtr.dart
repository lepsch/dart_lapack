      void serrtr(PATH, NUNIT ) {

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
      double               RCOND, SCALE, SCALES(0);
      // ..
      // .. Local Arrays ..
      int                IW( NMAX );
      double               A( NMAX, NMAX ), B( NMAX ), R1( NMAX ), R2( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SLATBS, SLATPS, SLATRS, SLATRS3, STBCON, STBRFS, STBTRS, STPCON, STPRFS, STPTRI, STPTRS, STRCON, STRRFS, STRTI2, STRTRI, STRTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A[1, 1] = 1.;
      A[1, 2] = 2.;
      A[2, 2] = 3.;
      A[2, 1] = 4.;
      OK = true;

      if ( lsamen( 2, C2, 'TR' ) ) {

         // Test error exits for the general triangular routines.

         // STRTRI

        srnamc.SRNAMT = 'STRTRI';
         INFOT = 1;
         strtri('/', 'N', 0, A, 1, INFO );
         chkxer('STRTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strtri('U', '/', 0, A, 1, INFO );
         chkxer('STRTRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         strtri('U', 'N', -1, A, 1, INFO );
         chkxer('STRTRI', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         strtri('U', 'N', 2, A, 1, INFO );
         chkxer('STRTRI', INFOT, NOUT, LERR, OK );

         // STRTI2

        srnamc.SRNAMT = 'STRTI2';
         INFOT = 1;
         strti2('/', 'N', 0, A, 1, INFO );
         chkxer('STRTI2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strti2('U', '/', 0, A, 1, INFO );
         chkxer('STRTI2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         strti2('U', 'N', -1, A, 1, INFO );
         chkxer('STRTI2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         strti2('U', 'N', 2, A, 1, INFO );
         chkxer('STRTI2', INFOT, NOUT, LERR, OK );

         // STRTRS

        srnamc.SRNAMT = 'STRTRS';
         INFOT = 1;
         strtrs('/', 'N', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strtrs('U', '/', 'N', 0, 0, A, 1, X, 1, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         strtrs('U', 'N', '/', 0, 0, A, 1, X, 1, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         strtrs('U', 'N', 'N', -1, 0, A, 1, X, 1, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         strtrs('U', 'N', 'N', 0, -1, A, 1, X, 1, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         strtrs('U', 'N', 'N', 2, 1, A, 1, X, 2, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         strtrs('U', 'N', 'N', 2, 1, A, 2, X, 1, INFO );
         chkxer('STRTRS', INFOT, NOUT, LERR, OK );

         // STRRFS

        srnamc.SRNAMT = 'STRRFS';
         INFOT = 1;
         strrfs('/', 'N', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strrfs('U', '/', 'N', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         strrfs('U', 'N', '/', 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         strrfs('U', 'N', 'N', -1, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         strrfs('U', 'N', 'N', 0, -1, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         strrfs('U', 'N', 'N', 2, 1, A, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         strrfs('U', 'N', 'N', 2, 1, A, 2, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         strrfs('U', 'N', 'N', 2, 1, A, 2, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('STRRFS', INFOT, NOUT, LERR, OK );

         // STRCON

        srnamc.SRNAMT = 'STRCON';
         INFOT = 1;
         strcon('/', 'U', 'N', 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STRCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strcon('1', '/', 'N', 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STRCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         strcon('1', 'U', '/', 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STRCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         strcon('1', 'U', 'N', -1, A, 1, RCOND, W, IW, INFO );
         chkxer('STRCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         strcon('1', 'U', 'N', 2, A, 1, RCOND, W, IW, INFO );
         chkxer('STRCON', INFOT, NOUT, LERR, OK );

         // SLATRS

        srnamc.SRNAMT = 'SLATRS';
         INFOT = 1;
         slatrs('/', 'N', 'N', 'N', 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         slatrs('U', '/', 'N', 'N', 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         slatrs('U', 'N', '/', 'N', 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         slatrs('U', 'N', 'N', '/', 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         slatrs('U', 'N', 'N', 'N', -1, A, 1, X, SCALE, W, INFO );
         chkxer('SLATRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         slatrs('U', 'N', 'N', 'N', 2, A, 1, X, SCALE, W, INFO );
         chkxer('SLATRS', INFOT, NOUT, LERR, OK );

         // SLATRS3

        srnamc.SRNAMT = 'SLATRS3';
         INFOT = 1;
         slatrs3('/', 'N', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         slatrs3('U', '/', 'N', 'N', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         slatrs3('U', 'N', '/', 'N', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         slatrs3('U', 'N', 'N', '/', 0, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         slatrs3('U', 'N', 'N', 'N', -1, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         slatrs3('U', 'N', 'N', 'N', 0, -1, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         slatrs3('U', 'N', 'N', 'N', 2, 0, A, 1, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         slatrs3('U', 'N', 'N', 'N', 2, 0, A, 2, X, 1, SCALES, W, W( 2 ), 1, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         slatrs3('U', 'N', 'N', 'N', 1, 0, A, 1, X, 1, SCALES, W, W( 2 ), 0, INFO );
         chkxer('SLATRS3', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // Test error exits for the packed triangular routines.

         // STPTRI

        srnamc.SRNAMT = 'STPTRI';
         INFOT = 1;
         stptri('/', 'N', 0, A, INFO );
         chkxer('STPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stptri('U', '/', 0, A, INFO );
         chkxer('STPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stptri('U', 'N', -1, A, INFO );
         chkxer('STPTRI', INFOT, NOUT, LERR, OK );

         // STPTRS

        srnamc.SRNAMT = 'STPTRS';
         INFOT = 1;
         stptrs('/', 'N', 'N', 0, 0, A, X, 1, INFO );
         chkxer('STPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stptrs('U', '/', 'N', 0, 0, A, X, 1, INFO );
         chkxer('STPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stptrs('U', 'N', '/', 0, 0, A, X, 1, INFO );
         chkxer('STPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stptrs('U', 'N', 'N', -1, 0, A, X, 1, INFO );
         chkxer('STPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stptrs('U', 'N', 'N', 0, -1, A, X, 1, INFO );
         chkxer('STPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stptrs('U', 'N', 'N', 2, 1, A, X, 1, INFO );
         chkxer('STPTRS', INFOT, NOUT, LERR, OK );

         // STPRFS

        srnamc.SRNAMT = 'STPRFS';
         INFOT = 1;
         stprfs('/', 'N', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stprfs('U', '/', 'N', 0, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stprfs('U', 'N', '/', 0, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stprfs('U', 'N', 'N', -1, 0, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stprfs('U', 'N', 'N', 0, -1, A, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stprfs('U', 'N', 'N', 2, 1, A, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stprfs('U', 'N', 'N', 2, 1, A, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('STPRFS', INFOT, NOUT, LERR, OK );

         // STPCON

        srnamc.SRNAMT = 'STPCON';
         INFOT = 1;
         stpcon('/', 'U', 'N', 0, A, RCOND, W, IW, INFO );
         chkxer('STPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stpcon('1', '/', 'N', 0, A, RCOND, W, IW, INFO );
         chkxer('STPCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stpcon('1', 'U', '/', 0, A, RCOND, W, IW, INFO );
         chkxer('STPCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stpcon('1', 'U', 'N', -1, A, RCOND, W, IW, INFO );
         chkxer('STPCON', INFOT, NOUT, LERR, OK );

         // SLATPS

        srnamc.SRNAMT = 'SLATPS';
         INFOT = 1;
         slatps('/', 'N', 'N', 'N', 0, A, X, SCALE, W, INFO );
         chkxer('SLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         slatps('U', '/', 'N', 'N', 0, A, X, SCALE, W, INFO );
         chkxer('SLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         slatps('U', 'N', '/', 'N', 0, A, X, SCALE, W, INFO );
         chkxer('SLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         slatps('U', 'N', 'N', '/', 0, A, X, SCALE, W, INFO );
         chkxer('SLATPS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         slatps('U', 'N', 'N', 'N', -1, A, X, SCALE, W, INFO );
         chkxer('SLATPS', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // Test error exits for the banded triangular routines.

         // STBTRS

        srnamc.SRNAMT = 'STBTRS';
         INFOT = 1;
         stbtrs('/', 'N', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stbtrs('U', '/', 'N', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stbtrs('U', 'N', '/', 0, 0, 0, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stbtrs('U', 'N', 'N', -1, 0, 0, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stbtrs('U', 'N', 'N', 0, -1, 0, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         stbtrs('U', 'N', 'N', 0, 0, -1, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stbtrs('U', 'N', 'N', 2, 1, 1, A, 1, X, 2, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stbtrs('U', 'N', 'N', 2, 0, 1, A, 1, X, 1, INFO );
         chkxer('STBTRS', INFOT, NOUT, LERR, OK );

         // STBRFS

        srnamc.SRNAMT = 'STBRFS';
         INFOT = 1;
         stbrfs('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stbrfs('U', '/', 'N', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stbrfs('U', 'N', '/', 0, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stbrfs('U', 'N', 'N', -1, 0, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stbrfs('U', 'N', 'N', 0, -1, 0, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         stbrfs('U', 'N', 'N', 0, 0, -1, A, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stbrfs('U', 'N', 'N', 2, 1, 1, A, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         stbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('STBRFS', INFOT, NOUT, LERR, OK );

         // STBCON

        srnamc.SRNAMT = 'STBCON';
         INFOT = 1;
         stbcon('/', 'U', 'N', 0, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stbcon('1', '/', 'N', 0, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stbcon('1', 'U', '/', 0, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STBCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stbcon('1', 'U', 'N', -1, 0, A, 1, RCOND, W, IW, INFO );
         chkxer('STBCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stbcon('1', 'U', 'N', 0, -1, A, 1, RCOND, W, IW, INFO );
         chkxer('STBCON', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         stbcon('1', 'U', 'N', 2, 1, A, 1, RCOND, W, IW, INFO );
         chkxer('STBCON', INFOT, NOUT, LERR, OK );

         // SLATBS

        srnamc.SRNAMT = 'SLATBS';
         INFOT = 1;
         slatbs('/', 'N', 'N', 'N', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         slatbs('U', '/', 'N', 'N', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         slatbs('U', 'N', '/', 'N', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         slatbs('U', 'N', 'N', '/', 0, 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         slatbs('U', 'N', 'N', 'N', -1, 0, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         slatbs('U', 'N', 'N', 'N', 1, -1, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         slatbs('U', 'N', 'N', 'N', 2, 1, A, 1, X, SCALE, W, INFO );
         chkxer('SLATBS', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
