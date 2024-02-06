      void zerrtz(PATH, infoc.NUNIT ) {

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
      // ..
      // .. Local Arrays ..
      Complex         A( NMAX, NMAX ), TAU( NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZTZRZF
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
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      C2 = PATH( 2: 3 );
      A[1][1] = DCMPLX( 1.0, -1.0 );
      A[1][2] = DCMPLX( 2.0, -2.0 );
      A[2][2] = DCMPLX( 3.0, -3.0 );
      A[2][1] = DCMPLX( 4.0, -4.0 );
      W[1] = DCMPLX( 0.0, 0.0 );
      W[2] = DCMPLX( 0.0, 0.0 );
      infoc.OK = true;

      // Test error exits for the trapezoidal routines.
      WRITE( NOUT, FMT = * );
      if ( lsamen( 2, C2, 'TZ' ) ) {


         // ZTZRZF

        srnamc.SRNAMT = 'ZTZRZF';
         infoc.INFOT = 1;
         ztzrzf(-1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         ztzrzf(1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         ztzrzf(2, 2, A, 1, TAU, W, 1, INFO );
         chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         ztzrzf(2, 2, A, 2, TAU, W, 0, INFO );
         chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         ztzrzf(2, 3, A, 2, TAU, W, 1, INFO );
         chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
