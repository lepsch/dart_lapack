      void zerrhs(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX, LW;
      const              NMAX = 3, LW = NMAX*NMAX ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, M, NT;
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      int                IFAILL( NMAX ), IFAILR( NMAX );
      double             RW( NMAX ), S( NMAX );
      Complex         A( NMAX, NMAX ), C( NMAX, NMAX ), TAU( NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), W( LW ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEIN, ZHSEQR, ZUNGHR, ZUNMHR, ZTREVC, ZTREVC3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = 1.0 / (I+J).toDouble();
         } // 10
         SEL[J] = true;
      } // 20
      OK = true;
      NT = 0;

      // Test error exits of the nonsymmetric eigenvalue routines.

      if ( lsamen( 2, C2, 'HS' ) ) {

         // ZGEBAL

         SRNAMT = 'ZGEBAL';
         INFOT = 1;
         zgebal('/', 0, A, 1, ILO, IHI, S, INFO );
         chkxer('ZGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgebal('N', -1, A, 1, ILO, IHI, S, INFO );
         chkxer('ZGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgebal('N', 2, A, 1, ILO, IHI, S, INFO );
         chkxer('ZGEBAL', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // ZGEBAK

         SRNAMT = 'ZGEBAK';
         INFOT = 1;
         zgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO );
         chkxer('ZGEBAK', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // ZGEHRD

         SRNAMT = 'ZGEHRD';
         INFOT = 1;
         zgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO );
         chkxer('ZGEHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // ZGEHD2

         SRNAMT = 'ZGEHD2';
         INFOT = 1;
         zgehd2(-1, 1, 1, A, 1, TAU, W, INFO );
         chkxer('ZGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgehd2(0, 0, 0, A, 1, TAU, W, INFO );
         chkxer('ZGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgehd2(0, 2, 0, A, 1, TAU, W, INFO );
         chkxer('ZGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgehd2(1, 1, 0, A, 1, TAU, W, INFO );
         chkxer('ZGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgehd2(0, 1, 1, A, 1, TAU, W, INFO );
         chkxer('ZGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgehd2(2, 1, 1, A, 1, TAU, W, INFO );
         chkxer('ZGEHD2', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // ZUNGHR

         SRNAMT = 'ZUNGHR';
         INFOT = 1;
         zunghr(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zunghr(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zunghr(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zunghr(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zunghr(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zunghr(2, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zunghr(3, 1, 3, A, 3, TAU, W, 1, INFO );
         chkxer('ZUNGHR', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // ZUNMHR

         SRNAMT = 'ZUNMHR';
         INFOT = 1;
         zunmhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zunmhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zunmhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zunmhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zunmhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zunmhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zunmhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zunmhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zunmhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zunmhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zunmhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zunmhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zunmhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zunmhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zunmhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zunmhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('ZUNMHR', INFOT, NOUT, LERR, OK );
         NT = NT + 16;

         // ZHSEQR

         SRNAMT = 'ZHSEQR';
         INFOT = 1;
         zhseqr('/', 'N', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhseqr('E', '/', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhseqr('E', 'N', -1, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhseqr('E', 'N', 0, 0, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhseqr('E', 'N', 0, 2, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhseqr('E', 'N', 1, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhseqr('E', 'N', 1, 1, 2, A, 1, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhseqr('E', 'N', 2, 1, 2, A, 1, X, C, 2, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zhseqr('E', 'V', 2, 1, 2, A, 2, X, C, 1, W, 1, INFO );
         chkxer('ZHSEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // ZHSEIN

         SRNAMT = 'ZHSEIN';
         INFOT = 1;
         zhsein('/', 'N', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhsein('R', '/', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhsein('R', 'N', '/', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhsein('R', 'N', 'N', SEL, -1, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhsein('R', 'N', 'N', SEL, 2, A, 1, X, VL, 1, VR, 2, 4, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zhsein('L', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zhsein('R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhsein('R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 2, 1, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('ZHSEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // ZTREVC

         SRNAMT = 'ZTREVC';
         INFOT = 1;
         ztrevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ztrevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ztrevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ztrevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ztrevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ztrevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ztrevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, RW, INFO );
         chkxer('ZTREVC', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // ZTREVC3

         SRNAMT = 'ZTREVC3';
         INFOT = 1;
         ztrevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ztrevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ztrevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ztrevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, RW, 2, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ztrevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ztrevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ztrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, RW, 2, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         ztrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, RW, 2, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         ztrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, LW, RW, 1, INFO );
         chkxer('ZTREVC3', INFOT, NOUT, LERR, OK );
         NT = NT + 9;
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits', ' (', I3, ' tests done)' );
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' );

      return;
      }
