      void derrhs(PATH, NUNIT ) {

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
      const              NMAX = 3, LW = ( NMAX+2 )*( NMAX+2 )+NMAX ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, M, NT;
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      int                IFAILL( NMAX ), IFAILR( NMAX );
      double             A( NMAX, NMAX ), C( NMAX, NMAX ), S( NMAX ), TAU( NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), W( LW ), WI( NMAX ), WR( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGEBAK, DGEBAL, DGEHRD, DHSEIN, DHSEQR, DORGHR, DORMHR, DTREVC, DTREVC3, DGEHD2
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
            A( I, J ) = 1.0 / DBLE( I+J );
         } // 10
         WI( J ) = DBLE( J );
         SEL( J ) = true;
      } // 20
      OK = true;
      NT = 0;

      // Test error exits of the nonsymmetric eigenvalue routines.

      if ( LSAMEN( 2, C2, 'HS' ) ) {

         // DGEBAL

         SRNAMT = 'DGEBAL';
         INFOT = 1;
         dgebal('/', 0, A, 1, ILO, IHI, S, INFO );
         chkxer('DGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgebal('N', -1, A, 1, ILO, IHI, S, INFO );
         chkxer('DGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgebal('N', 2, A, 1, ILO, IHI, S, INFO );
         chkxer('DGEBAL', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // DGEBAK

         SRNAMT = 'DGEBAK';
         INFOT = 1;
         dgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO );
         chkxer('DGEBAK', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DGEHRD

         SRNAMT = 'DGEHRD';
         INFOT = 1;
         dgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO );
         chkxer('DGEHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // DGEHD2

         SRNAMT = 'DGEHD2';
         INFOT = 1;
         dgehd2(-1, 1, 1, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgehd2(0, 0, 0, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgehd2(0, 2, 0, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgehd2(1, 1, 0, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgehd2(0, 1, 1, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgehd2(2, 1, 1, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DORGHR

         SRNAMT = 'DORGHR';
         INFOT = 1;
         dorghr(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dorghr(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dorghr(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dorghr(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dorghr(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dorghr(2, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dorghr(3, 1, 3, A, 3, TAU, W, 1, INFO );
         chkxer('DORGHR', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // DORMHR

         SRNAMT = 'DORMHR';
         INFOT = 1;
         dormhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dormhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dormhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dormhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dormhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dormhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dormhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dormhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dormhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dormhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dormhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dormhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dormhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dormhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dormhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dormhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('DORMHR', INFOT, NOUT, LERR, OK );
         NT = NT + 16;

         // DHSEQR

         SRNAMT = 'DHSEQR';
         INFOT = 1;
         dhseqr('/', 'N', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dhseqr('E', '/', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dhseqr('E', 'N', -1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dhseqr('E', 'N', 0, 0, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dhseqr('E', 'N', 0, 2, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dhseqr('E', 'N', 1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dhseqr('E', 'N', 1, 1, 2, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dhseqr('E', 'N', 2, 1, 2, A, 1, WR, WI, C, 2, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dhseqr('E', 'V', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dhseqr('E', 'N', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DHSEIN

         SRNAMT = 'DHSEIN';
         INFOT = 1;
         dhsein('/', 'N', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dhsein('R', '/', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dhsein('R', 'N', '/', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dhsein('R', 'N', 'N', SEL, -1, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dhsein('R', 'N', 'N', SEL, 2, A, 1, WR, WI, VL, 1, VR, 2, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dhsein('L', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dhsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dhsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 2, 1, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // DTREVC

         SRNAMT = 'DTREVC';
         INFOT = 1;
         dtrevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dtrevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dtrevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dtrevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dtrevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dtrevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dtrevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, INFO );
         chkxer('DTREVC', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // DTREVC3

         SRNAMT = 'DTREVC3';
         INFOT = 1;
         dtrevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dtrevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dtrevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dtrevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dtrevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dtrevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dtrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dtrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, INFO );
         chkxer('DTREVC3', INFOT, NOUT, LERR, OK );
         NT = NT + 8;
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

      // End of DERRHS

      }
