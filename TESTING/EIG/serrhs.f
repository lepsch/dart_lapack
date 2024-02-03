      SUBROUTINE SERRHS( PATH, NUNIT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX, LW;
      const              NMAX = 3, LW = ( NMAX+2 )*( NMAX+2 )+NMAX ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, ILO, IHI, INFO, J, M, NT;
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      int                IFAILL( NMAX ), IFAILR( NMAX );
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ), TAU( NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), W( LW ), WI( NMAX ), WR( NMAX ), S( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SGEBAK, SGEBAL, SGEHRD, SHSEIN, SHSEQR, SORGHR, SORMHR, STREVC, STREVC3, SGEHD2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
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
            A( I, J ) = 1. / REAL( I+J );
         } // 10
         WI( J ) = REAL( J );
         SEL( J ) = true;
      } // 20
      OK = true;
      NT = 0;

      // Test error exits of the nonsymmetric eigenvalue routines.

      if ( LSAMEN( 2, C2, 'HS' ) ) {

         // SGEBAL

         SRNAMT = 'SGEBAL';
         INFOT = 1;
         sgebal('/', 0, A, 1, ILO, IHI, S, INFO );
         chkxer('SGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgebal('N', -1, A, 1, ILO, IHI, S, INFO );
         chkxer('SGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgebal('N', 2, A, 1, ILO, IHI, S, INFO );
         chkxer('SGEBAL', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // SGEBAK

         SRNAMT = 'SGEBAK';
         INFOT = 1;
         sgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO );
         chkxer('SGEBAK', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SGEHRD

         SRNAMT = 'SGEHRD';
         INFOT = 1;
         sgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO );
         chkxer('SGEHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // SGEHD2

         SRNAMT = 'SGEHD2';
         INFOT = 1;
         sgehd2(-1, 1, 1, A, 1, TAU, W, INFO );
         chkxer('SGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgehd2(0, 0, 0, A, 1, TAU, W, INFO );
         chkxer('SGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgehd2(0, 2, 0, A, 1, TAU, W, INFO );
         chkxer('SGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgehd2(1, 1, 0, A, 1, TAU, W, INFO );
         chkxer('SGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgehd2(0, 1, 1, A, 1, TAU, W, INFO );
         chkxer('SGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgehd2(2, 1, 1, A, 1, TAU, W, INFO );
         chkxer('SGEHD2', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SORGHR

         SRNAMT = 'SORGHR';
         INFOT = 1;
         sorghr(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sorghr(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sorghr(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sorghr(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sorghr(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sorghr(2, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sorghr(3, 1, 3, A, 3, TAU, W, 1, INFO );
         chkxer('SORGHR', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // SORMHR

         SRNAMT = 'SORMHR';
         INFOT = 1;
         sormhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sormhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sormhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sormhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sormhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sormhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sormhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sormhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sormhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sormhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sormhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sormhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sormhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sormhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sormhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sormhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('SORMHR', INFOT, NOUT, LERR, OK );
         NT = NT + 16;

         // SHSEQR

         SRNAMT = 'SHSEQR';
         INFOT = 1;
         shseqr('/', 'N', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         shseqr('E', '/', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         shseqr('E', 'N', -1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         shseqr('E', 'N', 0, 0, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         shseqr('E', 'N', 0, 2, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         shseqr('E', 'N', 1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         shseqr('E', 'N', 1, 1, 2, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         shseqr('E', 'N', 2, 1, 2, A, 1, WR, WI, C, 2, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         shseqr('E', 'V', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         shseqr('E', 'N', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO );
         chkxer('SHSEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // SHSEIN

         SRNAMT = 'SHSEIN';
         INFOT = 1;
         shsein('/', 'N', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         shsein('R', '/', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         shsein('R', 'N', '/', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         shsein('R', 'N', 'N', SEL, -1, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         shsein('R', 'N', 'N', SEL, 2, A, 1, WR, WI, VL, 1, VR, 2, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         shsein('L', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         shsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         shsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 2, 1, M, W, IFAILL, IFAILR, INFO );
         chkxer('SHSEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // STREVC

         SRNAMT = 'STREVC';
         INFOT = 1;
         strevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         strevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         strevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         strevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         strevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         strevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, INFO );
         chkxer('STREVC', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // STREVC3

         SRNAMT = 'STREVC3';
         INFOT = 1;
         strevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         strevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         strevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         strevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         strevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         strevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         strevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         strevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, INFO );
         chkxer('STREVC3', INFOT, NOUT, LERR, OK );
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

      // End of SERRHS

      }
