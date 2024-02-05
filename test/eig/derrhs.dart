import 'common.dart';

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
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGEBAK, DGEBAL, DGEHRD, DHSEIN, DHSEQR, DORGHR, DORMHR, DTREVC, DTREVC3, DGEHD2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / infoc / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / srnamc / srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = 1.0 / (I+J).toDouble();
         } // 10
         WI[J] = J.toDouble();
         SEL[J] = true;
      } // 20
      infoc.OK = true;
      NT = 0;

      // Test error exits of the nonsymmetric eigenvalue routines.

      if ( lsamen( 2, C2, 'HS' ) ) {

         // DGEBAL

         srnamc.SRNAMT = 'DGEBAL';
         infoc.INFOT = 1;
         dgebal('/', 0, A, 1, ILO, IHI, S, INFO );
         chkxer('DGEBAL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgebal('N', -1, A, 1, ILO, IHI, S, INFO );
         chkxer('DGEBAL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgebal('N', 2, A, 1, ILO, IHI, S, INFO );
         chkxer('DGEBAL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 3;

         // DGEBAK

         srnamc.SRNAMT = 'DGEBAK';
         infoc.INFOT = 1;
         dgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO );
         chkxer('DGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 9;

         // DGEHRD

         srnamc.SRNAMT = 'DGEHRD';
         infoc.INFOT = 1;
         dgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO );
         chkxer('DGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 7;

         // DGEHD2

         srnamc.SRNAMT = 'DGEHD2';
         infoc.INFOT = 1;
         dgehd2(-1, 1, 1, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgehd2(0, 0, 0, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgehd2(0, 2, 0, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgehd2(1, 1, 0, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgehd2(0, 1, 1, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgehd2(2, 1, 1, A, 1, TAU, W, INFO );
         chkxer('DGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 6;

         // DORGHR

         srnamc.SRNAMT = 'DORGHR';
         infoc.INFOT = 1;
         dorghr(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dorghr(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dorghr(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dorghr(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dorghr(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dorghr(2, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dorghr(3, 1, 3, A, 3, TAU, W, 1, INFO );
         chkxer('DORGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 7;

         // DORMHR

         srnamc.SRNAMT = 'DORMHR';
         infoc.INFOT = 1;
         dormhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dormhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dormhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dormhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dormhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dormhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dormhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dormhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dormhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dormhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dormhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dormhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dormhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dormhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dormhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dormhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('DORMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 16;

         // DHSEQR

         srnamc.SRNAMT = 'DHSEQR';
         infoc.INFOT = 1;
         dhseqr('/', 'N', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dhseqr('E', '/', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dhseqr('E', 'N', -1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dhseqr('E', 'N', 0, 0, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dhseqr('E', 'N', 0, 2, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dhseqr('E', 'N', 1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dhseqr('E', 'N', 1, 1, 2, A, 1, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dhseqr('E', 'N', 2, 1, 2, A, 1, WR, WI, C, 2, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dhseqr('E', 'V', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dhseqr('E', 'N', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO );
         chkxer('DHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 10;

         // DHSEIN

         srnamc.SRNAMT = 'DHSEIN';
         infoc.INFOT = 1;
         dhsein('/', 'N', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dhsein('R', '/', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dhsein('R', 'N', '/', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dhsein('R', 'N', 'N', SEL, -1, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dhsein('R', 'N', 'N', SEL, 2, A, 1, WR, WI, VL, 1, VR, 2, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dhsein('L', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dhsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         dhsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 2, 1, M, W, IFAILL, IFAILR, INFO );
         chkxer('DHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 8;

         // DTREVC

         srnamc.SRNAMT = 'DTREVC';
         infoc.INFOT = 1;
         dtrevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtrevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dtrevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dtrevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dtrevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dtrevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, INFO );
         chkxer('DTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 7;

         // DTREVC3

         srnamc.SRNAMT = 'DTREVC3';
         infoc.INFOT = 1;
         dtrevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dtrevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dtrevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dtrevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dtrevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dtrevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dtrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         dtrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, INFO );
         chkxer('DTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 8;
      }

      // Print a summary line.

      if ( infoc.OK ) {
         WRITE( infoc.NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( infoc.NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits', ' (', I3, ' tests done)' );
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' );

      return;
      }
