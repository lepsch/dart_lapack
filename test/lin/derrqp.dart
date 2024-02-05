import 'common.dart';
      void derrqp(PATH, NUNIT ) {

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
      const              NMAX = 3 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                INFO, LW;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             A( NMAX, NMAX ), TAU( NMAX ), W( 3*NMAX+1 );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQP3
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      LW = 3*NMAX + 1;
      A[1, 1] = 1.0;
      A[1, 2] = 2.0;
      A[2, 2] = 3.0;
      A[2, 1] = 4.0;
      infoc.OK = true;

      if ( lsamen( 2, C2, 'QP' ) ) {

         // Test error exits for QR factorization with pivoting

         // DGEQP3

         srnamc.SRNAMT = 'DGEQP3';
         infoc.INFOT = 1;
         dgeqp3(-1, 0, A, 1, IP, TAU, W, LW, INFO );
         chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgeqp3(1, -1, A, 1, IP, TAU, W, LW, INFO );
         chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgeqp3(2, 3, A, 1, IP, TAU, W, LW, INFO );
         chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgeqp3(2, 2, A, 2, IP, TAU, W, LW-10, INFO );
         chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      return;
      }
