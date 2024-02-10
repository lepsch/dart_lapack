import 'package:lapack/src/complex.dart';

import 'common.dart';

      void zerrbd(final int PATH, final int NUNIT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX, LW;
      const              NMAX = 4, LW = NMAX ;
      String             C2;
      int                I, INFO, J, NT;
      double             D( NMAX ), E( NMAX ), RW( 4*NMAX );
      Complex         A( NMAX, NMAX ), TP( NMAX ), TQ( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZBDSQR, ZGEBD2, ZGEBRD, ZUNGBR, ZUNMBR
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
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1.0 / (I+J).toDouble();
         } // 10
      } // 20
      infoc.OK = true;
      NT = 0;

      // Test error exits of the SVD routines.

      if ( lsamen( 2, C2, 'BD' ) ) {

         // ZGEBRD

         srnamc.SRNAMT = 'ZGEBRD';
         infoc.INFOT = 1;
         zgebrd(-1, 0, A, 1, D, E, TQ, TP, W, 1, INFO );
         chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgebrd(0, -1, A, 1, D, E, TQ, TP, W, 1, INFO );
         chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgebrd(2, 1, A, 1, D, E, TQ, TP, W, 2, INFO );
         chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zgebrd(2, 1, A, 2, D, E, TQ, TP, W, 1, INFO );
         chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 4;

         // ZGEBD2

         srnamc.SRNAMT = 'ZGEBD2';
         infoc.INFOT = 1;
         zgebd2(-1, 0, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('ZGEBD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgebd2(0, -1, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('ZGEBD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgebd2(2, 1, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('ZGEBD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 3;

         // ZUNGBR

         srnamc.SRNAMT = 'ZUNGBR';
         infoc.INFOT = 1;
         zungbr('/', 0, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zungbr('Q', -1, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zungbr('Q', 0, -1, 0, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zungbr('Q', 0, 1, 0, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zungbr('Q', 1, 0, 1, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zungbr('P', 1, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zungbr('P', 0, 1, 1, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zungbr('Q', 0, 0, -1, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zungbr('Q', 2, 1, 1, A, 1, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zungbr('Q', 2, 2, 1, A, 2, TQ, W, 1, INFO );
         chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 10;

         // ZUNMBR

         srnamc.SRNAMT = 'ZUNMBR';
         infoc.INFOT = 1;
         zunmbr('/', 'L', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zunmbr('Q', '/', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zunmbr('Q', 'L', '/', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zunmbr('Q', 'L', 'C', -1, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zunmbr('Q', 'L', 'C', 0, -1, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zunmbr('Q', 'L', 'C', 0, 0, -1, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zunmbr('Q', 'L', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zunmbr('Q', 'R', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zunmbr('P', 'L', 'C', 2, 0, 2, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zunmbr('P', 'R', 'C', 0, 2, 2, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         zunmbr('Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         zunmbr('Q', 'L', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 0, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         zunmbr('Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 0, INFO );
         chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 13;

         // ZBDSQR

         srnamc.SRNAMT = 'ZBDSQR';
         infoc.INFOT = 1;
         zbdsqr('/', 0, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zbdsqr('U', -1, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zbdsqr('U', 0, -1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zbdsqr('U', 0, 0, -1, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zbdsqr('U', 0, 0, 0, -1, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zbdsqr('U', 2, 1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         zbdsqr('U', 0, 0, 2, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         zbdsqr('U', 2, 0, 0, 1, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         NT = NT + 8;
      }

      // Print a summary line.

      if ( infoc.OK ) {
         WRITE( infoc.NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( infoc.NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT(' ${.a3} routines passed the tests of the error exits (${.i3} tests done)' );
 9998 FORMAT( ' *** ${.a3} routines failed the tests of the error exits ***' );

      }
