import 'common.dart';

      void derrorhr_col(PATH, NUNIT ) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String   (LEN=3)   PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), D(NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DORHR_COL
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String   (LEN=32)  srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / (I+J).toDouble();
            T[I][J] = 1.0 / (I+J).toDouble();
         }
         D[J] = 0.0;
      }
      infoc.OK = true;

      // Error exits for Householder reconstruction

      // DORHR_COL

      srnamc.SRNAMT = 'DORHR_COL';

      infoc.INFOT = 1;
      dorhr_col(-1, 0, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 2;
      dorhr_col(0, -1, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      dorhr_col(1, 2, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 3;
      dorhr_col(0, 0, -1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      dorhr_col(0, 0, 0, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 5;
      dorhr_col(0, 0, 1, A, -1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      dorhr_col(0, 0, 1, A, 0, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      dorhr_col(2, 0, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 7;
      dorhr_col(0, 0, 1, A, 1, T, -1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      dorhr_col(0, 0, 1, A, 1, T, 0, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      dorhr_col(4, 3, 2, A, 4, T, 1, D, INFO );
      chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      return;
      }
