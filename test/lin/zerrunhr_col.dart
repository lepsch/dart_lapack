      void zerrunhr_col(PATH, infoc.NUNIT ) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String   (LEN=3)   PATH;
      int                infoc.NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      Complex         A( NMAX, NMAX ), T( NMAX, NMAX ), D(NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZUNHR_COL
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String   (LEN=32) srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = DCMPLX( 1.0 / (I+J).toDouble() );
            T[I][J] = DCMPLX( 1.0 / (I+J).toDouble() );
         }
         D[J] = ( 0.0, 0.0 );
      }
      infoc.OK = true;

      // Error exits for Householder reconstruction

      // ZUNHR_COL

     srnamc.SRNAMT = 'ZUNHR_COL';

      infoc.INFOT = 1;
      zunhr_col(-1, 0, 1, A, 1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 2;
      zunhr_col(0, -1, 1, A, 1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      zunhr_col(1, 2, 1, A, 1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 3;
      zunhr_col(0, 0, -1, A, 1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      zunhr_col(0, 0, 0, A, 1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 5;
      zunhr_col(0, 0, 1, A, -1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      zunhr_col(0, 0, 1, A, 0, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      zunhr_col(2, 0, 1, A, 1, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      infoc.INFOT = 7;
      zunhr_col(0, 0, 1, A, 1, T, -1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      zunhr_col(0, 0, 1, A, 1, T, 0, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      zunhr_col(4, 3, 2, A, 4, T, 1, D, INFO );
      chkxer('ZUNHR_COL', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
