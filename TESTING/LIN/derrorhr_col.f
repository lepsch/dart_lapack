      SUBROUTINE DERRORHR_COL( PATH, NUNIT );
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String   (LEN=3)   PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), D(NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DORHR_COL
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String   (LEN=32)  SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / DBLE( I+J );
            T( I, J ) = 1.0 / DBLE( I+J );
         }
         D( J ) = 0.0;
      }
      OK = true;

      // Error exits for Householder reconstruction

      // DORHR_COL

      SRNAMT = 'DORHR_COL';

      INFOT = 1;
      dorhr_col(-1, 0, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      INFOT = 2;
      dorhr_col(0, -1, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );
      dorhr_col(1, 2, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      INFOT = 3;
      dorhr_col(0, 0, -1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      dorhr_col(0, 0, 0, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      INFOT = 5;
      dorhr_col(0, 0, 1, A, -1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      dorhr_col(0, 0, 1, A, 0, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      dorhr_col(2, 0, 1, A, 1, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      INFOT = 7;
      dorhr_col(0, 0, 1, A, 1, T, -1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      dorhr_col(0, 0, 1, A, 1, T, 0, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      dorhr_col(4, 3, 2, A, 4, T, 1, D, INFO );
      chkxer('DORHR_COL', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN;

      // End of DERRORHR_COL

      }
