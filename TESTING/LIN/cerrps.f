      SUBROUTINE CERRPS( PATH, NUNIT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NUNIT;
      String             PATH;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, RANK;
      // ..
      // .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX );
      REAL               RWORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CPSTF2, CPSTRF
      // ..
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 110
         for (I = 1; I <= NMAX; I++) { // 100
            A( I, J ) = 1.0 / REAL( I+J );

         } // 100
         PIV( J ) = J;
         RWORK( J ) = 0.;
         RWORK( NMAX+J ) = 0.;

      } // 110
      OK = true;


         // Test error exits of the routines that use the Cholesky
         // decomposition of an Hermitian positive semidefinite matrix.

         // CPSTRF

      SRNAMT = 'CPSTRF';
      INFOT = 1;
      cpstrf('/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('CPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cpstrf('U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('CPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cpstrf('U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('CPSTRF', INFOT, NOUT, LERR, OK );

         // CPSTF2

      SRNAMT = 'CPSTF2';
      INFOT = 1;
      cpstf2('/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('CPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cpstf2('U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('CPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cpstf2('U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('CPSTF2', INFOT, NOUT, LERR, OK );


      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of CERRPS

      }
