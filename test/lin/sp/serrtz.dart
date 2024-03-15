      void serrtz(final int PATH, final int NUNIT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      String             C2;
      int                INFO;
      double               A( NMAX, NMAX ), TAU( NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, STZRZF
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A[1][1] = 1.0;
      A[1][2] = 2.0;
      A[2][2] = 3.0;
      A[2][1] = 4.0;
      W[1] = 0.0;
      W[2] = 0.0;
      OK = true;

      if ( lsamen( 2, C2, 'TZ' ) ) {

         // Test error exits for the trapezoidal routines.

         // STZRZF

        srnamc.SRNAMT = 'STZRZF';
         INFOT = 1;
         stzrzf(-1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('STZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stzrzf(1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('STZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stzrzf(2, 2, A, 1, TAU, W, 1, INFO );
         chkxer('STZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         stzrzf(2, 2, A, 2, TAU, W, 0, INFO );
         chkxer('STZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         stzrzf(2, 3, A, 2, TAU, W, 1, INFO );
         chkxer('STZRZF', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }