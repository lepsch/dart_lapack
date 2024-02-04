      void cerrtz(PATH, NUNIT ) {

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
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                INFO;
      // ..
      // .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), TAU( NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CTZRZF
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
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      C2 = PATH( 2: 3 );
      A[1, 1] = CMPLX( 1.0, -1.0 );
      A[1, 2] = CMPLX( 2.0, -2.0 );
      A[2, 2] = CMPLX( 3.0, -3.0 );
      A[2, 1] = CMPLX( 4.0, -4.0 );
      W[1] = CMPLX( 0.0, 0.0 );
      W[2] = CMPLX( 0.0, 0.0 );
      OK = true;

      // Test error exits for the trapezoidal routines.

      WRITE( NOUT, FMT = * );
      if ( LSAMEN( 2, C2, 'TZ' ) ) {

         // CTZRZF

         SRNAMT = 'CTZRZF';
         INFOT = 1;
         ctzrzf(-1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctzrzf(1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctzrzf(2, 2, A, 1, TAU, W, 1, INFO );
         chkxer('CTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ctzrzf(2, 2, A, 2, TAU, W, 0, INFO );
         chkxer('CTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ctzrzf(2, 3, A, 2, TAU, W, 1, INFO );
         chkxer('CTZRZF', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
