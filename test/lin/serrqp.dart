      void serrqp(PATH, NUNIT ) {

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
      double               A( NMAX, NMAX ), TAU( NMAX ), W( 3*NMAX+1 );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGEQP3
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
      LW = 3*NMAX + 1;
      A[1, 1] = 1.0;
      A[1, 2] = 2.0;
      A[2, 2] = 3.0;
      A[2, 1] = 4.0;
      OK = true;

      if ( LSAMEN( 2, C2, 'QP' ) ) {

         // Test error exits for QR factorization with pivoting

         // SGEQP3

         SRNAMT = 'SGEQP3';
         INFOT = 1;
         sgeqp3(-1, 0, A, 1, IP, TAU, W, LW, INFO );
         chkxer('SGEQP3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgeqp3(1, -1, A, 1, IP, TAU, W, LW, INFO );
         chkxer('SGEQP3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgeqp3(2, 3, A, 1, IP, TAU, W, LW, INFO );
         chkxer('SGEQP3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgeqp3(2, 2, A, 2, IP, TAU, W, LW-10, INFO );
         chkxer('SGEQP3', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }