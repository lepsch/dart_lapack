      SUBROUTINE SERRBD( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX, LW;
      const              NMAX = 4, LW = NMAX ;
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J, NS, NT;
      // ..
      // .. Local Arrays ..
      int                IQ( NMAX, NMAX ), IW( NMAX );
      REAL               A( NMAX, NMAX ), D( NMAX ), E( NMAX ), Q( NMAX, NMAX ), S( NMAX ), TP( NMAX ), TQ( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SBDSDC, SBDSQR, SBDSVDX, SGEBD2, SGEBRD, SORGBR, SORMBR
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
   20 CONTINUE
      OK = .TRUE.
      NT = 0

      // Test error exits of the SVD routines.

      if ( LSAMEN( 2, C2, 'BD' ) ) {

         // SGEBRD

         SRNAMT = 'SGEBRD'
         INFOT = 1
         sgebrd(-1, 0, A, 1, D, E, TQ, TP, W, 1, INFO );
         chkxer('SGEBRD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgebrd(0, -1, A, 1, D, E, TQ, TP, W, 1, INFO );
         chkxer('SGEBRD', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgebrd(2, 1, A, 1, D, E, TQ, TP, W, 2, INFO );
         chkxer('SGEBRD', INFOT, NOUT, LERR, OK );
         INFOT = 10
         sgebrd(2, 1, A, 2, D, E, TQ, TP, W, 1, INFO );
         chkxer('SGEBRD', INFOT, NOUT, LERR, OK );
         NT = NT + 4

         // SGEBD2

         SRNAMT = 'SGEBD2'
         INFOT = 1
         sgebd2(-1, 0, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('SGEBD2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgebd2(0, -1, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('SGEBD2', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgebd2(2, 1, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('SGEBD2', INFOT, NOUT, LERR, OK );
         NT = NT + 3

         // SORGBR

         SRNAMT = 'SORGBR'
         INFOT = 1
         sorgbr('/', 0, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sorgbr('Q', -1, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sorgbr('Q', 0, -1, 0, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sorgbr('Q', 0, 1, 0, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sorgbr('Q', 1, 0, 1, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sorgbr('P', 1, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sorgbr('P', 0, 1, 1, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sorgbr('Q', 0, 0, -1, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sorgbr('Q', 2, 1, 1, A, 1, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sorgbr('Q', 2, 2, 1, A, 2, TQ, W, 1, INFO );
         chkxer('SORGBR', INFOT, NOUT, LERR, OK );
         NT = NT + 10

         // SORMBR

         SRNAMT = 'SORMBR'
         INFOT = 1
         sormbr('/', 'L', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sormbr('Q', '/', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sormbr('Q', 'L', '/', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sormbr('Q', 'L', 'T', -1, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         sormbr('Q', 'L', 'T', 0, -1, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sormbr('Q', 'L', 'T', 0, 0, -1, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sormbr('Q', 'L', 'T', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sormbr('Q', 'R', 'T', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sormbr('P', 'L', 'T', 2, 0, 2, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sormbr('P', 'R', 'T', 0, 2, 2, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 11
         sormbr('Q', 'R', 'T', 2, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         sormbr('Q', 'L', 'T', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         sormbr('Q', 'R', 'T', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('SORMBR', INFOT, NOUT, LERR, OK );
         NT = NT + 13

         // SBDSQR

         SRNAMT = 'SBDSQR'
         INFOT = 1
         sbdsqr('/', 0, 0, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sbdsqr('U', -1, 0, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sbdsqr('U', 0, -1, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sbdsqr('U', 0, 0, -1, 0, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         sbdsqr('U', 0, 0, 0, -1, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sbdsqr('U', 2, 1, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 11
         sbdsqr('U', 0, 0, 2, 0, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         sbdsqr('U', 2, 0, 0, 1, D, E, V, 1, U, 1, A, 1, W, INFO );
         chkxer('SBDSQR', INFOT, NOUT, LERR, OK );
         NT = NT + 8

         // SBDSDC

         SRNAMT = 'SBDSDC'
         INFOT = 1
         sbdsdc('/', 'N', 0, D, E, U, 1, V, 1, Q, IQ, W, IW, INFO );
         chkxer('SBDSDC', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sbdsdc('U', '/', 0, D, E, U, 1, V, 1, Q, IQ, W, IW, INFO );
         chkxer('SBDSDC', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sbdsdc('U', 'N', -1, D, E, U, 1, V, 1, Q, IQ, W, IW, INFO );
         chkxer('SBDSDC', INFOT, NOUT, LERR, OK );
         INFOT = 7
         sbdsdc('U', 'I', 2, D, E, U, 1, V, 1, Q, IQ, W, IW, INFO );
         chkxer('SBDSDC', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sbdsdc('U', 'I', 2, D, E, U, 2, V, 1, Q, IQ, W, IW, INFO );
         chkxer('SBDSDC', INFOT, NOUT, LERR, OK );
         NT = NT + 5

         // SBDSVDX

         SRNAMT = 'SBDSVDX'
         INFOT = 1
         sbdsvdx('X', 'N', 'A', 1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sbdsvdx('U', 'X', 'A', 1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sbdsvdx('U', 'V', 'X', 1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sbdsvdx('U', 'V', 'A', -1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         sbdsvdx('U', 'V', 'V', 2, D, E, -ONE, ZERO, 0, 0, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sbdsvdx('U', 'V', 'V', 2, D, E, ONE, ZERO, 0, 0, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sbdsvdx('L', 'V', 'I', 2, D, E, ZERO, ZERO, 0, 2, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sbdsvdx('L', 'V', 'I', 4, D, E, ZERO, ZERO, 5, 2, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         sbdsvdx('L', 'V', 'I', 4, D, E, ZERO, ZERO, 3, 2, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         sbdsvdx('L', 'V', 'I', 4, D, E, ZERO, ZERO, 3, 5, NS, S, Q, 1, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         sbdsvdx('L', 'V', 'A', 4, D, E, ZERO, ZERO, 0, 0, NS, S, Q, 0, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         sbdsvdx('L', 'V', 'A', 4, D, E, ZERO, ZERO, 0, 0, NS, S, Q, 2, W, IW, INFO);
         chkxer('SBDSVDX', INFOT, NOUT, LERR, OK );
         NT = NT + 12
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits', ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' )

      RETURN

      // End of SERRBD

      }
