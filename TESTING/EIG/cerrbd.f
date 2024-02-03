      SUBROUTINE CERRBD( PATH, NUNIT )

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
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J, NT;
      // ..
      // .. Local Arrays ..
      REAL               D( NMAX ), E( NMAX ), RW( 4*NMAX )
      COMPLEX            A( NMAX, NMAX ), TP( NMAX ), TQ( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, CBDSQR, CGEBD2, CGEBRD, CUNGBR, CUNMBR
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
         } // 10
      } // 20
      OK = true;
      NT = 0

      // Test error exits of the SVD routines.

      if ( LSAMEN( 2, C2, 'BD' ) ) {

         // CGEBRD

         SRNAMT = 'CGEBRD'
         INFOT = 1
         cgebrd(-1, 0, A, 1, D, E, TQ, TP, W, 1, INFO );
         chkxer('CGEBRD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgebrd(0, -1, A, 1, D, E, TQ, TP, W, 1, INFO );
         chkxer('CGEBRD', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cgebrd(2, 1, A, 1, D, E, TQ, TP, W, 2, INFO );
         chkxer('CGEBRD', INFOT, NOUT, LERR, OK );
         INFOT = 10
         cgebrd(2, 1, A, 2, D, E, TQ, TP, W, 1, INFO );
         chkxer('CGEBRD', INFOT, NOUT, LERR, OK );
         NT = NT + 4

         // CGEBD2

         SRNAMT = 'CGEBD2'
         INFOT = 1
         cgebd2(-1, 0, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('CGEBD2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgebd2(0, -1, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('CGEBD2', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cgebd2(2, 1, A, 1, D, E, TQ, TP, W, INFO );
         chkxer('CGEBD2', INFOT, NOUT, LERR, OK );
         NT = NT + 3

         // CUNGBR

         SRNAMT = 'CUNGBR'
         INFOT = 1
         cungbr('/', 0, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cungbr('Q', -1, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cungbr('Q', 0, -1, 0, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cungbr('Q', 0, 1, 0, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cungbr('Q', 1, 0, 1, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cungbr('P', 1, 0, 0, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cungbr('P', 0, 1, 1, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cungbr('Q', 0, 0, -1, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         cungbr('Q', 2, 1, 1, A, 1, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         INFOT = 9
         cungbr('Q', 2, 2, 1, A, 2, TQ, W, 1, INFO );
         chkxer('CUNGBR', INFOT, NOUT, LERR, OK );
         NT = NT + 10

         // CUNMBR

         SRNAMT = 'CUNMBR'
         INFOT = 1
         cunmbr('/', 'L', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cunmbr('Q', '/', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cunmbr('Q', 'L', '/', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cunmbr('Q', 'L', 'C', -1, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cunmbr('Q', 'L', 'C', 0, -1, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         cunmbr('Q', 'L', 'C', 0, 0, -1, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunmbr('Q', 'L', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunmbr('Q', 'R', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunmbr('P', 'L', 'C', 2, 0, 2, A, 1, TQ, U, 2, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunmbr('P', 'R', 'C', 0, 2, 2, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 11
         cunmbr('Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 1, W, 1, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         cunmbr('Q', 'L', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 0, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         cunmbr('Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 0, INFO );
         chkxer('CUNMBR', INFOT, NOUT, LERR, OK );
         NT = NT + 13

         // CBDSQR

         SRNAMT = 'CBDSQR'
         INFOT = 1
         cbdsqr('/', 0, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cbdsqr('U', -1, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cbdsqr('U', 0, -1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cbdsqr('U', 0, 0, -1, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cbdsqr('U', 0, 0, 0, -1, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 9
         cbdsqr('U', 2, 1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 11
         cbdsqr('U', 0, 0, 2, 0, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         cbdsqr('U', 2, 0, 0, 1, D, E, V, 1, U, 1, A, 1, RW, INFO );
         chkxer('CBDSQR', INFOT, NOUT, LERR, OK );
         NT = NT + 8
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' )

      RETURN

      // End of CERRBD

      }
