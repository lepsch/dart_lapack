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
      PARAMETER          ( NMAX = 4, LW = NMAX )
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

      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
   20 CONTINUE
      OK = .TRUE.
      NT = 0

      // Test error exits of the SVD routines.

      IF( LSAMEN( 2, C2, 'BD' ) ) THEN

         // CGEBRD

         SRNAMT = 'CGEBRD'
         INFOT = 1
         CALL CGEBRD( -1, 0, A, 1, D, E, TQ, TP, W, 1, INFO )
         CALL CHKXER( 'CGEBRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEBRD( 0, -1, A, 1, D, E, TQ, TP, W, 1, INFO )
         CALL CHKXER( 'CGEBRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEBRD( 2, 1, A, 1, D, E, TQ, TP, W, 2, INFO )
         CALL CHKXER( 'CGEBRD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGEBRD( 2, 1, A, 2, D, E, TQ, TP, W, 1, INFO )
         CALL CHKXER( 'CGEBRD', INFOT, NOUT, LERR, OK )
         NT = NT + 4

         // CGEBD2

         SRNAMT = 'CGEBD2'
         INFOT = 1
         CALL CGEBD2( -1, 0, A, 1, D, E, TQ, TP, W, INFO )
         CALL CHKXER( 'CGEBD2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEBD2( 0, -1, A, 1, D, E, TQ, TP, W, INFO )
         CALL CHKXER( 'CGEBD2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEBD2( 2, 1, A, 1, D, E, TQ, TP, W, INFO )
         CALL CHKXER( 'CGEBD2', INFOT, NOUT, LERR, OK )
         NT = NT + 3

         // CUNGBR

         SRNAMT = 'CUNGBR'
         INFOT = 1
         CALL CUNGBR( '/', 0, 0, 0, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CUNGBR( 'Q', -1, 0, 0, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNGBR( 'Q', 0, -1, 0, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNGBR( 'Q', 0, 1, 0, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNGBR( 'Q', 1, 0, 1, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNGBR( 'P', 1, 0, 0, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNGBR( 'P', 0, 1, 1, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CUNGBR( 'Q', 0, 0, -1, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CUNGBR( 'Q', 2, 1, 1, A, 1, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CUNGBR( 'Q', 2, 2, 1, A, 2, TQ, W, 1, INFO )
         CALL CHKXER( 'CUNGBR', INFOT, NOUT, LERR, OK )
         NT = NT + 10

         // CUNMBR

         SRNAMT = 'CUNMBR'
         INFOT = 1
         CALL CUNMBR( '/', 'L', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CUNMBR( 'Q', '/', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNMBR( 'Q', 'L', '/', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CUNMBR( 'Q', 'L', 'C', -1, 0, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CUNMBR( 'Q', 'L', 'C', 0, -1, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CUNMBR( 'Q', 'L', 'C', 0, 0, -1, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CUNMBR( 'Q', 'L', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CUNMBR( 'Q', 'R', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CUNMBR( 'P', 'L', 'C', 2, 0, 2, A, 1, TQ, U, 2, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CUNMBR( 'P', 'R', 'C', 0, 2, 2, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CUNMBR( 'Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CUNMBR( 'Q', 'L', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 0, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CUNMBR( 'Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 0, INFO )
         CALL CHKXER( 'CUNMBR', INFOT, NOUT, LERR, OK )
         NT = NT + 13

         // CBDSQR

         SRNAMT = 'CBDSQR'
         INFOT = 1
         CALL CBDSQR( '/', 0, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CBDSQR( 'U', -1, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CBDSQR( 'U', 0, -1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CBDSQR( 'U', 0, 0, -1, 0, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CBDSQR( 'U', 0, 0, 0, -1, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CBDSQR( 'U', 2, 1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CBDSQR( 'U', 0, 0, 2, 0, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CBDSQR( 'U', 2, 0, 0, 1, D, E, V, 1, U, 1, A, 1, RW, INFO )
         CALL CHKXER( 'CBDSQR', INFOT, NOUT, LERR, OK )
         NT = NT + 8
      END IF

      // Print a summary line.

      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (',
     $        I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $        'exits ***' )

      RETURN

      // End of CERRBD

      END
