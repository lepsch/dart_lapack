      SUBROUTINE CERRHS( PATH, NUNIT )

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
      const              NMAX = 3, LW = NMAX*NMAX ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, M, NT;
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      int                IFAILL( NMAX ), IFAILR( NMAX );
      REAL               RW( NMAX ), S( NMAX )
      COMPLEX            A( NMAX, NMAX ), C( NMAX, NMAX ), TAU( NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), W( LW ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, CGEBAK, CGEBAL, CGEHRD, CHSEIN, CHSEQR, CUNGHR, CUNMHR, CTREVC, CTREVC3, CGEHD2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
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
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
         SEL( J ) = .TRUE.
   20 CONTINUE
      OK = .TRUE.
      NT = 0

      // Test error exits of the nonsymmetric eigenvalue routines.

      if ( LSAMEN( 2, C2, 'HS' ) ) {

         // CGEBAL

         SRNAMT = 'CGEBAL'
         INFOT = 1
         cgebal('/', 0, A, 1, ILO, IHI, S, INFO );
         chkxer('CGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgebal('N', -1, A, 1, ILO, IHI, S, INFO );
         chkxer('CGEBAL', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cgebal('N', 2, A, 1, ILO, IHI, S, INFO );
         chkxer('CGEBAL', INFOT, NOUT, LERR, OK );
         NT = NT + 3

         // CGEBAK

         SRNAMT = 'CGEBAK'
         INFOT = 1
         cgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 7
         cgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         INFOT = 9
         cgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO );
         chkxer('CGEBAK', INFOT, NOUT, LERR, OK );
         NT = NT + 9

         // CGEHRD

         SRNAMT = 'CGEHRD'
         INFOT = 1
         cgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO );
         chkxer('CGEHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 7

         // CGEHD2

         SRNAMT = 'CGEHD2'
         INFOT = 1
         cgehd2(-1, 1, 1, A, 1, TAU, W, INFO );
         chkxer('CGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgehd2(0, 0, 0, A, 1, TAU, W, INFO );
         chkxer('CGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgehd2(0, 2, 0, A, 1, TAU, W, INFO );
         chkxer('CGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgehd2(1, 1, 0, A, 1, TAU, W, INFO );
         chkxer('CGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgehd2(0, 1, 1, A, 1, TAU, W, INFO );
         chkxer('CGEHD2', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cgehd2(2, 1, 1, A, 1, TAU, W, INFO );
         chkxer('CGEHD2', INFOT, NOUT, LERR, OK );
         NT = NT + 6

         // CUNGHR

         SRNAMT = 'CUNGHR'
         INFOT = 1
         cunghr(-1, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cunghr(0, 0, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cunghr(0, 2, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cunghr(1, 1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cunghr(0, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cunghr(2, 1, 1, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunghr(3, 1, 3, A, 3, TAU, W, 1, INFO );
         chkxer('CUNGHR', INFOT, NOUT, LERR, OK );
         NT = NT + 7

         // CUNMHR

         SRNAMT = 'CUNMHR'
         INFOT = 1
         cunmhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cunmhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cunmhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cunmhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cunmhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cunmhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cunmhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         cunmhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         cunmhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         cunmhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 6
         cunmhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunmhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cunmhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 11
         cunmhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         cunmhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         INFOT = 13
         cunmhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('CUNMHR', INFOT, NOUT, LERR, OK );
         NT = NT + 16

         // CHSEQR

         SRNAMT = 'CHSEQR'
         INFOT = 1
         chseqr('/', 'N', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2
         chseqr('E', '/', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 3
         chseqr('E', 'N', -1, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         chseqr('E', 'N', 0, 0, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 4
         chseqr('E', 'N', 0, 2, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         chseqr('E', 'N', 1, 1, 0, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 5
         chseqr('E', 'N', 1, 1, 2, A, 1, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 7
         chseqr('E', 'N', 2, 1, 2, A, 1, X, C, 2, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         INFOT = 10
         chseqr('E', 'V', 2, 1, 2, A, 2, X, C, 1, W, 1, INFO );
         chkxer('CHSEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 9

         // CHSEIN

         SRNAMT = 'CHSEIN'
         INFOT = 1
         chsein('/', 'N', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 2
         chsein('R', '/', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 3
         chsein('R', 'N', '/', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 5
         chsein('R', 'N', 'N', SEL, -1, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 7
         chsein('R', 'N', 'N', SEL, 2, A, 1, X, VL, 1, VR, 2, 4, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 10
         chsein('L', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 12
         chsein('R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         INFOT = 13
         chsein('R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 2, 1, M, W, RW, IFAILL, IFAILR, INFO );
         chkxer('CHSEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 8

         // CTREVC

         SRNAMT = 'CTREVC'
         INFOT = 1
         ctrevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ctrevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ctrevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ctrevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ctrevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ctrevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ctrevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, RW, INFO );
         chkxer('CTREVC', INFOT, NOUT, LERR, OK );
         NT = NT + 7

         // CTREVC3

         SRNAMT = 'CTREVC3'
         INFOT = 1
         ctrevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ctrevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ctrevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ctrevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, RW, 2, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ctrevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ctrevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ctrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, RW, 2, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 14
         ctrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, RW, 2, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         INFOT = 16
         ctrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, LW, RW, 1, INFO );
         chkxer('CTREVC3', INFOT, NOUT, LERR, OK );
         NT = NT + 9
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

      // End of CERRHS

      }
