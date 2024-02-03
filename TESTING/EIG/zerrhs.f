      SUBROUTINE ZERRHS( PATH, NUNIT )

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
      double             RW( NMAX ), S( NMAX );
      COMPLEX*16         A( NMAX, NMAX ), C( NMAX, NMAX ), TAU( NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), W( LW ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEIN, ZHSEQR, ZUNGHR, ZUNMHR, ZTREVC, ZTREVC3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
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

      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
         SEL( J ) = .TRUE.
   20 CONTINUE
      OK = .TRUE.
      NT = 0

      // Test error exits of the nonsymmetric eigenvalue routines.

      if ( LSAMEN( 2, C2, 'HS' ) ) {

         // ZGEBAL

         SRNAMT = 'ZGEBAL'
         INFOT = 1
         CALL ZGEBAL( '/', 0, A, 1, ILO, IHI, S, INFO )
         CALL CHKXER( 'ZGEBAL', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEBAL( 'N', -1, A, 1, ILO, IHI, S, INFO )
         CALL CHKXER( 'ZGEBAL', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGEBAL( 'N', 2, A, 1, ILO, IHI, S, INFO )
         CALL CHKXER( 'ZGEBAL', INFOT, NOUT, LERR, OK )
         NT = NT + 3

         // ZGEBAK

         SRNAMT = 'ZGEBAK'
         INFOT = 1
         CALL ZGEBAK( '/', 'R', 0, 1, 0, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEBAK( 'N', '/', 0, 1, 0, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEBAK( 'N', 'R', -1, 1, 0, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGEBAK( 'N', 'R', 0, 0, 0, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGEBAK( 'N', 'R', 0, 2, 0, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGEBAK( 'N', 'R', 2, 2, 1, S, 0, A, 2, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGEBAK( 'N', 'R', 0, 1, 1, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGEBAK( 'N', 'R', 0, 1, 0, S, -1, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZGEBAK( 'N', 'R', 2, 1, 2, S, 0, A, 1, INFO )
         CALL CHKXER( 'ZGEBAK', INFOT, NOUT, LERR, OK )
         NT = NT + 9

         // ZGEHRD

         SRNAMT = 'ZGEHRD'
         INFOT = 1
         CALL ZGEHRD( -1, 1, 1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEHRD( 0, 0, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEHRD( 0, 2, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEHRD( 1, 1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEHRD( 0, 1, 1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGEHRD( 2, 1, 1, A, 1, TAU, W, 2, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGEHRD( 2, 1, 2, A, 2, TAU, W, 1, INFO )
         CALL CHKXER( 'ZGEHRD', INFOT, NOUT, LERR, OK )
         NT = NT + 7

         // ZGEHD2

         SRNAMT = 'ZGEHD2'
         INFOT = 1
         CALL ZGEHD2( -1, 1, 1, A, 1, TAU, W, INFO )
         CALL CHKXER( 'ZGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEHD2( 0, 0, 0, A, 1, TAU, W, INFO )
         CALL CHKXER( 'ZGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEHD2( 0, 2, 0, A, 1, TAU, W, INFO )
         CALL CHKXER( 'ZGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEHD2( 1, 1, 0, A, 1, TAU, W, INFO )
         CALL CHKXER( 'ZGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGEHD2( 0, 1, 1, A, 1, TAU, W, INFO )
         CALL CHKXER( 'ZGEHD2', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGEHD2( 2, 1, 1, A, 1, TAU, W, INFO )
         CALL CHKXER( 'ZGEHD2', INFOT, NOUT, LERR, OK )
         NT = NT + 6

         // ZUNGHR

         SRNAMT = 'ZUNGHR'
         INFOT = 1
         CALL ZUNGHR( -1, 1, 1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZUNGHR( 0, 0, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZUNGHR( 0, 2, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZUNGHR( 1, 1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZUNGHR( 0, 1, 1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZUNGHR( 2, 1, 1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZUNGHR( 3, 1, 3, A, 3, TAU, W, 1, INFO )
         CALL CHKXER( 'ZUNGHR', INFOT, NOUT, LERR, OK )
         NT = NT + 7

         // ZUNMHR

         SRNAMT = 'ZUNMHR'
         INFOT = 1
         CALL ZUNMHR( '/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZUNMHR( 'L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZUNMHR( 'L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZUNMHR( 'L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZUNMHR( 'L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZUNMHR( 'L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZUNMHR( 'L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZUNMHR( 'R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZUNMHR( 'L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZUNMHR( 'L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZUNMHR( 'R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZUNMHR( 'L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZUNMHR( 'R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZUNMHR( 'L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZUNMHR( 'L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZUNMHR( 'R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'ZUNMHR', INFOT, NOUT, LERR, OK )
         NT = NT + 16

         // ZHSEQR

         SRNAMT = 'ZHSEQR'
         INFOT = 1
         CALL ZHSEQR( '/', 'N', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHSEQR( 'E', '/', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHSEQR( 'E', 'N', -1, 1, 0, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZHSEQR( 'E', 'N', 0, 0, 0, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZHSEQR( 'E', 'N', 0, 2, 0, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHSEQR( 'E', 'N', 1, 1, 0, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHSEQR( 'E', 'N', 1, 1, 2, A, 1, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZHSEQR( 'E', 'N', 2, 1, 2, A, 1, X, C, 2, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHSEQR( 'E', 'V', 2, 1, 2, A, 2, X, C, 1, W, 1, INFO )
         CALL CHKXER( 'ZHSEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 9

         // ZHSEIN

         SRNAMT = 'ZHSEIN'
         INFOT = 1
         CALL ZHSEIN( '/', 'N', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHSEIN( 'R', '/', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHSEIN( 'R', 'N', '/', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHSEIN( 'R', 'N', 'N', SEL, -1, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZHSEIN( 'R', 'N', 'N', SEL, 2, A, 1, X, VL, 1, VR, 2, 4, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHSEIN( 'L', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZHSEIN( 'R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZHSEIN( 'R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 2, 1, M, W, RW, IFAILL, IFAILR, INFO )
         CALL CHKXER( 'ZHSEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 8

         // ZTREVC

         SRNAMT = 'ZTREVC'
         INFOT = 1
         CALL ZTREVC( '/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZTREVC( 'L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZTREVC( 'L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZTREVC( 'L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZTREVC( 'L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZTREVC( 'R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZTREVC( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, RW, INFO )
         CALL CHKXER( 'ZTREVC', INFOT, NOUT, LERR, OK )
         NT = NT + 7

         // ZTREVC3

         SRNAMT = 'ZTREVC3'
         INFOT = 1
         CALL ZTREVC3( '/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZTREVC3( 'L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZTREVC3( 'L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZTREVC3( 'L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, RW, 2, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZTREVC3( 'R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, RW, 2, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL ZTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, RW, 2, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL ZTREVC3( 'L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, LW, RW, 1, INFO )
         CALL CHKXER( 'ZTREVC3', INFOT, NOUT, LERR, OK )
         NT = NT + 9
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits',
     $      ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )

      RETURN

      // End of ZERRHS

      }
