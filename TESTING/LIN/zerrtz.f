      SUBROUTINE ZERRTZ( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                INFO;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), TAU( NMAX ), W( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZTZRZF
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
      // INTRINSIC DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = DCMPLX( 1.D+0, -1.D+0 )
      A( 1, 2 ) = DCMPLX( 2.D+0, -2.D+0 )
      A( 2, 2 ) = DCMPLX( 3.D+0, -3.D+0 )
      A( 2, 1 ) = DCMPLX( 4.D+0, -4.D+0 )
      W( 1 ) = DCMPLX( 0.D+0, 0.D+0 )
      W( 2 ) = DCMPLX( 0.D+0, 0.D+0 )
      OK = .TRUE.

      // Test error exits for the trapezoidal routines.
      WRITE( NOUT, FMT = * )
      if ( LSAMEN( 2, C2, 'TZ' ) ) {


         // ZTZRZF

         SRNAMT = 'ZTZRZF'
         INFOT = 1
         ztzrzf(-1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ztzrzf(1, 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ztzrzf(2, 2, A, 1, TAU, W, 1, INFO );
         chkxer('ZTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 7
         ztzrzf(2, 2, A, 2, TAU, W, 0, INFO );
         chkxer('ZTZRZF', INFOT, NOUT, LERR, OK );
         INFOT = 7
         ztzrzf(2, 3, A, 2, TAU, W, 1, INFO );
         chkxer('ZTZRZF', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRTZ

      }
