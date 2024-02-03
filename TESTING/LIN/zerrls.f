      SUBROUTINE ZERRLS( PATH, NUNIT )

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
      int                INFO, IRNK;
      double             RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             RW( NMAX ), S( NMAX );
      COMPLEX*16         A( NMAX, NMAX ), B( NMAX, NMAX ), W( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGELS, ZGELSD, ZGELSS, ZGELST, ZGELSY, ZGETSLS
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

      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
      A( 1, 2 ) = ( 2.0D+0, 0.0D+0 )
      A( 2, 2 ) = ( 3.0D+0, 0.0D+0 )
      A( 2, 1 ) = ( 4.0D+0, 0.0D+0 )
      OK = .TRUE.
      WRITE( NOUT, FMT = * )

      // Test error exits for the least squares driver routines.

      if ( LSAMEN( 2, C2, 'LS' ) ) {

         // ZGELS

         SRNAMT = 'ZGELS '
         INFOT = 1
         zgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('ZGELS', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', INFOT, NOUT, LERR, OK );

         // ZGELST

         SRNAMT = 'ZGELST'
         INFOT = 1
         zgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', INFOT, NOUT, LERR, OK );

         // ZGETSLS

         SRNAMT = 'ZGETSLS'
         INFOT = 1
         zgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('ZGETSLS', INFOT, NOUT, LERR, OK );

         // ZGELSS

         SRNAMT = 'ZGELSS'
         INFOT = 1
         zgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW, INFO );
         chkxer('ZGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW, INFO );
         chkxer('ZGELSS', INFOT, NOUT, LERR, OK );

         // ZGELSY

         SRNAMT = 'ZGELSY'
         INFOT = 1
         zgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zgelsy(0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSY', INFOT, NOUT, LERR, OK );

         // ZGELSD

         SRNAMT = 'ZGELSD'
         INFOT = 1
         zgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP, INFO );
         chkxer('ZGELSD', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRLS

      }
