      SUBROUTINE SERRVX( PATH, NUNIT )

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
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J;
      REAL               RCOND
      // ..
      // .. Local Arrays ..
      int                IP( NMAX ), IW( NMAX );
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), E( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), W( 2*NMAX ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SGBSV, SGBSVX, SGESV, SGESVX, SGTSV, SGTSVX, SPBSV, SPBSVX, SPOSV, SPOSVX, SPPSV, SPPSVX, SPTSV, SPTSVX, SSPSV, SSPSVX, SSYSV, SSYSV_AA, SSYSV_RK, SSYSV_ROOK, SSYSVX, SSYSV_AA_2STAGE
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
            AF( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
         B( J ) = 0.E+0
         E( J ) = 0.E+0
         R1( J ) = 0.E+0
         R2( J ) = 0.E+0
         W( J ) = 0.E+0
         X( J ) = 0.E+0
         C( J ) = 0.E+0
         R( J ) = 0.E+0
         IP( J ) = J
   20 CONTINUE
      EQ = ' '
      OK = .TRUE.

      if ( LSAMEN( 2, C2, 'GE' ) ) {

         // SGESV

         SRNAMT = 'SGESV '
         INFOT = 1
         sgesv(-1, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgesv(0, -1, A, 1, IP, B, 1, INFO );
         chkxer('SGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgesv(2, 1, A, 1, IP, B, 2, INFO );
         chkxer('SGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         sgesv(2, 1, A, 2, IP, B, 1, INFO );
         chkxer('SGESV ', INFOT, NOUT, LERR, OK );

         // SGESVX

         SRNAMT = 'SGESVX'
         INFOT = 1
         sgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = '/'
         sgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ = 'R'
         sgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         EQ = 'C'
         sgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         sgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         sgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGESVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // SGBSV

         SRNAMT = 'SGBSV '
         INFOT = 1
         sgbsv(-1, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgbsv(1, -1, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sgbsv(1, 0, -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgbsv(0, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('SGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sgbsv(1, 1, 1, 0, A, 3, IP, B, 1, INFO );
         chkxer('SGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sgbsv(2, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGBSV ', INFOT, NOUT, LERR, OK );

         // SGBSVX

         SRNAMT = 'SGBSVX'
         INFOT = 1
         sgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         sgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         sgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         EQ = '/'
         sgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         EQ = 'R'
         sgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         EQ = 'C'
         sgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         sgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18
         sgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGBSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // SGTSV

         SRNAMT = 'SGTSV '
         INFOT = 1
         sgtsv(-1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('SGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgtsv(0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('SGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         sgtsv(2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('SGTSV ', INFOT, NOUT, LERR, OK );

         // SGTSVX

         SRNAMT = 'SGTSVX'
         INFOT = 1
         sgtsvx('/', 'N', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sgtsvx('N', '/', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sgtsvx('N', 'N', -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sgtsvx('N', 'N', 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         sgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         sgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SGTSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // SPOSV

         SRNAMT = 'SPOSV '
         INFOT = 1
         sposv('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('SPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sposv('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('SPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sposv('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('SPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         sposv('U', 2, 0, A, 1, B, 2, INFO );
         chkxer('SPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         sposv('U', 2, 0, A, 2, B, 1, INFO );
         chkxer('SPOSV ', INFOT, NOUT, LERR, OK );

         // SPOSVX

         SRNAMT = 'SPOSVX'
         INFOT = 1
         sposvx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sposvx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sposvx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sposvx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sposvx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         sposvx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         EQ = '/'
         sposvx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = 'Y'
         sposvx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         sposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         sposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPOSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // SPPSV

         SRNAMT = 'SPPSV '
         INFOT = 1
         sppsv('/', 0, 0, A, B, 1, INFO );
         chkxer('SPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sppsv('U', -1, 0, A, B, 1, INFO );
         chkxer('SPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sppsv('U', 0, -1, A, B, 1, INFO );
         chkxer('SPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sppsv('U', 2, 0, A, B, 1, INFO );
         chkxer('SPPSV ', INFOT, NOUT, LERR, OK );

         // SPPSVX

         SRNAMT = 'SPPSVX'
         INFOT = 1
         sppsvx('/', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sppsvx('N', '/', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sppsvx('N', 'U', -1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sppsvx('N', 'U', 0, -1, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         EQ = '/'
         sppsvx('F', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         EQ = 'Y'
         sppsvx('F', 'U', 1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         sppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         sppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPPSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // SPBSV

         SRNAMT = 'SPBSV '
         INFOT = 1
         spbsv('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('SPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         spbsv('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('SPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         spbsv('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('SPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         spbsv('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('SPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         spbsv('U', 1, 1, 0, A, 1, B, 2, INFO );
         chkxer('SPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         spbsv('U', 2, 0, 0, A, 1, B, 1, INFO );
         chkxer('SPBSV ', INFOT, NOUT, LERR, OK );

         // SPBSVX

         SRNAMT = 'SPBSVX'
         INFOT = 1
         spbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         spbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         spbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         spbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         spbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         spbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         spbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = '/'
         spbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ = 'Y'
         spbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         spbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         spbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SPBSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // SPTSV

         SRNAMT = 'SPTSV '
         INFOT = 1
         sptsv(-1, 0, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('SPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sptsv(0, -1, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('SPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         sptsv(2, 0, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('SPTSV ', INFOT, NOUT, LERR, OK );

         // SPTSVX

         SRNAMT = 'SPTSVX'
         INFOT = 1
         sptsvx('/', 0, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('SPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sptsvx('N', -1, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('SPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sptsvx('N', 0, -1, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('SPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sptsvx('N', 2, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 2, RCOND, R1, R2, W, INFO );
         chkxer('SPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         sptsvx('N', 2, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 2, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('SPTSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // SSYSV

         SRNAMT = 'SSYSV '
         INFOT = 1
         ssysv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ssysv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ssysv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ssysv('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ssysv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ssysv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ssysv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('SSYSV ', INFOT, NOUT, LERR, OK );

         // SSYSVX

         SRNAMT = 'SSYSVX'
         INFOT = 1
         ssysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ssysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ssysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ssysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ssysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ssysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ssysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         ssysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18
         ssysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, IW, INFO );
         chkxer('SSYSVX', INFOT, NOUT, LERR, OK );


      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // SSYSV_ROOK

         SRNAMT = 'SSYSV_ROOK'
         INFOT = 1
         ssysv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ssysv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ssysv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ssysv_rook('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ssysv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ssysv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ssysv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('SSYSV_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // SSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         SRNAMT = 'SSYSV_RK'
         INFOT = 1
         ssysv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ssysv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ssysv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ssysv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 9
         ssysv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ssysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ssysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('SSYSV_RK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SA' ) ) {

         // SSYSV_AASEN

         SRNAMT = 'SSYSV_AA'
         INFOT = 1
         ssysv_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ssysv_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ssysv_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ssysv_aa('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('SSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ssysv_aa('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ssysv_aa('U', 3, 1, A, 3, IP, B, 3, W, 6, INFO );
         chkxer('SSYSV_AA', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // SSYSV_AASEN_2STAGE

         SRNAMT = 'SSYSV_AA_2STAGE'
         INFOT = 1
         ssysv_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ssysv_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ssysv_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ssysv_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7
         ssysv_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 2, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ssysv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13
         ssysv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 2, W, 1, INFO );
         chkxer('SSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // SSPSV

         SRNAMT = 'SSPSV '
         INFOT = 1
         sspsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('SSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sspsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('SSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sspsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('SSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         sspsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('SSPSV ', INFOT, NOUT, LERR, OK );

         // SSPSVX

         SRNAMT = 'SSPSVX'
         INFOT = 1
         sspsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         sspsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         sspsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         sspsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         sspsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('SSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         sspsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('SSPSVX', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

 9999 FORMAT( 1X, A3, ' drivers passed the tests of the error exits' )
 9998 FORMAT( ' *** ', A3, ' drivers failed the tests of the error ', 'exits ***' )

      RETURN

      // End of SERRVX

      }
