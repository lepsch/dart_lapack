      SUBROUTINE ZERRVX( PATH, NUNIT )

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
      REAL               ONE
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J, N_ERR_BNDS, NPARAMS;
      double             RCOND, RPVGRW, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), RF( NMAX ), RW( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGBSV, ZGBSVX, ZGESV, ZGESVX, ZGTSV, ZGTSVX, ZHESV, ZHESV_RK, ZHESV_ROOK, ZHESVX, ZHPSV, ZHPSVX, ZPBSV, ZPBSVX, ZPOSV, ZPOSVX, ZPPSV, ZPPSVX, ZPTSV, ZPTSVX, ZSPSV, ZSPSVX, ZSYSV, ZSYSV_RK, ZSYSV_ROOK, ZSYSVX, ZGESVXX, ZSYSVXX, ZPOSVXX, ZHESVXX, ZGBSVXX
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )
         } // 10
         B( J ) = 0.D0
         E( J ) = 0.D0
         R1( J ) = 0.D0
         R2( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
         C( J ) = 0.D0
         R( J ) = 0.D0
         IP( J ) = J
      } // 20
      EQ = ' '
      OK = .TRUE.

      if ( LSAMEN( 2, C2, 'GE' ) ) {

         // ZGESV

         SRNAMT = 'ZGESV '
         INFOT = 1
         zgesv(-1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgesv(0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgesv(2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgesv(2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZGESV ', INFOT, NOUT, LERR, OK );

         // ZGESVX

         SRNAMT = 'ZGESVX'
         INFOT = 1
         zgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = '/'
         zgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ = 'R'
         zgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         EQ = 'C'
         zgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGESVX', INFOT, NOUT, LERR, OK );

         // ZGESVXX

         N_ERR_BNDS = 3
         NPARAMS = 1
         SRNAMT = 'ZGESVXX'
         INFOT = 1
         zgesvxx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgesvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgesvxx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgesvxx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgesvxx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgesvxx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = '/'
         zgesvxx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ = 'R'
         zgesvxx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         EQ = 'C'
         zgesvxx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgesvxx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgesvxx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGESVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // ZGBSV

         SRNAMT = 'ZGBSV '
         INFOT = 1
         zgbsv(-1, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgbsv(1, -1, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgbsv(1, 0, -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgbsv(0, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgbsv(1, 1, 1, 0, A, 3, IP, B, 1, INFO );
         chkxer('ZGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zgbsv(2, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGBSV ', INFOT, NOUT, LERR, OK );

         // ZGBSVX

         SRNAMT = 'ZGBSVX'
         INFOT = 1
         zgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         EQ = '/'
         zgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         EQ = 'R'
         zgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         EQ = 'C'
         zgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGBSVX', INFOT, NOUT, LERR, OK );

         // ZGBSVXX

         N_ERR_BNDS = 3
         NPARAMS = 1
         SRNAMT = 'ZGBSVXX'
         INFOT = 1
         zgbsvxx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgbsvxx('N', '/', 0, 1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgbsvxx('N', 'N', -1, 1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgbsvxx('N', 'N', 2, -1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgbsvxx('N', 'N', 2, 1, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zgbsvxx('N', 'N', 0, 1, 1, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgbsvxx('N', 'N', 2, 1, 1, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 3, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         EQ = '/'
         zgbsvxx('F', 'N', 0, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         EQ = 'R'
         zgbsvxx('F', 'N', 1, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         EQ = 'C'
         zgbsvxx('F', 'N', 1, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 4, IP, EQ, R, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZGBSVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // ZGTSV

         SRNAMT = 'ZGTSV '
         INFOT = 1
         zgtsv(-1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('ZGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgtsv(0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('ZGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgtsv(2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('ZGTSV ', INFOT, NOUT, LERR, OK );

         // ZGTSVX

         SRNAMT = 'ZGTSVX'
         INFOT = 1
         zgtsvx('/', 'N', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgtsvx('N', '/', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgtsvx('N', 'N', -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgtsvx('N', 'N', 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZGTSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // ZHESV_ROOK

         SRNAMT = 'ZHESV_ROOK'
         INFOT = 1
         zhesv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhesv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhesv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhesv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // ZPOSV

         SRNAMT = 'ZPOSV '
         INFOT = 1
         zposv('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zposv('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('ZPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zposv('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('ZPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zposv('U', 2, 0, A, 1, B, 2, INFO );
         chkxer('ZPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zposv('U', 2, 0, A, 2, B, 1, INFO );
         chkxer('ZPOSV ', INFOT, NOUT, LERR, OK );

         // ZPOSVX

         SRNAMT = 'ZPOSVX'
         INFOT = 1
         zposvx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zposvx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zposvx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zposvx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zposvx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zposvx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         EQ = '/'
         zposvx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = 'Y'
         zposvx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPOSVX', INFOT, NOUT, LERR, OK );

         // ZPOSVXX

         N_ERR_BNDS = 3
         NPARAMS = 1
         SRNAMT = 'ZPOSVXX'
         INFOT = 1
         zposvxx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zposvxx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zposvxx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zposvxx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zposvxx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zposvxx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         EQ = '/'
         zposvxx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = 'Y'
         zposvxx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zposvxx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zposvxx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZPOSVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // ZPPSV

         SRNAMT = 'ZPPSV '
         INFOT = 1
         zppsv('/', 0, 0, A, B, 1, INFO );
         chkxer('ZPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zppsv('U', -1, 0, A, B, 1, INFO );
         chkxer('ZPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zppsv('U', 0, -1, A, B, 1, INFO );
         chkxer('ZPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zppsv('U', 2, 0, A, B, 1, INFO );
         chkxer('ZPPSV ', INFOT, NOUT, LERR, OK );

         // ZPPSVX

         SRNAMT = 'ZPPSVX'
         INFOT = 1
         zppsvx('/', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zppsvx('N', '/', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zppsvx('N', 'U', -1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zppsvx('N', 'U', 0, -1, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         EQ = '/'
         zppsvx('F', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         EQ = 'Y'
         zppsvx('F', 'U', 1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPPSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // ZPBSV

         SRNAMT = 'ZPBSV '
         INFOT = 1
         zpbsv('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbsv('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbsv('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('ZPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpbsv('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('ZPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zpbsv('U', 1, 1, 0, A, 1, B, 2, INFO );
         chkxer('ZPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zpbsv('U', 2, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBSV ', INFOT, NOUT, LERR, OK );

         // ZPBSVX

         SRNAMT = 'ZPBSVX'
         INFOT = 1
         zpbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zpbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zpbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = '/'
         zpbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ = 'Y'
         zpbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPBSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // ZPTSV

         SRNAMT = 'ZPTSV '
         INFOT = 1
         zptsv(-1, 0, R, A( 1, 1 ), B, 1, INFO );
         chkxer('ZPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zptsv(0, -1, R, A( 1, 1 ), B, 1, INFO );
         chkxer('ZPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zptsv(2, 0, R, A( 1, 1 ), B, 1, INFO );
         chkxer('ZPTSV ', INFOT, NOUT, LERR, OK );

         // ZPTSVX

         SRNAMT = 'ZPTSVX'
         INFOT = 1
         zptsvx('/', 0, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zptsvx('N', -1, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zptsvx('N', 0, -1, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zptsvx('N', 2, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zptsvx('N', 2, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZPTSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HE' ) ) {

         // ZHESV

         SRNAMT = 'ZHESV '
         INFOT = 1
         zhesv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhesv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhesv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhesv('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhesv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhesv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhesv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('ZHESV ', INFOT, NOUT, LERR, OK );

         // ZHESVX

         SRNAMT = 'ZHESVX'
         INFOT = 1
         zhesvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhesvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhesvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhesvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zhesvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhesvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zhesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zhesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zhesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, RW, INFO );
         chkxer('ZHESVX', INFOT, NOUT, LERR, OK );

         // ZHESVXX

         N_ERR_BNDS = 3
         NPARAMS = 1
         SRNAMT = 'ZHESVXX'
         INFOT = 1
         zhesvxx('/', 'U', 0, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhesvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhesvxx('N', 'U', -1, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhesvxx('N', 'U', 0, -1, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zhesvxx('N', 'U', 2, 0, A, 1, AF, 2, IP, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhesvxx('N', 'U', 2, 0, A, 2, AF, 1, IP, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C,  NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         EQ = '/'
         zhesvxx('F', 'U', 0, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         EQ = 'Y'
         zhesvxx('F', 'U', 1, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zhesvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zhesvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZHESVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // ZHESV_ROOK

         SRNAMT = 'ZHESV_ROOK'
         INFOT = 1
         zhesv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhesv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhesv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhesv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhesv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhesv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('ZHESV_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HK' ) ) {

         // ZSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a Hermitian indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         SRNAMT = 'ZHESV_RK'
         INFOT = 1
         zhesv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhesv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhesv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhesv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zhesv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zhesv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zhesv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('ZHESV_RK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HP' ) ) {

         // ZHPSV

         SRNAMT = 'ZHPSV '
         INFOT = 1
         zhpsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('ZHPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhpsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('ZHPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhpsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('ZHPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhpsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('ZHPSV ', INFOT, NOUT, LERR, OK );

         // ZHPSVX

         SRNAMT = 'ZHPSVX'
         INFOT = 1
         zhpsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhpsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhpsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhpsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zhpsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zhpsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZHPSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // ZSYSV

         SRNAMT = 'ZSYSV '
         INFOT = 1
         zsysv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zsysv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zsysv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zsysv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zsysv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('ZSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zsysv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('ZSYSV ', INFOT, NOUT, LERR, OK );

         // ZSYSVX

         SRNAMT = 'ZSYSVX'
         INFOT = 1
         zsysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zsysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zsysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zsysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zsysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zsysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, RW, INFO );
         chkxer('ZSYSVX', INFOT, NOUT, LERR, OK );

         // ZSYSVXX

         N_ERR_BNDS = 3
         NPARAMS = 1
         SRNAMT = 'ZSYSVXX'
         INFOT = 1
         EQ = 'N'
         zsysvxx('/', 'U', 0, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zsysvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zsysvxx('N', 'U', -1, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         EQ = '/'
         zsysvxx('N', 'U', 0, -1, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         EQ = 'Y'
         INFOT = 6
         zsysvxx('N', 'U', 2, 0, A, 1, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zsysvxx('N', 'U', 2, 0, A, 2, AF, 1, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zsysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, 'A', R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ='Y'
         zsysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         EQ='Y'
         R(1) = -ONE
         zsysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         EQ = 'N'
         zsysvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zsysvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('ZSYSVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // ZSYSV_ROOK

         SRNAMT = 'ZSYSV_ROOK'
         INFOT = 1
         zsysv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zsysv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zsysv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zsysv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zsysv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('ZSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zsysv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // ZSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         SRNAMT = 'ZSYSV_RK'
         INFOT = 1
         zsysv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zsysv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zsysv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zsysv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zsysv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zsysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zsysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('ZSYSV_RK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // ZSPSV

         SRNAMT = 'ZSPSV '
         INFOT = 1
         zspsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('ZSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zspsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('ZSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zspsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('ZSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zspsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('ZSPSV ', INFOT, NOUT, LERR, OK );

         // ZSPSVX

         SRNAMT = 'ZSPSVX'
         INFOT = 1
         zspsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zspsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zspsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zspsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zspsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zspsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('ZSPSVX', INFOT, NOUT, LERR, OK );
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

      // End of ZERRVXX

      }
