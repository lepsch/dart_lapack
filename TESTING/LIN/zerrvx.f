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
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J;
      double             RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), RF( NMAX ), RW( NMAX );
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGBSV, ZGBSVX, ZGESV, ZGESVX, ZGTSV, ZGTSVX, ZHESV, ZHESV_RK, ZHESV_ROOK, ZHESVX, ZHPSV, ZHPSVX, ZPBSV, ZPBSVX, ZPOSV, ZPOSVX, ZPPSV, ZPPSVX, ZPTSV, ZPTSVX, ZSPSV, ZSPSVX, ZSYSV, ZSYSV_AA, ZSYSV_RK, ZSYSV_ROOK, ZSYSVX, ZHESV_AA_2STAGE
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      // Set the variables to innocuous values.

      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )
   10    CONTINUE
         B( J ) = 0.D0
         E( J ) = 0.D0
         R1( J ) = 0.D0
         R2( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
         C( J ) = 0.D0
         R( J ) = 0.D0
         IP( J ) = J
   20 CONTINUE
      EQ = ' '
      OK = .TRUE.

      if ( LSAMEN( 2, C2, 'GE' ) ) {

         // ZGESV

         SRNAMT = 'ZGESV '
         INFOT = 1
         CALL ZGESV( -1, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGESV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGESV( 0, -1, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGESV ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGESV( 2, 1, A, 1, IP, B, 2, INFO )
         CALL CHKXER( 'ZGESV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGESV( 2, 1, A, 2, IP, B, 1, INFO )
         CALL CHKXER( 'ZGESV ', INFOT, NOUT, LERR, OK )

         // ZGESVX

         SRNAMT = 'ZGESVX'
         INFOT = 1
         CALL ZGESVX( '/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGESVX( 'N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGESVX( 'N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGESVX( 'N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGESVX( 'N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGESVX( 'N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         EQ = '/'
         CALL ZGESVX( 'F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         EQ = 'R'
         CALL ZGESVX( 'F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
         EQ = 'C'
         CALL ZGESVX( 'F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL ZGESVX( 'N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL ZGESVX( 'N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGESVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // ZGBSV

         SRNAMT = 'ZGBSV '
         INFOT = 1
         CALL ZGBSV( -1, 0, 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGBSV( 1, -1, 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGBSV( 1, 0, -1, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGBSV( 0, 0, 0, -1, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGBSV( 1, 1, 1, 0, A, 3, IP, B, 1, INFO )
         CALL CHKXER( 'ZGBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZGBSV( 2, 0, 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'ZGBSV ', INFOT, NOUT, LERR, OK )

         // ZGBSVX

         SRNAMT = 'ZGBSVX'
         INFOT = 1
         CALL ZGBSVX( '/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGBSVX( 'N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGBSVX( 'N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGBSVX( 'N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGBSVX( 'N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGBSVX( 'N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGBSVX( 'N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGBSVX( 'N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
         EQ = '/'
         CALL ZGBSVX( 'F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         EQ = 'R'
         CALL ZGBSVX( 'F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
         EQ = 'C'
         CALL ZGBSVX( 'F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL ZGBSVX( 'N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL ZGBSVX( 'N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGBSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // ZGTSV

         SRNAMT = 'ZGTSV '
         INFOT = 1
         CALL ZGTSV( -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO )
         CALL CHKXER( 'ZGTSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGTSV( 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO )
         CALL CHKXER( 'ZGTSV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGTSV( 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO )
         CALL CHKXER( 'ZGTSV ', INFOT, NOUT, LERR, OK )

         // ZGTSVX

         SRNAMT = 'ZGTSVX'
         INFOT = 1
         CALL ZGTSVX( '/', 'N', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGTSVX( 'N', '/', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGTSVX( 'N', 'N', -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGTSVX( 'N', 'N', 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL ZGTSVX( 'N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL ZGTSVX( 'N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZGTSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // ZPOSV

         SRNAMT = 'ZPOSV '
         INFOT = 1
         CALL ZPOSV( '/', 0, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPOSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPOSV( 'U', -1, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPOSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPOSV( 'U', 0, -1, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPOSV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZPOSV( 'U', 2, 0, A, 1, B, 2, INFO )
         CALL CHKXER( 'ZPOSV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZPOSV( 'U', 2, 0, A, 2, B, 1, INFO )
         CALL CHKXER( 'ZPOSV ', INFOT, NOUT, LERR, OK )

         // ZPOSVX

         SRNAMT = 'ZPOSVX'
         INFOT = 1
         CALL ZPOSVX( '/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPOSVX( 'N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPOSVX( 'N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZPOSVX( 'N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZPOSVX( 'N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZPOSVX( 'N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         EQ = '/'
         CALL ZPOSVX( 'F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         EQ = 'Y'
         CALL ZPOSVX( 'F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZPOSVX( 'N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL ZPOSVX( 'N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPOSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // ZPPSV

         SRNAMT = 'ZPPSV '
         INFOT = 1
         CALL ZPPSV( '/', 0, 0, A, B, 1, INFO )
         CALL CHKXER( 'ZPPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPPSV( 'U', -1, 0, A, B, 1, INFO )
         CALL CHKXER( 'ZPPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPPSV( 'U', 0, -1, A, B, 1, INFO )
         CALL CHKXER( 'ZPPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZPPSV( 'U', 2, 0, A, B, 1, INFO )
         CALL CHKXER( 'ZPPSV ', INFOT, NOUT, LERR, OK )

         // ZPPSVX

         SRNAMT = 'ZPPSVX'
         INFOT = 1
         CALL ZPPSVX( '/', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPPSVX( 'N', '/', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPPSVX( 'N', 'U', -1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZPPSVX( 'N', 'U', 0, -1, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         EQ = '/'
         CALL ZPPSVX( 'F', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         EQ = 'Y'
         CALL ZPPSVX( 'F', 'U', 1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZPPSVX( 'N', 'U', 2, 0, A, AF, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZPPSVX( 'N', 'U', 2, 0, A, AF, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPPSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // ZPBSV

         SRNAMT = 'ZPBSV '
         INFOT = 1
         CALL ZPBSV( '/', 0, 0, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPBSV( 'U', -1, 0, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPBSV( 'U', 1, -1, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZPBSV( 'U', 0, 0, -1, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZPBSV( 'U', 1, 1, 0, A, 1, B, 2, INFO )
         CALL CHKXER( 'ZPBSV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZPBSV( 'U', 2, 0, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'ZPBSV ', INFOT, NOUT, LERR, OK )

         // ZPBSVX

         SRNAMT = 'ZPBSVX'
         INFOT = 1
         CALL ZPBSVX( '/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPBSVX( 'N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPBSVX( 'N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZPBSVX( 'N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZPBSVX( 'N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZPBSVX( 'N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZPBSVX( 'N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         EQ = '/'
         CALL ZPBSVX( 'F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         EQ = 'Y'
         CALL ZPBSVX( 'F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZPBSVX( 'N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL ZPBSVX( 'N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPBSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // ZPTSV

         SRNAMT = 'ZPTSV '
         INFOT = 1
         CALL ZPTSV( -1, 0, R, A( 1, 1 ), B, 1, INFO )
         CALL CHKXER( 'ZPTSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPTSV( 0, -1, R, A( 1, 1 ), B, 1, INFO )
         CALL CHKXER( 'ZPTSV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZPTSV( 2, 0, R, A( 1, 1 ), B, 1, INFO )
         CALL CHKXER( 'ZPTSV ', INFOT, NOUT, LERR, OK )

         // ZPTSVX

         SRNAMT = 'ZPTSVX'
         INFOT = 1
         CALL ZPTSVX( '/', 0, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZPTSVX( 'N', -1, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZPTSVX( 'N', 0, -1, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZPTSVX( 'N', 2, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPTSVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZPTSVX( 'N', 2, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZPTSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'HE' ) ) {

         // ZHESV

         SRNAMT = 'ZHESV '
         INFOT = 1
         CALL ZHESV( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHESV( 'U', -1, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHESV( 'U', 0, -1, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHESV( 'U', 2, 0, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZHESV( 'U', 2, 0, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHESV( 'U', 0, 0, A, 1, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHESV( 'U', 0, 0, A, 1, IP, B, 1, W, -2, INFO )
         CALL CHKXER( 'ZHESV ', INFOT, NOUT, LERR, OK )

         // ZHESVX

         SRNAMT = 'ZHESVX'
         INFOT = 1
         CALL ZHESVX( '/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHESVX( 'N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHESVX( 'N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZHESVX( 'N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZHESVX( 'N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZHESVX( 'N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZHESVX( 'N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZHESVX( 'N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL ZHESVX( 'N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, RW, INFO )
         CALL CHKXER( 'ZHESVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // ZHESV_ROOK

         SRNAMT = 'ZHESV_ROOK'
         INFOT = 1
         CALL ZHESV_ROOK( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHESV_ROOK( 'U', -1, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHESV_ROOK( 'U', 0, -1, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHESV_ROOK( 'U', 2, 0, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZHESV_ROOK( 'U', 2, 0, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHESV_ROOK( 'U', 0, 0, A, 1, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHESV_ROOK( 'U', 0, 0, A, 1, IP, B, 1, W, -2, INFO )
         CALL CHKXER( 'ZHESV_ROOK', INFOT, NOUT, LERR, OK )

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
         CALL ZHESV_RK( '/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHESV_RK( 'U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHESV_RK( 'U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHESV_RK( 'U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZHESV_RK( 'U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZHESV_RK( 'U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZHESV_RK( 'U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO )
         CALL CHKXER( 'ZHESV_RK', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'HA' ) ) {

         // ZHESV_AASEN

         SRNAMT = 'ZHESV_AA'
         INFOT = 1
         CALL ZHESV_AA( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHESV_AA( 'U', -1, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHESV_AA( 'U', 0, -1, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHESV_AA( 'U', 2, 0, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZHESV_AA( 'U', 2, 0, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHESV_AA( 'U', 3, 1, A, 3, IP, B, 3, W, 6, INFO )
         CALL CHKXER( 'ZHESV_AA', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'H2' ) ) {

         // ZHESV_AASEN_2STAGE

         SRNAMT = 'ZHESV_AA_2STAGE'
         INFOT = 1
         CALL ZHESV_AA_2STAGE( '/', 0, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHESV_AA_2STAGE( 'U', -1, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHESV_AA_2STAGE( 'U', 0, -1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHESV_AA_2STAGE( 'U', 2, 1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZHESV_AA_2STAGE( 'U', 2, 1, A, 2, A, 1, IP, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZHESV_AA_2STAGE( 'U', 2, 1, A, 2, A, 8, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZHESV_AA_2STAGE( 'U', 2, 1, A, 2, A, 8, IP, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZHESV_AA_2STAGE', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'SA' ) ) {

         // ZSYSV_AASEN

         SRNAMT = 'ZSYSV_AA'
         INFOT = 1
         CALL ZSYSV_AA( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSYSV_AA( 'U', -1, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSYSV_AA( 'U', 0, -1, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZSYSV_AA( 'U', 2, 0, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZSYSV_AA( 'U', 2, 0, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZSYSV_AA( 'U', 3, 1, A, 3, IP, B, 3, W, 6, INFO )
         CALL CHKXER( 'ZSYSV_AA', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // ZSYSV_AASEN_2STAGE

         SRNAMT = 'ZSYSV_AA_2STAGE'
         INFOT = 1
         CALL ZSYSV_AA_2STAGE( '/', 0, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSYSV_AA_2STAGE( 'U', -1, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSYSV_AA_2STAGE( 'U', 0, -1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZSYSV_AA_2STAGE( 'U', 2, 1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZSYSV_AA_2STAGE( 'U', 2, 1, A, 2, A, 1, IP, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZSYSV_AA_2STAGE( 'U', 2, 1, A, 2, A, 8, IP, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZSYSV_AA_2STAGE( 'U', 2, 1, A, 2, A, 8, IP, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'HP' ) ) {

         // ZHPSV

         SRNAMT = 'ZHPSV '
         INFOT = 1
         CALL ZHPSV( '/', 0, 0, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZHPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHPSV( 'U', -1, 0, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZHPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHPSV( 'U', 0, -1, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZHPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZHPSV( 'U', 2, 0, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZHPSV ', INFOT, NOUT, LERR, OK )

         // ZHPSVX

         SRNAMT = 'ZHPSVX'
         INFOT = 1
         CALL ZHPSVX( '/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZHPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHPSVX( 'N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZHPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHPSVX( 'N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZHPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZHPSVX( 'N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZHPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZHPSVX( 'N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZHPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZHPSVX( 'N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZHPSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // ZSYSV

         SRNAMT = 'ZSYSV '
         INFOT = 1
         CALL ZSYSV( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSYSV( 'U', -1, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSYSV( 'U', 0, -1, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZSYSV( 'U', 2, 0, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZSYSV( 'U', 2, 0, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZSYSV( 'U', 0, 0, A, 1, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZSYSV( 'U', 0, 0, A, 1, IP, B, 1, W, -2, INFO )
         CALL CHKXER( 'ZSYSV ', INFOT, NOUT, LERR, OK )

         // ZSYSVX

         SRNAMT = 'ZSYSVX'
         INFOT = 1
         CALL ZSYSVX( '/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSYSVX( 'N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSYSVX( 'N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZSYSVX( 'N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZSYSVX( 'N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZSYSVX( 'N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZSYSVX( 'N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZSYSVX( 'N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL ZSYSVX( 'N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, RW, INFO )
         CALL CHKXER( 'ZSYSVX', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // ZSYSV_ROOK

         SRNAMT = 'ZSYSV_ROOK'
         INFOT = 1
         CALL ZSYSV_ROOK( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSYSV_ROOK( 'U', -1, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSYSV_ROOK( 'U', 0, -1, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZSYSV_ROOK( 'U', 2, 0, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZSYSV_ROOK( 'U', 2, 0, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZSYSV_ROOK( 'U', 0, 0, A, 1, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'ZSYSV_ROOK', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZSYSV_ROOK( 'U', 0, 0, A, 1, IP, B, 1, W, -2, INFO )

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
         CALL ZSYSV_RK( '/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSYSV_RK( 'U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSYSV_RK( 'U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZSYSV_RK( 'U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZSYSV_RK( 'U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZSYSV_RK( 'U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZSYSV_RK( 'U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO )
         CALL CHKXER( 'ZSYSV_RK', INFOT, NOUT, LERR, OK )

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // ZSPSV

         SRNAMT = 'ZSPSV '
         INFOT = 1
         CALL ZSPSV( '/', 0, 0, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZSPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSPSV( 'U', -1, 0, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZSPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSPSV( 'U', 0, -1, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZSPSV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZSPSV( 'U', 2, 0, A, IP, B, 1, INFO )
         CALL CHKXER( 'ZSPSV ', INFOT, NOUT, LERR, OK )

         // ZSPSVX

         SRNAMT = 'ZSPSVX'
         INFOT = 1
         CALL ZSPSVX( '/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZSPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZSPSVX( 'N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZSPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZSPSVX( 'N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZSPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZSPSVX( 'N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZSPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZSPSVX( 'N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZSPSVX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZSPSVX( 'N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO )
         CALL CHKXER( 'ZSPSVX', INFOT, NOUT, LERR, OK )
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

 9999 FORMAT( 1X, A3, ' drivers passed the tests of the error exits' )
 9998 FORMAT( ' *** ', A3, ' drivers failed the tests of the error ',
     $      'exits ***' )

      RETURN

      // End of ZERRVX

      }
