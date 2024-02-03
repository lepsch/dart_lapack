      SUBROUTINE ZERRPO( PATH, NUNIT )

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
      String             C2;
      int                I, INFO, J;
      double             ANRM, RCOND;
      // ..
      // .. Local Arrays ..
      double             R( NMAX ), R1( NMAX ), R2( NMAX );
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZPBCON, ZPBEQU, ZPBRFS, ZPBTF2, ZPBTRF, ZPBTRS, ZPOCON, ZPOEQU, ZPORFS, ZPOTF2, ZPOTRF, ZPOTRI, ZPOTRS, ZPPCON, ZPPEQU, ZPPRFS, ZPPTRF, ZPPTRI, ZPPTRS
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

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )
   10    CONTINUE
         B( J ) = 0.D0
         R1( J ) = 0.D0
         R2( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
   20 CONTINUE
      ANRM = 1.D0
      OK = .TRUE.

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite matrix.

      if ( LSAMEN( 2, C2, 'PO' ) ) {

         // ZPOTRF

         SRNAMT = 'ZPOTRF'
         INFOT = 1
         zpotrf('/', 0, A, 1, INFO );
         chkxer('ZPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpotrf('U', -1, A, 1, INFO );
         chkxer('ZPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpotrf('U', 2, A, 1, INFO );
         chkxer('ZPOTRF', INFOT, NOUT, LERR, OK );

         // ZPOTF2

         SRNAMT = 'ZPOTF2'
         INFOT = 1
         zpotf2('/', 0, A, 1, INFO );
         chkxer('ZPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpotf2('U', -1, A, 1, INFO );
         chkxer('ZPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpotf2('U', 2, A, 1, INFO );
         chkxer('ZPOTF2', INFOT, NOUT, LERR, OK );

         // ZPOTRI

         SRNAMT = 'ZPOTRI'
         INFOT = 1
         zpotri('/', 0, A, 1, INFO );
         chkxer('ZPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpotri('U', -1, A, 1, INFO );
         chkxer('ZPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpotri('U', 2, A, 1, INFO );
         chkxer('ZPOTRI', INFOT, NOUT, LERR, OK );

         // ZPOTRS

         SRNAMT = 'ZPOTRS'
         INFOT = 1
         zpotrs('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpotrs('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('ZPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpotrs('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('ZPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpotrs('U', 2, 1, A, 1, B, 2, INFO );
         chkxer('ZPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zpotrs('U', 2, 1, A, 2, B, 1, INFO );
         chkxer('ZPOTRS', INFOT, NOUT, LERR, OK );

         // ZPORFS

         SRNAMT = 'ZPORFS'
         INFOT = 1
         zporfs('/', 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zporfs('U', -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zporfs('U', 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zporfs('U', 2, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zporfs('U', 2, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zporfs('U', 2, 1, A, 2, AF, 2, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zporfs('U', 2, 1, A, 2, AF, 2, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', INFOT, NOUT, LERR, OK );

         // ZPOCON

         SRNAMT = 'ZPOCON'
         INFOT = 1
         zpocon('/', 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpocon('U', -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpocon('U', 2, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpocon('U', 1, A, 1, -ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', INFOT, NOUT, LERR, OK );

         // ZPOEQU

         SRNAMT = 'ZPOEQU'
         INFOT = 1
         zpoequ(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPOEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpoequ(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPOEQU', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite packed matrix.

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // ZPPTRF

         SRNAMT = 'ZPPTRF'
         INFOT = 1
         zpptrf('/', 0, A, INFO );
         chkxer('ZPPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpptrf('U', -1, A, INFO );
         chkxer('ZPPTRF', INFOT, NOUT, LERR, OK );

         // ZPPTRI

         SRNAMT = 'ZPPTRI'
         INFOT = 1
         zpptri('/', 0, A, INFO );
         chkxer('ZPPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpptri('U', -1, A, INFO );
         chkxer('ZPPTRI', INFOT, NOUT, LERR, OK );

         // ZPPTRS

         SRNAMT = 'ZPPTRS'
         INFOT = 1
         zpptrs('/', 0, 0, A, B, 1, INFO );
         chkxer('ZPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpptrs('U', -1, 0, A, B, 1, INFO );
         chkxer('ZPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpptrs('U', 0, -1, A, B, 1, INFO );
         chkxer('ZPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zpptrs('U', 2, 1, A, B, 1, INFO );
         chkxer('ZPPTRS', INFOT, NOUT, LERR, OK );

         // ZPPRFS

         SRNAMT = 'ZPPRFS'
         INFOT = 1
         zpprfs('/', 0, 0, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpprfs('U', -1, 0, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpprfs('U', 0, -1, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zpprfs('U', 2, 1, A, AF, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zpprfs('U', 2, 1, A, AF, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', INFOT, NOUT, LERR, OK );

         // ZPPCON

         SRNAMT = 'ZPPCON'
         INFOT = 1
         zppcon('/', 0, A, ANRM, RCOND, W, R, INFO );
         chkxer('ZPPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zppcon('U', -1, A, ANRM, RCOND, W, R, INFO );
         chkxer('ZPPCON', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zppcon('U', 1, A, -ANRM, RCOND, W, R, INFO );
         chkxer('ZPPCON', INFOT, NOUT, LERR, OK );

         // ZPPEQU

         SRNAMT = 'ZPPEQU'
         INFOT = 1
         zppequ('/', 0, A, R1, RCOND, ANRM, INFO );
         chkxer('ZPPEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zppequ('U', -1, A, R1, RCOND, ANRM, INFO );
         chkxer('ZPPEQU', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite band matrix.

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // ZPBTRF

         SRNAMT = 'ZPBTRF'
         INFOT = 1
         zpbtrf('/', 0, 0, A, 1, INFO );
         chkxer('ZPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbtrf('U', -1, 0, A, 1, INFO );
         chkxer('ZPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbtrf('U', 1, -1, A, 1, INFO );
         chkxer('ZPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpbtrf('U', 2, 1, A, 1, INFO );
         chkxer('ZPBTRF', INFOT, NOUT, LERR, OK );

         // ZPBTF2

         SRNAMT = 'ZPBTF2'
         INFOT = 1
         zpbtf2('/', 0, 0, A, 1, INFO );
         chkxer('ZPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbtf2('U', -1, 0, A, 1, INFO );
         chkxer('ZPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbtf2('U', 1, -1, A, 1, INFO );
         chkxer('ZPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpbtf2('U', 2, 1, A, 1, INFO );
         chkxer('ZPBTF2', INFOT, NOUT, LERR, OK );

         // ZPBTRS

         SRNAMT = 'ZPBTRS'
         INFOT = 1
         zpbtrs('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbtrs('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbtrs('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpbtrs('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zpbtrs('U', 2, 1, 1, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zpbtrs('U', 2, 0, 1, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', INFOT, NOUT, LERR, OK );

         // ZPBRFS

         SRNAMT = 'ZPBRFS'
         INFOT = 1
         zpbrfs('/', 0, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbrfs('U', -1, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbrfs('U', 1, -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zpbrfs('U', 0, 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zpbrfs('U', 2, 1, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zpbrfs('U', 2, 1, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', INFOT, NOUT, LERR, OK );

         // ZPBCON

         SRNAMT = 'ZPBCON'
         INFOT = 1
         zpbcon('/', 0, 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbcon('U', -1, 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbcon('U', 1, -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpbcon('U', 2, 1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zpbcon('U', 1, 0, A, 1, -ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', INFOT, NOUT, LERR, OK );

         // ZPBEQU

         SRNAMT = 'ZPBEQU'
         INFOT = 1
         zpbequ('/', 0, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpbequ('U', -1, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpbequ('U', 1, -1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zpbequ('U', 2, 1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRPO

      }
