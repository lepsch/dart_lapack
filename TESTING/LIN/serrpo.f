      SUBROUTINE SERRPO( PATH, NUNIT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J;
      REAL               ANRM, RCOND;
      // ..
      // .. Local Arrays ..
      int                IW( NMAX );
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SPBCON, SPBEQU, SPBRFS, SPBTF2, SPBTRF, SPBTRS, SPOCON, SPOEQU, SPORFS, SPOTF2, SPOTRF, SPOTRI, SPOTRS, SPPCON, SPPEQU, SPPRFS, SPPTRF, SPPTRI, SPPTRS
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

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1. / REAL( I+J );
            AF( I, J ) = 1. / REAL( I+J );
         } // 10
         B( J ) = 0.;
         R1( J ) = 0.;
         R2( J ) = 0.;
         W( J ) = 0.;
         X( J ) = 0.;
         IW( J ) = J;
      } // 20
      OK = true;

      if ( LSAMEN( 2, C2, 'PO' ) ) {

         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive definite matrix.

         // SPOTRF

         SRNAMT = 'SPOTRF';
         INFOT = 1;
         spotrf('/', 0, A, 1, INFO );
         chkxer('SPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spotrf('U', -1, A, 1, INFO );
         chkxer('SPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         spotrf('U', 2, A, 1, INFO );
         chkxer('SPOTRF', INFOT, NOUT, LERR, OK );

         // SPOTF2

         SRNAMT = 'SPOTF2';
         INFOT = 1;
         spotf2('/', 0, A, 1, INFO );
         chkxer('SPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spotf2('U', -1, A, 1, INFO );
         chkxer('SPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         spotf2('U', 2, A, 1, INFO );
         chkxer('SPOTF2', INFOT, NOUT, LERR, OK );

         // SPOTRI

         SRNAMT = 'SPOTRI';
         INFOT = 1;
         spotri('/', 0, A, 1, INFO );
         chkxer('SPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spotri('U', -1, A, 1, INFO );
         chkxer('SPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         spotri('U', 2, A, 1, INFO );
         chkxer('SPOTRI', INFOT, NOUT, LERR, OK );

         // SPOTRS

         SRNAMT = 'SPOTRS';
         INFOT = 1;
         spotrs('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('SPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spotrs('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('SPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spotrs('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('SPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         spotrs('U', 2, 1, A, 1, B, 2, INFO );
         chkxer('SPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         spotrs('U', 2, 1, A, 2, B, 1, INFO );
         chkxer('SPOTRS', INFOT, NOUT, LERR, OK );

         // SPORFS

         SRNAMT = 'SPORFS';
         INFOT = 1;
         sporfs('/', 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sporfs('U', -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sporfs('U', 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sporfs('U', 2, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sporfs('U', 2, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sporfs('U', 2, 1, A, 2, AF, 2, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sporfs('U', 2, 1, A, 2, AF, 2, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPORFS', INFOT, NOUT, LERR, OK );

         // SPOCON

         SRNAMT = 'SPOCON';
         INFOT = 1;
         spocon('/', 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spocon('U', -1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         spocon('U', 2, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPOCON', INFOT, NOUT, LERR, OK );

         // SPOEQU

         SRNAMT = 'SPOEQU';
         INFOT = 1;
         spoequ(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('SPOEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spoequ(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('SPOEQU', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive definite packed matrix.

         // SPPTRF

         SRNAMT = 'SPPTRF';
         INFOT = 1;
         spptrf('/', 0, A, INFO );
         chkxer('SPPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spptrf('U', -1, A, INFO );
         chkxer('SPPTRF', INFOT, NOUT, LERR, OK );

         // SPPTRI

         SRNAMT = 'SPPTRI';
         INFOT = 1;
         spptri('/', 0, A, INFO );
         chkxer('SPPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spptri('U', -1, A, INFO );
         chkxer('SPPTRI', INFOT, NOUT, LERR, OK );

         // SPPTRS

         SRNAMT = 'SPPTRS';
         INFOT = 1;
         spptrs('/', 0, 0, A, B, 1, INFO );
         chkxer('SPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spptrs('U', -1, 0, A, B, 1, INFO );
         chkxer('SPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spptrs('U', 0, -1, A, B, 1, INFO );
         chkxer('SPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         spptrs('U', 2, 1, A, B, 1, INFO );
         chkxer('SPPTRS', INFOT, NOUT, LERR, OK );

         // SPPRFS

         SRNAMT = 'SPPRFS';
         INFOT = 1;
         spprfs('/', 0, 0, A, AF, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spprfs('U', -1, 0, A, AF, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spprfs('U', 0, -1, A, AF, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         spprfs('U', 2, 1, A, AF, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         spprfs('U', 2, 1, A, AF, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPPRFS', INFOT, NOUT, LERR, OK );

         // SPPCON

         SRNAMT = 'SPPCON';
         INFOT = 1;
         sppcon('/', 0, A, ANRM, RCOND, W, IW, INFO );
         chkxer('SPPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sppcon('U', -1, A, ANRM, RCOND, W, IW, INFO );
         chkxer('SPPCON', INFOT, NOUT, LERR, OK );

         // SPPEQU

         SRNAMT = 'SPPEQU';
         INFOT = 1;
         sppequ('/', 0, A, R1, RCOND, ANRM, INFO );
         chkxer('SPPEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sppequ('U', -1, A, R1, RCOND, ANRM, INFO );
         chkxer('SPPEQU', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive definite band matrix.

         // SPBTRF

         SRNAMT = 'SPBTRF';
         INFOT = 1;
         spbtrf('/', 0, 0, A, 1, INFO );
         chkxer('SPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spbtrf('U', -1, 0, A, 1, INFO );
         chkxer('SPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spbtrf('U', 1, -1, A, 1, INFO );
         chkxer('SPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         spbtrf('U', 2, 1, A, 1, INFO );
         chkxer('SPBTRF', INFOT, NOUT, LERR, OK );

         // SPBTF2

         SRNAMT = 'SPBTF2';
         INFOT = 1;
         spbtf2('/', 0, 0, A, 1, INFO );
         chkxer('SPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spbtf2('U', -1, 0, A, 1, INFO );
         chkxer('SPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spbtf2('U', 1, -1, A, 1, INFO );
         chkxer('SPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         spbtf2('U', 2, 1, A, 1, INFO );
         chkxer('SPBTF2', INFOT, NOUT, LERR, OK );

         // SPBTRS

         SRNAMT = 'SPBTRS';
         INFOT = 1;
         spbtrs('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('SPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spbtrs('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('SPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spbtrs('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('SPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         spbtrs('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('SPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         spbtrs('U', 2, 1, 1, A, 1, B, 1, INFO );
         chkxer('SPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         spbtrs('U', 2, 0, 1, A, 1, B, 1, INFO );
         chkxer('SPBTRS', INFOT, NOUT, LERR, OK );

         // SPBRFS

         SRNAMT = 'SPBRFS';
         INFOT = 1;
         spbrfs('/', 0, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spbrfs('U', -1, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spbrfs('U', 1, -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         spbrfs('U', 0, 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         spbrfs('U', 2, 1, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         spbrfs('U', 2, 1, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         spbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         spbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SPBRFS', INFOT, NOUT, LERR, OK );

         // SPBCON

         SRNAMT = 'SPBCON';
         INFOT = 1;
         spbcon('/', 0, 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spbcon('U', -1, 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spbcon('U', 1, -1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         spbcon('U', 2, 1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SPBCON', INFOT, NOUT, LERR, OK );

         // SPBEQU

         SRNAMT = 'SPBEQU';
         INFOT = 1;
         spbequ('/', 0, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('SPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spbequ('U', -1, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('SPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         spbequ('U', 1, -1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('SPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         spbequ('U', 2, 1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('SPBEQU', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of SERRPO

      }
