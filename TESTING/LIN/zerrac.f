      SUBROUTINE ZERRAC( NUNIT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, ITER, J;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), W( 2*NMAX ), X( NMAX );
      double             RWORK( NMAX );
      COMPLEX*16         WORK(NMAX*NMAX);
      COMPLEX            SWORK(NMAX*NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZCPOSV
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1.0 / DBLE( I+J );
            AF( I, J ) = 1.0 / DBLE( I+J );
         } // 10
         B( J ) = 0.0;
         R1( J ) = 0.0;
         R2( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
         C( J ) = 0.0;
         R( J ) = 0.0;
      } // 20
      OK = true;

      SRNAMT = 'ZCPOSV';
      INFOT = 1;
      zcposv('/',0,0,A,1,B,1,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCPOSV', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zcposv('U',-1,0,A,1,B,1,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCPOSV', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zcposv('U',0,-1,A,1,B,1,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCPOSV', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zcposv('U',2,1,A,1,B,2,X,2,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCPOSV', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zcposv('U',2,1,A,2,B,1,X,2,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCPOSV', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      zcposv('U',2,1,A,2,B,2,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCPOSV', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )'ZCPOSV';
      } else {
         WRITE( NOUT, FMT = 9998 )'ZCPOSV';
      }

 9999 FORMAT( 1X, A6, ' drivers passed the tests of the error exits' );
 9998 FORMAT( ' *** ', A6, ' drivers failed the tests of the error ', 'exits ***' );

      return;

      // End of ZERRAC

      }
