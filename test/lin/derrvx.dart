      void derrvx(final int PATH, final int NUNIT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      String             EQ;
      String             C2;
      int                I, INFO, J;
      double             RCOND;
      int                IP( NMAX ), IW( NMAX );
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), E( NMAX ),  R( NMAX ), R1( NMAX ), R2( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGBSV, DGBSVX, DGESV, DGESVX, DGTSV, DGTSVX, DPBSV, DPBSVX, DPOSV, DPOSVX, DPPSV, DPPSVX, DPTSV, DPTSVX, DSPSV, DSPSVX, DSYSV, DSYSV_AA, DSYSV_RK, DSYSV_ROOK, DSYSVX, DSYSV_AA_2STAGE
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1.0 / (I+J).toDouble();
            AF[I][J] = 1.0 / (I+J).toDouble();
         } // 10
         B[J] = 0.0;
         E[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         C[J] = 0.0;
         R[J] = 0.0;
         IP[J] = J;
      } // 20
      EQ = ' ';
      infoc.OK = true;

      if ( lsamen( 2, C2, 'GE' ) ) {

         // DGESV

        srnamc.SRNAMT = 'DGESV ';
         infoc.INFOT = 1;
         dgesv(-1, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGESV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgesv(0, -1, A, 1, IP, B, 1, INFO );
         chkxer('DGESV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgesv(2, 1, A, 1, IP, B, 2, INFO );
         chkxer('DGESV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dgesv(2, 1, A, 2, IP, B, 1, INFO );
         chkxer('DGESV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGESVX

        srnamc.SRNAMT = 'DGESVX';
         infoc.INFOT = 1;
         dgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         EQ = '/';
         dgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         EQ = 'R';
         dgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         EQ = 'C';
         dgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         dgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 16;
         dgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // DGBSV

        srnamc.SRNAMT = 'DGBSV ';
         infoc.INFOT = 1;
         dgbsv(-1, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgbsv(1, -1, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgbsv(1, 0, -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgbsv(0, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dgbsv(1, 1, 1, 0, A, 3, IP, B, 1, INFO );
         chkxer('DGBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dgbsv(2, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGBSVX

        srnamc.SRNAMT = 'DGBSVX';
         infoc.INFOT = 1;
         dgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         EQ = '/';
         dgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         EQ = 'R';
         dgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         EQ = 'C';
         dgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 16;
         dgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 18;
         dgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'GT' ) ) {

         // DGTSV

        srnamc.SRNAMT = 'DGTSV ';
         infoc.INFOT = 1;
         dgtsv(-1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('DGTSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgtsv(0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('DGTSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dgtsv(2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('DGTSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGTSVX

        srnamc.SRNAMT = 'DGTSVX';
         infoc.INFOT = 1;
         dgtsvx('/', 'N', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgtsvx('N', '/', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgtsvx('N', 'N', -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgtsvx('N', 'N', 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         dgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 16;
         dgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'PO' ) ) {

         // DPOSV

        srnamc.SRNAMT = 'DPOSV ';
         infoc.INFOT = 1;
         dposv('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('DPOSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dposv('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('DPOSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dposv('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('DPOSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dposv('U', 2, 0, A, 1, B, 2, INFO );
         chkxer('DPOSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dposv('U', 2, 0, A, 2, B, 1, INFO );
         chkxer('DPOSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPOSVX

        srnamc.SRNAMT = 'DPOSVX';
         infoc.INFOT = 1;
         dposvx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dposvx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dposvx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dposvx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dposvx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dposvx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         EQ = '/';
         dposvx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         EQ = 'Y';
         dposvx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         dposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         dposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'PP' ) ) {

         // DPPSV

        srnamc.SRNAMT = 'DPPSV ';
         infoc.INFOT = 1;
         dppsv('/', 0, 0, A, B, 1, INFO );
         chkxer('DPPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dppsv('U', -1, 0, A, B, 1, INFO );
         chkxer('DPPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dppsv('U', 0, -1, A, B, 1, INFO );
         chkxer('DPPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dppsv('U', 2, 0, A, B, 1, INFO );
         chkxer('DPPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPPSVX

        srnamc.SRNAMT = 'DPPSVX';
         infoc.INFOT = 1;
         dppsvx('/', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dppsvx('N', '/', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dppsvx('N', 'U', -1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dppsvx('N', 'U', 0, -1, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         EQ = '/';
         dppsvx('F', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         EQ = 'Y';
         dppsvx('F', 'U', 1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         dppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'PB' ) ) {

         // DPBSV

        srnamc.SRNAMT = 'DPBSV ';
         infoc.INFOT = 1;
         dpbsv('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dpbsv('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dpbsv('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dpbsv('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('DPBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dpbsv('U', 1, 1, 0, A, 1, B, 2, INFO );
         chkxer('DPBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dpbsv('U', 2, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPBSVX

        srnamc.SRNAMT = 'DPBSVX';
         infoc.INFOT = 1;
         dpbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dpbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dpbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dpbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dpbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dpbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dpbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         EQ = '/';
         dpbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         EQ = 'Y';
         dpbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 15;
         dpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'PT' ) ) {

         // DPTSV

        srnamc.SRNAMT = 'DPTSV ';
         infoc.INFOT = 1;
         dptsv(-1, 0, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('DPTSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dptsv(0, -1, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('DPTSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dptsv(2, 0, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('DPTSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPTSVX

        srnamc.SRNAMT = 'DPTSVX';
         infoc.INFOT = 1;
         dptsvx('/', 0, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dptsvx('N', -1, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dptsvx('N', 0, -1, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dptsvx('N', 2, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 2, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dptsvx('N', 2, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 2, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SY' ) ) {

         // DSYSV

        srnamc.SRNAMT = 'DSYSV ';
         infoc.INFOT = 1;
         dsysv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dsysv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dsysv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dsysv('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dsysv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dsysv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('DSYSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dsysv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('DSYSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DSYSVX

        srnamc.SRNAMT = 'DSYSVX';
         infoc.INFOT = 1;
         dsysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dsysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dsysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dsysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dsysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 18;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, IW, INFO );
         chkxer('DSYSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SR' ) ) {

         // DSYSV_ROOK

        srnamc.SRNAMT = 'DSYSV_ROOK';
         infoc.INFOT = 1;
         dsysv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dsysv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dsysv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dsysv_rook('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dsysv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dsysv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dsysv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('DSYSV_ROOK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SK' ) ) {

         // DSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

        srnamc.SRNAMT = 'DSYSV_RK';
         infoc.INFOT = 1;
         dsysv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dsysv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dsysv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dsysv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dsysv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dsysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dsysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('DSYSV_RK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SA' ) ) {

         // DSYSV_AASEN

        srnamc.SRNAMT = 'DSYSV_AA';
         infoc.INFOT = 1;
         dsysv_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_AA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dsysv_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_AA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dsysv_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_AA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dsysv_aa('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_AA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dsysv_aa('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_AA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dsysv_aa('U', 3, 1, A, 3, IP, B, 3, W, 6, INFO );
         chkxer('DSYSV_AA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'S2' ) ) {

         // DSYSV_AASEN_2STAGE

        srnamc.SRNAMT = 'DSYSV_AA_2STAGE';
         infoc.INFOT = 1;
         dsysv_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1,  W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dsysv_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1,  W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dsysv_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1,  W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dsysv_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dsysv_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dsysv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dsysv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_AA_2STAGE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SP' ) ) {

         // DSPSV

        srnamc.SRNAMT = 'DSPSV ';
         infoc.INFOT = 1;
         dspsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('DSPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dspsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('DSPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dspsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('DSPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dspsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('DSPSV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DSPSVX

        srnamc.SRNAMT = 'DSPSVX';
         infoc.INFOT = 1;
         dspsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dspsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dspsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dspsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         dspsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         dspsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      if ( infoc.OK ) {
         WRITE( infoc.NOUT, FMT = 9999 )PATH;
      } else {
         WRITE( infoc.NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT(' ${.a3} drivers passed the tests of the error exits' );
 9998 FORMAT( ' *** ${.a3} drivers failed the tests of the error exits ***' );

      }
