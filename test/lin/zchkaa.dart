      void main() {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                NMAX;
      const              NMAX = 132 ;
      int                MAXIN;
      const              MAXIN = 12 ;
      int                MAXRHS;
      const              MAXRHS = 16 ;
      int                MATMAX;
      const              MATMAX = 30 ;
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      int                KDMAX;
      const              KDMAX = NMAX+( NMAX+1 ) / 4 ;
      bool               FATAL, TSTCHK, TSTDRV, TSTERR;
      String             C1;
      String             C2;
      String             PATH;
      String             INTSTR;
      String             ALINE;
      int                I, IC, J, K, LA, LAFAC, LDA, NB, NM, NMATS, NN, NNB, NNB2, NNS, NRHS, NTYPES, NRANK, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      double             EPS, S1, S2, THREQ, THRESH;
      bool               DOTYPE( MATMAX );
      int                IWORK( 25*NMAX ), MVAL( MAXIN ), NBVAL( MAXIN ), NBVAL2( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), RANKVAL( MAXIN ), PIV( NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      double          , DIMENSION(:), ALLOCATABLE::  RWORK, S;
      Complex, DIMENSION(:), ALLOCATABLE :: E;
      Complex, DIMENSION(:,:), ALLOCATABLE::  A, B, WORK;
      // ..
      // .. External Functions ..
      //- bool               lsame, LSAMEN;
      //- double             DLAMCH, DSECND;
      // EXTERNAL lsame, LSAMEN, DLAMCH, DSECND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, ZCHKEQ, ZCHKGB, ZCHKGE, ZCHKGT, ZCHKHE, ZCHKHE_ROOK, ZCHKHE_RK, ZCHKHE_AA, ZCHKHP, ZCHKLQ, ZCHKUNHR_COL, ZCHKPB, ZCHKPO, ZCHKPS, ZCHKPP, ZCHKPT, ZCHKQ3, ZCHKQP3RK, ZCHKQL, ZCHKQR, ZCHKRQ, ZCHKSP, ZCHKSY, ZCHKSY_ROOK, ZCHKSY_RK, ZCHKSY_AA, ZCHKTB, ZCHKTP, ZCHKTR, ZCHKTZ, ZDRVGB, ZDRVGE, ZDRVGT, ZDRVHE, ZDRVHE_ROOK, ZDRVHE_RK, ZDRVHE_AA, ZDRVHE_AA_2STAGE, ZDRVHP, ZDRVLS, ZDRVPB, ZDRVPO, ZDRVPP, ZDRVPT, ZDRVSP, ZDRVSY, ZDRVSY_ROOK, ZDRVSY_RK, ZDRVSY_AA, ZDRVSY_AA_2STAGE, ILAVER, ZCHKQRT, ZCHKQRTP, ZCHKLQT, ZCHKLQTP, ZCHKTSQR
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Arrays in Common ..
      // int                claenv.IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // COMMON / CLAENV / claenv.IPARMS
      // ..
      // .. Data statements ..
      const THREQ = 2.0, INTSTR = '0123456789';

      // .. Allocate memory dynamically ..

      ALLOCATE ( A ( (KDMAX+1) * NMAX, 7 ), STAT = AllocateStatus);
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( B ( NMAX * MAXRHS, 4 ), STAT = AllocateStatus);
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( WORK ( NMAX, NMAX+MAXRHS+10 ), STAT = AllocateStatus);
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( E( NMAX ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( S( 2*NMAX ), STAT = AllocateStatus);
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( RWORK( 150*NMAX+2*MAXRHS ), STAT = AllocateStatus);
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";

      S1 = DSECND( );
      LDA = NMAX;
      FATAL = false;

      // Read a dummy line.

      READ( NIN, FMT = * );

      // Report values of parameters.

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH;

      // Read the values of M

      READ( NIN, FMT = * )NM;
      if ( NM < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NM ', NM, 1;
         NM = 0;
         FATAL = true;
      } else if ( NM > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN;
         NM = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM );
      for (I = 1; I <= NM; I++) { // 10
         if ( MVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0;
            FATAL = true;
         } else if ( MVAL( I ) > NMAX ) {
            WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX;
            FATAL = true;
         }
      } // 10
      if (NM > 0) WRITE( NOUT, FMT = 9993 )'M   ', ( MVAL( I ), I = 1, NM );

      // Read the values of N

      READ( NIN, FMT = * )NN;
      if ( NN < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NN ', NN, 1;
         NN = 0;
         FATAL = true;
      } else if ( NN > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN;
         NN = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN );
      for (I = 1; I <= NN; I++) { // 20
         if ( NVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )' N  ', NVAL( I ), 0;
            FATAL = true;
         } else if ( NVAL( I ) > NMAX ) {
            WRITE( NOUT, FMT = 9995 )' N  ', NVAL( I ), NMAX;
            FATAL = true;
         }
      } // 20
      if (NN > 0) WRITE( NOUT, FMT = 9993 )'N   ', ( NVAL( I ), I = 1, NN );

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS;
      if ( NNS < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1;
         NNS = 0;
         FATAL = true;
      } else if ( NNS > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN;
         NNS = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS );
      for (I = 1; I <= NNS; I++) { // 30
         if ( NSVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0;
            FATAL = true;
         } else if ( NSVAL( I ) > MAXRHS ) {
            WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS;
            FATAL = true;
         }
      } // 30
      if (NNS > 0) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS );

      // Read the values of NB

      READ( NIN, FMT = * )NNB;
      if ( NNB < 1 ) {
         WRITE( NOUT, FMT = 9996 )'NNB ', NNB, 1;
         NNB = 0;
         FATAL = true;
      } else if ( NNB > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )'NNB ', NNB, MAXIN;
         NNB = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB );
      for (I = 1; I <= NNB; I++) { // 40
         if ( NBVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )' NB ', NBVAL( I ), 0;
            FATAL = true;
         }
      } // 40
      if (NNB > 0) WRITE( NOUT, FMT = 9993 )'NB  ', ( NBVAL( I ), I = 1, NNB );

      // Set NBVAL2 to be the set of unique values of NB

      NNB2 = 0;
      for (I = 1; I <= NNB; I++) { // 60
         NB = NBVAL( I );
         for (J = 1; J <= NNB2; J++) { // 50
            if( NB == NBVAL2( J ) ) GO TO 60;
         } // 50
         NNB2 = NNB2 + 1;
         NBVAL2[NNB2] = NB;
      } // 60

      // Read the values of NX

      READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB );
      for (I = 1; I <= NNB; I++) { // 70
         if ( NXVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )' NX ', NXVAL( I ), 0;
            FATAL = true;
         }
      } // 70
      if (NNB > 0) WRITE( NOUT, FMT = 9993 )'NX  ', ( NXVAL( I ), I = 1, NNB );

      // Read the values of RANKVAL

      READ( NIN, FMT = * )NRANK;
      if ( NN < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NRANK ', NRANK, 1;
         NRANK = 0;
         FATAL = true;
      } else if ( NN > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NRANK ', NRANK, MAXIN;
         NRANK = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( RANKVAL( I ), I = 1, NRANK );
      for (I = 1; I <= NRANK; I++) {
         if ( RANKVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )' RANK  ', RANKVAL( I ), 0;
            FATAL = true;
         } else if ( RANKVAL( I ) > 100 ) {
            WRITE( NOUT, FMT = 9995 )' RANK  ', RANKVAL( I ), 100;
            FATAL = true;
         }
      }
      if (NRANK > 0) WRITE( NOUT, FMT = 9993 )'RANK % OF N', ( RANKVAL( I ), I = 1, NRANK );

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH;
      WRITE( NOUT, FMT = 9992 )THRESH;

      // Read the flag that indicates whether to test the LAPACK routines.

      READ( NIN, FMT = * )TSTCHK;

      // Read the flag that indicates whether to test the driver routines.

      READ( NIN, FMT = * )TSTDRV;

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR;

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 );
         STOP;
      }

      // Calculate and print the machine dependent constants.

      EPS = dlamch( 'Underflow threshold' );
      WRITE( NOUT, FMT = 9991 )'underflow', EPS;
      EPS = dlamch( 'Overflow threshold' );
      WRITE( NOUT, FMT = 9991 )'overflow ', EPS;
      EPS = dlamch( 'Epsilon' );
      WRITE( NOUT, FMT = 9991 )'precision', EPS;
      WRITE( NOUT, FMT = * );
      NRHS = NSVAL( 1 );

      } // 80

      // Read a test path and the number of matrix types to use.

      READ( NIN, FMT = '(A72)', END = 140 )ALINE;
      PATH = ALINE( 1: 3 );
      NMATS = MATMAX;
      I = 3;
      } // 90
      I = I + 1;
      if (I > 72) GO TO 130;
      IF( ALINE( I: I ) == ' ' ) GO TO 90;
      NMATS = 0;
      } // 100
      C1 = ALINE( I: I );
      for (K = 1; K <= 10; K++) { // 110
         if ( C1 == INTSTR( K: K ) ) {
            IC = K - 1;
            GO TO 120;
         }
      } // 110
      GO TO 130;
      } // 120
      NMATS = NMATS*10 + IC;
      I = I + 1;
      if (I > 72) GO TO 130;
      GO TO 100;
      } // 130
      C1 = PATH( 1: 1 );
      C2 = PATH( 2: 3 );

      // Check first character for correct precision.

      if ( !lsame( C1, 'Zomplex precision' ) ) {
         WRITE( NOUT, FMT = 9990 )PATH;

      } else if ( NMATS <= 0 ) {

         // Check for a positive number of tests requested.

         WRITE( NOUT, FMT = 9989 )PATH;

      } else if ( lsamen( 2, C2, 'GE' ) ) {

         // GE:  general matrices

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkge(DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvge(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // GB:  general banded matrices

         LA = ( 2*KDMAX+1 )*NMAX;
         LAFAC = ( 3*KDMAX+1 )*NMAX;
         NTYPES = 8;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkgb(DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), LA, A( 1, 3 ), LAFAC, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvgb(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), LA, A( 1, 3 ), LAFAC, A( 1, 6 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'GT' ) ) {

         // GT:  general tridiagonal matrices

         NTYPES = 12;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkgt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvgt(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'PO' ) ) {

         // PO:  positive definite matrices

         NTYPES = 9;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkpo(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvpo(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'PS' ) ) {

         // PS:  positive semi-definite matrices

         NTYPES = 9;

         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkps(DOTYPE, NN, NVAL, NNB2, NBVAL2, NRANK, RANKVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), PIV, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'PP' ) ) {

         // PP:  positive definite packed matrices

         NTYPES = 9;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkpp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvpp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'PB' ) ) {

         // PB:  positive definite banded matrices

         NTYPES = 8;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkpb(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvpb(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'PT' ) ) {

         // PT:  positive definite tridiagonal matrices

         NTYPES = 12;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkpt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvpt(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'HE' ) ) {

         // HE:  Hermitian indefinite matrices

         NTYPES = 10;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkhe(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvhe(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'HR' ) ) {

         // HR:  Hermitian indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm,

         NTYPES = 10;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkhe_rook(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvhe_rook(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'HK' ) ) {

         // HK:  Hermitian indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm,
              // different matrix storage format than HR path version.

         NTYPES = 10;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkhe_rk(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvhe_rk(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'HA' ) ) {

         // HA:  Hermitian matrices,
              // Aasen Algorithm

         NTYPES = 10;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkhe_aa(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvhe_aa(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'H2' ) ) {

         // H2:  Hermitian matrices,
              // with partial (Aasen's) pivoting algorithm

         NTYPES = 10;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkhe_aa_2stage(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvhe_aa_2stage(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }


      } else if ( lsamen( 2, C2, 'HP' ) ) {

         // HP:  Hermitian indefinite packed matrices

         NTYPES = 10;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkhp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvhp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'SY' ) ) {

         // SY:  symmetric indefinite matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchksy(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvsy(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'SR' ) ) {

         // SR:  symmetric indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchksy_rook(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvsy_rook(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'SK' ) ) {

         // SK:  symmetric indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm,
              // different matrix storage format than SR path version.

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchksy_rk(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvsy_rk(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'SA' ) ) {

         // SA:  symmetric indefinite matrices with Aasen's algorithm,

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchksy_aa(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvsy_aa(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'S2' ) ) {

         // S2:  symmetric indefinite matrices with Aasen's algorithm
              // 2 stage

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchksy_aa_2stage(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvsy_aa_2stage(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'SP' ) ) {

         // SP:  symmetric indefinite packed matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 11;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchksp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

         if ( TSTDRV ) {
            zdrvsp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TR' ) ) {

         // TR:  triangular matrices

         NTYPES = 18;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchktr(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // TP:  triangular packed matrices

         NTYPES = 18;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchktp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // TB:  triangular banded matrices

         NTYPES = 17;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchktb(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'QR' ) ) {

         // QR:  QR factorization

         NTYPES = 8;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkqr(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'LQ' ) ) {

         // LQ:  LQ factorization

         NTYPES = 8;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchklq(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'QL' ) ) {

         // QL:  QL factorization

         NTYPES = 8;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkql(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'RQ' ) ) {

         // RQ:  RQ factorization

         NTYPES = 8;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkrq(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'EQ' ) ) {

         // EQ:  Equilibration routines for general and positive definite
              // matrices (THREQ should be between 2 and 10)

         if ( TSTCHK ) {
            zchkeq(THREQ, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TZ' ) ) {

         // TZ:  Trapezoidal matrix

         NTYPES = 3;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchktz(DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), S( 1 ), B( 1, 1 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'QP' ) ) {

         // QP:  QR factorization with pivoting

         NTYPES = 6;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkq3(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, THRESH, A( 1, 1 ), A( 1, 2 ), S( 1 ), B( 1, 1 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'QK' ) ) {

         // QK: truncated QR factorization with pivoting

         NTYPES = 19;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            zchkqp3rk(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), S( 1 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'LS' ) ) {

         // LS:  Least squares drivers

         NTYPES = 6;
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTDRV ) {
            zdrvls(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), S( 1 ), S( NMAX+1 ), NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }


      } else if ( lsamen( 2, C2, 'QT' ) ) {

         // QT:  QRT routines for general matrices

         if ( TSTCHK ) {
            zchkqrt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'QX' ) ) {

         // QX:  QRT routines for triangular-pentagonal matrices

         if ( TSTCHK ) {
            zchkqrtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TQ' ) ) {

         // TQ:  LQT routines for general matrices

         if ( TSTCHK ) {
            zchklqt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'XQ' ) ) {

         // XQ:  LQT routines for triangular-pentagonal matrices

         if ( TSTCHK ) {
            zchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TS' ) ) {

         // TS:  QR routines for tall-skinny matrices

         if ( TSTCHK ) {
            zchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TQ' ) ) {

         // TQ:  LQT routines for general matrices

         if ( TSTCHK ) {
            zchklqt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'XQ' ) ) {

         // XQ:  LQT routines for triangular-pentagonal matrices

         if ( TSTCHK ) {
            zchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'TS' ) ) {

         // TS:  QR routines for tall-skinny matrices

         if ( TSTCHK ) {
            zchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH;
         }

      } else if ( lsamen( 2, C2, 'HH' ) ) {

         // HH:  Householder reconstruction for tall-skinny matrices

         if ( TSTCHK ) {
            zchkunhr_col(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 ) PATH;
         }

      } else {

         WRITE( NOUT, FMT = 9990 )PATH;
      }

      // Go back to get another input line.

      GO TO 80;

      // Branch to this line when the last record is read.

      } // 140
      CLOSE ( NIN );
      S2 = DSECND( );
      WRITE( NOUT, FMT = 9998 );
      WRITE( NOUT, FMT = 9997 )S2 - S1;

      DEALLOCATE (A, STAT = AllocateStatus);
      DEALLOCATE (B, STAT = AllocateStatus);
      DEALLOCATE (E, STAT = AllocateStatus);
      DEALLOCATE (S, STAT = AllocateStatus);
      DEALLOCATE (RWORK, STAT = AllocateStatus);
      DEALLOCATE (WORK,  STAT = AllocateStatus);

 9999 FORMAT('\n Execution not attempted due to input errors' );
 9998 FORMAT('\n End of tests' );
 9997 FORMAT( ' Total time used = ${.f12_2} seconds', / );
 9996 FORMAT( ' Invalid input value: ${.a4}=${.i6}; must be >=', I6 )
 9995 FORMAT( ' Invalid input value: ${.a4}=${.i6}; must be <=', I6 )
 9994 FORMAT( ' Tests of the Complex LAPACK routines \n LAPACK VERSION ${.i1}.${.i1}.', I1, / / ' The following parameter values will be used:' );
 9993 FORMAT('    ${.a4}:  ', 10I6, / 11X, 10I6 );
 9992 FORMAT('\n Routines pass computational tests if test ratio is less than', F8.2, / );
 9991 FORMAT( ' Relative machine ${} is taken to be', D16.6 );
 9990 FORMAT( / 1X, '${.a3}:  Unrecognized path name' );
 9989 FORMAT( / 1X, '${.a3} routines were not tested' );
 9988 FORMAT( / 1X, '${.a3} driver routines were not tested' );
      }
