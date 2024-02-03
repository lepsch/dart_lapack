      PROGRAM CCHKAA

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
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
      // ..
      // .. Local Scalars ..
      bool               FATAL, TSTCHK, TSTDRV, TSTERR;
      String             C1;
      String             C2;
      String             PATH;
      String             INTSTR;
      String             ALINE;
      int                I, IC, J, K, LA, LAFAC, LDA, NB, NM, NMATS, NN, NNB, NNB2, NNS, NRHS, NTYPES, NRANK, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      REAL               EPS, S1, S2, THREQ, THRESH
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MATMAX );
      int                IWORK( 25*NMAX ), MVAL( MAXIN ), NBVAL( MAXIN ), NBVAL2( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), RANKVAL( MAXIN ), PIV( NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      REAL, DIMENSION(:), ALLOCATABLE :: RWORK, S
      COMPLEX, DIMENSION(:), ALLOCATABLE :: E
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: A, B, WORK
      // ..
      // .. External Functions ..
      bool               LSAME, LSAMEN;
      REAL               SECOND, SLAMCH
      // EXTERNAL LSAME, LSAMEN, SECOND, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, CCHKEQ, CCHKGB, CCHKGE, CCHKGT, CCHKHE, CCHKHE_ROOK, CCHKHE_RK, CCHKHE_AA, CCHKHP, CCHKLQ, CCHKUNHR_COL, CCHKPB, CCHKPO, CCHKPS, CCHKPP, CCHKPT, CCHKQ3, CCHKQP3RK, CCHKQL, CCHKQR, CCHKRQ, CCHKSP, CCHKSY, CCHKSY_ROOK, CCHKSY_RK, CCHKSY_AA, CCHKTB, CCHKTP, CCHKTR, CCHKTZ, CDRVGB, CDRVGE, CDRVGT, CDRVHE, CDRVHE_ROOK, CDRVHE_RK, CDRVHE_AA, CDRVHP, CDRVLS, CDRVPB, CDRVPO, CDRVPP, CDRVPT, CDRVSP, CDRVSY, CDRVSY_ROOK, CDRVSY_RK, CDRVSY_AA, ILAVER, CCHKQRT, CCHKQRTP
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Arrays in Common ..
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / IPARMS
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               THREQ / 2.0 / , INTSTR / '0123456789' /
      // ..
      // .. Allocate memory dynamically ..

      ALLOCATE ( A( ( KDMAX+1 )*NMAX, 7 ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( B( NMAX*MAXRHS, 4 ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( WORK( NMAX, NMAX+MAXRHS+10 ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( E( NMAX ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( S( 2*NMAX ), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( RWORK( 150*NMAX+2*MAXRHS ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      // ..
      // .. Executable Statements ..

      S1 = SECOND( )
      LDA = NMAX
      FATAL = .FALSE.

      // Read a dummy line.

      READ( NIN, FMT = * )

      // Report values of parameters.

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH

      // Read the values of M

      READ( NIN, FMT = * )NM
      if ( NM.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NM ', NM, 1
         NM = 0
         FATAL = .TRUE.
      } else if ( NM.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN
         NM = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      for (I = 1; I <= NM; I++) { // 10
         if ( MVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0
            FATAL = .TRUE.
         } else if ( MVAL( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX
            FATAL = .TRUE.
         }
      } // 10
      IF( NM.GT.0 ) WRITE( NOUT, FMT = 9993 )'M   ', ( MVAL( I ), I = 1, NM )

      // Read the values of N

      READ( NIN, FMT = * )NN
      if ( NN.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      } else if ( NN.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      for (I = 1; I <= NN; I++) { // 20
         if ( NVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' N  ', NVAL( I ), 0
            FATAL = .TRUE.
         } else if ( NVAL( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9995 )' N  ', NVAL( I ), NMAX
            FATAL = .TRUE.
         }
      } // 20
      IF( NN.GT.0 ) WRITE( NOUT, FMT = 9993 )'N   ', ( NVAL( I ), I = 1, NN )

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS
      if ( NNS.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1
         NNS = 0
         FATAL = .TRUE.
      } else if ( NNS.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN
         NNS = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS )
      for (I = 1; I <= NNS; I++) { // 30
         if ( NSVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0
            FATAL = .TRUE.
         } else if ( NSVAL( I ).GT.MAXRHS ) {
            WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS
            FATAL = .TRUE.
         }
      } // 30
      IF( NNS.GT.0 ) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS )

      // Read the values of NB

      READ( NIN, FMT = * )NNB
      if ( NNB.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )'NNB ', NNB, 1
         NNB = 0
         FATAL = .TRUE.
      } else if ( NNB.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )'NNB ', NNB, MAXIN
         NNB = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB )
      for (I = 1; I <= NNB; I++) { // 40
         if ( NBVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' NB ', NBVAL( I ), 0
            FATAL = .TRUE.
         }
      } // 40
      IF( NNB.GT.0 ) WRITE( NOUT, FMT = 9993 )'NB  ', ( NBVAL( I ), I = 1, NNB )

      // Set NBVAL2 to be the set of unique values of NB

      NNB2 = 0
      for (I = 1; I <= NNB; I++) { // 60
         NB = NBVAL( I )
         for (J = 1; J <= NNB2; J++) { // 50
            IF( NB.EQ.NBVAL2( J ) ) GO TO 60
         } // 50
         NNB2 = NNB2 + 1
         NBVAL2( NNB2 ) = NB
      } // 60

      // Read the values of NX

      READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB )
      for (I = 1; I <= NNB; I++) { // 70
         if ( NXVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' NX ', NXVAL( I ), 0
            FATAL = .TRUE.
         }
      } // 70
      IF( NNB.GT.0 ) WRITE( NOUT, FMT = 9993 )'NX  ', ( NXVAL( I ), I = 1, NNB )

      // Read the values of RANKVAL

      READ( NIN, FMT = * )NRANK
      if ( NN.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NRANK ', NRANK, 1
         NRANK = 0
         FATAL = .TRUE.
      } else if ( NN.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NRANK ', NRANK, MAXIN
         NRANK = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( RANKVAL( I ), I = 1, NRANK )
      for (I = 1; I <= NRANK; I++) {
         if ( RANKVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' RANK  ', RANKVAL( I ), 0
            FATAL = .TRUE.
         } else if ( RANKVAL( I ).GT.100 ) {
            WRITE( NOUT, FMT = 9995 )' RANK  ', RANKVAL( I ), 100
            FATAL = .TRUE.
         }
      }
      IF( NRANK.GT.0 ) WRITE( NOUT, FMT = 9993 )'RANK % OF N', ( RANKVAL( I ), I = 1, NRANK )

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9992 )THRESH

      // Read the flag that indicates whether to test the LAPACK routines.

      READ( NIN, FMT = * )TSTCHK

      // Read the flag that indicates whether to test the driver routines.

      READ( NIN, FMT = * )TSTDRV

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 )
         STOP
      }

      // Calculate and print the machine dependent constants.

      EPS = SLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9991 )'underflow', EPS
      EPS = SLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9991 )'overflow ', EPS
      EPS = SLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9991 )'precision', EPS
      WRITE( NOUT, FMT = * )
      NRHS = NSVAL( 1 )

      } // 80

      // Read a test path and the number of matrix types to use.

      READ( NIN, FMT = '(A72)', END = 140 )ALINE
      PATH = ALINE( 1: 3 )
      NMATS = MATMAX
      I = 3
      } // 90
      I = I + 1
      IF( I.GT.72 ) GO TO 130       IF( ALINE( I: I ).EQ.' ' ) GO TO 90
      NMATS = 0
      } // 100
      C1 = ALINE( I: I )
      for (K = 1; K <= 10; K++) { // 110
         if ( C1.EQ.INTSTR( K: K ) ) {
            IC = K - 1
            GO TO 120
         }
      } // 110
      GO TO 130
      } // 120
      NMATS = NMATS*10 + IC
      I = I + 1
      IF( I.GT.72 ) GO TO 130
      GO TO 100
      } // 130
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )

      // Check first character for correct precision.

      if ( .NOT.LSAME( C1, 'Complex precision' ) ) {
         WRITE( NOUT, FMT = 9990 )PATH

      } else if ( NMATS.LE.0 ) {

         // Check for a positive number of tests requested.

         WRITE( NOUT, FMT = 9989 )PATH

      } else if ( LSAMEN( 2, C2, 'GE' ) ) {

         // GE:  general matrices

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkge(DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvge(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // GB:  general banded matrices

         LA = ( 2*KDMAX+1 )*NMAX
         LAFAC = ( 3*KDMAX+1 )*NMAX
         NTYPES = 8
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkgb(DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), LA, A( 1, 3 ), LAFAC, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvgb(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), LA, A( 1, 3 ), LAFAC, A( 1, 6 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // GT:  general tridiagonal matrices

         NTYPES = 12
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkgt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvgt(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // PO:  positive definite matrices

         NTYPES = 9
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkpo(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvpo(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'PS' ) ) {

         // PS:  positive semi-definite matrices

         NTYPES = 9

         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkps(DOTYPE, NN, NVAL, NNB2, NBVAL2, NRANK, RANKVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), PIV, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // PP:  positive definite packed matrices

         NTYPES = 9
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkpp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvpp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // PB:  positive definite banded matrices

         NTYPES = 8
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkpb(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvpb(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // PT:  positive definite tridiagonal matrices

         NTYPES = 12
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkpt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvpt(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), S, A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'HE' ) ) {

         // HE:  Hermitian indefinite matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 10
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkhe(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvhe(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // HR:  Hermitian indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm

         NTYPES = 10
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkhe_rook(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvhe_rook(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'HK' ) ) {

         // HK:  Hermitian indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm,
              // different matrix storage format than HR path version.

         NTYPES = 10
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkhe_rk(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvhe_rk(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'HA' ) ) {

         // HA:  Hermitian matrices,
              // Aasen Algorithm

         NTYPES = 10
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkhe_aa(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvhe_aa(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'H2' ) ) {

         // H2:  Hermitian matrices,
              // with partial (Aasen's) pivoting algorithm

         NTYPES = 10
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkhe_aa_2stage(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvhe_aa_2stage(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'HP' ) ) {

         // HP:  Hermitian indefinite packed matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 10
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkhp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvhp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // SY:  symmetric indefinite matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchksy(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvsy(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // SR:  symmetric indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchksy_rook(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvsy_rook(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // SK:  symmetric indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm,
              // different matrix storage format than SR path version.

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchksy_rk(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvsy_rk(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'SA' ) ) {

         // SA:  symmetric indefinite matrices with Aasen's algorithm,

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchksy_aa(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvsy_aa(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // S2:  symmetric indefinite matrices with Aasen's algorithm
              // 2 stage

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchksy_aa_2stage(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvsy_aa_2stage(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // SP:  symmetric indefinite packed matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 11
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchksp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

         if ( TSTDRV ) {
            cdrvsp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9988 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'TR' ) ) {

         // TR:  triangular matrices

         NTYPES = 18
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchktr(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'TP' ) ) {

         // TP:  triangular packed matrices

         NTYPES = 18
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchktp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'TB' ) ) {

         // TB:  triangular banded matrices

         NTYPES = 17
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchktb(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'QR' ) ) {

         // QR:  QR factorization

         NTYPES = 8
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkqr(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'LQ' ) ) {

         // LQ:  LQ factorization

         NTYPES = 8
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchklq(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'QL' ) ) {

         // QL:  QL factorization

         NTYPES = 8
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkql(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'RQ' ) ) {

         // RQ:  RQ factorization

         NTYPES = 8
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkrq(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'EQ' ) ) {

         // EQ:  Equilibration routines for general and positive definite
              // matrices (THREQ should be between 2 and 10)

         if ( TSTCHK ) {
            cchkeq(THREQ, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'TZ' ) ) {

         // TZ:  Trapezoidal matrix

         NTYPES = 3
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchktz(DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), S( 1 ), B( 1, 1 ), WORK, RWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'QP' ) ) {

         // QP:  QR factorization with pivoting

         NTYPES = 6
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkq3(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, THRESH, A( 1, 1 ), A( 1, 2 ), S( 1 ), B( 1, 1 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'QK' ) ) {

         // QK: truncated QR factorization with pivoting

         NTYPES = 19
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTCHK ) {
            cchkqp3rk(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), S( 1 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'LS' ) ) {

         // LS:  Least squares drivers

         NTYPES = 6
         alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if ( TSTDRV ) {
            cdrvls(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), S( 1 ), S( NMAX+1 ), NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'QT' ) ) {

         // QT:  QRT routines for general matrices

         if ( TSTCHK ) {
            cchkqrt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'QX' ) ) {

         // QX:  QRT routines for triangular-pentagonal matrices

         if ( TSTCHK ) {
            cchkqrtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'TQ' ) ) {

         // TQ:  LQT routines for general matrices

         if ( TSTCHK ) {
            cchklqt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'XQ' ) ) {

         // XQ:  LQT routines for triangular-pentagonal matrices

         if ( TSTCHK ) {
            cchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'TS' ) ) {

         // TS:  QR routines for tall-skinny matrices

         if ( TSTCHK ) {
            cchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }

      } else if ( LSAMEN( 2, C2, 'HH' ) ) {

         // HH:  Householder reconstruction for tall-skinny matrices

         if ( TSTCHK ) {
            cchkunhr_col(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 ) PATH
         }

      } else {

         WRITE( NOUT, FMT = 9990 )PATH
      }

      // Go back to get another input line.

      GO TO 80

      // Branch to this line when the last record is read.

      } // 140
      CLOSE ( NIN )
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9998 )
      WRITE( NOUT, FMT = 9997 )S2 - S1

      DEALLOCATE (A, STAT = AllocateStatus)
      DEALLOCATE (B, STAT = AllocateStatus)
      DEALLOCATE (E, STAT = AllocateStatus)
      DEALLOCATE (S, STAT = AllocateStatus)
      DEALLOCATE (WORK, STAT = AllocateStatus)
      DEALLOCATE (RWORK,  STAT = AllocateStatus)

 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=', I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=', I6 )
 9994 FORMAT( ' Tests of the COMPLEX LAPACK routines ', / ' LAPACK VERSION ', I1, '.', I1, '.', I1, / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ', 'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 )
 9990 FORMAT( / 1X, A3, ':  Unrecognized path name' )
 9989 FORMAT( / 1X, A3, ' routines were not tested' )
 9988 FORMAT( / 1X, A3, ' driver routines were not tested' )

      // End of CCHKAA

      }
