      PROGRAM CCHKEE

#if defined(_OPENMP)
      use omp_lib
#endif

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 132 ;
      int                NCMAX;
      const              NCMAX = 20 ;
      int                NEED;
      const              NEED = 14 ;
      int                LWORK;
      const              LWORK = NMAX*( 5*NMAX+20 ) ;
      int                LIWORK;
      const              LIWORK = NMAX*( NMAX+20 ) ;
      int                MAXIN;
      const              MAXIN = 20 ;
      int                MAXT;
      const              MAXT = 30 ;
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      // ..
      // .. Local Scalars ..
      bool               CBB, CBK, CBL, CES, CEV, CGG, CGK, CGL, CGS, CGV, CGX, CHB, CSD, CSX, CVX, CXV, FATAL, GLM, GQR, GSV, LSE, NEP, SEP, SVD, TSTCHK, TSTDIF, TSTDRV, TSTERR;
      String             C1;
      String             C3, PATH;
      String             VNAME;
      String             INTSTR;
      String             LINE;
      int                I, I1, IC, INFO, ITMP, K, LENP, MAXTYP, NEWSD, NK, NN, NPARMS, NRHS, NTYPES, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      int    *4          N_THREADS, ONE_THREAD;
      REAL               EPS, S1, S2, THRESH, THRSHN
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MAXT ), LOGWRK( NMAX );
      int                IOLDSD( 4 ), ISEED( 4 ), IWORK( LIWORK ), KVAL( MAXIN ), MVAL( MAXIN ), MXBVAL( MAXIN ), NBCOL( MAXIN ), NBMIN( MAXIN ), NBVAL( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), PVAL( MAXIN );
      int                INMIN( MAXIN ), INWIN( MAXIN ), INIBL( MAXIN ), ISHFTS( MAXIN ), IACC22( MAXIN );
      REAL               ALPHA( NMAX ), BETA( NMAX ), DR( NMAX, 12 ), RESULT( 500 );
      COMPLEX            DC( NMAX, 6 ), TAUA( NMAX ), TAUB( NMAX ), X( 5*NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      REAL, DIMENSION(:), ALLOCATABLE :: RWORK, S
      COMPLEX, DIMENSION(:), ALLOCATABLE :: WORK
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: A, B, C
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      REAL               SECOND, SLAMCH
      // EXTERNAL LSAMEN, SECOND, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, CCHKBB, CCHKBD, CCHKBK, CCHKBL, CCHKEC, CCHKGG, CCHKGK, CCHKGL, CCHKHB, CCHKHS, CCHKST, CCKCSD, CCKGLM, CCKGQR, CCKGSV, CCKLSE, CDRGES, CDRGEV, CDRGSX, CDRGVX, CDRVBD, CDRVES, CDRVEV, CDRVSG, CDRVST, CDRVSX, CDRVVX, CERRBD, CERRED, CERRGG, CERRHS, CERRST, ILAVER, XLAENV, CDRGES3, CDRGEV3, CCHKST2STG, CDRVST2STG, CCHKHB2STG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, MAXB, NPROC, NSHIFT, NUNIT, SELDIM, SELOPT;
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      int                IPARMS( 100 );
      REAL               SELWI( 20 ), SELWR( 20 )
      // ..
      // .. Common blocks ..
      COMMON             / CENVIR / NPROC, NSHIFT, MAXB
      COMMON             / CLAENV / IPARMS
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Data statements ..
      DATA               INTSTR / '0123456789' /
      DATA               IOLDSD / 0, 0, 0, 1 /
      // ..
      // .. Allocate memory dynamically ..

      ALLOCATE ( S(NMAX*NMAX), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( A(NMAX*NMAX,NEED), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( B(NMAX*NMAX,5), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( C(NCMAX*NCMAX,NCMAX*NCMAX), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( RWORK(LWORK), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( WORK(LWORK), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      // ..
      // .. Executable Statements ..

      A = 0.0
      B = 0.0
      C = 0.0
      DC = 0.0
      S1 = SECOND( )
      FATAL = .FALSE.
      NUNIT = NOUT

      // Return to here to read multiple sets of data

      } // 10

      // Read the first line and set the 3-character test path

      READ( NIN, FMT = '(A80)', END = 380 )LINE
      PATH = LINE( 1: 3 )
      NEP = LSAMEN( 3, PATH, 'NEP' ) .OR. LSAMEN( 3, PATH, 'CHS' )
      SEP = LSAMEN( 3, PATH, 'SEP' ) .OR. LSAMEN( 3, PATH, 'CST' ) .OR. LSAMEN( 3, PATH, 'CSG' ) .OR. LSAMEN( 3, PATH, 'SE2' )
      SVD = LSAMEN( 3, PATH, 'SVD' ) .OR. LSAMEN( 3, PATH, 'CBD' )
      CEV = LSAMEN( 3, PATH, 'CEV' )
      CES = LSAMEN( 3, PATH, 'CES' )
      CVX = LSAMEN( 3, PATH, 'CVX' )
      CSX = LSAMEN( 3, PATH, 'CSX' )
      CGG = LSAMEN( 3, PATH, 'CGG' )
      CGS = LSAMEN( 3, PATH, 'CGS' )
      CGX = LSAMEN( 3, PATH, 'CGX' )
      CGV = LSAMEN( 3, PATH, 'CGV' )
      CXV = LSAMEN( 3, PATH, 'CXV' )
      CHB = LSAMEN( 3, PATH, 'CHB' )
      CBB = LSAMEN( 3, PATH, 'CBB' )
      GLM = LSAMEN( 3, PATH, 'GLM' )
      GQR = LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' )
      GSV = LSAMEN( 3, PATH, 'GSV' )
      CSD = LSAMEN( 3, PATH, 'CSD' )
      LSE = LSAMEN( 3, PATH, 'LSE' )
      CBL = LSAMEN( 3, PATH, 'CBL' )
      CBK = LSAMEN( 3, PATH, 'CBK' )
      CGL = LSAMEN( 3, PATH, 'CGL' )
      CGK = LSAMEN( 3, PATH, 'CGK' )

      // Report values of parameters.

      if ( PATH.EQ.'   ' ) {
         GO TO 10
      } else if ( NEP ) {
         WRITE( NOUT, FMT = 9987 )
      } else if ( SEP ) {
         WRITE( NOUT, FMT = 9986 )
      } else if ( SVD ) {
         WRITE( NOUT, FMT = 9985 )
      } else if ( CEV ) {
         WRITE( NOUT, FMT = 9979 )
      } else if ( CES ) {
         WRITE( NOUT, FMT = 9978 )
      } else if ( CVX ) {
         WRITE( NOUT, FMT = 9977 )
      } else if ( CSX ) {
         WRITE( NOUT, FMT = 9976 )
      } else if ( CGG ) {
         WRITE( NOUT, FMT = 9975 )
      } else if ( CGS ) {
         WRITE( NOUT, FMT = 9964 )
      } else if ( CGX ) {
         WRITE( NOUT, FMT = 9965 )
      } else if ( CGV ) {
         WRITE( NOUT, FMT = 9963 )
      } else if ( CXV ) {
         WRITE( NOUT, FMT = 9962 )
      } else if ( CHB ) {
         WRITE( NOUT, FMT = 9974 )
      } else if ( CBB ) {
         WRITE( NOUT, FMT = 9967 )
      } else if ( GLM ) {
         WRITE( NOUT, FMT = 9971 )
      } else if ( GQR ) {
         WRITE( NOUT, FMT = 9970 )
      } else if ( GSV ) {
         WRITE( NOUT, FMT = 9969 )
      } else if ( CSD ) {
         WRITE( NOUT, FMT = 9960 )
      } else if ( LSE ) {
         WRITE( NOUT, FMT = 9968 )
      } else if ( CBL ) {

         // CGEBAL:  Balancing

         cchkbl(NIN, NOUT );
         GO TO 380
      } else if ( CBK ) {

         // CGEBAK:  Back transformation

         cchkbk(NIN, NOUT );
         GO TO 380
      } else if ( CGL ) {

         // CGGBAL:  Balancing

         cchkgl(NIN, NOUT );
         GO TO 380
      } else if ( CGK ) {

         // CGGBAK:  Back transformation

         cchkgk(NIN, NOUT );
         GO TO 380
      } else if ( LSAMEN( 3, PATH, 'CEC' ) ) {

         // CEC:  Eigencondition estimation

         READ( NIN, FMT = * )THRESH
         xlaenv(1, 1 );
         xlaenv(12, 1 );
         TSTERR = .TRUE.
         cchkec(THRESH, TSTERR, NIN, NOUT );
         GO TO 380
      } else {
         WRITE( NOUT, FMT = 9992 )PATH
         GO TO 380
      }
      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9972 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH
      WRITE( NOUT, FMT = 9984 )

      // Read the number of values of M, P, and N.

      READ( NIN, FMT = * )NN
      if ( NN.LT.0 ) {
         WRITE( NOUT, FMT = 9989 )'   NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      } else if ( NN.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9988 )'   NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      }

      // Read the values of M

      if ( .NOT.( CGX .OR. CXV ) ) {
         READ( NIN, FMT = * )( MVAL( I ), I = 1, NN )
         if ( SVD ) {
            VNAME = '    M '
         } else {
            VNAME = '    N '
         }
         for (I = 1; I <= NN; I++) { // 20
            if ( MVAL( I ).LT.0 ) {
               WRITE( NOUT, FMT = 9989 )VNAME, MVAL( I ), 0
               FATAL = .TRUE.
            } else if ( MVAL( I ).GT.NMAX ) {
               WRITE( NOUT, FMT = 9988 )VNAME, MVAL( I ), NMAX
               FATAL = .TRUE.
            }
         } // 20
         WRITE( NOUT, FMT = 9983 )'M:    ', ( MVAL( I ), I = 1, NN )
      }

      // Read the values of P

      if ( GLM .OR. GQR .OR. GSV .OR. CSD .OR. LSE ) {
         READ( NIN, FMT = * )( PVAL( I ), I = 1, NN )
         for (I = 1; I <= NN; I++) { // 30
            if ( PVAL( I ).LT.0 ) {
               WRITE( NOUT, FMT = 9989 )' P  ', PVAL( I ), 0
               FATAL = .TRUE.
            } else if ( PVAL( I ).GT.NMAX ) {
               WRITE( NOUT, FMT = 9988 )' P  ', PVAL( I ), NMAX
               FATAL = .TRUE.
            }
         } // 30
         WRITE( NOUT, FMT = 9983 )'P:    ', ( PVAL( I ), I = 1, NN )
      }

      // Read the values of N

      if ( SVD .OR. CBB .OR. GLM .OR. GQR .OR. GSV .OR. CSD .OR. LSE ) {
         READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
         for (I = 1; I <= NN; I++) { // 40
            if ( NVAL( I ).LT.0 ) {
               WRITE( NOUT, FMT = 9989 )'    N ', NVAL( I ), 0
               FATAL = .TRUE.
            } else if ( NVAL( I ).GT.NMAX ) {
               WRITE( NOUT, FMT = 9988 )'    N ', NVAL( I ), NMAX
               FATAL = .TRUE.
            }
         } // 40
      } else {
         for (I = 1; I <= NN; I++) { // 50
            NVAL( I ) = MVAL( I )
         } // 50
      }
      if ( .NOT.( CGX .OR. CXV ) ) {
         WRITE( NOUT, FMT = 9983 )'N:    ', ( NVAL( I ), I = 1, NN )
      } else {
         WRITE( NOUT, FMT = 9983 )'N:    ', NN
      }

      // Read the number of values of K, followed by the values of K

      if ( CHB .OR. CBB ) {
         READ( NIN, FMT = * )NK
         READ( NIN, FMT = * )( KVAL( I ), I = 1, NK )
         for (I = 1; I <= NK; I++) { // 60
            if ( KVAL( I ).LT.0 ) {
               WRITE( NOUT, FMT = 9989 )'    K ', KVAL( I ), 0
               FATAL = .TRUE.
            } else if ( KVAL( I ).GT.NMAX ) {
               WRITE( NOUT, FMT = 9988 )'    K ', KVAL( I ), NMAX
               FATAL = .TRUE.
            }
         } // 60
         WRITE( NOUT, FMT = 9983 )'K:    ', ( KVAL( I ), I = 1, NK )
      }

      if ( CEV .OR. CES .OR. CVX .OR. CSX ) {

         // For the nonsymmetric QR driver routines, only one set of
         // parameters is allowed.

         READ( NIN, FMT = * )NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), INMIN( 1 ), INWIN( 1 ), INIBL(1), ISHFTS(1), IACC22(1)
         if ( NBVAL( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1
            FATAL = .TRUE.
         } else if ( NBMIN( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1
            FATAL = .TRUE.
         } else if ( NXVAL( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1
            FATAL = .TRUE.
         } else if ( INMIN( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   INMIN ', INMIN( 1 ), 1
            FATAL = .TRUE.
         } else if ( INWIN( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   INWIN ', INWIN( 1 ), 1
            FATAL = .TRUE.
         } else if ( INIBL( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   INIBL ', INIBL( 1 ), 1
            FATAL = .TRUE.
         } else if ( ISHFTS( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   ISHFTS ', ISHFTS( 1 ), 1
            FATAL = .TRUE.
         } else if ( IACC22( 1 ).LT.0 ) {
            WRITE( NOUT, FMT = 9989 )'   IACC22 ', IACC22( 1 ), 0
            FATAL = .TRUE.
         }
         xlaenv(1, NBVAL( 1 ) );
         xlaenv(2, NBMIN( 1 ) );
         xlaenv(3, NXVAL( 1 ) );
         xlaenv(12, MAX( 11, INMIN( 1 ) ) );
         xlaenv(13, INWIN( 1 ) );
         xlaenv(14, INIBL( 1 ) );
         xlaenv(15, ISHFTS( 1 ) );
         xlaenv(16, IACC22( 1 ) );
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 )
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'INMIN:   ', INMIN( 1 )
         WRITE( NOUT, FMT = 9983 )'INWIN: ', INWIN( 1 )
         WRITE( NOUT, FMT = 9983 )'INIBL: ', INIBL( 1 )
         WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ISHFTS( 1 )
         WRITE( NOUT, FMT = 9983 )'IACC22: ', IACC22( 1 )

      } else if ( CGS .OR. CGX .OR. CGV .OR. CXV ) {

         // For the nonsymmetric generalized driver routines, only one set of
         // parameters is allowed.

         READ( NIN, FMT = * )NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), NSVAL( 1 ), MXBVAL( 1 )
         if ( NBVAL( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1
            FATAL = .TRUE.
         } else if ( NBMIN( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1
            FATAL = .TRUE.
         } else if ( NXVAL( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1
            FATAL = .TRUE.
         } else if ( NSVAL( 1 ).LT.2 ) {
            WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( 1 ), 2
            FATAL = .TRUE.
         } else if ( MXBVAL( 1 ).LT.1 ) {
            WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( 1 ), 1
            FATAL = .TRUE.
         }
         xlaenv(1, NBVAL( 1 ) );
         xlaenv(2, NBMIN( 1 ) );
         xlaenv(3, NXVAL( 1 ) );
         xlaenv(4, NSVAL( 1 ) );
         xlaenv(8, MXBVAL( 1 ) );
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 )
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'NS:   ', NSVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'MAXB: ', MXBVAL( 1 )
      } else if ( .NOT.CHB .AND. .NOT.GLM .AND. .NOT.GQR .AND. .NOT. GSV .AND. .NOT.CSD .AND. .NOT.LSE ) {

         // For the other paths, the number of parameters can be varied
         // from the input file.  Read the number of parameter values.

         READ( NIN, FMT = * )NPARMS
         if ( NPARMS.LT.1 ) {
            WRITE( NOUT, FMT = 9989 )'NPARMS', NPARMS, 1
            NPARMS = 0
            FATAL = .TRUE.
         } else if ( NPARMS.GT.MAXIN ) {
            WRITE( NOUT, FMT = 9988 )'NPARMS', NPARMS, MAXIN
            NPARMS = 0
            FATAL = .TRUE.
         }

         // Read the values of NB

         if ( .NOT.CBB ) {
            READ( NIN, FMT = * )( NBVAL( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 70
               if ( NBVAL( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( I ), 0
                  FATAL = .TRUE.
               } else if ( NBVAL( I ).GT.NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'   NB ', NBVAL( I ), NMAX
                  FATAL = .TRUE.
               }
            } // 70
            WRITE( NOUT, FMT = 9983 )'NB:   ', ( NBVAL( I ), I = 1, NPARMS )
         }

         // Read the values of NBMIN

         if ( NEP .OR. SEP .OR. SVD .OR. CGG ) {
            READ( NIN, FMT = * )( NBMIN( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 80
               if ( NBMIN( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( I ), 0
                  FATAL = .TRUE.
               } else if ( NBMIN( I ).GT.NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'NBMIN ', NBMIN( I ), NMAX
                  FATAL = .TRUE.
               }
            } // 80
            WRITE( NOUT, FMT = 9983 )'NBMIN:', ( NBMIN( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 90
               NBMIN( I ) = 1
            } // 90
         }

         // Read the values of NX

         if ( NEP .OR. SEP .OR. SVD ) {
            READ( NIN, FMT = * )( NXVAL( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 100
               if ( NXVAL( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( I ), 0
                  FATAL = .TRUE.
               } else if ( NXVAL( I ).GT.NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'   NX ', NXVAL( I ), NMAX
                  FATAL = .TRUE.
               }
            } // 100
            WRITE( NOUT, FMT = 9983 )'NX:   ', ( NXVAL( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 110
               NXVAL( I ) = 1
            } // 110
         }

         // Read the values of NSHIFT (if CGG) or NRHS (if SVD
         // or CBB).

         if ( SVD .OR. CBB .OR. CGG ) {
            READ( NIN, FMT = * )( NSVAL( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 120
               if ( NSVAL( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( I ), 0
                  FATAL = .TRUE.
               } else if ( NSVAL( I ).GT.NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'   NS ', NSVAL( I ), NMAX
                  FATAL = .TRUE.
               }
            } // 120
            WRITE( NOUT, FMT = 9983 )'NS:   ', ( NSVAL( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 130
               NSVAL( I ) = 1
            } // 130
         }

         // Read the values for MAXB.

         if ( CGG ) {
            READ( NIN, FMT = * )( MXBVAL( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 140
               if ( MXBVAL( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( I ), 0
                  FATAL = .TRUE.
               } else if ( MXBVAL( I ).GT.NMAX ) {
                  WRITE( NOUT, FMT = 9988 )' MAXB ', MXBVAL( I ), NMAX
                  FATAL = .TRUE.
               }
            } // 140
            WRITE( NOUT, FMT = 9983 )'MAXB: ', ( MXBVAL( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 150
               MXBVAL( I ) = 1
            } // 150
         }

         // Read the values for INMIN.

         if ( NEP ) {
            READ( NIN, FMT = * )( INMIN( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 540
               if ( INMIN( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )' INMIN ', INMIN( I ), 0
                  FATAL = .TRUE.
               }
            } // 540
            WRITE( NOUT, FMT = 9983 )'INMIN: ', ( INMIN( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 550
               INMIN( I ) = 1
            } // 550
         }

         // Read the values for INWIN.

         if ( NEP ) {
            READ( NIN, FMT = * )( INWIN( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 560
               if ( INWIN( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )' INWIN ', INWIN( I ), 0
                  FATAL = .TRUE.
               }
            } // 560
            WRITE( NOUT, FMT = 9983 )'INWIN: ', ( INWIN( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 570
               INWIN( I ) = 1
            } // 570
         }

         // Read the values for INIBL.

         if ( NEP ) {
            READ( NIN, FMT = * )( INIBL( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 580
               if ( INIBL( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )' INIBL ', INIBL( I ), 0
                  FATAL = .TRUE.
               }
            } // 580
            WRITE( NOUT, FMT = 9983 )'INIBL: ', ( INIBL( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 590
               INIBL( I ) = 1
            } // 590
         }

         // Read the values for ISHFTS.

         if ( NEP ) {
            READ( NIN, FMT = * )( ISHFTS( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 600
               if ( ISHFTS( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )' ISHFTS ', ISHFTS( I ), 0
                  FATAL = .TRUE.
               }
            } // 600
            WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ( ISHFTS( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 610
               ISHFTS( I ) = 1
            } // 610
         }

         // Read the values for IACC22.

         if ( NEP .OR. CGG ) {
            READ( NIN, FMT = * )( IACC22( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 620
               if ( IACC22( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )' IACC22 ', IACC22( I ), 0
                  FATAL = .TRUE.
               }
            } // 620
            WRITE( NOUT, FMT = 9983 )'IACC22: ', ( IACC22( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 630
               IACC22( I ) = 1
            } // 630
         }

         // Read the values for NBCOL.

         if ( CGG ) {
            READ( NIN, FMT = * )( NBCOL( I ), I = 1, NPARMS )
            for (I = 1; I <= NPARMS; I++) { // 160
               if ( NBCOL( I ).LT.0 ) {
                  WRITE( NOUT, FMT = 9989 )'NBCOL ', NBCOL( I ), 0
                  FATAL = .TRUE.
               } else if ( NBCOL( I ).GT.NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'NBCOL ', NBCOL( I ), NMAX
                  FATAL = .TRUE.
               }
            } // 160
            WRITE( NOUT, FMT = 9983 )'NBCOL:', ( NBCOL( I ), I = 1, NPARMS )
         } else {
            for (I = 1; I <= NPARMS; I++) { // 170
               NBCOL( I ) = 1
            } // 170
         }
      }

      // Calculate and print the machine dependent constants.

      WRITE( NOUT, FMT = * )
      EPS = SLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9981 )'underflow', EPS
      EPS = SLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9981 )'overflow ', EPS
      EPS = SLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9981 )'precision', EPS

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9982 )THRESH
      if ( SEP .OR. SVD .OR. CGG ) {

         // Read the flag that indicates whether to test LAPACK routines.

         READ( NIN, FMT = * )TSTCHK

         // Read the flag that indicates whether to test driver routines.

         READ( NIN, FMT = * )TSTDRV
      }

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR

      // Read the code describing how to set the random number seed.

      READ( NIN, FMT = * )NEWSD

      // If NEWSD = 2, read another line with 4 integers for the seed.

      IF( NEWSD.EQ.2 ) READ( NIN, FMT = * )( IOLDSD( I ), I = 1, 4 )

      for (I = 1; I <= 4; I++) { // 180
         ISEED( I ) = IOLDSD( I )
      } // 180

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 )
         STOP
      }

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.

      } // 190

      if ( .NOT.( CGX .OR. CXV ) ) {

         } // 200
         READ( NIN, FMT = '(A80)', END = 380 )LINE
         C3 = LINE( 1: 3 )
         LENP = LEN( LINE )
         I = 3
         ITMP = 0
         I1 = 0
         } // 210
         I = I + 1
         if ( I.GT.LENP ) {
            if ( I1.GT.0 ) {
               GO TO 240
            } else {
               NTYPES = MAXT
               GO TO 240
            }
         }
         if ( LINE( I: I ).NE.' ' .AND. LINE( I: I ).NE.',' ) {
            I1 = I
            C1 = LINE( I1: I1 )

         // Check that a valid integer was read

            for (K = 1; K <= 10; K++) { // 220
               if ( C1.EQ.INTSTR( K: K ) ) {
                  IC = K - 1
                  GO TO 230
               }
            } // 220
            WRITE( NOUT, FMT = 9991 )I, LINE
            GO TO 200
            } // 230
            ITMP = 10*ITMP + IC
            GO TO 210
         } else if ( I1.GT.0 ) {
            GO TO 240
         } else {
            GO TO 210
         }
         } // 240
         NTYPES = ITMP

      // Skip the tests if NTYPES is <= 0.

         if ( .NOT.( CEV .OR. CES .OR. CVX .OR. CSX .OR. CGV .OR. CGS ) .AND. NTYPES.LE.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
            GO TO 200
         }

      } else {
         IF( CGX ) C3 = 'CGX'          IF( CXV ) C3 = 'CXV'
      }

      // Reset the random number seed.

      if ( NEWSD.EQ.0 ) {
         for (K = 1; K <= 4; K++) { // 250
            ISEED( K ) = IOLDSD( K )
         } // 250
      }

      if ( LSAMEN( 3, C3, 'CHS' ) .OR. LSAMEN( 3, C3, 'NEP' ) ) {

         // -------------------------------------
         // NEP:  Nonsymmetric Eigenvalue Problem
         // -------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point
            // NS    = number of shifts
            // MAXB  = minimum submatrix size

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         IF( TSTERR ) CALL CERRHS( 'CHSEQR', NOUT )
         for (I = 1; I <= NPARMS; I++) { // 270
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );
            xlaenv(12, MAX( 11, INMIN( I ) ) );
            xlaenv(13, INWIN( I ) );
            xlaenv(14, INIBL( I ) );
            xlaenv(15, ISHFTS( I ) );
            xlaenv(16, IACC22( I ) );

            if ( NEWSD.EQ.0 ) {
               for (K = 1; K <= 4; K++) { // 260
                  ISEED( K ) = IOLDSD( K )
               } // 260
            }
            WRITE( NOUT, FMT = 9961 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), MAX( 11, INMIN(I)), INWIN( I ), INIBL( I ), ISHFTS( I ), IACC22( I )             CALL CCHKHS( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), A( 1, 7 ), DC( 1, 1 ), DC( 1, 2 ), A( 1, 8 ), A( 1, 9 ), A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), DC( 1, 3 ), WORK, LWORK, RWORK, IWORK, LOGWRK, RESULT, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCHKHS', INFO
         } // 270

      } else if ( LSAMEN( 3, C3, 'CST' ) .OR. LSAMEN( 3, C3, 'SEP' ) .OR. LSAMEN( 3, C3, 'SE2' ) ) {

         // ----------------------------------
         // SEP:  Symmetric Eigenvalue Problem
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         xlaenv(9, 25 );
         if ( TSTERR ) {
#if defined(_OPENMP)
            N_THREADS = OMP_GET_MAX_THREADS()
            ONE_THREAD = 1
            omp_set_num_threads(ONE_THREAD);
#endif
            cerrst('CST', NOUT );
#if defined(_OPENMP)
            omp_set_num_threads(N_THREADS);
#endif
         }
         for (I = 1; I <= NPARMS; I++) { // 290
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );

            if ( NEWSD.EQ.0 ) {
               for (K = 1; K <= 4; K++) { // 280
                  ISEED( K ) = IOLDSD( K )
               } // 280
            }
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I )
            if ( TSTCHK ) {
               if ( LSAMEN( 3, C3, 'SE2' ) ) {
               cchkst2stg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), DR( 1, 7 ), DR( 1, 8 ), DR( 1, 9 ), DR( 1, 10 ), DR( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), DC( 1, 1 ), A( 1, 6 ), WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               } else {
               cchkst(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), DR( 1, 7 ), DR( 1, 8 ), DR( 1, 9 ), DR( 1, 10 ), DR( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), DC( 1, 1 ), A( 1, 6 ), WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               ENDIF
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCHKST', INFO
            }
            if ( TSTDRV ) {
               if ( LSAMEN( 3, C3, 'SE2' ) ) {
               cdrvst2stg(NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 8 ), DR( 1, 9 ), DR( 1, 10 ), A( 1, 2 ), NMAX, A( 1, 3 ), DC( 1, 1 ), A( 1, 4 ), WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               } else {
               cdrvst(NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 8 ), DR( 1, 9 ), DR( 1, 10 ), A( 1, 2 ), NMAX, A( 1, 3 ), DC( 1, 1 ), A( 1, 4 ), WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT, INFO );
           ENDIF
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRVST', INFO
            }
         } // 290

      } else if ( LSAMEN( 3, C3, 'CSG' ) ) {

         // ----------------------------------------------
         // CSG:  Hermitian Generalized Eigenvalue Problem
         // ----------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(9, 25 );
         for (I = 1; I <= NPARMS; I++) { // 310
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );

            if ( NEWSD.EQ.0 ) {
               for (K = 1; K <= 4; K++) { // 300
                  ISEED( K ) = IOLDSD( K )
               } // 300
            }
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I )
            if ( TSTCHK ) {
                // CALL CDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
      // $                      DR( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
      // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
      // $                      LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT,
      // $                      INFO )
               cdrvsg2stg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, DR( 1, 3 ), DR( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK, LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRVSG', INFO
            }
         } // 310

      } else if ( LSAMEN( 3, C3, 'CBD' ) .OR. LSAMEN( 3, C3, 'SVD' ) ) {

         // ----------------------------------
         // SVD:  Singular Value Decomposition
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point
            // NRHS  = number of right hand sides

         MAXTYP = 16
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(9, 25 );

         // Test the error exits

         xlaenv(1, 1 );
         IF( TSTERR .AND. TSTCHK ) CALL CERRBD( 'CBD', NOUT )          IF( TSTERR .AND. TSTDRV ) CALL CERRED( 'CBD', NOUT )

         for (I = 1; I <= NPARMS; I++) { // 330
            NRHS = NSVAL( I )
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );
            if ( NEWSD.EQ.0 ) {
               for (K = 1; K <= 4; K++) { // 320
                  ISEED( K ) = IOLDSD( K )
               } // 320
            }
            WRITE( NOUT, FMT = 9995 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), NRHS
            if ( TSTCHK ) {
               cchkbd(NN, MVAL, NVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), A( 1, 2 ), NMAX, A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), A( 1, 8 ), WORK, LWORK, RWORK, NOUT, INFO );
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCHKBD', INFO
            }
            IF( TSTDRV ) CALL CDRVBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), WORK, LWORK, RWORK, IWORK, NOUT, INFO )
         } // 330

      } else if ( LSAMEN( 3, C3, 'CEV' ) ) {

         // --------------------------------------------
         // CEV:  Nonsymmetric Eigenvalue Problem Driver
               // CGEEV (eigenvalues and eigenvectors)
         // --------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         if ( NTYPES.LE.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRED( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrvev(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), DC( 1, 1 ), DC( 1, 2 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, RESULT, WORK, LWORK, RWORK, IWORK, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CGEEV', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( LSAMEN( 3, C3, 'CES' ) ) {

         // --------------------------------------------
         // CES:  Nonsymmetric Eigenvalue Problem Driver
               // CGEES (Schur form)
         // --------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         if ( NTYPES.LE.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRED( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrves(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), DC( 1, 1 ), DC( 1, 2 ), A( 1, 4 ), NMAX, RESULT, WORK, LWORK, RWORK, IWORK, LOGWRK, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CGEES', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( LSAMEN( 3, C3, 'CVX' ) ) {

         // --------------------------------------------------------------
         // CVX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // CGEEVX (eigenvalues, eigenvectors and condition numbers)
         // --------------------------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         if ( NTYPES.LT.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRED( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrvvx(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), DC( 1, 1 ), DC( 1, 2 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), DR( 1, 7 ), DR( 1, 8 ), RESULT, WORK, LWORK, RWORK, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CGEEVX', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( LSAMEN( 3, C3, 'CSX' ) ) {

         // ---------------------------------------------------
         // CSX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // CGEESX (Schur form and condition numbers)
         // ---------------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         if ( NTYPES.LT.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRED( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrvsx(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), DC( 1, 1 ), DC( 1, 2 ), DC( 1, 3 ), A( 1, 4 ), NMAX, A( 1, 5 ), RESULT, WORK, LWORK, RWORK, LOGWRK, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CGEESX', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( LSAMEN( 3, C3, 'CGG' ) ) {

         // -------------------------------------------------
         // CGG:  Generalized Nonsymmetric Eigenvalue Problem
         // -------------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NS    = number of shifts
            // MAXB  = minimum submatrix size
            // IACC22: structured matrix multiply
            // NBCOL = minimum column dimension for blocks

         MAXTYP = 26
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1,1);
         IF( TSTCHK .AND. TSTERR ) CALL CERRGG( C3, NOUT )
         for (I = 1; I <= NPARMS; I++) { // 350
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(4, NSVAL( I ) );
            xlaenv(8, MXBVAL( I ) );
            xlaenv(16, IACC22( I ) );
            xlaenv(5, NBCOL( I ) );

            if ( NEWSD.EQ.0 ) {
               for (K = 1; K <= 4; K++) { // 340
                  ISEED( K ) = IOLDSD( K )
               } // 340
            }
            WRITE( NOUT, FMT = 9996 )C3, NBVAL( I ), NBMIN( I ), NSVAL( I ), MXBVAL( I ), IACC22( I ), NBCOL( I )
            TSTDIF = .FALSE.
            THRSHN = 10.
            if ( TSTCHK ) {
               cchkgg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, TSTDIF, THRSHN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), A( 1, 8 ), A( 1, 9 ), NMAX, A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), DC( 1, 1 ), DC( 1, 2 ), DC( 1, 3 ), DC( 1, 4 ), A( 1, 13 ), A( 1, 14 ), WORK, LWORK, RWORK, LOGWRK, RESULT, INFO );
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCHKGG', INFO
            }
         } // 350

      } else if ( LSAMEN( 3, C3, 'CGS' ) ) {

         // -------------------------------------------------
         // CGS:  Generalized Nonsymmetric Eigenvalue Problem
               // CGGES (Schur form)
         // -------------------------------------------------

         MAXTYP = 26
         NTYPES = MIN( MAXTYP, NTYPES )
         if ( NTYPES.LE.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRGG( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrges(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), DC( 1, 1 ), DC( 1, 2 ), WORK, LWORK, RWORK, RESULT, LOGWRK, INFO );

            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRGES', INFO

* Blocked version

            xlaenv(16,2);
            cdrges3(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), DC( 1, 1 ), DC( 1, 2 ), WORK, LWORK, RWORK, RESULT, LOGWRK, INFO );

            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRGES3', INFO
         }
         WRITE( NOUT, FMT = 9973 )

         GO TO 10

      } else if ( CGX ) {

         // -------------------------------------------------
         // CGX  Generalized Nonsymmetric Eigenvalue Problem
               // CGGESX (Schur form and condition numbers)
         // -------------------------------------------------

         MAXTYP = 5
         NTYPES = MAXTYP
         if ( NN.LT.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRGG( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            xlaenv(5, 2 );
            cdrgsx(NN, NCMAX, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), DC( 1, 1 ), DC( 1, 2 ), C, NCMAX*NCMAX, S, WORK, LWORK, RWORK, IWORK, LIWORK, LOGWRK, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRGSX', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( LSAMEN( 3, C3, 'CGV' ) ) {

         // -------------------------------------------------
         // CGV:  Generalized Nonsymmetric Eigenvalue Problem
               // CGGEV (Eigenvalue/vector form)
         // -------------------------------------------------

         MAXTYP = 26
         NTYPES = MIN( MAXTYP, NTYPES )
         if ( NTYPES.LE.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRGG( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrgev(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, DC( 1, 1 ), DC( 1, 2 ), DC( 1, 3 ), DC( 1, 4 ), WORK, LWORK, RWORK, RESULT, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRGEV', INFO

* Blocked version

            xlaenv(16,2);
            cdrgev3(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, DC( 1, 1 ), DC( 1, 2 ), DC( 1, 3 ), DC( 1, 4 ), WORK, LWORK, RWORK, RESULT, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRGEV3', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( CXV ) {

         // -------------------------------------------------
         // CXV:  Generalized Nonsymmetric Eigenvalue Problem
               // CGGEVX (eigenvalue/vector with condition numbers)
         // -------------------------------------------------

         MAXTYP = 2
         NTYPES = MAXTYP
         if ( NN.LT.0 ) {
            WRITE( NOUT, FMT = 9990 )C3
         } else {
            IF( TSTERR ) CALL CERRGG( C3, NOUT )
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            cdrgvx(NN, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), DC( 1, 1 ), DC( 1, 2 ), A( 1, 5 ), A( 1, 6 ), IWORK( 1 ), IWORK( 2 ), DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), DR( 1, 6 ), WORK, LWORK, RWORK, IWORK( 3 ), LIWORK-2, RESULT, LOGWRK, INFO );

            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CDRGVX', INFO
         }
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      } else if ( LSAMEN( 3, C3, 'CHB' ) ) {

         // ------------------------------
         // CHB:  Hermitian Band Reduction
         // ------------------------------

         MAXTYP = 15
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         if ( TSTERR ) {
#if defined(_OPENMP)
            N_THREADS = OMP_GET_MAX_THREADS()
            ONE_THREAD = 1
            omp_set_num_threads(ONE_THREAD);
#endif
            cerrst('CHB', NOUT );
#if defined(_OPENMP)
            omp_set_num_threads(N_THREADS);
#endif
         }
          // CALL CCHKHB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ),
      // $                A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT,
      // $                INFO )
         cchkhb2stg(NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ), DR( 1, 3 ), DR( 1, 4 ), DR( 1, 5 ), A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT, INFO );
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCHKHB', INFO

      } else if ( LSAMEN( 3, C3, 'CBB' ) ) {

         // ------------------------------
         // CBB:  General Band Reduction
         // ------------------------------

         MAXTYP = 15
         NTYPES = MIN( MAXTYP, NTYPES )
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         for (I = 1; I <= NPARMS; I++) { // 370
            NRHS = NSVAL( I )

            if ( NEWSD.EQ.0 ) {
               for (K = 1; K <= 4; K++) { // 360
                  ISEED( K ) = IOLDSD( K )
               } // 360
            }
            WRITE( NOUT, FMT = 9966 )C3, NRHS
            cchkbb(NN, MVAL, NVAL, NK, KVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), 2*NMAX, DR( 1, 1 ), DR( 1, 2 ), A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), WORK, LWORK, RWORK, RESULT, INFO );
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCHKBB', INFO
         } // 370

      } else if ( LSAMEN( 3, C3, 'GLM' ) ) {

         // -----------------------------------------
         // GLM:  Generalized Linear Regression Model
         // -----------------------------------------

         xlaenv(1, 1 );
         IF( TSTERR ) CALL CERRGG( 'GLM', NOUT )          CALL CCKGLM( NN, NVAL, MVAL, PVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, DR( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCKGLM', INFO

      } else if ( LSAMEN( 3, C3, 'GQR' ) ) {

         // ------------------------------------------
         // GQR:  Generalized QR and RQ factorizations
         // ------------------------------------------

         xlaenv(1, 1 );
         IF( TSTERR ) CALL CERRGG( 'GQR', NOUT )          CALL CCKGQR( NN, MVAL, NN, PVAL, NN, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), TAUA, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), B( 1, 5 ), TAUB, WORK, DR( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCKGQR', INFO

      } else if ( LSAMEN( 3, C3, 'GSV' ) ) {

         // ----------------------------------------------
         // GSV:  Generalized Singular Value Decomposition
         // ----------------------------------------------

         xlaenv(1,1);
         IF( TSTERR ) CALL CERRGG( 'GSV', NOUT )          CALL CCKGSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), ALPHA, BETA, B( 1, 4 ), IWORK, WORK, DR( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCKGSV', INFO

      } else if ( LSAMEN( 3, C3, 'CSD' ) ) {

         // ----------------------------------------------
         // CSD:  CS Decomposition
         // ----------------------------------------------

         xlaenv(1,1);
         IF( TSTERR ) CALL CERRGG( 'CSD', NOUT )          CALL CCKCSD( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), RWORK, IWORK, WORK, DR( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCKCSD', INFO

      } else if ( LSAMEN( 3, C3, 'LSE' ) ) {

         // --------------------------------------
         // LSE:  Constrained Linear Least Squares
         // --------------------------------------

         xlaenv(1, 1 );
         IF( TSTERR ) CALL CERRGG( 'LSE', NOUT )          CALL CCKLSE( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, DR( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'CCKLSE', INFO
      } else {
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9992 )C3
      }
      IF( .NOT.( CGX .OR. CXV ) ) GO TO 190
      } // 380
      WRITE( NOUT, FMT = 9994 )
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9993 )S2 - S1

      DEALLOCATE (S, STAT = AllocateStatus)
      DEALLOCATE (A, STAT = AllocateStatus)
      DEALLOCATE (B, STAT = AllocateStatus)
      DEALLOCATE (C, STAT = AllocateStatus)
      DEALLOCATE (RWORK, STAT = AllocateStatus)
      DEALLOCATE (WORK,  STAT = AllocateStatus)

 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9997 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4 )
 9996 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NS =', I4, ', MAXB =', I4, ', IACC22 =', I4, ', NBCOL =', I4 )
 9995 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, ', NRHS =', I4 )
 9994 FORMAT( / / ' End of tests' )
 9993 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9992 FORMAT( 1X, A3, ':  Unrecognized path name' )
 9991 FORMAT( / / ' *** Invalid int     value in column ', I2,; ' of input', ' line:', / A79 )
 9990 FORMAT( / / 1X, A3, ' routines were not tested' )
 9989 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be >=', I6 )
 9988 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be <=', I6 )
 9987 FORMAT( ' Tests of the Nonsymmetric Eigenvalue Problem routines' )
 9986 FORMAT( ' Tests of the Hermitian Eigenvalue Problem routines' )
 9985 FORMAT( ' Tests of the Singular Value Decomposition routines' )
 9984 FORMAT( / ' The following parameter values will be used:' )
 9983 FORMAT( 4X, A, 10I6, / 10X, 10I6 )
 9982 FORMAT( / ' Routines pass computational tests if test ratio is ', 'less than', F8.2, / )
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 )
 9980 FORMAT( ' *** Error code from ', A, ' = ', I4 )
 9979 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', / '    CGEEV (eigenvalues and eigevectors)' )
 9978 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', / '    CGEES (Schur form)' )
 9977 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', ' Driver', / '    CGEEVX (eigenvalues, eigenvectors and', ' condition numbers)' )
 9976 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', ' Driver', / '    CGEESX (Schur form and condition', ' numbers)' )
 9975 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem routines' )
 9974 FORMAT( ' Tests of CHBTRD', / ' (reduction of a Hermitian band ', 'matrix to real tridiagonal form)' )
 9973 FORMAT( / 1X, 71( '-' ) )
 9972 FORMAT( / ' LAPACK VERSION ', I1, '.', I1, '.', I1 )
 9971 FORMAT( / ' Tests of the Generalized Linear Regression Model ', 'routines' )
 9970 FORMAT( / ' Tests of the Generalized QR and RQ routines' )
 9969 FORMAT( / ' Tests of the Generalized Singular Value', ' Decomposition routines' )
 9968 FORMAT( / ' Tests of the Linear Least Squares routines' )
 9967 FORMAT( ' Tests of CGBBRD', / ' (reduction of a general band ', 'matrix to real bidiagonal form)' )
 9966 FORMAT( / / 1X, A3, ':  NRHS =', I4 )
 9965 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Expert Driver CGGESX' )
 9964 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Driver CGGES' )
 9963 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Driver CGGEV' )
 9962 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Expert Driver CGGEVX' )
 9961 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, ', INMIN=', I4, ', INWIN =', I4, ', INIBL =', I4, ', ISHFTS =', I4, ', IACC22 =', I4)
 9960 FORMAT( / ' Tests of the CS Decomposition routines' )

      // End of CCHKEE

      }
