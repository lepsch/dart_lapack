      void main() {
// #if defined(_OPENMP)
      use omp_lib;
// #endif

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 132 ;
      int                NCMAX;
      const              NCMAX = 20 ;
      int                NEED;
      const              NEED = 14 ;
      int                LWORK;
      const              LWORK = NMAX*( 5*NMAX+5 )+1 ;
      int                LIWORK;
      const              LIWORK = NMAX*( 5*NMAX+20 ) ;
      int                MAXIN;
      const              MAXIN = 20 ;
      int                MAXT;
      const              MAXT = 30 ;
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      // ..
      // .. Local Scalars ..
      bool               CSD, FATAL, GLM, GQR, GSV, LSE, NEP, SBB, SBK, SBL, SEP, SES, SEV, SGG, SGK, SGL, SGS, SGV, SGX, SSB, SSX, SVD, SVX, SXV, TSTCHK, TSTDIF, TSTDRV, TSTERR;
      String             C1;
      String             C3, PATH;
      String             VNAME;
      String             INTSTR;
      String             LINE;
      int                I, I1, IC, INFO, ITMP, K, LENP, MAXTYP, NEWSD, NK, NN, NPARMS, NRHS, NTYPES, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      int    *4          N_THREADS, ONE_THREAD;
      REAL               EPS, S1, S2, THRESH, THRSHN;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MAXT ), LOGWRK( NMAX );
      int                IOLDSD( 4 ), ISEED( 4 ), IWORK( LIWORK ), KVAL( MAXIN ), MVAL( MAXIN ), MXBVAL( MAXIN ), NBCOL( MAXIN ), NBMIN( MAXIN ), NBVAL( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), PVAL( MAXIN );
      int                INMIN( MAXIN ), INWIN( MAXIN ), INIBL( MAXIN ), ISHFTS( MAXIN ), IACC22( MAXIN );
      REAL               D( NMAX, 12 ), RESULT( 500 ), TAUA( NMAX ), TAUB( NMAX ), X( 5*NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      REAL, DIMENSION(:), ALLOCATABLE :: WORK;
      REAL, DIMENSION(:,:), ALLOCATABLE :: A, B, C;
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      REAL               SECOND, SLAMCH;
      // EXTERNAL LSAMEN, SECOND, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, SCHKBB, SCHKBD, SCHKBK, SCHKBL, SCHKEC, SCHKGG, SCHKGK, SCHKGL, SCHKHS, SCHKSB, SCHKST, SCKCSD, SCKGLM, SCKGQR, SCKGSV, SCKLSE, SDRGES, SDRGEV, SDRGSX, SDRGVX, SDRVBD, SDRVES, SDRVEV, SDRVSG, SDRVST, SDRVSX, SDRVVX, SERRBD, SERRED, SERRGG, SERRHS, SERRST, ILAVER, XLAENV, SDRGES3, SDRGEV3, SCHKST2STG, SDRVST2STG, SCHKSB2STG, SDRVSG2STG
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
      REAL               SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Common blocks ..
      // COMMON / CENVIR / NPROC, NSHIFT, MAXB
      // COMMON / CLAENV / IPARMS
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Data statements ..
      const INTSTR = '0123456789';
      const IOLDSD = [ 0, 0, 0, 1 ];
      // ..
      // .. Allocate memory dynamically ..

      ALLOCATE ( A(NMAX*NMAX,NEED), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( B(NMAX*NMAX,5), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( C(NCMAX*NCMAX,NCMAX*NCMAX), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( WORK(LWORK), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      // ..
      // .. Executable Statements ..

      A = 0.0;
      B = 0.0;
      C = 0.0;
      D = 0.0;
      S1 = SECOND( );
      FATAL = false;
      NUNIT = NOUT;

      // Return to here to read multiple sets of data

      } // 10

      // Read the first line and set the 3-character test path

      READ( NIN, FMT = '(A80)', END = 380 )LINE;
      PATH = LINE( 1: 3 );
      NEP = LSAMEN( 3, PATH, 'NEP' ) || LSAMEN( 3, PATH, 'SHS' );
      SEP = LSAMEN( 3, PATH, 'SEP' ) || LSAMEN( 3, PATH, 'SST' ) || LSAMEN( 3, PATH, 'SSG' ) || LSAMEN( 3, PATH, 'SE2' );
      SVD = LSAMEN( 3, PATH, 'SVD' ) || LSAMEN( 3, PATH, 'DBD' );
      SVD = LSAMEN( 3, PATH, 'SVD' ) || LSAMEN( 3, PATH, 'SBD' );
      SEV = LSAMEN( 3, PATH, 'SEV' );
      SES = LSAMEN( 3, PATH, 'SES' );
      SVX = LSAMEN( 3, PATH, 'SVX' );
      SSX = LSAMEN( 3, PATH, 'SSX' );
      SGG = LSAMEN( 3, PATH, 'SGG' );
      SGS = LSAMEN( 3, PATH, 'SGS' );
      SGX = LSAMEN( 3, PATH, 'SGX' );
      SGV = LSAMEN( 3, PATH, 'SGV' );
      SXV = LSAMEN( 3, PATH, 'SXV' );
      SSB = LSAMEN( 3, PATH, 'SSB' );
      SBB = LSAMEN( 3, PATH, 'SBB' );
      GLM = LSAMEN( 3, PATH, 'GLM' );
      GQR = LSAMEN( 3, PATH, 'GQR' ) || LSAMEN( 3, PATH, 'GRQ' );
      GSV = LSAMEN( 3, PATH, 'GSV' );
      CSD = LSAMEN( 3, PATH, 'CSD' );
      LSE = LSAMEN( 3, PATH, 'LSE' );
      SBL = LSAMEN( 3, PATH, 'SBL' );
      SBK = LSAMEN( 3, PATH, 'SBK' );
      SGL = LSAMEN( 3, PATH, 'SGL' );
      SGK = LSAMEN( 3, PATH, 'SGK' );

      // Report values of parameters.

      if ( PATH == '   ' ) {
         GO TO 10;
      } else if ( NEP ) {
         WRITE( NOUT, FMT = 9987 );
      } else if ( SEP ) {
         WRITE( NOUT, FMT = 9986 );
      } else if ( SVD ) {
         WRITE( NOUT, FMT = 9985 );
      } else if ( SEV ) {
         WRITE( NOUT, FMT = 9979 );
      } else if ( SES ) {
         WRITE( NOUT, FMT = 9978 );
      } else if ( SVX ) {
         WRITE( NOUT, FMT = 9977 );
      } else if ( SSX ) {
         WRITE( NOUT, FMT = 9976 );
      } else if ( SGG ) {
         WRITE( NOUT, FMT = 9975 );
      } else if ( SGS ) {
         WRITE( NOUT, FMT = 9964 );
      } else if ( SGX ) {
         WRITE( NOUT, FMT = 9965 );
      } else if ( SGV ) {
         WRITE( NOUT, FMT = 9963 );
      } else if ( SXV ) {
         WRITE( NOUT, FMT = 9962 );
      } else if ( SSB ) {
         WRITE( NOUT, FMT = 9974 );
      } else if ( SBB ) {
         WRITE( NOUT, FMT = 9967 );
      } else if ( GLM ) {
         WRITE( NOUT, FMT = 9971 );
      } else if ( GQR ) {
         WRITE( NOUT, FMT = 9970 );
      } else if ( GSV ) {
         WRITE( NOUT, FMT = 9969 );
      } else if ( CSD ) {
         WRITE( NOUT, FMT = 9960 );
      } else if ( LSE ) {
         WRITE( NOUT, FMT = 9968 );
      } else if ( SBL ) {

         // SGEBAL:  Balancing

         schkbl(NIN, NOUT );
         GO TO 10;
      } else if ( SBK ) {

         // SGEBAK:  Back transformation

         schkbk(NIN, NOUT );
         GO TO 10;
      } else if ( SGL ) {

         // SGGBAL:  Balancing

         schkgl(NIN, NOUT );
         GO TO 10;
      } else if ( SGK ) {

         // SGGBAK:  Back transformation

         schkgk(NIN, NOUT );
         GO TO 10;
      } else if ( LSAMEN( 3, PATH, 'SEC' ) ) {

         // SEC:  Eigencondition estimation

         READ( NIN, FMT = * )THRESH;
         xlaenv(1, 1 );
         xlaenv(12, 11 );
         xlaenv(13, 2 );
         xlaenv(14, 0 );
         xlaenv(15, 2 );
         xlaenv(16, 2 );
         TSTERR = true;
         schkec(THRESH, TSTERR, NIN, NOUT );
         GO TO 10;
      } else {
         WRITE( NOUT, FMT = 9992 )PATH;
         GO TO 10;
      }
      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9972 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      WRITE( NOUT, FMT = 9984 );

      // Read the number of values of M, P, and N.

      READ( NIN, FMT = * )NN;
      if ( NN < 0 ) {
         WRITE( NOUT, FMT = 9989 )'   NN ', NN, 1;
         NN = 0;
         FATAL = true;
      } else if ( NN > MAXIN ) {
         WRITE( NOUT, FMT = 9988 )'   NN ', NN, MAXIN;
         NN = 0;
         FATAL = true;
      }

      // Read the values of M

      if ( !( SGX || SXV ) ) {
         READ( NIN, FMT = * )( MVAL( I ), I = 1, NN );
         if ( SVD ) {
            VNAME = '    M ';
         } else {
            VNAME = '    N ';
         }
         for (I = 1; I <= NN; I++) { // 20
            if ( MVAL( I ) < 0 ) {
               WRITE( NOUT, FMT = 9989 )VNAME, MVAL( I ), 0;
               FATAL = true;
            } else if ( MVAL( I ) > NMAX ) {
               WRITE( NOUT, FMT = 9988 )VNAME, MVAL( I ), NMAX;
               FATAL = true;
            }
         } // 20
         WRITE( NOUT, FMT = 9983 )'M:    ', ( MVAL( I ), I = 1, NN );
      }

      // Read the values of P

      if ( GLM || GQR || GSV || CSD || LSE ) {
         READ( NIN, FMT = * )( PVAL( I ), I = 1, NN );
         for (I = 1; I <= NN; I++) { // 30
            if ( PVAL( I ) < 0 ) {
               WRITE( NOUT, FMT = 9989 )' P  ', PVAL( I ), 0;
               FATAL = true;
            } else if ( PVAL( I ) > NMAX ) {
               WRITE( NOUT, FMT = 9988 )' P  ', PVAL( I ), NMAX;
               FATAL = true;
            }
         } // 30
         WRITE( NOUT, FMT = 9983 )'P:    ', ( PVAL( I ), I = 1, NN );
      }

      // Read the values of N

      if ( SVD || SBB || GLM || GQR || GSV || CSD || LSE ) {
         READ( NIN, FMT = * )( NVAL( I ), I = 1, NN );
         for (I = 1; I <= NN; I++) { // 40
            if ( NVAL( I ) < 0 ) {
               WRITE( NOUT, FMT = 9989 )'    N ', NVAL( I ), 0;
               FATAL = true;
            } else if ( NVAL( I ) > NMAX ) {
               WRITE( NOUT, FMT = 9988 )'    N ', NVAL( I ), NMAX;
               FATAL = true;
            }
         } // 40
      } else {
         for (I = 1; I <= NN; I++) { // 50
            NVAL( I ) = MVAL( I );
         } // 50
      }
      if ( !( SGX || SXV ) ) {
         WRITE( NOUT, FMT = 9983 )'N:    ', ( NVAL( I ), I = 1, NN );
      } else {
         WRITE( NOUT, FMT = 9983 )'N:    ', NN;
      }

      // Read the number of values of K, followed by the values of K

      if ( SSB || SBB ) {
         READ( NIN, FMT = * )NK;
         READ( NIN, FMT = * )( KVAL( I ), I = 1, NK );
         for (I = 1; I <= NK; I++) { // 60
            if ( KVAL( I ) < 0 ) {
               WRITE( NOUT, FMT = 9989 )'    K ', KVAL( I ), 0;
               FATAL = true;
            } else if ( KVAL( I ) > NMAX ) {
               WRITE( NOUT, FMT = 9988 )'    K ', KVAL( I ), NMAX;
               FATAL = true;
            }
         } // 60
         WRITE( NOUT, FMT = 9983 )'K:    ', ( KVAL( I ), I = 1, NK );
      }

      if ( SEV || SES || SVX || SSX ) {

         // For the nonsymmetric QR driver routines, only one set of
         // parameters is allowed.

         READ( NIN, FMT = * )NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), INMIN( 1 ), INWIN( 1 ), INIBL(1), ISHFTS(1), IACC22(1);
         if ( NBVAL( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1;
            FATAL = true;
         } else if ( NBMIN( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1;
            FATAL = true;
         } else if ( NXVAL( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1;
            FATAL = true;
         } else if ( INMIN( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   INMIN ', INMIN( 1 ), 1;
            FATAL = true;
         } else if ( INWIN( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   INWIN ', INWIN( 1 ), 1;
            FATAL = true;
         } else if ( INIBL( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   INIBL ', INIBL( 1 ), 1;
            FATAL = true;
         } else if ( ISHFTS( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   ISHFTS ', ISHFTS( 1 ), 1;
            FATAL = true;
         } else if ( IACC22( 1 ) < 0 ) {
            WRITE( NOUT, FMT = 9989 )'   IACC22 ', IACC22( 1 ), 0;
            FATAL = true;
         }
         xlaenv(1, NBVAL( 1 ) );
         xlaenv(2, NBMIN( 1 ) );
         xlaenv(3, NXVAL( 1 ) );
         xlaenv(12, MAX( 11, INMIN( 1 ) ) );
         xlaenv(13, INWIN( 1 ) );
         xlaenv(14, INIBL( 1 ) );
         xlaenv(15, ISHFTS( 1 ) );
         xlaenv(16, IACC22( 1 ) );
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 );
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 );
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 );
         WRITE( NOUT, FMT = 9983 )'INMIN:   ', INMIN( 1 );
         WRITE( NOUT, FMT = 9983 )'INWIN: ', INWIN( 1 );
         WRITE( NOUT, FMT = 9983 )'INIBL: ', INIBL( 1 );
         WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ISHFTS( 1 );
         WRITE( NOUT, FMT = 9983 )'IACC22: ', IACC22( 1 );

      } else if ( SGS || SGX || SGV || SXV ) {

         // For the nonsymmetric generalized driver routines, only one set
         // of parameters is allowed.

         READ( NIN, FMT = * )NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), NSVAL( 1 ), MXBVAL( 1 );
         if ( NBVAL( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1;
            FATAL = true;
         } else if ( NBMIN( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1;
            FATAL = true;
         } else if ( NXVAL( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1;
            FATAL = true;
         } else if ( NSVAL( 1 ) < 2 ) {
            WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( 1 ), 2;
            FATAL = true;
         } else if ( MXBVAL( 1 ) < 1 ) {
            WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( 1 ), 1;
            FATAL = true;
         }
         xlaenv(1, NBVAL( 1 ) );
         xlaenv(2, NBMIN( 1 ) );
         xlaenv(3, NXVAL( 1 ) );
         xlaenv(4, NSVAL( 1 ) );
         xlaenv(8, MXBVAL( 1 ) );
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 );
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 );
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 );
         WRITE( NOUT, FMT = 9983 )'NS:   ', NSVAL( 1 );
         WRITE( NOUT, FMT = 9983 )'MAXB: ', MXBVAL( 1 );

      } else if ( !SSB && !GLM && !GQR && !GSV && !CSD && !LSE ) {

         // For the other paths, the number of parameters can be varied
         // from the input file.  Read the number of parameter values.

         READ( NIN, FMT = * )NPARMS;
         if ( NPARMS < 1 ) {
            WRITE( NOUT, FMT = 9989 )'NPARMS', NPARMS, 1;
            NPARMS = 0;
            FATAL = true;
         } else if ( NPARMS > MAXIN ) {
            WRITE( NOUT, FMT = 9988 )'NPARMS', NPARMS, MAXIN;
            NPARMS = 0;
            FATAL = true;
         }

         // Read the values of NB

         if ( !SBB ) {
            READ( NIN, FMT = * )( NBVAL( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 70
               if ( NBVAL( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( I ), 0;
                  FATAL = true;
               } else if ( NBVAL( I ) > NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'   NB ', NBVAL( I ), NMAX;
                  FATAL = true;
               }
            } // 70
            WRITE( NOUT, FMT = 9983 )'NB:   ', ( NBVAL( I ), I = 1, NPARMS );
         }

         // Read the values of NBMIN

         if ( NEP || SEP || SVD || SGG ) {
            READ( NIN, FMT = * )( NBMIN( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 80
               if ( NBMIN( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( I ), 0;
                  FATAL = true;
               } else if ( NBMIN( I ) > NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'NBMIN ', NBMIN( I ), NMAX;
                  FATAL = true;
               }
            } // 80
            WRITE( NOUT, FMT = 9983 )'NBMIN:', ( NBMIN( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 90
               NBMIN( I ) = 1;
            } // 90
         }

         // Read the values of NX

         if ( NEP || SEP || SVD ) {
            READ( NIN, FMT = * )( NXVAL( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 100
               if ( NXVAL( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( I ), 0;
                  FATAL = true;
               } else if ( NXVAL( I ) > NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'   NX ', NXVAL( I ), NMAX;
                  FATAL = true;
               }
            } // 100
            WRITE( NOUT, FMT = 9983 )'NX:   ', ( NXVAL( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 110
               NXVAL( I ) = 1;
            } // 110
         }

         // Read the values of NSHIFT (if SGG) or NRHS (if SVD
         // or SBB).

         if ( SVD || SBB || SGG ) {
            READ( NIN, FMT = * )( NSVAL( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 120
               if ( NSVAL( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( I ), 0;
                  FATAL = true;
               } else if ( NSVAL( I ) > NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'   NS ', NSVAL( I ), NMAX;
                  FATAL = true;
               }
            } // 120
            WRITE( NOUT, FMT = 9983 )'NS:   ', ( NSVAL( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 130
               NSVAL( I ) = 1;
            } // 130
         }

         // Read the values for MAXB.

         if ( SGG ) {
            READ( NIN, FMT = * )( MXBVAL( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 140
               if ( MXBVAL( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( I ), 0;
                  FATAL = true;
               } else if ( MXBVAL( I ) > NMAX ) {
                  WRITE( NOUT, FMT = 9988 )' MAXB ', MXBVAL( I ), NMAX;
                  FATAL = true;
               }
            } // 140
            WRITE( NOUT, FMT = 9983 )'MAXB: ', ( MXBVAL( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 150
               MXBVAL( I ) = 1;
            } // 150
         }

         // Read the values for INMIN.

         if ( NEP ) {
            READ( NIN, FMT = * )( INMIN( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 540
               if ( INMIN( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )' INMIN ', INMIN( I ), 0;
                  FATAL = true;
               }
            } // 540
            WRITE( NOUT, FMT = 9983 )'INMIN: ', ( INMIN( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 550
               INMIN( I ) = 1;
            } // 550
         }

         // Read the values for INWIN.

         if ( NEP ) {
            READ( NIN, FMT = * )( INWIN( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 560
               if ( INWIN( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )' INWIN ', INWIN( I ), 0;
                  FATAL = true;
               }
            } // 560
            WRITE( NOUT, FMT = 9983 )'INWIN: ', ( INWIN( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 570
               INWIN( I ) = 1;
            } // 570
         }

         // Read the values for INIBL.

         if ( NEP ) {
            READ( NIN, FMT = * )( INIBL( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 580
               if ( INIBL( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )' INIBL ', INIBL( I ), 0;
                  FATAL = true;
               }
            } // 580
            WRITE( NOUT, FMT = 9983 )'INIBL: ', ( INIBL( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 590
               INIBL( I ) = 1;
            } // 590
         }

         // Read the values for ISHFTS.

         if ( NEP ) {
            READ( NIN, FMT = * )( ISHFTS( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 600
               if ( ISHFTS( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )' ISHFTS ', ISHFTS( I ), 0;
                  FATAL = true;
               }
            } // 600
            WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ( ISHFTS( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 610
               ISHFTS( I ) = 1;
            } // 610
         }

         // Read the values for IACC22.

         if ( NEP || SGG ) {
            READ( NIN, FMT = * )( IACC22( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 620
               if ( IACC22( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )' IACC22 ', IACC22( I ), 0;
                  FATAL = true;
               }
            } // 620
            WRITE( NOUT, FMT = 9983 )'IACC22: ', ( IACC22( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 630
               IACC22( I ) = 1;
            } // 630
         }

         // Read the values for NBCOL.

         if ( SGG ) {
            READ( NIN, FMT = * )( NBCOL( I ), I = 1, NPARMS );
            for (I = 1; I <= NPARMS; I++) { // 160
               if ( NBCOL( I ) < 0 ) {
                  WRITE( NOUT, FMT = 9989 )'NBCOL ', NBCOL( I ), 0;
                  FATAL = true;
               } else if ( NBCOL( I ) > NMAX ) {
                  WRITE( NOUT, FMT = 9988 )'NBCOL ', NBCOL( I ), NMAX;
                  FATAL = true;
               }
            } // 160
            WRITE( NOUT, FMT = 9983 )'NBCOL:', ( NBCOL( I ), I = 1, NPARMS );
         } else {
            for (I = 1; I <= NPARMS; I++) { // 170
               NBCOL( I ) = 1;
            } // 170
         }
      }

      // Calculate and print the machine dependent constants.

      WRITE( NOUT, FMT = * );
      EPS = SLAMCH( 'Underflow threshold' );
      WRITE( NOUT, FMT = 9981 )'underflow', EPS;
      EPS = SLAMCH( 'Overflow threshold' );
      WRITE( NOUT, FMT = 9981 )'overflow ', EPS;
      EPS = SLAMCH( 'Epsilon' );
      WRITE( NOUT, FMT = 9981 )'precision', EPS;

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH;
      WRITE( NOUT, FMT = 9982 )THRESH;
      if ( SEP || SVD || SGG ) {

         // Read the flag that indicates whether to test LAPACK routines.

         READ( NIN, FMT = * )TSTCHK;

         // Read the flag that indicates whether to test driver routines.

         READ( NIN, FMT = * )TSTDRV;
      }

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR;

      // Read the code describing how to set the random number seed.

      READ( NIN, FMT = * )NEWSD;

      // If NEWSD = 2, read another line with 4 integers for the seed.

      if (NEWSD == 2) READ( NIN, FMT = * )( IOLDSD( I ), I = 1, 4 );

      for (I = 1; I <= 4; I++) { // 180
         ISEED( I ) = IOLDSD( I );
      } // 180

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 );
         STOP;
      }

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.

      } // 190

      if ( !( SGX || SXV ) ) {

         } // 200
         READ( NIN, FMT = '(A80)', END = 380 )LINE;
         C3 = LINE( 1: 3 );
         LENP = LEN( LINE );
         I = 3;
         ITMP = 0;
         I1 = 0;
         } // 210
         I = I + 1;
         if ( I > LENP ) {
            if ( I1 > 0 ) {
               GO TO 240;
            } else {
               NTYPES = MAXT;
               GO TO 240;
            }
         }
         if ( LINE( I: I ) != ' ' && LINE( I: I ) != ',' ) {
            I1 = I;
            C1 = LINE( I1: I1 );

         // Check that a valid integer was read

            for (K = 1; K <= 10; K++) { // 220
               if ( C1 == INTSTR( K: K ) ) {
                  IC = K - 1;
                  GO TO 230;
               }
            } // 220
            WRITE( NOUT, FMT = 9991 )I, LINE;
            GO TO 200;
            } // 230
            ITMP = 10*ITMP + IC;
            GO TO 210;
         } else if ( I1 > 0 ) {
            GO TO 240;
         } else {
            GO TO 210;
         }
         } // 240
         NTYPES = ITMP;

      // Skip the tests if NTYPES is <= 0.

         if ( !( SEV || SES || SVX || SSX || SGV || SGS ) && NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
            GO TO 200;
         }

      } else {
         if (SXV) C3 = 'SXV'          IF( SGX ) C3 = 'SGX';
      }

      // Reset the random number seed.

      if ( NEWSD == 0 ) {
         for (K = 1; K <= 4; K++) { // 250
            ISEED( K ) = IOLDSD( K );
         } // 250
      }

      if ( LSAMEN( 3, C3, 'SHS' ) || LSAMEN( 3, C3, 'NEP' ) ) {

         // -------------------------------------
         // NEP:  Nonsymmetric Eigenvalue Problem
         // -------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point
            // NS    = number of shifts
            // MAXB  = minimum submatrix size

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         if (TSTERR) CALL SERRHS( 'SHSEQR', NOUT );
         for (I = 1; I <= NPARMS; I++) { // 270
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );
            xlaenv(12, MAX( 11, INMIN( I ) ) );
            xlaenv(13, INWIN( I ) );
            xlaenv(14, INIBL( I ) );
            xlaenv(15, ISHFTS( I ) );
            xlaenv(16, IACC22( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 260
                  ISEED( K ) = IOLDSD( K );
               } // 260
            }
            WRITE( NOUT, FMT = 9961 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), MAX( 11, INMIN(I)), INWIN( I ), INIBL( I ), ISHFTS( I ), IACC22( I );
            schkhs(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), A( 1, 7 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 8 ), A( 1, 9 ), A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), D( 1, 7 ), WORK, LWORK, IWORK, LOGWRK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCHKHS', INFO;
         } // 270

      } else if ( LSAMEN( 3, C3, 'SST' ) || LSAMEN( 3, C3, 'SEP' ) || LSAMEN( 3, C3, 'SE2' ) ) {

         // ----------------------------------
         // SEP:  Symmetric Eigenvalue Problem
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         xlaenv(9, 25 );
         if ( TSTERR ) {
// #if defined(_OPENMP)
            N_THREADS = OMP_GET_MAX_THREADS();
            ONE_THREAD = 1;
            omp_set_num_threads(ONE_THREAD);
// #endif
            serrst('SST', NOUT );
// #if defined(_OPENMP)
            omp_set_num_threads(N_THREADS);
// #endif
         }
         for (I = 1; I <= NPARMS; I++) { // 290
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 280
                  ISEED( K ) = IOLDSD( K );
               } // 280
            }
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I );
            if ( TSTCHK ) {
               if ( LSAMEN( 3, C3, 'SE2' ) ) {
               schkst2stg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               } else {
               schkst(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               }
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCHKST', INFO;
            }
            if ( TSTDRV ) {
               if ( LSAMEN( 3, C3, 'SE2' ) ) {
               sdrvst2stg(NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11), A( 1, 2 ), NMAX, A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               } else {
               sdrvst(NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11), A( 1, 2 ), NMAX, A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               }
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRVST', INFO;
            }
         } // 290

      } else if ( LSAMEN( 3, C3, 'SSG' ) ) {

         // ----------------------------------------------
         // SSG:  Symmetric Generalized Eigenvalue Problem
         // ----------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(9, 25 );
         for (I = 1; I <= NPARMS; I++) { // 310
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 300
                  ISEED( K ) = IOLDSD( K );
               } // 300
            }
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I );
            if ( TSTCHK ) {
                // CALL SDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
      // $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
      // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
      // $                      LWORK, IWORK, LIWORK, RESULT, INFO )
               sdrvsg2stg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, D( 1, 3 ), D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRVSG', INFO;
            }
         } // 310

      } else if ( LSAMEN( 3, C3, 'SBD' ) || LSAMEN( 3, C3, 'SVD' ) ) {

         // ----------------------------------
         // SVD:  Singular Value Decomposition
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point
            // NRHS  = number of right hand sides

         MAXTYP = 16;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         xlaenv(9, 25 );

         // Test the error exits

         if (TSTERR && TSTCHK) CALL SERRBD( 'SBD', NOUT )          IF( TSTERR && TSTDRV ) CALL SERRED( 'SBD', NOUT );

         for (I = 1; I <= NPARMS; I++) { // 330
            NRHS = NSVAL( I );
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );
            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 320
                  ISEED( K ) = IOLDSD( K );
               } // 320
            }
            WRITE( NOUT, FMT = 9995 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), NRHS;
            if ( TSTCHK ) {
               schkbd(NN, MVAL, NVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 2 ), NMAX, A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), A( 1, 8 ), WORK, LWORK, IWORK, NOUT, INFO );
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCHKBD', INFO;
            }
            if (TSTDRV) CALL SDRVBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, IWORK, NOUT, INFO );
         } // 330

      } else if ( LSAMEN( 3, C3, 'SEV' ) ) {

         // --------------------------------------------
         // SEV:  Nonsymmetric Eigenvalue Problem Driver
               // SGEEV (eigenvalues and eigenvectors)
         // --------------------------------------------

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRED( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrvev(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, RESULT, WORK, LWORK, IWORK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SGEEV', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'SES' ) ) {

         // --------------------------------------------
         // SES:  Nonsymmetric Eigenvalue Problem Driver
               // SGEES (Schur form)
         // --------------------------------------------

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRED( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrves(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 4 ), NMAX, RESULT, WORK, LWORK, IWORK, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SGEES', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

         // --------------------------------------------------------------
         // SVX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // SGEEVX (eigenvalues, eigenvectors and condition numbers)
         // --------------------------------------------------------------

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         if ( NTYPES < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRED( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrvvx(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), D( 1, 12 ), RESULT, WORK, LWORK, IWORK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SGEEVX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'SSX' ) ) {

         // ---------------------------------------------------
         // SSX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // SGEESX (Schur form and condition numbers)
         // ---------------------------------------------------

         MAXTYP = 21;
         NTYPES = MIN( MAXTYP, NTYPES );
         if ( NTYPES < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRED( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrvsx(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 4 ), NMAX, A( 1, 5 ), RESULT, WORK, LWORK, IWORK, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SGEESX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'SGG' ) ) {

         // -------------------------------------------------
         // SGG:  Generalized Nonsymmetric Eigenvalue Problem
         // -------------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NS    = number of shifts
            // MAXB  = minimum submatrix size
            // IACC22: structured matrix multiply
            // NBCOL = minimum column dimension for blocks

         MAXTYP = 26;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1,1);
         if (TSTCHK && TSTERR) CALL SERRGG( C3, NOUT );
         for (I = 1; I <= NPARMS; I++) { // 350
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(4, NSVAL( I ) );
            xlaenv(8, MXBVAL( I ) );
            xlaenv(16, IACC22( I ) );
            xlaenv(5, NBCOL( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 340
                  ISEED( K ) = IOLDSD( K );
               } // 340
            }
            WRITE( NOUT, FMT = 9996 )C3, NBVAL( I ), NBMIN( I ), NSVAL( I ), MXBVAL( I ), IACC22( I ), NBCOL( I );
            TSTDIF = false;
            THRSHN = 10.;
            if ( TSTCHK ) {
               schkgg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, TSTDIF, THRSHN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), A( 1, 8 ), A( 1, 9 ), NMAX, A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 13 ), A( 1, 14 ), WORK, LWORK, LOGWRK, RESULT, INFO );
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCHKGG', INFO;
            }
         } // 350

      } else if ( LSAMEN( 3, C3, 'SGS' ) ) {

         // -------------------------------------------------
         // SGS:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGES (Schur form)
         // -------------------------------------------------

         MAXTYP = 26;
         NTYPES = MIN( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRGG( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrges(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, RESULT, LOGWRK, INFO );

            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRGES', INFO;

      // Blocked version

            xlaenv(16,1);
            sdrges3(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, RESULT, LOGWRK, INFO );

            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRGES3', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( SGX ) {

         // -------------------------------------------------
         // SGX:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGESX (Schur form and condition numbers)
         // -------------------------------------------------

         MAXTYP = 5;
         NTYPES = MAXTYP;
         if ( NN < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRGG( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            xlaenv(5, 2 );
            sdrgsx(NN, NCMAX, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), C( 1, 1 ), NCMAX*NCMAX, A( 1, 12 ), WORK, LWORK, IWORK, LIWORK, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRGSX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'SGV' ) ) {

         // -------------------------------------------------
         // SGV:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGEV (Eigenvalue/vector form)
         // -------------------------------------------------

         MAXTYP = 26;
         NTYPES = MIN( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRGG( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrgev(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), WORK, LWORK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRGEV', INFO;

// Blocked version

            sdrgev3(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), WORK, LWORK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRGEV3', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( SXV ) {

         // -------------------------------------------------
         // SXV:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGEVX (eigenvalue/vector with condition numbers)
         // -------------------------------------------------

         MAXTYP = 2;
         NTYPES = MAXTYP;
         if ( NN < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) CALL SERRGG( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            sdrgvx(NN, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), A( 1, 5 ), A( 1, 6 ), IWORK( 1 ), IWORK( 2 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), WORK, LWORK, IWORK( 3 ), LIWORK-2, RESULT, LOGWRK, INFO );

            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SDRGVX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'SSB' ) ) {

         // ------------------------------
         // SSB:  Symmetric Band Reduction
         // ------------------------------

         MAXTYP = 15;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         if (TSTERR) CALL SERRST( 'SSB', NOUT );
          // CALL SCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
      // $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
         schksb2stg(NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCHKSB', INFO;

      } else if ( LSAMEN( 3, C3, 'SBB' ) ) {

         // ------------------------------
         // SBB:  General Band Reduction
         // ------------------------------

         MAXTYP = 15;
         NTYPES = MIN( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         for (I = 1; I <= NPARMS; I++) { // 370
            NRHS = NSVAL( I );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 360
                  ISEED( K ) = IOLDSD( K );
               } // 360
            }
            WRITE( NOUT, FMT = 9966 )C3, NRHS;
            schkbb(NN, MVAL, NVAL, NK, KVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), 2*NMAX, D( 1, 1 ), D( 1, 2 ), A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), WORK, LWORK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCHKBB', INFO;
         } // 370

      } else if ( LSAMEN( 3, C3, 'GLM' ) ) {

         // -----------------------------------------
         // GLM:  Generalized Linear Regression Model
         // -----------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) CALL SERRGG( 'GLM', NOUT );
         sckglm(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCKGLM', INFO;

      } else if ( LSAMEN( 3, C3, 'GQR' ) ) {

         // ------------------------------------------
         // GQR:  Generalized QR and RQ factorizations
         // ------------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) CALL SERRGG( 'GQR', NOUT );
         sckgqr(NN, MVAL, NN, PVAL, NN, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), TAUA, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), B( 1, 5 ), TAUB, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCKGQR', INFO;

      } else if ( LSAMEN( 3, C3, 'GSV' ) ) {

         // ----------------------------------------------
         // GSV:  Generalized Singular Value Decomposition
         // ----------------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) CALL SERRGG( 'GSV', NOUT );
         sckgsv(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB, B( 1, 4 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCKGSV', INFO;

      } else if ( LSAMEN( 3, C3, 'CSD' ) ) {

         // ----------------------------------------------
         // CSD:  CS Decomposition
         // ----------------------------------------------

         xlaenv(1,1);
         if (TSTERR) CALL SERRGG( 'CSD', NOUT );
         sckcsd(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCKCSD', INFO;

      } else if ( LSAMEN( 3, C3, 'LSE' ) ) {

         // --------------------------------------
         // LSE:  Constrained Linear Least Squares
         // --------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) CALL SERRGG( 'LSE', NOUT );
         scklse(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'SCKLSE', INFO;

      } else {
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = 9992 )C3;
      }
      IF( !( SGX || SXV ) ) GO TO 190;
      } // 380
      WRITE( NOUT, FMT = 9994 );
      S2 = SECOND( );
      WRITE( NOUT, FMT = 9993 )S2 - S1;

      DEALLOCATE (A, STAT = AllocateStatus);
      DEALLOCATE (B, STAT = AllocateStatus);
      DEALLOCATE (C, STAT = AllocateStatus);
      DEALLOCATE (WORK,  STAT = AllocateStatus);

 9999 FORMAT( / ' Execution not attempted due to input errors' );
 9997 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4 );
 9996 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NS =', I4, ', MAXB =', I4, ', IACC22 =', I4, ', NBCOL =', I4 );
 9995 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, ', NRHS =', I4 );
 9994 FORMAT( / / ' End of tests' );
 9993 FORMAT( ' Total time used = ', F12.2, ' seconds', / );
 9992 FORMAT( 1X, A3, ':  Unrecognized path name' );
 9991 FORMAT( / / ' *** Invalid int     value in column ', I2,; ' of input', ' line:', / A79 )
 9990 FORMAT( / / 1X, A3, ' routines were not tested' );
 9989 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be >=', I6 )
 9988 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be <=', I6 )
 9987 FORMAT( ' Tests of the Nonsymmetric Eigenvalue Problem routines' );
 9986 FORMAT( ' Tests of the Symmetric Eigenvalue Problem routines' );
 9985 FORMAT( ' Tests of the Singular Value Decomposition routines' );
 9984 FORMAT( / ' The following parameter values will be used:' );
 9983 FORMAT( 4X, A, 10I6, / 10X, 10I6 );
 9982 FORMAT( / ' Routines pass computational tests if test ratio is ', 'less than', F8.2, / );
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 );
 9980 FORMAT( ' *** Error code from ', A, ' = ', I4 );
 9979 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', / '    SGEEV (eigenvalues and eigevectors)' );
 9978 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', / '    SGEES (Schur form)' );
 9977 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', ' Driver', / '    SGEEVX (eigenvalues, eigenvectors and', ' condition numbers)' );
 9976 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', ' Driver', / '    SGEESX (Schur form and condition', ' numbers)' );
 9975 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem routines' );
 9974 FORMAT( ' Tests of SSBTRD', / ' (reduction of a symmetric band ', 'matrix to tridiagonal form)' );
 9973 FORMAT( / 1X, 71( '-' ) );
 9972 FORMAT( / ' LAPACK VERSION ', I1, '.', I1, '.', I1 );
 9971 FORMAT( / ' Tests of the Generalized Linear Regression Model ', 'routines' );
 9970 FORMAT( / ' Tests of the Generalized QR and RQ routines' );
 9969 FORMAT( / ' Tests of the Generalized Singular Value', ' Decomposition routines' );
 9968 FORMAT( / ' Tests of the Linear Least Squares routines' );
 9967 FORMAT( ' Tests of SGBBRD', / ' (reduction of a general band ', 'matrix to real bidiagonal form)' );
 9966 FORMAT( / / 1X, A3, ':  NRHS =', I4 );
 9965 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Expert Driver SGGESX' );
 9964 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Driver SGGES' );
 9963 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Driver SGGEV' );
 9962 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Expert Driver SGGEVX' );
 9961 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, ', INMIN=', I4, ', INWIN =', I4, ', INIBL =', I4, ', ISHFTS =', I4, ', IACC22 =', I4);
 9960 FORMAT( / ' Tests of the CS Decomposition routines' );

      // End of SCHKEE

      }
