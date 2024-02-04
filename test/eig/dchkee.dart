      void main() {
// #if defined(_OPENMP)
      use omp_lib;
// #endif

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Parameters ..
      const              NMAX = 132 ;
      const              NCMAX = 20 ;
      const              NEED = 14 ;
      const              LWORK = NMAX*( 5*NMAX+5 )+1 ;
      const              LIWORK = NMAX*( 5*NMAX+20 ) ;
      const              MAXIN = 20 ;
      const              MAXT = 30 ;
      const              NIN = 5, NOUT = 6 ;
      // ..
      // .. Local Scalars ..
      bool               CSD, DBB, DGG, DSB, FATAL, GLM, GQR, GSV, LSE, NEP, DBK, DBL, SEP, DES, DEV, DGK, DGL, DGS, DGV, DGX, DSX, SVD, DVX, DXV, TSTCHK, TSTDIF, TSTDRV, TSTERR;
      String             C1;
      String             C3, PATH;
      String             VNAME;
      String             INTSTR;
      String             LINE;
      int                I, I1, IC, INFO, ITMP, K, LENP, MAXTYP, NEWSD, NK, NN, NPARMS, NRHS, NTYPES, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      int    *4          N_THREADS, ONE_THREAD;
      double             EPS, S1, S2, THRESH, THRSHN;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MAXT ), LOGWRK( NMAX );
      int                IOLDSD( 4 ), ISEED( 4 ), IWORK( LIWORK ), KVAL( MAXIN ), MVAL( MAXIN ), MXBVAL( MAXIN ), NBCOL( MAXIN ), NBMIN( MAXIN ), NBVAL( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), PVAL( MAXIN );
      int                INMIN( MAXIN ), INWIN( MAXIN ), INIBL( MAXIN ), ISHFTS( MAXIN ), IACC22( MAXIN );
      double             D( NMAX, 12 ), RESULT( 500 ), TAUA( NMAX ), TAUB( NMAX ), X( 5*NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      double          , DIMENSION(:), ALLOCATABLE :: WORK;
      double          , DIMENSION(:,:), ALLOCATABLE :: A, B, C;
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      //- double             DLAMCH, DSECND;
      // EXTERNAL LSAMEN, DLAMCH, DSECND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, DCHKBB, DCHKBD, DCHKBK, DCHKBL, DCHKEC, DCHKGG, DCHKGK, DCHKGL, DCHKHS, DCHKSB, DCHKST, DCKCSD, DCKGLM, DCKGQR, DCKGSV, DCKLSE, DDRGES, DDRGEV, DDRGSX, DDRGVX, DDRVBD, DDRVES, DDRVEV, DDRVSG, DDRVST, DDRVSX, DDRVVX, DERRBD, DERRED, DERRGG, DERRHS, DERRST, ILAVER, XLAENV, DDRGES3, DDRGEV3, DCHKST2STG, DDRVST2STG, DCHKSB2STG, DDRVSG2STG
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
      double             SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Common blocks ..
      // COMMON / CENVIR / NPROC, NSHIFT, MAXB
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // COMMON / CLAENV / IPARMS
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
      S1 = DSECND( );
      FATAL = false;
      NUNIT = NOUT;

      // Return to here to read multiple sets of data

      } // 10

      // Read the first line and set the 3-character test path

      READ( NIN, FMT = '(A80)', END = 380 )LINE;
      PATH = LINE( 1: 3 );
      NEP = LSAMEN( 3, PATH, 'NEP' ) || LSAMEN( 3, PATH, 'DHS' );
      SEP = LSAMEN( 3, PATH, 'SEP' ) || LSAMEN( 3, PATH, 'DST' ) || LSAMEN( 3, PATH, 'DSG' ) || LSAMEN( 3, PATH, 'SE2' );
      SVD = LSAMEN( 3, PATH, 'SVD' ) || LSAMEN( 3, PATH, 'DBD' );
      DEV = LSAMEN( 3, PATH, 'DEV' );
      DES = LSAMEN( 3, PATH, 'DES' );
      DVX = LSAMEN( 3, PATH, 'DVX' );
      DSX = LSAMEN( 3, PATH, 'DSX' );
      DGG = LSAMEN( 3, PATH, 'DGG' );
      DGS = LSAMEN( 3, PATH, 'DGS' );
      DGX = LSAMEN( 3, PATH, 'DGX' );
      DGV = LSAMEN( 3, PATH, 'DGV' );
      DXV = LSAMEN( 3, PATH, 'DXV' );
      DSB = LSAMEN( 3, PATH, 'DSB' );
      DBB = LSAMEN( 3, PATH, 'DBB' );
      GLM = LSAMEN( 3, PATH, 'GLM' );
      GQR = LSAMEN( 3, PATH, 'GQR' ) || LSAMEN( 3, PATH, 'GRQ' );
      GSV = LSAMEN( 3, PATH, 'GSV' );
      CSD = LSAMEN( 3, PATH, 'CSD' );
      LSE = LSAMEN( 3, PATH, 'LSE' );
      DBL = LSAMEN( 3, PATH, 'DBL' );
      DBK = LSAMEN( 3, PATH, 'DBK' );
      DGL = LSAMEN( 3, PATH, 'DGL' );
      DGK = LSAMEN( 3, PATH, 'DGK' );

      // Report values of parameters.

      if ( PATH == '   ' ) {
         GO TO 10;
      } else if ( NEP ) {
         WRITE( NOUT, FMT = 9987 );
      } else if ( SEP ) {
         WRITE( NOUT, FMT = 9986 );
      } else if ( SVD ) {
         WRITE( NOUT, FMT = 9985 );
      } else if ( DEV ) {
         WRITE( NOUT, FMT = 9979 );
      } else if ( DES ) {
         WRITE( NOUT, FMT = 9978 );
      } else if ( DVX ) {
         WRITE( NOUT, FMT = 9977 );
      } else if ( DSX ) {
         WRITE( NOUT, FMT = 9976 );
      } else if ( DGG ) {
         WRITE( NOUT, FMT = 9975 );
      } else if ( DGS ) {
         WRITE( NOUT, FMT = 9964 );
      } else if ( DGX ) {
         WRITE( NOUT, FMT = 9965 );
      } else if ( DGV ) {
         WRITE( NOUT, FMT = 9963 );
      } else if ( DXV ) {
         WRITE( NOUT, FMT = 9962 );
      } else if ( DSB ) {
         WRITE( NOUT, FMT = 9974 );
      } else if ( DBB ) {
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
      } else if ( DBL ) {

         // DGEBAL:  Balancing

         dchkbl(NIN, NOUT );
         GO TO 10;
      } else if ( DBK ) {

         // DGEBAK:  Back transformation

         dchkbk(NIN, NOUT );
         GO TO 10;
      } else if ( DGL ) {

         // DGGBAL:  Balancing

         dchkgl(NIN, NOUT );
         GO TO 10;
      } else if ( DGK ) {

         // DGGBAK:  Back transformation

         dchkgk(NIN, NOUT );
         GO TO 10;
      } else if ( LSAMEN( 3, PATH, 'DEC' ) ) {

         // DEC:  Eigencondition estimation

         READ( NIN, FMT = * )THRESH;
         xlaenv(1, 1 );
         xlaenv(12, 11 );
         xlaenv(13, 2 );
         xlaenv(14, 0 );
         xlaenv(15, 2 );
         xlaenv(16, 2 );
         TSTERR = true;
         dchkec(THRESH, TSTERR, NIN, NOUT );
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

      if ( !( DGX || DXV ) ) {
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

      if ( SVD || DBB || GLM || GQR || GSV || CSD || LSE ) {
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
            NVAL[I] = MVAL( I );
         } // 50
      }
      if ( !( DGX || DXV ) ) {
         WRITE( NOUT, FMT = 9983 )'N:    ', ( NVAL( I ), I = 1, NN );
      } else {
         WRITE( NOUT, FMT = 9983 )'N:    ', NN;
      }

      // Read the number of values of K, followed by the values of K

      if ( DSB || DBB ) {
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

      if ( DEV || DES || DVX || DSX ) {

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
         xlaenv(12, max( 11, INMIN( 1 ) ) );
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

      } else if ( DGS || DGX || DGV || DXV ) {

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

      } else if ( !DSB && !GLM && !GQR && !GSV && !CSD && !LSE ) {

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

         if ( !DBB ) {
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

         if ( NEP || SEP || SVD || DGG ) {
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
               NBMIN[I] = 1;
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
               NXVAL[I] = 1;
            } // 110
         }

         // Read the values of NSHIFT (if DGG) or NRHS (if SVD
         // or DBB).

         if ( SVD || DBB || DGG ) {
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
               NSVAL[I] = 1;
            } // 130
         }

         // Read the values for MAXB.

         if ( DGG ) {
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
               MXBVAL[I] = 1;
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
               INMIN[I] = 1;
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
               INWIN[I] = 1;
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
               INIBL[I] = 1;
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
               ISHFTS[I] = 1;
            } // 610
         }

         // Read the values for IACC22.

         if ( NEP || DGG ) {
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
               IACC22[I] = 1;
            } // 630
         }

         // Read the values for NBCOL.

         if ( DGG ) {
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
               NBCOL[I] = 1;
            } // 170
         }
      }

      // Calculate and print the machine dependent constants.

      WRITE( NOUT, FMT = * );
      EPS = DLAMCH( 'Underflow threshold' );
      WRITE( NOUT, FMT = 9981 )'underflow', EPS;
      EPS = DLAMCH( 'Overflow threshold' );
      WRITE( NOUT, FMT = 9981 )'overflow ', EPS;
      EPS = DLAMCH( 'Epsilon' );
      WRITE( NOUT, FMT = 9981 )'precision', EPS;

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH;
      WRITE( NOUT, FMT = 9982 )THRESH;
      if ( SEP || SVD || DGG ) {

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
         ISEED[I] = IOLDSD( I );
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

      if ( !( DGX || DXV ) ) {

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

         if ( !( DEV || DES || DVX || DSX || DGV || DGS ) && NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
            GO TO 200;
         }

      } else {
         if (DXV) C3 = 'DXV';
         IF( DGX ) C3 = 'DGX';
      }

      // Reset the random number seed.

      if ( NEWSD == 0 ) {
         for (K = 1; K <= 4; K++) { // 250
            ISEED[K] = IOLDSD( K );
         } // 250
      }

      if ( LSAMEN( 3, C3, 'DHS' ) || LSAMEN( 3, C3, 'NEP' ) ) {

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
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         if (TSTERR) derrhs( 'DHSEQR', NOUT );
         for (I = 1; I <= NPARMS; I++) { // 270
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );
            xlaenv(12, max( 11, INMIN( I ) ) );
            xlaenv(13, INWIN( I ) );
            xlaenv(14, INIBL( I ) );
            xlaenv(15, ISHFTS( I ) );
            xlaenv(16, IACC22( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 260
                  ISEED[K] = IOLDSD( K );
               } // 260
            }
            WRITE( NOUT, FMT = 9961 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), max( 11, INMIN(I)), INWIN( I ), INIBL( I ), ISHFTS( I ), IACC22( I );
            dchkhs(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), A( 1, 7 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 8 ), A( 1, 9 ), A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), D( 1, 7 ), WORK, LWORK, IWORK, LOGWRK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCHKHS', INFO;
         } // 270

      } else if ( LSAMEN( 3, C3, 'DST' ) || LSAMEN( 3, C3, 'SEP' ) || LSAMEN( 3, C3, 'SE2' ) ) {

         // ----------------------------------
         // SEP:  Symmetric Eigenvalue Problem
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21;
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         xlaenv(9, 25 );
         if ( TSTERR ) {
// #if defined(_OPENMP)
            N_THREADS = OMP_GET_MAX_THREADS();
            ONE_THREAD = 1;
            omp_set_num_threads(ONE_THREAD);
// #endif
            derrst('DST', NOUT );
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
                  ISEED[K] = IOLDSD( K );
               } // 280
            }
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I );
            if ( TSTCHK ) {
               if ( LSAMEN( 3, C3, 'SE2' ) ) {
               dchkst2stg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               } else {
               dchkst(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               }
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCHKST', INFO;
            }
            if ( TSTDRV ) {
               if ( LSAMEN( 3, C3, 'SE2' ) ) {
               ddrvst2stg(NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 2 ), NMAX, A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               } else {
               ddrvst(NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 2 ), NMAX, A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               }
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRVST', INFO;
            }
         } // 290

      } else if ( LSAMEN( 3, C3, 'DSG' ) ) {

         // ----------------------------------------------
         // DSG:  Symmetric Generalized Eigenvalue Problem
         // ----------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21;
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(9, 25 );
         for (I = 1; I <= NPARMS; I++) { // 310
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 300
                  ISEED[K] = IOLDSD( K );
               } // 300
            }
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I );
            if ( TSTCHK ) {
                // CALL DDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
      // $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
      // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
      // $                      LWORK, IWORK, LIWORK, RESULT, INFO )
               ddrvsg2stg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, D( 1, 3 ), D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO );
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRVSG', INFO;
            }
         } // 310

      } else if ( LSAMEN( 3, C3, 'DBD' ) || LSAMEN( 3, C3, 'SVD' ) ) {

         // ----------------------------------
         // SVD:  Singular Value Decomposition
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point
            // NRHS  = number of right hand sides

         MAXTYP = 16;
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1, 1 );
         xlaenv(9, 25 );

         // Test the error exits

         if (TSTERR && TSTCHK) derrbd( 'DBD', NOUT );
         IF( TSTERR && TSTDRV ) derred( 'DBD', NOUT );

         for (I = 1; I <= NPARMS; I++) { // 330
            NRHS = NSVAL( I );
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(3, NXVAL( I ) );
            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 320
                  ISEED[K] = IOLDSD( K );
               } // 320
            }
            WRITE( NOUT, FMT = 9995 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), NRHS;
            if ( TSTCHK ) {
               dchkbd(NN, MVAL, NVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 2 ), NMAX, A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), A( 1, 8 ), WORK, LWORK, IWORK, NOUT, INFO );
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCHKBD', INFO;
            }
            if (TSTDRV) ddrvbd( NN, MVAL, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, IWORK, NOUT, INFO );
         } // 330

      } else if ( LSAMEN( 3, C3, 'DEV' ) ) {

         // --------------------------------------------
         // DEV:  Nonsymmetric Eigenvalue Problem Driver
               // DGEEV (eigenvalues and eigenvectors)
         // --------------------------------------------

         MAXTYP = 21;
         NTYPES = min( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derred( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrvev(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, RESULT, WORK, LWORK, IWORK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DGEEV', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'DES' ) ) {

         // --------------------------------------------
         // DES:  Nonsymmetric Eigenvalue Problem Driver
               // DGEES (Schur form)
         // --------------------------------------------

         MAXTYP = 21;
         NTYPES = min( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derred( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrves(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 4 ), NMAX, RESULT, WORK, LWORK, IWORK, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DGEES', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'DVX' ) ) {

         // --------------------------------------------------------------
         // DVX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // DGEEVX (eigenvalues, eigenvectors and condition numbers)
         // --------------------------------------------------------------

         MAXTYP = 21;
         NTYPES = min( MAXTYP, NTYPES );
         if ( NTYPES < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derred( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrvvx(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), D( 1, 12 ), RESULT, WORK, LWORK, IWORK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DGEEVX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'DSX' ) ) {

         // ---------------------------------------------------
         // DSX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // DGEESX (Schur form and condition numbers)
         // ---------------------------------------------------

         MAXTYP = 21;
         NTYPES = min( MAXTYP, NTYPES );
         if ( NTYPES < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derred( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrvsx(NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 4 ), NMAX, A( 1, 5 ), RESULT, WORK, LWORK, IWORK, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DGEESX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'DGG' ) ) {

         // -------------------------------------------------
         // DGG:  Generalized Nonsymmetric Eigenvalue Problem
         // -------------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NS    = number of shifts
            // MAXB  = minimum submatrix size
            // IACC22: structured matrix multiply
            // NBCOL = minimum column dimension for blocks

         MAXTYP = 26;
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         xlaenv(1,1);
         if (TSTCHK && TSTERR) derrgg( C3, NOUT );
         for (I = 1; I <= NPARMS; I++) { // 350
            xlaenv(1, NBVAL( I ) );
            xlaenv(2, NBMIN( I ) );
            xlaenv(4, NSVAL( I ) );
            xlaenv(8, MXBVAL( I ) );
            xlaenv(16, IACC22( I ) );
            xlaenv(5, NBCOL( I ) );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 340
                  ISEED[K] = IOLDSD( K );
               } // 340
            }
            WRITE( NOUT, FMT = 9996 )C3, NBVAL( I ), NBMIN( I ), NSVAL( I ), MXBVAL( I ), IACC22( I ), NBCOL( I );
            TSTDIF = false;
            THRSHN = 10.0;
            if ( TSTCHK ) {
               dchkgg(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, TSTDIF, THRSHN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), A( 1, 8 ), A( 1, 9 ), NMAX, A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 13 ), A( 1, 14 ), WORK, LWORK, LOGWRK, RESULT, INFO );
               if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCHKGG', INFO;
            }
         } // 350

      } else if ( LSAMEN( 3, C3, 'DGS' ) ) {

         // -------------------------------------------------
         // DGS:  Generalized Nonsymmetric Eigenvalue Problem
               // DGGES (Schur form)
         // -------------------------------------------------

         MAXTYP = 26;
         NTYPES = min( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derrgg( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrges(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, RESULT, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRGES', INFO;

      // Blocked version

            xlaenv(16, 2);
            ddrges3(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, RESULT, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRGES3', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( DGX ) {

         // -------------------------------------------------
         // DGX:  Generalized Nonsymmetric Eigenvalue Problem
               // DGGESX (Schur form and condition numbers)
         // -------------------------------------------------

         MAXTYP = 5;
         NTYPES = MAXTYP;
         if ( NN < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derrgg( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            xlaenv(5, 2 );
            ddrgsx(NN, NCMAX, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), C( 1, 1 ), NCMAX*NCMAX, A( 1, 12 ), WORK, LWORK, IWORK, LIWORK, LOGWRK, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRGSX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'DGV' ) ) {

         // -------------------------------------------------
         // DGV:  Generalized Nonsymmetric Eigenvalue Problem
               // DGGEV (Eigenvalue/vector form)
         // -------------------------------------------------

         MAXTYP = 26;
         NTYPES = min( MAXTYP, NTYPES );
         if ( NTYPES <= 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derrgg( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrgev(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), WORK, LWORK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRGEV', INFO;

      // Blocked version

            ddrgev3(NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), WORK, LWORK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRGEV3', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( DXV ) {

         // -------------------------------------------------
         // DXV:  Generalized Nonsymmetric Eigenvalue Problem
               // DGGEVX (eigenvalue/vector with condition numbers)
         // -------------------------------------------------

         MAXTYP = 2;
         NTYPES = MAXTYP;
         if ( NN < 0 ) {
            WRITE( NOUT, FMT = 9990 )C3;
         } else {
            if (TSTERR) derrgg( C3, NOUT );
            alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
            ddrgvx(NN, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), A( 1, 5 ), A( 1, 6 ), IWORK( 1 ), IWORK( 2 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), WORK, LWORK, IWORK( 3 ), LIWORK-2, RESULT, LOGWRK, INFO );

            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DDRGVX', INFO;
         }
         WRITE( NOUT, FMT = 9973 );
         GO TO 10;

      } else if ( LSAMEN( 3, C3, 'DSB' ) ) {

         // ------------------------------
         // DSB:  Symmetric Band Reduction
         // ------------------------------

         MAXTYP = 15;
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         if (TSTERR) derrst( 'DSB', NOUT );
          // CALL DCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
      // $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
         dchksb2stg(NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCHKSB', INFO;

      } else if ( LSAMEN( 3, C3, 'DBB' ) ) {

         // ------------------------------
         // DBB:  General Band Reduction
         // ------------------------------

         MAXTYP = 15;
         NTYPES = min( MAXTYP, NTYPES );
         alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT );
         for (I = 1; I <= NPARMS; I++) { // 370
            NRHS = NSVAL( I );

            if ( NEWSD == 0 ) {
               for (K = 1; K <= 4; K++) { // 360
                  ISEED[K] = IOLDSD( K );
               } // 360
            }
            WRITE( NOUT, FMT = 9966 )C3, NRHS;
            dchkbb(NN, MVAL, NVAL, NK, KVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), 2*NMAX, D( 1, 1 ), D( 1, 2 ), A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), WORK, LWORK, RESULT, INFO );
            if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCHKBB', INFO;
         } // 370

      } else if ( LSAMEN( 3, C3, 'GLM' ) ) {

         // -----------------------------------------
         // GLM:  Generalized Linear Regression Model
         // -----------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) derrgg( 'GLM', NOUT );
         dckglm(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCKGLM', INFO;

      } else if ( LSAMEN( 3, C3, 'GQR' ) ) {

         // ------------------------------------------
         // GQR:  Generalized QR and RQ factorizations
         // ------------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) derrgg( 'GQR', NOUT );
         dckgqr(NN, MVAL, NN, PVAL, NN, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), TAUA, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), B( 1, 5 ), TAUB, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCKGQR', INFO;

      } else if ( LSAMEN( 3, C3, 'GSV' ) ) {

         // ----------------------------------------------
         // GSV:  Generalized Singular Value Decomposition
         // ----------------------------------------------

         xlaenv(1,1);
         if (TSTERR) derrgg( 'GSV', NOUT );
         dckgsv(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB, B( 1, 4 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCKGSV', INFO;

      } else if ( LSAMEN( 3, C3, 'CSD' ) ) {

         // ----------------------------------------------
         // CSD:  CS Decomposition
         // ----------------------------------------------

         xlaenv(1,1);
         if (TSTERR) derrgg( 'CSD', NOUT );
         dckcsd(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCKCSD', INFO;

      } else if ( LSAMEN( 3, C3, 'LSE' ) ) {

         // --------------------------------------
         // LSE:  Constrained Linear Least Squares
         // --------------------------------------

         xlaenv(1, 1 );
         if (TSTERR) derrgg( 'LSE', NOUT );
         dcklse(NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, D( 1, 1 ), NIN, NOUT, INFO );
         if (INFO != 0) WRITE( NOUT, FMT = 9980 )'DCKLSE', INFO;

      } else {
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = 9992 )C3;
      }
      if( !( DGX || DXV ) ) GO TO 190;
      } // 380
      WRITE( NOUT, FMT = 9994 );
      S2 = DSECND( );
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
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 );
 9980 FORMAT( ' *** Error code from ', A, ' = ', I4 );
 9979 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', / '    DGEEV (eigenvalues and eigevectors)' );
 9978 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver', / '    DGEES (Schur form)' );
 9977 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', ' Driver', / '    DGEEVX (eigenvalues, eigenvectors and', ' condition numbers)' );
 9976 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert', ' Driver', / '    DGEESX (Schur form and condition', ' numbers)' );
 9975 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem routines' );
 9974 FORMAT( ' Tests of DSBTRD', / ' (reduction of a symmetric band ', 'matrix to tridiagonal form)' );
 9973 FORMAT( / 1X, 71( '-' ) );
 9972 FORMAT( / ' LAPACK VERSION ', I1, '.', I1, '.', I1 );
 9971 FORMAT( / ' Tests of the Generalized Linear Regression Model ', 'routines' );
 9970 FORMAT( / ' Tests of the Generalized QR and RQ routines' );
 9969 FORMAT( / ' Tests of the Generalized Singular Value', ' Decomposition routines' );
 9968 FORMAT( / ' Tests of the Linear Least Squares routines' );
 9967 FORMAT( ' Tests of DGBBRD', / ' (reduction of a general band ', 'matrix to real bidiagonal form)' );
 9966 FORMAT( / / 1X, A3, ':  NRHS =', I4 );
 9965 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Expert Driver DGGESX' );
 9964 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Driver DGGES' );
 9963 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Driver DGGEV' );
 9962 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ', 'Problem Expert Driver DGGEVX' );
 9961 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4, ', INMIN=', I4, ', INWIN =', I4, ', INIBL =', I4, ', ISHFTS =', I4, ', IACC22 =', I4);
 9960 FORMAT( / ' Tests of the CS Decomposition routines' );
      }
