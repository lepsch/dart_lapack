      PROGRAM SCHKEE

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
      REAL               EPS, S1, S2, THRESH, THRSHN
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MAXT ), LOGWRK( NMAX );
      int                IOLDSD( 4 ), ISEED( 4 ), IWORK( LIWORK ), KVAL( MAXIN ), MVAL( MAXIN ), MXBVAL( MAXIN ), NBCOL( MAXIN ), NBMIN( MAXIN ), NBVAL( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), PVAL( MAXIN );
      int                INMIN( MAXIN ), INWIN( MAXIN ), INIBL( MAXIN ), ISHFTS( MAXIN ), IACC22( MAXIN );
      REAL               D( NMAX, 12 ), RESULT( 500 ), TAUA( NMAX ), TAUB( NMAX ), X( 5*NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      REAL, DIMENSION(:), ALLOCATABLE :: WORK
      REAL, DIMENSION(:,:), ALLOCATABLE :: A, B, C
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      REAL               SECOND, SLAMCH
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

      ALLOCATE ( A(NMAX*NMAX,NEED), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( B(NMAX*NMAX,5), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( C(NCMAX*NCMAX,NCMAX*NCMAX), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( WORK(LWORK), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      // ..
      // .. Executable Statements ..

      A = 0.0
      B = 0.0
      C = 0.0
      D = 0.0
      S1 = SECOND( )
      FATAL = .FALSE.
      NUNIT = NOUT

      // Return to here to read multiple sets of data

   10 CONTINUE

      // Read the first line and set the 3-character test path

      READ( NIN, FMT = '(A80)', END = 380 )LINE
      PATH = LINE( 1: 3 )
      NEP = LSAMEN( 3, PATH, 'NEP' ) .OR. LSAMEN( 3, PATH, 'SHS' )
      SEP = LSAMEN( 3, PATH, 'SEP' ) .OR. LSAMEN( 3, PATH, 'SST' ) .OR. LSAMEN( 3, PATH, 'SSG' ) .OR. LSAMEN( 3, PATH, 'SE2' )
      SVD = LSAMEN( 3, PATH, 'SVD' ) .OR. LSAMEN( 3, PATH, 'DBD' )
      SVD = LSAMEN( 3, PATH, 'SVD' ) .OR. LSAMEN( 3, PATH, 'SBD' )
      SEV = LSAMEN( 3, PATH, 'SEV' )
      SES = LSAMEN( 3, PATH, 'SES' )
      SVX = LSAMEN( 3, PATH, 'SVX' )
      SSX = LSAMEN( 3, PATH, 'SSX' )
      SGG = LSAMEN( 3, PATH, 'SGG' )
      SGS = LSAMEN( 3, PATH, 'SGS' )
      SGX = LSAMEN( 3, PATH, 'SGX' )
      SGV = LSAMEN( 3, PATH, 'SGV' )
      SXV = LSAMEN( 3, PATH, 'SXV' )
      SSB = LSAMEN( 3, PATH, 'SSB' )
      SBB = LSAMEN( 3, PATH, 'SBB' )
      GLM = LSAMEN( 3, PATH, 'GLM' )
      GQR = LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' )
      GSV = LSAMEN( 3, PATH, 'GSV' )
      CSD = LSAMEN( 3, PATH, 'CSD' )
      LSE = LSAMEN( 3, PATH, 'LSE' )
      SBL = LSAMEN( 3, PATH, 'SBL' )
      SBK = LSAMEN( 3, PATH, 'SBK' )
      SGL = LSAMEN( 3, PATH, 'SGL' )
      SGK = LSAMEN( 3, PATH, 'SGK' )

      // Report values of parameters.

      IF( PATH.EQ.'   ' ) THEN
         GO TO 10
      ELSE IF( NEP ) THEN
         WRITE( NOUT, FMT = 9987 )
      ELSE IF( SEP ) THEN
         WRITE( NOUT, FMT = 9986 )
      ELSE IF( SVD ) THEN
         WRITE( NOUT, FMT = 9985 )
      ELSE IF( SEV ) THEN
         WRITE( NOUT, FMT = 9979 )
      ELSE IF( SES ) THEN
         WRITE( NOUT, FMT = 9978 )
      ELSE IF( SVX ) THEN
         WRITE( NOUT, FMT = 9977 )
      ELSE IF( SSX ) THEN
         WRITE( NOUT, FMT = 9976 )
      ELSE IF( SGG ) THEN
         WRITE( NOUT, FMT = 9975 )
      ELSE IF( SGS ) THEN
         WRITE( NOUT, FMT = 9964 )
      ELSE IF( SGX ) THEN
         WRITE( NOUT, FMT = 9965 )
      ELSE IF( SGV ) THEN
         WRITE( NOUT, FMT = 9963 )
      ELSE IF( SXV ) THEN
         WRITE( NOUT, FMT = 9962 )
      ELSE IF( SSB ) THEN
         WRITE( NOUT, FMT = 9974 )
      ELSE IF( SBB ) THEN
         WRITE( NOUT, FMT = 9967 )
      ELSE IF( GLM ) THEN
         WRITE( NOUT, FMT = 9971 )
      ELSE IF( GQR ) THEN
         WRITE( NOUT, FMT = 9970 )
      ELSE IF( GSV ) THEN
         WRITE( NOUT, FMT = 9969 )
      ELSE IF( CSD ) THEN
         WRITE( NOUT, FMT = 9960 )
      ELSE IF( LSE ) THEN
         WRITE( NOUT, FMT = 9968 )
      ELSE IF( SBL ) THEN

         // SGEBAL:  Balancing

         CALL SCHKBL( NIN, NOUT )
         GO TO 10
      ELSE IF( SBK ) THEN

         // SGEBAK:  Back transformation

         CALL SCHKBK( NIN, NOUT )
         GO TO 10
      ELSE IF( SGL ) THEN

         // SGGBAL:  Balancing

         CALL SCHKGL( NIN, NOUT )
         GO TO 10
      ELSE IF( SGK ) THEN

         // SGGBAK:  Back transformation

         CALL SCHKGK( NIN, NOUT )
         GO TO 10
      ELSE IF( LSAMEN( 3, PATH, 'SEC' ) ) THEN

         // SEC:  Eigencondition estimation

         READ( NIN, FMT = * )THRESH
         CALL XLAENV( 1, 1 )
         CALL XLAENV( 12, 11 )
         CALL XLAENV( 13, 2 )
         CALL XLAENV( 14, 0 )
         CALL XLAENV( 15, 2 )
         CALL XLAENV( 16, 2 )
         TSTERR = .TRUE.
         CALL SCHKEC( THRESH, TSTERR, NIN, NOUT )
         GO TO 10
      ELSE
         WRITE( NOUT, FMT = 9992 )PATH
         GO TO 10
      END IF
      CALL ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
      WRITE( NOUT, FMT = 9972 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH
      WRITE( NOUT, FMT = 9984 )

      // Read the number of values of M, P, and N.

      READ( NIN, FMT = * )NN
      IF( NN.LT.0 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      ELSE IF( NN.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9988 )'   NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      END IF

      // Read the values of M

      IF( .NOT.( SGX .OR. SXV ) ) THEN
         READ( NIN, FMT = * )( MVAL( I ), I = 1, NN )
         IF( SVD ) THEN
            VNAME = '    M '
         ELSE
            VNAME = '    N '
         END IF
         DO 20 I = 1, NN
            IF( MVAL( I ).LT.0 ) THEN
               WRITE( NOUT, FMT = 9989 )VNAME, MVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( MVAL( I ).GT.NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )VNAME, MVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
   20    CONTINUE
         WRITE( NOUT, FMT = 9983 )'M:    ', ( MVAL( I ), I = 1, NN )
      END IF

      // Read the values of P

      IF( GLM .OR. GQR .OR. GSV .OR. CSD .OR. LSE ) THEN
         READ( NIN, FMT = * )( PVAL( I ), I = 1, NN )
         DO 30 I = 1, NN
            IF( PVAL( I ).LT.0 ) THEN
               WRITE( NOUT, FMT = 9989 )' P  ', PVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( PVAL( I ).GT.NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )' P  ', PVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
   30    CONTINUE
         WRITE( NOUT, FMT = 9983 )'P:    ', ( PVAL( I ), I = 1, NN )
      END IF

      // Read the values of N

      IF( SVD .OR. SBB .OR. GLM .OR. GQR .OR. GSV .OR. CSD .OR. LSE ) THEN
         READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
         DO 40 I = 1, NN
            IF( NVAL( I ).LT.0 ) THEN
               WRITE( NOUT, FMT = 9989 )'    N ', NVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( NVAL( I ).GT.NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'    N ', NVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
   40    CONTINUE
      ELSE
         DO 50 I = 1, NN
            NVAL( I ) = MVAL( I )
   50    CONTINUE
      END IF
      IF( .NOT.( SGX .OR. SXV ) ) THEN
         WRITE( NOUT, FMT = 9983 )'N:    ', ( NVAL( I ), I = 1, NN )
      ELSE
         WRITE( NOUT, FMT = 9983 )'N:    ', NN
      END IF

      // Read the number of values of K, followed by the values of K

      IF( SSB .OR. SBB ) THEN
         READ( NIN, FMT = * )NK
         READ( NIN, FMT = * )( KVAL( I ), I = 1, NK )
         DO 60 I = 1, NK
            IF( KVAL( I ).LT.0 ) THEN
               WRITE( NOUT, FMT = 9989 )'    K ', KVAL( I ), 0
               FATAL = .TRUE.
            ELSE IF( KVAL( I ).GT.NMAX ) THEN
               WRITE( NOUT, FMT = 9988 )'    K ', KVAL( I ), NMAX
               FATAL = .TRUE.
            END IF
   60    CONTINUE
         WRITE( NOUT, FMT = 9983 )'K:    ', ( KVAL( I ), I = 1, NK )
      END IF

      IF( SEV .OR. SES .OR. SVX .OR. SSX ) THEN

         // For the nonsymmetric QR driver routines, only one set of
         // parameters is allowed.

         READ( NIN, FMT = * )NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), INMIN( 1 ), INWIN( 1 ), INIBL(1), ISHFTS(1), IACC22(1)
         IF( NBVAL( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( NBMIN( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( NXVAL( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( INMIN( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   INMIN ', INMIN( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( INWIN( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   INWIN ', INWIN( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( INIBL( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   INIBL ', INIBL( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( ISHFTS( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   ISHFTS ', ISHFTS( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( IACC22( 1 ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9989 )'   IACC22 ', IACC22( 1 ), 0
            FATAL = .TRUE.
         END IF
         CALL XLAENV( 1, NBVAL( 1 ) )
         CALL XLAENV( 2, NBMIN( 1 ) )
         CALL XLAENV( 3, NXVAL( 1 ) )
         CALL XLAENV(12, MAX( 11, INMIN( 1 ) ) )
         CALL XLAENV(13, INWIN( 1 ) )
         CALL XLAENV(14, INIBL( 1 ) )
         CALL XLAENV(15, ISHFTS( 1 ) )
         CALL XLAENV(16, IACC22( 1 ) )
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 )
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'INMIN:   ', INMIN( 1 )
         WRITE( NOUT, FMT = 9983 )'INWIN: ', INWIN( 1 )
         WRITE( NOUT, FMT = 9983 )'INIBL: ', INIBL( 1 )
         WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ISHFTS( 1 )
         WRITE( NOUT, FMT = 9983 )'IACC22: ', IACC22( 1 )

      ELSE IF( SGS .OR. SGX .OR. SGV .OR. SXV ) THEN

         // For the nonsymmetric generalized driver routines, only one set
         // of parameters is allowed.

         READ( NIN, FMT = * )NBVAL( 1 ), NBMIN( 1 ), NXVAL( 1 ), NSVAL( 1 ), MXBVAL( 1 )
         IF( NBVAL( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( NBMIN( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( NXVAL( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( 1 ), 1
            FATAL = .TRUE.
         ELSE IF( NSVAL( 1 ).LT.2 ) THEN
            WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( 1 ), 2
            FATAL = .TRUE.
         ELSE IF( MXBVAL( 1 ).LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( 1 ), 1
            FATAL = .TRUE.
         END IF
         CALL XLAENV( 1, NBVAL( 1 ) )
         CALL XLAENV( 2, NBMIN( 1 ) )
         CALL XLAENV( 3, NXVAL( 1 ) )
         CALL XLAENV( 4, NSVAL( 1 ) )
         CALL XLAENV( 8, MXBVAL( 1 ) )
         WRITE( NOUT, FMT = 9983 )'NB:   ', NBVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'NBMIN:', NBMIN( 1 )
         WRITE( NOUT, FMT = 9983 )'NX:   ', NXVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'NS:   ', NSVAL( 1 )
         WRITE( NOUT, FMT = 9983 )'MAXB: ', MXBVAL( 1 )

      ELSE IF( .NOT.SSB .AND. .NOT.GLM .AND. .NOT.GQR .AND. .NOT. GSV .AND. .NOT.CSD .AND. .NOT.LSE ) THEN

         // For the other paths, the number of parameters can be varied
         // from the input file.  Read the number of parameter values.

         READ( NIN, FMT = * )NPARMS
         IF( NPARMS.LT.1 ) THEN
            WRITE( NOUT, FMT = 9989 )'NPARMS', NPARMS, 1
            NPARMS = 0
            FATAL = .TRUE.
         ELSE IF( NPARMS.GT.MAXIN ) THEN
            WRITE( NOUT, FMT = 9988 )'NPARMS', NPARMS, MAXIN
            NPARMS = 0
            FATAL = .TRUE.
         END IF

         // Read the values of NB

         IF( .NOT.SBB ) THEN
            READ( NIN, FMT = * )( NBVAL( I ), I = 1, NPARMS )
            DO 70 I = 1, NPARMS
               IF( NBVAL( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )'   NB ', NBVAL( I ), 0
                  FATAL = .TRUE.
               ELSE IF( NBVAL( I ).GT.NMAX ) THEN
                  WRITE( NOUT, FMT = 9988 )'   NB ', NBVAL( I ), NMAX
                  FATAL = .TRUE.
               END IF
   70       CONTINUE
            WRITE( NOUT, FMT = 9983 )'NB:   ', ( NBVAL( I ), I = 1, NPARMS )
         END IF

         // Read the values of NBMIN

         IF( NEP .OR. SEP .OR. SVD .OR. SGG ) THEN
            READ( NIN, FMT = * )( NBMIN( I ), I = 1, NPARMS )
            DO 80 I = 1, NPARMS
               IF( NBMIN( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )'NBMIN ', NBMIN( I ), 0
                  FATAL = .TRUE.
               ELSE IF( NBMIN( I ).GT.NMAX ) THEN
                  WRITE( NOUT, FMT = 9988 )'NBMIN ', NBMIN( I ), NMAX
                  FATAL = .TRUE.
               END IF
   80       CONTINUE
            WRITE( NOUT, FMT = 9983 )'NBMIN:', ( NBMIN( I ), I = 1, NPARMS )
         ELSE
            DO 90 I = 1, NPARMS
               NBMIN( I ) = 1
   90       CONTINUE
         END IF

         // Read the values of NX

         IF( NEP .OR. SEP .OR. SVD ) THEN
            READ( NIN, FMT = * )( NXVAL( I ), I = 1, NPARMS )
            DO 100 I = 1, NPARMS
               IF( NXVAL( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )'   NX ', NXVAL( I ), 0
                  FATAL = .TRUE.
               ELSE IF( NXVAL( I ).GT.NMAX ) THEN
                  WRITE( NOUT, FMT = 9988 )'   NX ', NXVAL( I ), NMAX
                  FATAL = .TRUE.
               END IF
  100       CONTINUE
            WRITE( NOUT, FMT = 9983 )'NX:   ', ( NXVAL( I ), I = 1, NPARMS )
         ELSE
            DO 110 I = 1, NPARMS
               NXVAL( I ) = 1
  110       CONTINUE
         END IF

         // Read the values of NSHIFT (if SGG) or NRHS (if SVD
         // or SBB).

         IF( SVD .OR. SBB .OR. SGG ) THEN
            READ( NIN, FMT = * )( NSVAL( I ), I = 1, NPARMS )
            DO 120 I = 1, NPARMS
               IF( NSVAL( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )'   NS ', NSVAL( I ), 0
                  FATAL = .TRUE.
               ELSE IF( NSVAL( I ).GT.NMAX ) THEN
                  WRITE( NOUT, FMT = 9988 )'   NS ', NSVAL( I ), NMAX
                  FATAL = .TRUE.
               END IF
  120       CONTINUE
            WRITE( NOUT, FMT = 9983 )'NS:   ', ( NSVAL( I ), I = 1, NPARMS )
         ELSE
            DO 130 I = 1, NPARMS
               NSVAL( I ) = 1
  130       CONTINUE
         END IF

         // Read the values for MAXB.

         IF( SGG ) THEN
            READ( NIN, FMT = * )( MXBVAL( I ), I = 1, NPARMS )
            DO 140 I = 1, NPARMS
               IF( MXBVAL( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )' MAXB ', MXBVAL( I ), 0
                  FATAL = .TRUE.
               ELSE IF( MXBVAL( I ).GT.NMAX ) THEN
                  WRITE( NOUT, FMT = 9988 )' MAXB ', MXBVAL( I ), NMAX
                  FATAL = .TRUE.
               END IF
  140       CONTINUE
            WRITE( NOUT, FMT = 9983 )'MAXB: ', ( MXBVAL( I ), I = 1, NPARMS )
         ELSE
            DO 150 I = 1, NPARMS
               MXBVAL( I ) = 1
  150       CONTINUE
         END IF

         // Read the values for INMIN.

         IF( NEP ) THEN
            READ( NIN, FMT = * )( INMIN( I ), I = 1, NPARMS )
            DO 540 I = 1, NPARMS
               IF( INMIN( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )' INMIN ', INMIN( I ), 0
                  FATAL = .TRUE.
               END IF
  540       CONTINUE
            WRITE( NOUT, FMT = 9983 )'INMIN: ', ( INMIN( I ), I = 1, NPARMS )
         ELSE
            DO 550 I = 1, NPARMS
               INMIN( I ) = 1
  550       CONTINUE
         END IF

         // Read the values for INWIN.

         IF( NEP ) THEN
            READ( NIN, FMT = * )( INWIN( I ), I = 1, NPARMS )
            DO 560 I = 1, NPARMS
               IF( INWIN( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )' INWIN ', INWIN( I ), 0
                  FATAL = .TRUE.
               END IF
  560       CONTINUE
            WRITE( NOUT, FMT = 9983 )'INWIN: ', ( INWIN( I ), I = 1, NPARMS )
         ELSE
            DO 570 I = 1, NPARMS
               INWIN( I ) = 1
  570       CONTINUE
         END IF

         // Read the values for INIBL.

         IF( NEP ) THEN
            READ( NIN, FMT = * )( INIBL( I ), I = 1, NPARMS )
            DO 580 I = 1, NPARMS
               IF( INIBL( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )' INIBL ', INIBL( I ), 0
                  FATAL = .TRUE.
               END IF
  580       CONTINUE
            WRITE( NOUT, FMT = 9983 )'INIBL: ', ( INIBL( I ), I = 1, NPARMS )
         ELSE
            DO 590 I = 1, NPARMS
               INIBL( I ) = 1
  590       CONTINUE
         END IF

         // Read the values for ISHFTS.

         IF( NEP ) THEN
            READ( NIN, FMT = * )( ISHFTS( I ), I = 1, NPARMS )
            DO 600 I = 1, NPARMS
               IF( ISHFTS( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )' ISHFTS ', ISHFTS( I ), 0
                  FATAL = .TRUE.
               END IF
  600       CONTINUE
            WRITE( NOUT, FMT = 9983 )'ISHFTS: ', ( ISHFTS( I ), I = 1, NPARMS )
         ELSE
            DO 610 I = 1, NPARMS
               ISHFTS( I ) = 1
  610       CONTINUE
         END IF

         // Read the values for IACC22.

         IF( NEP .OR. SGG ) THEN
            READ( NIN, FMT = * )( IACC22( I ), I = 1, NPARMS )
            DO 620 I = 1, NPARMS
               IF( IACC22( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )' IACC22 ', IACC22( I ), 0
                  FATAL = .TRUE.
               END IF
  620       CONTINUE
            WRITE( NOUT, FMT = 9983 )'IACC22: ', ( IACC22( I ), I = 1, NPARMS )
         ELSE
            DO 630 I = 1, NPARMS
               IACC22( I ) = 1
  630       CONTINUE
         END IF

         // Read the values for NBCOL.

         IF( SGG ) THEN
            READ( NIN, FMT = * )( NBCOL( I ), I = 1, NPARMS )
            DO 160 I = 1, NPARMS
               IF( NBCOL( I ).LT.0 ) THEN
                  WRITE( NOUT, FMT = 9989 )'NBCOL ', NBCOL( I ), 0
                  FATAL = .TRUE.
               ELSE IF( NBCOL( I ).GT.NMAX ) THEN
                  WRITE( NOUT, FMT = 9988 )'NBCOL ', NBCOL( I ), NMAX
                  FATAL = .TRUE.
               END IF
  160       CONTINUE
            WRITE( NOUT, FMT = 9983 )'NBCOL:', ( NBCOL( I ), I = 1, NPARMS )
         ELSE
            DO 170 I = 1, NPARMS
               NBCOL( I ) = 1
  170       CONTINUE
         END IF
      END IF

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
      IF( SEP .OR. SVD .OR. SGG ) THEN

         // Read the flag that indicates whether to test LAPACK routines.

         READ( NIN, FMT = * )TSTCHK

         // Read the flag that indicates whether to test driver routines.

         READ( NIN, FMT = * )TSTDRV
      END IF

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR

      // Read the code describing how to set the random number seed.

      READ( NIN, FMT = * )NEWSD

      // If NEWSD = 2, read another line with 4 integers for the seed.

      IF( NEWSD.EQ.2 ) READ( NIN, FMT = * )( IOLDSD( I ), I = 1, 4 )

      DO 180 I = 1, 4
         ISEED( I ) = IOLDSD( I )
  180 CONTINUE

      IF( FATAL ) THEN
         WRITE( NOUT, FMT = 9999 )
         STOP
      END IF

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.

  190 CONTINUE

      IF( .NOT.( SGX .OR. SXV ) ) THEN

  200    CONTINUE
         READ( NIN, FMT = '(A80)', END = 380 )LINE
         C3 = LINE( 1: 3 )
         LENP = LEN( LINE )
         I = 3
         ITMP = 0
         I1 = 0
  210    CONTINUE
         I = I + 1
         IF( I.GT.LENP ) THEN
            IF( I1.GT.0 ) THEN
               GO TO 240
            ELSE
               NTYPES = MAXT
               GO TO 240
            END IF
         END IF
         IF( LINE( I: I ).NE.' ' .AND. LINE( I: I ).NE.',' ) THEN
            I1 = I
            C1 = LINE( I1: I1 )

         // Check that a valid integer was read

            DO 220 K = 1, 10
               IF( C1.EQ.INTSTR( K: K ) ) THEN
                  IC = K - 1
                  GO TO 230
               END IF
  220       CONTINUE
            WRITE( NOUT, FMT = 9991 )I, LINE
            GO TO 200
  230       CONTINUE
            ITMP = 10*ITMP + IC
            GO TO 210
         ELSE IF( I1.GT.0 ) THEN
            GO TO 240
         ELSE
            GO TO 210
         END IF
  240    CONTINUE
         NTYPES = ITMP

      // Skip the tests if NTYPES is <= 0.

         IF( .NOT.( SEV .OR. SES .OR. SVX .OR. SSX .OR. SGV .OR. SGS ) .AND. NTYPES.LE.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
            GO TO 200
         END IF

      ELSE
         IF( SXV ) C3 = 'SXV'          IF( SGX ) C3 = 'SGX'
      END IF

      // Reset the random number seed.

      IF( NEWSD.EQ.0 ) THEN
         DO 250 K = 1, 4
            ISEED( K ) = IOLDSD( K )
  250    CONTINUE
      END IF

      IF( LSAMEN( 3, C3, 'SHS' ) .OR. LSAMEN( 3, C3, 'NEP' ) ) THEN

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
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV( 1, 1 )
         IF( TSTERR ) CALL SERRHS( 'SHSEQR', NOUT )
         DO 270 I = 1, NPARMS
            CALL XLAENV( 1, NBVAL( I ) )
            CALL XLAENV( 2, NBMIN( I ) )
            CALL XLAENV( 3, NXVAL( I ) )
            CALL XLAENV(12, MAX( 11, INMIN( I ) ) )
            CALL XLAENV(13, INWIN( I ) )
            CALL XLAENV(14, INIBL( I ) )
            CALL XLAENV(15, ISHFTS( I ) )
            CALL XLAENV(16, IACC22( I ) )

            IF( NEWSD.EQ.0 ) THEN
               DO 260 K = 1, 4
                  ISEED( K ) = IOLDSD( K )
  260          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9961 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), MAX( 11, INMIN(I)), INWIN( I ), INIBL( I ), ISHFTS( I ), IACC22( I )             CALL SCHKHS( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), A( 1, 7 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 8 ), A( 1, 9 ), A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), D( 1, 7 ), WORK, LWORK, IWORK, LOGWRK, RESULT, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCHKHS', INFO
  270    CONTINUE

      ELSE IF( LSAMEN( 3, C3, 'SST' ) .OR. LSAMEN( 3, C3, 'SEP' ) .OR. LSAMEN( 3, C3, 'SE2' ) ) THEN

         // ----------------------------------
         // SEP:  Symmetric Eigenvalue Problem
         // ----------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV( 1, 1 )
         CALL XLAENV( 9, 25 )
         IF( TSTERR ) THEN
#if defined(_OPENMP)
            N_THREADS = OMP_GET_MAX_THREADS()
            ONE_THREAD = 1
            CALL OMP_SET_NUM_THREADS(ONE_THREAD)
#endif
            CALL SERRST( 'SST', NOUT )
#if defined(_OPENMP)
            CALL OMP_SET_NUM_THREADS(N_THREADS)
#endif
         END IF
         DO 290 I = 1, NPARMS
            CALL XLAENV( 1, NBVAL( I ) )
            CALL XLAENV( 2, NBMIN( I ) )
            CALL XLAENV( 3, NXVAL( I ) )

            IF( NEWSD.EQ.0 ) THEN
               DO 280 K = 1, 4
                  ISEED( K ) = IOLDSD( K )
  280          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I )
            IF( TSTCHK ) THEN
               IF( LSAMEN( 3, C3, 'SE2' ) ) THEN
               CALL SCHKST2STG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
               ELSE
               CALL SCHKST( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), D( 1, 12 ), A( 1, 6 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
               ENDIF
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCHKST', INFO
            END IF
            IF( TSTDRV ) THEN
               IF( LSAMEN( 3, C3, 'SE2' ) ) THEN
               CALL SDRVST2STG( NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11), A( 1, 2 ), NMAX, A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
               ELSE
               CALL SDRVST( NN, NVAL, 18, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11), A( 1, 2 ), NMAX, A( 1, 3 ), D( 1, 12 ), A( 1, 4 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
               ENDIF
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRVST', INFO
            END IF
  290    CONTINUE

      ELSE IF( LSAMEN( 3, C3, 'SSG' ) ) THEN

         // ----------------------------------------------
         // SSG:  Symmetric Generalized Eigenvalue Problem
         // ----------------------------------------------
         // Vary the parameters
            // NB    = block size
            // NBMIN = minimum block size
            // NX    = crossover point

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV( 9, 25 )
         DO 310 I = 1, NPARMS
            CALL XLAENV( 1, NBVAL( I ) )
            CALL XLAENV( 2, NBMIN( I ) )
            CALL XLAENV( 3, NXVAL( I ) )

            IF( NEWSD.EQ.0 ) THEN
               DO 300 K = 1, 4
                  ISEED( K ) = IOLDSD( K )
  300          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9997 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I )
            IF( TSTCHK ) THEN
                // CALL SDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
      // $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
      // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
      // $                      LWORK, IWORK, LIWORK, RESULT, INFO )
               CALL SDRVSG2STG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, D( 1, 3 ), D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK, LWORK, IWORK, LIWORK, RESULT, INFO )
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRVSG', INFO
            END IF
  310    CONTINUE

      ELSE IF( LSAMEN( 3, C3, 'SBD' ) .OR. LSAMEN( 3, C3, 'SVD' ) ) THEN

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
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV( 1, 1 )
         CALL XLAENV( 9, 25 )

         // Test the error exits

         IF( TSTERR .AND. TSTCHK ) CALL SERRBD( 'SBD', NOUT )          IF( TSTERR .AND. TSTDRV ) CALL SERRED( 'SBD', NOUT )

         DO 330 I = 1, NPARMS
            NRHS = NSVAL( I )
            CALL XLAENV( 1, NBVAL( I ) )
            CALL XLAENV( 2, NBMIN( I ) )
            CALL XLAENV( 3, NXVAL( I ) )
            IF( NEWSD.EQ.0 ) THEN
               DO 320 K = 1, 4
                  ISEED( K ) = IOLDSD( K )
  320          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9995 )C3, NBVAL( I ), NBMIN( I ), NXVAL( I ), NRHS
            IF( TSTCHK ) THEN
               CALL SCHKBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 2 ), NMAX, A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), A( 1, 8 ), WORK, LWORK, IWORK, NOUT, INFO )
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCHKBD', INFO
            END IF
            IF( TSTDRV ) CALL SDRVBD( NN, MVAL, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX, A( 1, 3 ), NMAX, A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, IWORK, NOUT, INFO )
  330    CONTINUE

      ELSE IF( LSAMEN( 3, C3, 'SEV' ) ) THEN

         // --------------------------------------------
         // SEV:  Nonsymmetric Eigenvalue Problem Driver
               // SGEEV (eigenvalues and eigenvectors)
         // --------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.LE.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRED( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRVEV( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, RESULT, WORK, LWORK, IWORK, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SGEEV', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( LSAMEN( 3, C3, 'SES' ) ) THEN

         // --------------------------------------------
         // SES:  Nonsymmetric Eigenvalue Problem Driver
               // SGEES (Schur form)
         // --------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.LE.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRED( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRVES( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 4 ), NMAX, RESULT, WORK, LWORK, IWORK, LOGWRK, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SGEES', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( LSAMEN( 3, C3, 'SVX' ) ) THEN

         // --------------------------------------------------------------
         // SVX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // SGEEVX (eigenvalues, eigenvectors and condition numbers)
         // --------------------------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.LT.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRED( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRVVX( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), A( 1, 3 ), NMAX, A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), D( 1, 10 ), D( 1, 11 ), D( 1, 12 ), RESULT, WORK, LWORK, IWORK, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SGEEVX', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( LSAMEN( 3, C3, 'SSX' ) ) THEN

         // ---------------------------------------------------
         // SSX:  Nonsymmetric Eigenvalue Problem Expert Driver
               // SGEESX (Schur form and condition numbers)
         // ---------------------------------------------------

         MAXTYP = 21
         NTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.LT.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRED( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRVSX( NN, NVAL, NTYPES, DOTYPE, ISEED, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 4 ), NMAX, A( 1, 5 ), RESULT, WORK, LWORK, IWORK, LOGWRK, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SGEESX', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( LSAMEN( 3, C3, 'SGG' ) ) THEN

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

         MAXTYP = 26
         NTYPES = MIN( MAXTYP, NTYPES )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         CALL XLAENV(1,1)
         IF( TSTCHK .AND. TSTERR ) CALL SERRGG( C3, NOUT )
         DO 350 I = 1, NPARMS
            CALL XLAENV( 1, NBVAL( I ) )
            CALL XLAENV( 2, NBMIN( I ) )
            CALL XLAENV( 4, NSVAL( I ) )
            CALL XLAENV( 8, MXBVAL( I ) )
            CALL XLAENV( 16, IACC22( I ) )
            CALL XLAENV( 5, NBCOL( I ) )

            IF( NEWSD.EQ.0 ) THEN
               DO 340 K = 1, 4
                  ISEED( K ) = IOLDSD( K )
  340          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9996 )C3, NBVAL( I ), NBMIN( I ), NSVAL( I ), MXBVAL( I ), IACC22( I ), NBCOL( I )
            TSTDIF = .FALSE.
            THRSHN = 10.
            IF( TSTCHK ) THEN
               CALL SCHKGG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, TSTDIF, THRSHN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), A( 1, 8 ), A( 1, 9 ), NMAX, A( 1, 10 ), A( 1, 11 ), A( 1, 12 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), A( 1, 13 ), A( 1, 14 ), WORK, LWORK, LOGWRK, RESULT, INFO )
               IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCHKGG', INFO
            END IF
  350    CONTINUE

      ELSE IF( LSAMEN( 3, C3, 'SGS' ) ) THEN

         // -------------------------------------------------
         // SGS:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGES (Schur form)
         // -------------------------------------------------

         MAXTYP = 26
         NTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.LE.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRGG( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRGES( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, RESULT, LOGWRK, INFO )

            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRGES', INFO

      // Blocked version

            CALL XLAENV(16,1)
            CALL SDRGES3( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), WORK, LWORK, RESULT, LOGWRK, INFO )

            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRGES3', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( SGX ) THEN

         // -------------------------------------------------
         // SGX:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGESX (Schur form and condition numbers)
         // -------------------------------------------------

         MAXTYP = 5
         NTYPES = MAXTYP
         IF( NN.LT.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRGG( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL XLAENV( 5, 2 )
            CALL SDRGSX( NN, NCMAX, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), C( 1, 1 ), NCMAX*NCMAX, A( 1, 12 ), WORK, LWORK, IWORK, LIWORK, LOGWRK, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRGSX', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( LSAMEN( 3, C3, 'SGV' ) ) THEN

         // -------------------------------------------------
         // SGV:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGEV (Eigenvalue/vector form)
         // -------------------------------------------------

         MAXTYP = 26
         NTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.LE.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRGG( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRGEV( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), WORK, LWORK, RESULT, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRGEV', INFO

* Blocked version

            CALL SDRGEV3( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 7 ), NMAX, A( 1, 8 ), A( 1, 9 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), WORK, LWORK, RESULT, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRGEV3', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( SXV ) THEN

         // -------------------------------------------------
         // SXV:  Generalized Nonsymmetric Eigenvalue Problem
               // SGGEVX (eigenvalue/vector with condition numbers)
         // -------------------------------------------------

         MAXTYP = 2
         NTYPES = MAXTYP
         IF( NN.LT.0 ) THEN
            WRITE( NOUT, FMT = 9990 )C3
         ELSE
            IF( TSTERR ) CALL SERRGG( C3, NOUT )
            CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
            CALL SDRGVX( NN, THRESH, NIN, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), A( 1, 5 ), A( 1, 6 ), IWORK( 1 ), IWORK( 2 ), D( 1, 4 ), D( 1, 5 ), D( 1, 6 ), D( 1, 7 ), D( 1, 8 ), D( 1, 9 ), WORK, LWORK, IWORK( 3 ), LIWORK-2, RESULT, LOGWRK, INFO )

            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SDRGVX', INFO
         END IF
         WRITE( NOUT, FMT = 9973 )
         GO TO 10

      ELSE IF( LSAMEN( 3, C3, 'SSB' ) ) THEN

         // ------------------------------
         // SSB:  Symmetric Band Reduction
         // ------------------------------

         MAXTYP = 15
         NTYPES = MIN( MAXTYP, NTYPES )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         IF( TSTERR ) CALL SERRST( 'SSB', NOUT )
          // CALL SCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
      // $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
      // $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
         CALL SCHKSB2STG( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ), D( 1, 3 ), D( 1, 4 ), D( 1, 5 ), A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCHKSB', INFO

      ELSE IF( LSAMEN( 3, C3, 'SBB' ) ) THEN

         // ------------------------------
         // SBB:  General Band Reduction
         // ------------------------------

         MAXTYP = 15
         NTYPES = MIN( MAXTYP, NTYPES )
         CALL ALAREQ( C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT )
         DO 370 I = 1, NPARMS
            NRHS = NSVAL( I )

            IF( NEWSD.EQ.0 ) THEN
               DO 360 K = 1, 4
                  ISEED( K ) = IOLDSD( K )
  360          CONTINUE
            END IF
            WRITE( NOUT, FMT = 9966 )C3, NRHS
            CALL SCHKBB( NN, MVAL, NVAL, NK, KVAL, MAXTYP, DOTYPE, NRHS, ISEED, THRESH, NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), 2*NMAX, D( 1, 1 ), D( 1, 2 ), A( 1, 4 ), NMAX, A( 1, 5 ), NMAX, A( 1, 6 ), NMAX, A( 1, 7 ), WORK, LWORK, RESULT, INFO )
            IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCHKBB', INFO
  370    CONTINUE

      ELSE IF( LSAMEN( 3, C3, 'GLM' ) ) THEN

         // -----------------------------------------
         // GLM:  Generalized Linear Regression Model
         // -----------------------------------------

         CALL XLAENV( 1, 1 )
         IF( TSTERR ) CALL SERRGG( 'GLM', NOUT )          CALL SCKGLM( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, D( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCKGLM', INFO

      ELSE IF( LSAMEN( 3, C3, 'GQR' ) ) THEN

         // ------------------------------------------
         // GQR:  Generalized QR and RQ factorizations
         // ------------------------------------------

         CALL XLAENV( 1, 1 )
         IF( TSTERR ) CALL SERRGG( 'GQR', NOUT )          CALL SCKGQR( NN, MVAL, NN, PVAL, NN, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), TAUA, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), B( 1, 5 ), TAUB, WORK, D( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCKGQR', INFO

      ELSE IF( LSAMEN( 3, C3, 'GSV' ) ) THEN

         // ----------------------------------------------
         // GSV:  Generalized Singular Value Decomposition
         // ----------------------------------------------

         CALL XLAENV( 1, 1 )
         IF( TSTERR ) CALL SERRGG( 'GSV', NOUT )          CALL SCKGSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB, B( 1, 4 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCKGSV', INFO

      ELSE IF( LSAMEN( 3, C3, 'CSD' ) ) THEN

         // ----------------------------------------------
         // CSD:  CS Decomposition
         // ----------------------------------------------

         CALL XLAENV(1,1)
         IF( TSTERR ) CALL SERRGG( 'CSD', NOUT )          CALL SCKCSD( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCKCSD', INFO

      ELSE IF( LSAMEN( 3, C3, 'LSE' ) ) THEN

         // --------------------------------------
         // LSE:  Constrained Linear Least Squares
         // --------------------------------------

         CALL XLAENV( 1, 1 )
         IF( TSTERR ) CALL SERRGG( 'LSE', NOUT )          CALL SCKLSE( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), X, WORK, D( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 ) WRITE( NOUT, FMT = 9980 )'SCKLSE', INFO

      ELSE
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9992 )C3
      END IF
      IF( .NOT.( SGX .OR. SXV ) ) GO TO 190
  380 CONTINUE
      WRITE( NOUT, FMT = 9994 )
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9993 )S2 - S1

      DEALLOCATE (A, STAT = AllocateStatus)
      DEALLOCATE (B, STAT = AllocateStatus)
      DEALLOCATE (C, STAT = AllocateStatus)
      DEALLOCATE (WORK,  STAT = AllocateStatus)

 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9997 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4 )
 9996 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NS =', I4,
     $      ', MAXB =', I4, ', IACC22 =', I4, ', NBCOL =', I4 )
 9995 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4,
     $      ', NRHS =', I4 )
 9994 FORMAT( / / ' End of tests' )
 9993 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9992 FORMAT( 1X, A3, ':  Unrecognized path name' )
 9991 FORMAT( / / ' *** Invalid int     value in column ', I2,;
     $      ' of input', ' line:', / A79 )
 9990 FORMAT( / / 1X, A3, ' routines were not tested' )
 9989 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be >=',
     $      I6 )
 9988 FORMAT( ' Invalid input value: ', A, '=', I6, '; must be <=',
     $      I6 )
 9987 FORMAT( ' Tests of the Nonsymmetric Eigenvalue Problem routines' )
 9986 FORMAT( ' Tests of the Symmetric Eigenvalue Problem routines' )
 9985 FORMAT( ' Tests of the Singular Value Decomposition routines' )
 9984 FORMAT( / ' The following parameter values will be used:' )
 9983 FORMAT( 4X, A, 10I6, / 10X, 10I6 )
 9982 FORMAT( / ' Routines pass computational tests if test ratio is ',
     $      'less than', F8.2, / )
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 )
 9980 FORMAT( ' *** Error code from ', A, ' = ', I4 )
 9979 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver',
     $      / '    SGEEV (eigenvalues and eigevectors)' )
 9978 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Driver',
     $      / '    SGEES (Schur form)' )
 9977 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert',
     $      ' Driver', / '    SGEEVX (eigenvalues, eigenvectors and',
     $      ' condition numbers)' )
 9976 FORMAT( / ' Tests of the Nonsymmetric Eigenvalue Problem Expert',
     $      ' Driver', / '    SGEESX (Schur form and condition',
     $      ' numbers)' )
 9975 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ',
     $      'Problem routines' )
 9974 FORMAT( ' Tests of SSBTRD', / ' (reduction of a symmetric band ',
     $      'matrix to tridiagonal form)' )
 9973 FORMAT( / 1X, 71( '-' ) )
 9972 FORMAT( / ' LAPACK VERSION ', I1, '.', I1, '.', I1 )
 9971 FORMAT( / ' Tests of the Generalized Linear Regression Model ',
     $      'routines' )
 9970 FORMAT( / ' Tests of the Generalized QR and RQ routines' )
 9969 FORMAT( / ' Tests of the Generalized Singular Value',
     $      ' Decomposition routines' )
 9968 FORMAT( / ' Tests of the Linear Least Squares routines' )
 9967 FORMAT( ' Tests of SGBBRD', / ' (reduction of a general band ',
     $      'matrix to real bidiagonal form)' )
 9966 FORMAT( / / 1X, A3, ':  NRHS =', I4 )
 9965 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ',
     $      'Problem Expert Driver SGGESX' )
 9964 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ',
     $      'Problem Driver SGGES' )
 9963 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ',
     $      'Problem Driver SGGEV' )
 9962 FORMAT( / ' Tests of the Generalized Nonsymmetric Eigenvalue ',
     $      'Problem Expert Driver SGGEVX' )
 9961 FORMAT( / / 1X, A3, ':  NB =', I4, ', NBMIN =', I4, ', NX =', I4,
     $      ', INMIN=', I4,
     $      ', INWIN =', I4, ', INIBL =', I4, ', ISHFTS =', I4,
     $      ', IACC22 =', I4)
 9960 FORMAT( / ' Tests of the CS Decomposition routines' )

      // End of SCHKEE

      }
