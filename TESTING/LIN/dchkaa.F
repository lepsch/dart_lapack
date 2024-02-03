      PROGRAM DCHKAA

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      PARAMETER          ( NMAX = 132 )
      int                MAXIN;
      PARAMETER          ( MAXIN = 12 )
      int                MAXRHS;
      PARAMETER          ( MAXRHS = 16 )
      int                MATMAX;
      PARAMETER          ( MATMAX = 30 )
      int                NIN, NOUT;
      PARAMETER          ( NIN = 5, NOUT = 6 )
      int                KDMAX;
      PARAMETER          ( KDMAX = NMAX+( NMAX+1 ) / 4 )
      // ..
      // .. Local Scalars ..
      bool               FATAL, TSTCHK, TSTDRV, TSTERR;
      String             C1;
      String             C2;
      String             PATH;
      String             INTSTR;
      String             ALINE;
      int                I, IC, J, K, LA, LAFAC, LDA, NB, NM, NMATS, NN, NNB, NNB2, NNS, NRHS, NTYPES, NRANK, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      double             EPS, S1, S2, THREQ, THRESH;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MATMAX );
      int                IWORK( 25*NMAX ), MVAL( MAXIN ), NBVAL( MAXIN ), NBVAL2( MAXIN ), NSVAL( MAXIN ), NVAL( MAXIN ), NXVAL( MAXIN ), RANKVAL( MAXIN ), PIV( NMAX );
      // ..
      // .. Allocatable Arrays ..
      int     AllocateStatus;
      double          , DIMENSION(:), ALLOCATABLE :: RWORK, S;
      double          , DIMENSION(:), ALLOCATABLE :: E;
      double          , DIMENSION(:,:), ALLOCATABLE :: A, B, WORK;
      // ..
      // .. External Functions ..
      bool               LSAME, LSAMEN;
      double             DLAMCH, DSECND;
      // EXTERNAL LSAME, LSAMEN, DLAMCH, DSECND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, DCHKEQ, DCHKGB, DCHKGE, DCHKGT, DCHKLQ, DCHKORHR_COL, DCHKPB, DCHKPO, DCHKPS, DCHKPP, DCHKPT, DCHKQ3, DCHKQP3RK, DCHKQL, DCHKQR, DCHKRQ, DCHKSP, DCHKSY, DCHKSY_ROOK, DCHKSY_RK, DCHKSY_AA, DCHKTB, DCHKTP, DCHKTR, DCHKTZ, DDRVGB, DDRVGE, DDRVGT, DDRVLS, DDRVPB, DDRVPO, DDRVPP, DDRVPT, DDRVSP, DDRVSY, DDRVSY_ROOK, DDRVSY_RK, DDRVSY_AA, ILAVER, DCHKLQTP, DCHKQRT, DCHKQRTP, DCHKLQT,DCHKTSQR
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
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      COMMON             / CLAENV / IPARMS
      // ..
      // .. Data statements ..
      DATA               THREQ / 2.0D0 / , INTSTR / '0123456789' /
      // ..
      // ..
      // .. Allocate memory dynamically ..

      ALLOCATE ( A( ( KDMAX+1 )*NMAX, 7 ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( B( NMAX*MAXRHS, 4 ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( WORK( NMAX, 3*NMAX+MAXRHS+30 ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( E( NMAX ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( S( 2*NMAX ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ALLOCATE ( RWORK( 5*NMAX+2*MAXRHS ), STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      // .. Executable Statements ..

      S1 = DSECND( )
      LDA = NMAX
      FATAL = .FALSE.

      // Read a dummy line.

      READ( NIN, FMT = * )

      // Report values of parameters.

      CALL ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH

      // Read the values of M

      READ( NIN, FMT = * )NM
      IF( NM.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NM ', NM, 1
         NM = 0
         FATAL = .TRUE.
      ELSE IF( NM.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN
         NM = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      DO 10 I = 1, NM
         IF( MVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( MVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   10 CONTINUE
      IF( NM.GT.0 ) WRITE( NOUT, FMT = 9993 )'M   ', ( MVAL( I ), I = 1, NM )

      // Read the values of N

      READ( NIN, FMT = * )NN
      IF( NN.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      ELSE IF( NN.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      DO 20 I = 1, NN
         IF( NVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' N  ', NVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( NVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9995 )' N  ', NVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   20 CONTINUE
      IF( NN.GT.0 ) WRITE( NOUT, FMT = 9993 )'N   ', ( NVAL( I ), I = 1, NN )

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS
      IF( NNS.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1
         NNS = 0
         FATAL = .TRUE.
      ELSE IF( NNS.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN
         NNS = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS )
      DO 30 I = 1, NNS
         IF( NSVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( NSVAL( I ).GT.MAXRHS ) THEN
            WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS
            FATAL = .TRUE.
         END IF
   30 CONTINUE
      IF( NNS.GT.0 ) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS )

      // Read the values of NB

      READ( NIN, FMT = * )NNB
      IF( NNB.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )'NNB ', NNB, 1
         NNB = 0
         FATAL = .TRUE.
      ELSE IF( NNB.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )'NNB ', NNB, MAXIN
         NNB = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB )
      DO 40 I = 1, NNB
         IF( NBVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NB ', NBVAL( I ), 0
            FATAL = .TRUE.
         END IF
   40 CONTINUE
      IF( NNB.GT.0 ) WRITE( NOUT, FMT = 9993 )'NB  ', ( NBVAL( I ), I = 1, NNB )

      // Set NBVAL2 to be the set of unique values of NB

      NNB2 = 0
      DO 60 I = 1, NNB
         NB = NBVAL( I )
         DO 50 J = 1, NNB2
            IF( NB.EQ.NBVAL2( J ) ) GO TO 60
   50    CONTINUE
         NNB2 = NNB2 + 1
         NBVAL2( NNB2 ) = NB
   60 CONTINUE

      // Read the values of NX

      READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB )
      DO 70 I = 1, NNB
         IF( NXVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NX ', NXVAL( I ), 0
            FATAL = .TRUE.
         END IF
   70 CONTINUE
      IF( NNB.GT.0 ) WRITE( NOUT, FMT = 9993 )'NX  ', ( NXVAL( I ), I = 1, NNB )

      // Read the values of RANKVAL

      READ( NIN, FMT = * )NRANK
      IF( NN.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NRANK ', NRANK, 1
         NRANK = 0
         FATAL = .TRUE.
      ELSE IF( NN.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NRANK ', NRANK, MAXIN
         NRANK = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( RANKVAL( I ), I = 1, NRANK )
      DO I = 1, NRANK
         IF( RANKVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' RANK  ', RANKVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( RANKVAL( I ).GT.100 ) THEN
            WRITE( NOUT, FMT = 9995 )' RANK  ', RANKVAL( I ), 100
            FATAL = .TRUE.
         END IF
      END DO
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

      IF( FATAL ) THEN
         WRITE( NOUT, FMT = 9999 )
         STOP
      END IF

      // Calculate and print the machine dependent constants.

      EPS = DLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9991 )'underflow', EPS
      EPS = DLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9991 )'overflow ', EPS
      EPS = DLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9991 )'precision', EPS
      WRITE( NOUT, FMT = * )

   80 CONTINUE

      // Read a test path and the number of matrix types to use.

      READ( NIN, FMT = '(A72)', END = 140 )ALINE
      PATH = ALINE( 1: 3 )
      NMATS = MATMAX
      I = 3
   90 CONTINUE
      I = I + 1
      IF( I.GT.72 ) THEN
         NMATS = MATMAX
         GO TO 130
      END IF
      IF( ALINE( I: I ).EQ.' ' ) GO TO 90
      NMATS = 0
  100 CONTINUE
      C1 = ALINE( I: I )
      DO 110 K = 1, 10
         IF( C1.EQ.INTSTR( K: K ) ) THEN
            IC = K - 1
            GO TO 120
         END IF
  110 CONTINUE
      GO TO 130
  120 CONTINUE
      NMATS = NMATS*10 + IC
      I = I + 1
      IF( I.GT.72 ) GO TO 130
      GO TO 100
  130 CONTINUE
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )
      NRHS = NSVAL( 1 )

      // Check first character for correct precision.

      IF( .NOT.LSAME( C1, 'double          ' ) ) THEN;
         WRITE( NOUT, FMT = 9990 )PATH

      ELSE IF( NMATS.LE.0 ) THEN

         // Check for a positive number of tests requested.

         WRITE( NOUT, FMT = 9989 )PATH

      ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN

         // GE:  general matrices

         NTYPES = 11
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN

         // GB:  general banded matrices

         LA = ( 2*KDMAX+1 )*NMAX
         LAFAC = ( 3*KDMAX+1 )*NMAX
         NTYPES = 8
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKGB( DOTYPE, NM, MVAL, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), LA, A( 1, 3 ), LAFAC, B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), LA, A( 1, 3 ), LAFAC, A( 1, 6 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'GT' ) ) THEN

         // GT:  general tridiagonal matrices

         NTYPES = 12
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'PO' ) ) THEN

         // PO:  positive definite matrices

         NTYPES = 9
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKPO( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'PS' ) ) THEN

         // PS:  positive semi-definite matrices

         NTYPES = 9

         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKPS( DOTYPE, NN, NVAL, NNB2, NBVAL2, NRANK, RANKVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), PIV, WORK, RWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'PP' ) ) THEN

         // PP:  positive definite packed matrices

         NTYPES = 9
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN

         // PB:  positive definite banded matrices

         NTYPES = 8
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKPB( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVPB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), S, WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'PT' ) ) THEN

         // PT:  positive definite tridiagonal matrices

         NTYPES = 12
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN

         // SY:  symmetric indefinite matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 10
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKSY( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'SR' ) ) THEN

         // SR:  symmetric indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm

         NTYPES = 10
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKSY_ROOK(DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVSY_ROOK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'SK' ) ) THEN

         // SK:  symmetric indefinite matrices,
              // with bounded Bunch-Kaufman (rook) pivoting algorithm,
              // different matrix storage format than SR path version.

         NTYPES = 10
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKSY_RK( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVSY_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'SA' ) ) THEN

         // SA:  symmetric indefinite matrices,
              // with partial (Aasen's) pivoting algorithm

         NTYPES = 10
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKSY_AA( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVSY_AA( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF


      ELSE IF( LSAMEN( 2, C2, 'S2' ) ) THEN

         // SA:  symmetric indefinite matrices,
              // with partial (Aasen's) pivoting algorithm

         NTYPES = 10
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVSY_AA_2STAGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF


      ELSE IF( LSAMEN( 2, C2, 'SP' ) ) THEN

         // SP:  symmetric indefinite packed matrices,
              // with partial (Bunch-Kaufman) pivoting algorithm

         NTYPES = 10
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKSP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

         IF( TSTDRV ) THEN
            CALL DDRVSP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN

         // TR:  triangular matrices

         NTYPES = 18
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKTR( DOTYPE, NN, NVAL, NNB2, NBVAL2, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'TP' ) ) THEN

         // TP:  triangular packed matrices

         NTYPES = 18
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKTP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN

         // TB:  triangular banded matrices

         NTYPES = 17
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'QR' ) ) THEN

         // QR:  QR factorization

         NTYPES = 8
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKQR( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'LQ' ) ) THEN

         // LQ:  LQ factorization

         NTYPES = 8
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'QL' ) ) THEN

         // QL:  QL factorization

         NTYPES = 8
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKQL( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'RQ' ) ) THEN

         // RQ:  RQ factorization

         NTYPES = 8
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKRQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), A( 1, 4 ), A( 1, 5 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, RWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'QP' ) ) THEN

         // QP:  QR factorization with pivoting

         NTYPES = 6
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKQ3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, THRESH, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 3 ), WORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'QK' ) ) THEN

         // QK: truncated QR factorization with pivoting

         NTYPES = 19
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKQP3RK( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), B( 1, 4 ), WORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'TZ' ) ) THEN

         // TZ:  Trapezoidal matrix

         NTYPES = 3
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTCHK ) THEN
            CALL DCHKTZ( DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 3 ), WORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'LS' ) ) THEN

         // LS:  Least squares drivers

         NTYPES = 6
         CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

         IF( TSTDRV ) THEN
            CALL DDRVLS( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), B( 1, 3 ), RWORK, RWORK( NMAX+1 ), NOUT )
         ELSE
            WRITE( NOUT, FMT = 9988 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'EQ' ) ) THEN

         // EQ:  Equilibration routines for general and positive definite
              // matrices (THREQ should be between 2 and 10)

         IF( TSTCHK ) THEN
            CALL DCHKEQ( THREQ, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'QT' ) ) THEN

         // QT:  QRT routines for general matrices

         IF( TSTCHK ) THEN
            CALL DCHKQRT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'QX' ) ) THEN

         // QX:  QRT routines for triangular-pentagonal matrices

         IF( TSTCHK ) THEN
            CALL DCHKQRTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'TQ' ) ) THEN

         // TQ:  LQT routines for general matrices

         IF( TSTCHK ) THEN
            CALL DCHKLQT( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'XQ' ) ) THEN

         // XQ:  LQT routines for triangular-pentagonal matrices

         IF( TSTCHK ) THEN
            CALL DCHKLQTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'TS' ) ) THEN

         // TS:  QR routines for tall-skinny matrices

         IF( TSTCHK ) THEN
            CALL DCHKTSQR( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF

      ELSE IF( LSAMEN( 2, C2, 'HH' ) ) THEN

         // HH:  Householder reconstruction for tall-skinny matrices

         IF( TSTCHK ) THEN
            CALL DCHKORHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 ) PATH
         END IF

      ELSE


         WRITE( NOUT, FMT = 9990 )PATH
      END IF

      // Go back to get another input line.

      GO TO 80

      // Branch to this line when the last record is read.

  140 CONTINUE
      CLOSE ( NIN )
      S2 = DSECND( )
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
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=',
     $      I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=',
     $      I6 )
 9994 FORMAT( ' Tests of the double           LAPACK routines ',;
     $      / ' LAPACK VERSION ', I1, '.', I1, '.', I1,
     $      / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ',
     $      'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
 9990 FORMAT( / 1X, A3, ':  Unrecognized path name' )
 9989 FORMAT( / 1X, A3, ' routines were not tested' )
 9988 FORMAT( / 1X, A3, ' driver routines were not tested' )

      // End of DCHKAA

      }
