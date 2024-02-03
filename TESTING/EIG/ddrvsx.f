      SUBROUTINE DDRVSX( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NIUNIT, NOUNIT, A, LDA, H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS, LDVS, VS1, RESULT, WORK, LWORK, IWORK, BWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDVS, LWORK, NIUNIT, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             A( LDA, * ), H( LDA, * ), HT( LDA, * ), RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ), WI( * ), WIT( * ), WITMP( * ), WORK( * ), WR( * ), WRT( * ), WRTMP( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 21 )
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             PATH;
      int                I, IINFO, IMODE, ITYPE, IWK, J, JCOL, JSIZE, JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NSLCT, NTEST, NTESTF, NTESTT;
      double             ANORM, COND, CONDS, OVFL, RCDEIN, RCDVIN, RTULP, RTULPI, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      String             ADUMMA( 1 );
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISLCT( 20 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGET24, DLASET, DLASUM, DLATME, DLATMR, DLATMS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
      DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
      DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'SX'

      // Check for errors

      NTESTT = 0
      NTESTF = 0
      INFO = 0

      // Important constants

      BADNN = .FALSE.

      // 12 is the largest dimension in the input file of precomputed
      // problems

      NMAX = 12
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      // Check for errors

      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( NIUNIT.LE.0 ) THEN
         INFO = -7
      ELSE IF( NOUNIT.LE.0 ) THEN
         INFO = -8
      ELSE IF( LDA.LT.1 .OR. LDA.LT.NMAX ) THEN
         INFO = -10
      ELSE IF( LDVS.LT.1 .OR. LDVS.LT.NMAX ) THEN
         INFO = -20
      ELSE IF( MAX( 3*NMAX, 2*NMAX**2 ).GT.LWORK ) THEN
         INFO = -24
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DDRVSX', -INFO )
         RETURN
      END IF

      // If nothing to do check on NIUNIT

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) GO TO 150

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      RTULP = SQRT( ULP )
      RTULPI = ONE / RTULP

      // Loop over sizes, types

      NERRS = 0

      DO 140 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF

         DO 130 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 130

            // Save ISEED in case of an error.

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Compute "A"

            // Control parameters:

            // KMAGN  KCONDS  KMODE        KTYPE
        // =1  O(1)   1       clustered 1  zero
        // =2  large  large   clustered 2  identity
        // =3  small          exponential  Jordan
        // =4                 arithmetic   diagonal, (w/ eigenvalues)
        // =5                 random log   symmetric, w/ eigenvalues
        // =6                 random       general, w/ eigenvalues
        // =7                              random diagonal
        // =8                              random symmetric
        // =9                              random general
        // =10                             random triangular

            IF( MTYPES.GT.MAXTYP ) GO TO 90

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 30, 40, 50 )KMAGN( JTYPE )

   30       CONTINUE
            ANORM = ONE
            GO TO 60

   40       CONTINUE
            ANORM = OVFL*ULP
            GO TO 60

   50       CONTINUE
            ANORM = UNFL*ULPINV
            GO TO 60

   60       CONTINUE

            CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            COND = ULPINV

            // Special Matrices -- Identity & Jordan block

               // Zero

            IF( ITYPE.EQ.1 ) THEN
               IINFO = 0

            ELSE IF( ITYPE.EQ.2 ) THEN

               // Identity

               DO 70 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   70          CONTINUE

            ELSE IF( ITYPE.EQ.3 ) THEN

               // Jordan Block

               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
                  IF( JCOL.GT.1 ) A( JCOL, JCOL-1 ) = ONE
   80          CONTINUE

            ELSE IF( ITYPE.EQ.4 ) THEN

               // Diagonal Matrix, [Eigen]values Specified

               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO )

            ELSE IF( ITYPE.EQ.5 ) THEN

               // Symmetric, eigenvalues specified

               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO )

            ELSE IF( ITYPE.EQ.6 ) THEN

               // General, eigenvalues specified

               IF( KCONDS( JTYPE ).EQ.1 ) THEN
                  CONDS = ONE
               ELSE IF( KCONDS( JTYPE ).EQ.2 ) THEN
                  CONDS = RTULPI
               ELSE
                  CONDS = ZERO
               END IF

               ADUMMA( 1 ) = ' '
               CALL DLATME( N, 'S', ISEED, WORK, IMODE, COND, ONE, ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO )

            ELSE IF( ITYPE.EQ.7 ) THEN

               // Diagonal, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            ELSE IF( ITYPE.EQ.8 ) THEN

               // Symmetric, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            ELSE IF( ITYPE.EQ.9 ) THEN

               // General, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
               IF( N.GE.4 ) THEN
                  CALL DLASET( 'Full', 2, N, ZERO, ZERO, A, LDA )
                  CALL DLASET( 'Full', N-3, 1, ZERO, ZERO, A( 3, 1 ), LDA )                   CALL DLASET( 'Full', N-3, 2, ZERO, ZERO, A( 3, N-1 ), LDA )                   CALL DLASET( 'Full', 1, N, ZERO, ZERO, A( N, 1 ), LDA )
               END IF

            ELSE IF( ITYPE.EQ.10 ) THEN

               // Triangular, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            ELSE

               IINFO = 1
            END IF

            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9991 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF

   90       CONTINUE

            // Test for minimal and generous workspace

            DO 120 IWK = 1, 2
               IF( IWK.EQ.1 ) THEN
                  NNWORK = 3*N
               ELSE
                  NNWORK = MAX( 3*N, 2*N*N )
               END IF
               NNWORK = MAX( NNWORK, 1 )

               CALL DGET24( .FALSE., JTYPE, THRESH, IOLDSD, NOUNIT, N, A, LDA, H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, RESULT, WORK, NNWORK, IWORK, BWORK, INFO )

               // Check for RESULT(j) > THRESH

               NTEST = 0
               NFAIL = 0
               DO 100 J = 1, 15
                  IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1                   IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  100          CONTINUE

               IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
               IF( NTESTF.EQ.1 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )PATH
                  WRITE( NOUNIT, FMT = 9998 )
                  WRITE( NOUNIT, FMT = 9997 )
                  WRITE( NOUNIT, FMT = 9996 )
                  WRITE( NOUNIT, FMT = 9995 )THRESH
                  WRITE( NOUNIT, FMT = 9994 )
                  NTESTF = 2
               END IF

               DO 110 J = 1, 15
                  IF( RESULT( J ).GE.THRESH ) THEN
                     WRITE( NOUNIT, FMT = 9993 )N, IWK, IOLDSD, JTYPE, J, RESULT( J )
                  END IF
  110          CONTINUE

               NERRS = NERRS + NFAIL
               NTESTT = NTESTT + NTEST

  120       CONTINUE
  130    CONTINUE
  140 CONTINUE

  150 CONTINUE

      // Read in data from file to check accuracy of condition estimation
      // Read input data until N=0

      JTYPE = 0
  160 CONTINUE
      READ( NIUNIT, FMT = *, END = 200 )N, NSLCT
      IF( N.EQ.0 ) GO TO 200
      JTYPE = JTYPE + 1
      ISEED( 1 ) = JTYPE
      IF( NSLCT.GT.0 ) READ( NIUNIT, FMT = * )( ISLCT( I ), I = 1, NSLCT )
      DO 170 I = 1, N
         READ( NIUNIT, FMT = * )( A( I, J ), J = 1, N )
  170 CONTINUE
      READ( NIUNIT, FMT = * )RCDEIN, RCDVIN

      CALL DGET24( .TRUE., 22, THRESH, ISEED, NOUNIT, N, A, LDA, H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, RESULT, WORK, LWORK, IWORK, BWORK, INFO )

      // Check for RESULT(j) > THRESH

      NTEST = 0
      NFAIL = 0
      DO 180 J = 1, 17
         IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1          IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  180 CONTINUE

      IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
      IF( NTESTF.EQ.1 ) THEN
         WRITE( NOUNIT, FMT = 9999 )PATH
         WRITE( NOUNIT, FMT = 9998 )
         WRITE( NOUNIT, FMT = 9997 )
         WRITE( NOUNIT, FMT = 9996 )
         WRITE( NOUNIT, FMT = 9995 )THRESH
         WRITE( NOUNIT, FMT = 9994 )
         NTESTF = 2
      END IF
      DO 190 J = 1, 17
         IF( RESULT( J ).GE.THRESH ) THEN
            WRITE( NOUNIT, FMT = 9992 )N, JTYPE, J, RESULT( J )
         END IF
  190 CONTINUE

      NERRS = NERRS + NFAIL
      NTESTT = NTESTT + NTEST
      GO TO 160
  200 CONTINUE

      // Summary

      CALL DLASUM( PATH, NOUNIT, NERRS, NTESTT )

 9999 FORMAT( / 1X, A3, ' -- Real Schur Form Decomposition Expert ',
     $      'Driver', / ' Matrix types (see DDRVSX for details):' )

 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ',
     $      '           ', '  5=Diagonal: geometr. spaced entries.',
     $      / '  2=Identity matrix.                    ', '  6=Diagona',
     $      'l: clustered entries.', / '  3=Transposed Jordan block.  ',
     $      '          ', '  7=Diagonal: large, evenly spaced.', / '  ',
     $      '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s',
     $      'mall, evenly spaced.' )
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev',
     $      'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e',
     $      'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ',
     $      ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond',
     $      'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp',
     $      'lex ', / ' 12=Well-cond., random complex ', '         ',
     $      ' 17=Ill-cond., large rand. complx ', / ' 13=Ill-condi',
     $      'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.',
     $      ' complx ' )
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ',
     $      'with small random entries.', / ' 20=Matrix with large ran',
     $      'dom entries.   ', / )
 9995 FORMAT( ' Tests performed with test threshold =', F8.2,
     $      / ' ( A denotes A on input and T denotes A on output)',
     $      / / ' 1 = 0 if T in Schur form (no sort), ',
     $      '  1/ulp otherwise', /
     $      ' 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)',
     $      / ' 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) ', /
     $      ' 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),',
     $      '  1/ulp otherwise', /
     $      ' 5 = 0 if T same no matter if VS computed (no sort),',
     $      '  1/ulp otherwise', /
     $      ' 6 = 0 if WR, WI same no matter if VS computed (no sort)',
     $      ',  1/ulp otherwise' )
 9994 FORMAT( ' 7 = 0 if T in Schur form (sort), ', '  1/ulp otherwise',
     $      / ' 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)',
     $      / ' 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) ',
     $      / ' 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),',
     $      '  1/ulp otherwise', /
     $      ' 11 = 0 if T same no matter what else computed (sort),',
     $      '  1/ulp otherwise', /
     $      ' 12 = 0 if WR, WI same no matter what else computed ',
     $      '(sort), 1/ulp otherwise', /
     $      ' 13 = 0 if sorting successful, 1/ulp otherwise',
     $      / ' 14 = 0 if RCONDE same no matter what else computed,',
     $      ' 1/ulp otherwise', /
     $      ' 15 = 0 if RCONDv same no matter what else computed,',
     $      ' 1/ulp otherwise', /
     $      ' 16 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),',
     $      / ' 17 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),' )
 9993 FORMAT( ' N=', I5, ', IWK=', I2, ', seed=', 4( I4, ',' ),
     $      ' type ', I2, ', test(', I2, ')=', G10.3 )
 9992 FORMAT( ' N=', I5, ', input example =', I3, ',  test(', I2, ')=',
     $      G10.3 )
 9991 FORMAT( ' DDRVSX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of DDRVSX

      END
