      SUBROUTINE ZDRVSX( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NIUNIT, NOUNIT, A, LDA, H, HT, W, WT, WTMP, VS, LDVS, VS1, RESULT, WORK, LWORK, RWORK, BWORK, INFO );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDVS, LWORK, NIUNIT, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      double             RESULT( 17 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), H( LDA, * ), HT( LDA, * ), VS( LDVS, * ), VS1( LDVS, * ), W( * ), WORK( * ), WT( * ), WTMP( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      COMPLEX*16         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             PATH;
      int                I, IINFO, IMODE, ISRT, ITYPE, IWK, J, JCOL, JSIZE, JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NSLCT, NTEST, NTESTF, NTESTT;
      double             ANORM, COND, CONDS, OVFL, RCDEIN, RCDVIN, RTULP, RTULPI, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
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
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASUM, XERBLA, ZGET24, ZLASET, ZLATME, ZLATMR, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9,/;
      DATA               KMAGN / 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /;
      DATA               KMODE / 0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /;
      DATA               KCONDS / 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0,/;
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Zomplex precision';
      PATH( 2: 3 ) = 'SX';

      // Check for errors

      NTESTT = 0;
      NTESTF = 0;
      INFO = 0;

      // Important constants

      BADNN = false;

      // 8 is the largest dimension in the input file of precomputed
      // problems

      NMAX = 8;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = MAX( NMAX, NN( J ) );
         IF( NN( J ) < 0 ) BADNN = true;
      } // 10

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADNN ) {
         INFO = -2;
      } else if ( NTYPES < 0 ) {
         INFO = -3;
      } else if ( THRESH < ZERO ) {
         INFO = -6;
      } else if ( NIUNIT <= 0 ) {
         INFO = -7;
      } else if ( NOUNIT <= 0 ) {
         INFO = -8;
      } else if ( LDA < 1 || LDA < NMAX ) {
         INFO = -10;
      } else if ( LDVS < 1 || LDVS < NMAX ) {
         INFO = -20;
      } else if ( MAX( 3*NMAX, 2*NMAX**2 ) > LWORK ) {
         INFO = -24;
      }

      if ( INFO != 0 ) {
         xerbla('ZDRVSX', -INFO );
         RETURN;
      }

      // If nothing to do check on NIUNIT

      if (NSIZES == 0 || NTYPES == 0) GO TO 150;

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = DLAMCH( 'Precision' );
      ULPINV = ONE / ULP;
      RTULP = SQRT( ULP );
      RTULPI = ONE / RTULP;

      // Loop over sizes, types

      NERRS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 140
         N = NN( JSIZE );
         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES );
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 130
            IF( !DOTYPE( JTYPE ) ) GO TO 130;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD( J ) = ISEED( J );
            } // 20

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

            if (MTYPES > MAXTYP) GO TO 90;

            ITYPE = KTYPE( JTYPE );
            IMODE = KMODE( JTYPE );

            // Compute norm

            GO TO ( 30, 40, 50 )KMAGN( JTYPE );

            } // 30
            ANORM = ONE;
            GO TO 60;

            } // 40
            ANORM = OVFL*ULP;
            GO TO 60;

            } // 50
            ANORM = UNFL*ULPINV;
            GO TO 60;

            } // 60

            zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

            if ( ITYPE == 1 ) {

               // Zero

               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 70
                  A( JCOL, JCOL ) = ANORM;
               } // 70

            } else if ( ITYPE == 3 ) {

               // Jordan Block

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A( JCOL, JCOL ) = ANORM;
                  if (JCOL > 1) A( JCOL, JCOL-1 ) = CONE;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // Symmetric, eigenvalues specified

               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 6 ) {

               // General, eigenvalues specified

               if ( KCONDS( JTYPE ) == 1 ) {
                  CONDS = ONE;
               } else if ( KCONDS( JTYPE ) == 2 ) {
                  CONDS = RTULPI;
               } else {
                  CONDS = ZERO;
               }

               zlatme(N, 'D', ISEED, WORK, IMODE, COND, CONE, 'T', 'T', 'T', RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IDUMMA, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IDUMMA, IINFO );

            } else if ( ITYPE == 9 ) {

               // General, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IDUMMA, IINFO );
               if ( N >= 4 ) {
                  zlaset('Full', 2, N, CZERO, CZERO, A, LDA );
                  zlaset('Full', N-3, 1, CZERO, CZERO, A( 3, 1 ), LDA );
                  zlaset('Full', N-3, 2, CZERO, CZERO, A( 3, N-1 ), LDA );
                  zlaset('Full', 1, N, CZERO, CZERO, A( N, 1 ), LDA );
               }

            } else if ( ITYPE == 10 ) {

               // Triangular, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IDUMMA, IINFO );

            } else {

               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9991 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ABS( IINFO );
               RETURN;
            }

            } // 90

            // Test for minimal and generous workspace

            for (IWK = 1; IWK <= 2; IWK++) { // 120
               if ( IWK == 1 ) {
                  NNWORK = 2*N;
               } else {
                  NNWORK = MAX( 2*N, N*( N+1 ) / 2 );
               }
               NNWORK = MAX( NNWORK, 1 );

               zget24( false , JTYPE, THRESH, IOLDSD, NOUNIT, N, A, LDA, H, HT, W, WT, WTMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, 0, RESULT, WORK, NNWORK, RWORK, BWORK, INFO );

               // Check for RESULT(j) > THRESH

               NTEST = 0;
               NFAIL = 0;
               for (J = 1; J <= 15; J++) { // 100
                  IF( RESULT( J ) >= ZERO ) NTEST = NTEST + 1                   IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1;
               } // 100

               if (NFAIL > 0) NTESTF = NTESTF + 1;
               if ( NTESTF == 1 ) {
                  WRITE( NOUNIT, FMT = 9999 )PATH;
                  WRITE( NOUNIT, FMT = 9998 );
                  WRITE( NOUNIT, FMT = 9997 );
                  WRITE( NOUNIT, FMT = 9996 );
                  WRITE( NOUNIT, FMT = 9995 )THRESH;
                  WRITE( NOUNIT, FMT = 9994 );
                  NTESTF = 2;
               }

               for (J = 1; J <= 15; J++) { // 110
                  if ( RESULT( J ) >= THRESH ) {
                     WRITE( NOUNIT, FMT = 9993 )N, IWK, IOLDSD, JTYPE, J, RESULT( J );
                  }
               } // 110

               NERRS = NERRS + NFAIL;
               NTESTT = NTESTT + NTEST;

            } // 120
         } // 130
      } // 140

      } // 150

      // Read in data from file to check accuracy of condition estimation
      // Read input data until N=0

      JTYPE = 0;
      } // 160
      READ( NIUNIT, FMT = *, END = 200 )N, NSLCT, ISRT;
      if (N == 0) GO TO 200;
      JTYPE = JTYPE + 1;
      ISEED( 1 ) = JTYPE;
      READ( NIUNIT, FMT = * )( ISLCT( I ), I = 1, NSLCT );
      for (I = 1; I <= N; I++) { // 170
         READ( NIUNIT, FMT = * )( A( I, J ), J = 1, N );
      } // 170
      READ( NIUNIT, FMT = * )RCDEIN, RCDVIN;

      zget24( true , 22, THRESH, ISEED, NOUNIT, N, A, LDA, H, HT, W, WT, WTMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, ISRT, RESULT, WORK, LWORK, RWORK, BWORK, INFO );

      // Check for RESULT(j) > THRESH

      NTEST = 0;
      NFAIL = 0;
      for (J = 1; J <= 17; J++) { // 180
         IF( RESULT( J ) >= ZERO ) NTEST = NTEST + 1          IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1;
      } // 180

      if (NFAIL > 0) NTESTF = NTESTF + 1;
      if ( NTESTF == 1 ) {
         WRITE( NOUNIT, FMT = 9999 )PATH;
         WRITE( NOUNIT, FMT = 9998 );
         WRITE( NOUNIT, FMT = 9997 );
         WRITE( NOUNIT, FMT = 9996 );
         WRITE( NOUNIT, FMT = 9995 )THRESH;
         WRITE( NOUNIT, FMT = 9994 );
         NTESTF = 2;
      }
      for (J = 1; J <= 17; J++) { // 190
         if ( RESULT( J ) >= THRESH ) {
            WRITE( NOUNIT, FMT = 9992 )N, JTYPE, J, RESULT( J );
         }
      } // 190

      NERRS = NERRS + NFAIL;
      NTESTT = NTESTT + NTEST;
      GO TO 160;
      } // 200

      // Summary

      dlasum(PATH, NOUNIT, NERRS, NTESTT );

 9999 FORMAT( / 1X, A3, ' -- Complex Schur Form Decomposition Expert ', 'Driver', / ' Matrix types (see ZDRVSX for details): ' );

 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ', '           ', '  5=Diagonal: geometr. spaced entries.', / '  2=Identity matrix.                    ', '  6=Diagona', 'l: clustered entries.', / '  3=Transposed Jordan block.  ', '          ', '  7=Diagonal: large, evenly spaced.', / '  ', '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s', 'mall, evenly spaced.' );
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev', 'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e', 'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ', ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond', 'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp', 'lex ', / ' 12=Well-cond., random complex ', '         ', ' 17=Ill-cond., large rand. complx ', / ' 13=Ill-condi', 'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.', ' complx ' );
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ', 'with small random entries.', / ' 20=Matrix with large ran', 'dom entries.   ', / );
 9995 FORMAT( ' Tests performed with test threshold =', F8.2, / ' ( A denotes A on input and T denotes A on output)', / / ' 1 = 0 if T in Schur form (no sort), ', '  1/ulp otherwise', / ' 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)', / ' 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) ', / ' 4 = 0 if W are eigenvalues of T (no sort),', '  1/ulp otherwise', / ' 5 = 0 if T same no matter if VS computed (no sort),', '  1/ulp otherwise', / ' 6 = 0 if W same no matter if VS computed (no sort)', ',  1/ulp otherwise' );
 9994 FORMAT( ' 7 = 0 if T in Schur form (sort), ', '  1/ulp otherwise', / ' 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)', / ' 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) ', / ' 10 = 0 if W are eigenvalues of T (sort),', '  1/ulp otherwise', / ' 11 = 0 if T same no matter what else computed (sort),', '  1/ulp otherwise', / ' 12 = 0 if W same no matter what else computed ', '(sort), 1/ulp otherwise', / ' 13 = 0 if sorting successful, 1/ulp otherwise', / ' 14 = 0 if RCONDE same no matter what else computed,', ' 1/ulp otherwise', / ' 15 = 0 if RCONDv same no matter what else computed,', ' 1/ulp otherwise', / ' 16 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),', / ' 17 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),' );
 9993 FORMAT( ' N=', I5, ', IWK=', I2, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 );
 9992 FORMAT( ' N=', I5, ', input example =', I3, ',  test(', I2, ')=', G10.3 );
 9991 FORMAT( ' ZDRVSX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

      RETURN;

      // End of ZDRVSX

      }
