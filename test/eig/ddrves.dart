import 'common.dart';

      void ddrves(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, H, HT, WR, WI, WRT, WIT, VS, LDVS, RESULT, WORK, NWORK, IWORK, BWORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDVS, NOUNIT, NSIZES, NTYPES, NWORK;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             A( LDA, * ), H( LDA, * ), HT( LDA, * ), RESULT( 13 ), VS( LDVS, * ), WI( * ), WIT( * ), WORK( * ), WR( * ), WRT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             SORT;
      String             PATH;
      int                I, IINFO, IMODE, ISORT, ITYPE, IWK, J, JCOL, JSIZE, JTYPE, KNTEIG, LWORK, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NTEST, NTESTF, NTESTT, RSUB, SDIM;
      double             ANORM, COND, CONDS, OVFL, RTULP, RTULPI, TMP, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      String             ADUMMA( 1 );
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double             RES( 2 );
      // ..
      // .. Arrays in Common ..
      // bool               sslct.SELVAL( 20 );
      // double             sslct.SELWI( 20 ), sslct.SELWR( 20 );
      // // ..
      // // .. Scalars in Common ..
      // int                sslct.SELDIM, sslct.SELOPT;
      // ..
      // .. Common blocks ..
      // COMMON / sslct / sslct.SELOPT, sslct.SELDIM, sslct.SELVAL, sslct.SELWR, sslct.SELWI
      // ..
      // .. External Functions ..
      //- bool               DSLECT;
      //- double             DLAMCH;
      // EXTERNAL DSLECT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEES, DHST01, DLACPY, DLASET, DLASUM, DLATME, DLATMR, DLATMS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 ];
      const KCONDS = [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0,];
      // ..
      // .. Executable Statements ..

      PATH[1: 1] = 'double          ';
      PATH[2: 3] = 'ES';

      // Check for errors

      NTESTT = 0;
      NTESTF = 0;
      INFO = 0;
      sslct.SELOPT = 0;

      // Important constants

      BADNN = false;
      NMAX = 0;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
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
      } else if ( NOUNIT <= 0 ) {
         INFO = -7;
      } else if ( LDA < 1 || LDA < NMAX ) {
         INFO = -9;
      } else if ( LDVS < 1 || LDVS < NMAX ) {
         INFO = -17;
      } else if ( 5*NMAX+2*NMAX**2 > NWORK ) {
         INFO = -20;
      }

      if ( INFO != 0 ) {
         xerbla('DDRVES', -INFO );
         return;
      }

      // Quick return if nothing to do

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = dlamch( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = dlamch( 'Precision' );
      ULPINV = ONE / ULP;
      RTULP = sqrt( ULP );
      RTULPI = ONE / RTULP;

      // Loop over sizes, types

      NERRS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 270
         N = NN( JSIZE );
         MTYPES = MAXTYP;
         if (NSIZES == 1 && NTYPES == MAXTYP+1) MTYPES = MTYPES + 1;

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 260
            if( !DOTYPE( JTYPE ) ) GO TO 260;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
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

            dlaset('Full', LDA, N, ZERO, ZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 70
                  A[JCOL, JCOL] = ANORM;
               } // 70

            } else if ( ITYPE == 3 ) {

               // Jordan Block

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL, JCOL] = ANORM;
                  if (JCOL > 1) A( JCOL, JCOL-1 ) = ONE;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // Symmetric, eigenvalues specified

               dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 6 ) {

               // General, eigenvalues specified

               if ( KCONDS( JTYPE ) == 1 ) {
                  CONDS = ONE;
               } else if ( KCONDS( JTYPE ) == 2 ) {
                  CONDS = RTULPI;
               } else {
                  CONDS = ZERO;
               }

               ADUMMA[1] = ' ';
               dlatme(N, 'S', ISEED, WORK, IMODE, COND, ONE, ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // General, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );
               if ( N >= 4 ) {
                  dlaset('Full', 2, N, ZERO, ZERO, A, LDA );
                  dlaset('Full', N-3, 1, ZERO, ZERO, A( 3, 1 ), LDA );
                  dlaset('Full', N-3, 2, ZERO, ZERO, A( 3, N-1 ), LDA );
                  dlaset('Full', 1, N, ZERO, ZERO, A( N, 1 ), LDA );
               }

            } else if ( ITYPE == 10 ) {

               // Triangular, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else {

               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9992 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 90

            // Test for minimal and generous workspace

            for (IWK = 1; IWK <= 2; IWK++) { // 250
               if ( IWK == 1 ) {
                  NNWORK = 3*N;
               } else {
                  NNWORK = 5*N + 2*N**2;
               }
               NNWORK = max( NNWORK, 1 );

               // Initialize RESULT

               for (J = 1; J <= 13; J++) { // 100
                  RESULT[J] = -ONE;
               } // 100

               // Test with and without sorting of eigenvalues

               for (ISORT = 0; ISORT <= 1; ISORT++) { // 210
                  if ( ISORT == 0 ) {
                     SORT = 'N';
                     RSUB = 0;
                  } else {
                     SORT = 'S';
                     RSUB = 6;
                  }

                  // Compute Schur form and Schur vectors, and test them

                  dlacpy('F', N, N, A, LDA, H, LDA );
                  dgees('V', SORT, DSLECT, N, H, LDA, SDIM, WR, WI, VS, LDVS, WORK, NNWORK, BWORK, IINFO );
                  if ( IINFO != 0 && IINFO != N+2 ) {
                     RESULT[1+RSUB] = ULPINV;
                     WRITE( NOUNIT, FMT = 9992 )'DGEES1', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     GO TO 220;
                  }

                  // Do Test (1) or Test (7)

                  RESULT[1+RSUB] = ZERO;
                  for (J = 1; J <= N - 2; J++) { // 120
                     for (I = J + 2; I <= N; I++) { // 110
                        if( H( I, J ) != ZERO ) RESULT( 1+RSUB ) = ULPINV;
                     } // 110
                  } // 120
                  for (I = 1; I <= N - 2; I++) { // 130
                     if( H( I+1, I ) != ZERO && H( I+2, I+1 ) != ZERO )RESULT( 1+RSUB ) = ULPINV;
                  } // 130
                  for (I = 1; I <= N - 1; I++) { // 140
                     if ( H( I+1, I ) != ZERO ) {
                        if( H( I, I ) != H( I+1, I+1 ) || H( I, I+1 ) == ZERO || sign( ONE, H( I+1, I ) ) == sign( ONE, H( I, I+1 ) ) )RESULT( 1+RSUB ) = ULPINV;
                     }
                  } // 140

                  // Do Tests (2) and (3) or Tests (8) and (9)

                  LWORK = max( 1, 2*N*N );
                  dhst01(N, 1, N, A, LDA, H, LDA, VS, LDVS, WORK, LWORK, RES );
                  RESULT[2+RSUB] = RES( 1 );
                  RESULT[3+RSUB] = RES( 2 );

                  // Do Test (4) or Test (10)

                  RESULT[4+RSUB] = ZERO;
                  for (I = 1; I <= N; I++) { // 150
                     if( H( I, I ) != WR( I ) ) RESULT( 4+RSUB ) = ULPINV;
                  } // 150
                  if ( N > 1 ) {
                     if( H( 2, 1 ) == ZERO && WI( 1 ) != ZERO ) RESULT( 4+RSUB ) = ULPINV;
                     IF( H( N, N-1 ) == ZERO && WI( N ) != ZERO ) RESULT( 4+RSUB ) = ULPINV;
                  }
                  for (I = 1; I <= N - 1; I++) { // 160
                     if ( H( I+1, I ) != ZERO ) {
                        TMP = sqrt( ( H( I+1, I ) ) ).abs()* sqrt( ( H( I, I+1 ) ) ).abs()                         RESULT( 4+RSUB ) = max( RESULT( 4+RSUB ), ABS( WI( I )-TMP ) / max( ULP*TMP, UNFL ) )                         RESULT( 4+RSUB ) = max( RESULT( 4+RSUB ), ABS( WI( I+1 )+TMP ) / max( ULP*TMP, UNFL ) );
                     } else if ( I > 1 ) {
                        if( H( I+1, I ) == ZERO && H( I, I-1 ) == ZERO && WI( I ) != ZERO )RESULT( 4+RSUB ) = ULPINV;
                     }
                  } // 160

                  // Do Test (5) or Test (11)

                  dlacpy('F', N, N, A, LDA, HT, LDA );
                  dgees('N', SORT, DSLECT, N, HT, LDA, SDIM, WRT, WIT, VS, LDVS, WORK, NNWORK, BWORK, IINFO );
                  if ( IINFO != 0 && IINFO != N+2 ) {
                     RESULT[5+RSUB] = ULPINV;
                     WRITE( NOUNIT, FMT = 9992 )'DGEES2', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     GO TO 220;
                  }

                  RESULT[5+RSUB] = ZERO;
                  for (J = 1; J <= N; J++) { // 180
                     for (I = 1; I <= N; I++) { // 170
                        if( H( I, J ) != HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV;
                     } // 170
                  } // 180

                  // Do Test (6) or Test (12)

                  RESULT[6+RSUB] = ZERO;
                  for (I = 1; I <= N; I++) { // 190
                     if( WR( I ) != WRT( I ) || WI( I ) != WIT( I ) ) RESULT( 6+RSUB ) = ULPINV;
                  } // 190

                  // Do Test (13)

                  if ( ISORT == 1 ) {
                     RESULT[13] = ZERO;
                     KNTEIG = 0;
                     for (I = 1; I <= N; I++) { // 200
                        if( DSLECT( WR( I ), WI( I ) ) || DSLECT( WR( I ), -WI( I ) ) ) KNTEIG = KNTEIG + 1;
                        if ( I < N ) {
                           if( ( DSLECT( WR( I+1 ), WI( I+1 ) ) || DSLECT( WR( I+1 ), -WI( I+1 ) ) ) && ( !( DSLECT( WR( I ), WI( I ) ) || DSLECT( WR( I ), -WI( I ) ) ) ) && IINFO != N+2 ) RESULT( 13 ) = ULPINV;
                        }
                     } // 200
                     if ( SDIM != KNTEIG ) {
                        RESULT[13] = ULPINV;
                     }
                  }

               } // 210

               // End of Loop -- Check for RESULT(j) > THRESH

               } // 220

               NTEST = 0;
               NFAIL = 0;
               for (J = 1; J <= 13; J++) { // 230
                  if( RESULT( J ) >= ZERO ) NTEST = NTEST + 1;
                  IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1;
               } // 230

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

               for (J = 1; J <= 13; J++) { // 240
                  if ( RESULT( J ) >= THRESH ) {
                     WRITE( NOUNIT, FMT = 9993 )N, IWK, IOLDSD, JTYPE, J, RESULT( J );
                  }
               } // 240

               NERRS = NERRS + NFAIL;
               NTESTT = NTESTT + NTEST;

            } // 250
         } // 260
      } // 270

      // Summary

      dlasum(PATH, NOUNIT, NERRS, NTESTT );

 9999 FORMAT( / 1X, A3, ' -- Real Schur Form Decomposition Driver', / ' Matrix types (see DDRVES for details): ' );

 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ', '           ', '  5=Diagonal: geometr. spaced entries.', / '  2=Identity matrix.                    ', '  6=Diagona', 'l: clustered entries.', / '  3=Transposed Jordan block.  ', '          ', '  7=Diagonal: large, evenly spaced.', / '  ', '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s', 'mall, evenly spaced.' );
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev', 'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e', 'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ', ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond', 'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp', 'lex ', / ' 12=Well-cond., random complex ', 6X, '   ', ' 17=Ill-cond., large rand. complx ', / ' 13=Ill-condi', 'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.', ' complx ' );
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ', 'with small random entries.', / ' 20=Matrix with large ran', 'dom entries.   ', / );
 9995 FORMAT( ' Tests performed with test threshold =', F8.2, / ' ( A denotes A on input and T denotes A on output)', / / ' 1 = 0 if T in Schur form (no sort), ', '  1/ulp otherwise', / ' 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)', / ' 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) ', / ' 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),', '  1/ulp otherwise', / ' 5 = 0 if T same no matter if VS computed (no sort),', '  1/ulp otherwise', / ' 6 = 0 if WR, WI same no matter if VS computed (no sort)', ',  1/ulp otherwise' );
 9994 FORMAT( ' 7 = 0 if T in Schur form (sort), ', '  1/ulp otherwise', / ' 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)', / ' 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) ', / ' 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),', '  1/ulp otherwise', / ' 11 = 0 if T same no matter if VS computed (sort),', '  1/ulp otherwise', / ' 12 = 0 if WR, WI same no matter if VS computed (sort),', '  1/ulp otherwise', / ' 13 = 0 if sorting successful, 1/ulp otherwise', / );
 9993 FORMAT( ' N=', I5, ', IWK=', I2, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 );
 9992 FORMAT( ' DDRVES: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

      return;
      }
