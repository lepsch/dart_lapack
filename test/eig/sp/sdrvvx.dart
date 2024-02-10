      void sdrvvx(NSIZES, NN, NTYPES, final Array<bool> DOTYPE, final Array<int> ISEED, THRESH, NIUNIT, NOUNIT, final Matrix<double> A, final int LDA, H, WR, WI, WR1, WI1, final Matrix<double> VL, final int LDVL, final Matrix<double> VR, final int LDVR, final Matrix<double> LRE, final int LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, NWORK, IWORK, Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDLRE, LDVL, LDVR, NIUNIT, NOUNIT, NSIZES, NTYPES, NWORK;
      double               THRESH;
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double               A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), SCALE( * ), SCALE1( * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      bool               BADNN;
      String             BALANC;
      String             PATH;
      int                I, IBAL, IINFO, IMODE, ITYPE, IWK, J, JCOL, JSIZE, JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NTEST, NTESTF, NTESTT;
      double               ANORM, COND, CONDS, OVFL, RTULP, RTULPI, ULP, ULPINV, UNFL;
      String             ADUMMA( 1 ), BAL( 4 );
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGET23, SLASUM, SLATME, SLATMR, SLATMS, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 ];
      const KCONDS = [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0,];
      const BAL = [ 'N', 'P', 'S', 'B' ];

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'VX';

      // Check for errors

      NTESTT = 0;
      NTESTF = 0;
      INFO = 0;

      // Important constants

      BADNN = false;

      // 12 is the largest dimension in the input file of precomputed
      // problems

      NMAX = 12;
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
      } else if ( LDA < 1 || LDA < NMAX ) {
         INFO = -10;
      } else if ( LDVL < 1 || LDVL < NMAX ) {
         INFO = -17;
      } else if ( LDVR < 1 || LDVR < NMAX ) {
         INFO = -19;
      } else if ( LDLRE < 1 || LDLRE < NMAX ) {
         INFO = -21;
      } else if ( 6*NMAX+2*NMAX**2 > NWORK ) {
         INFO = -32;
      }

      if ( INFO != 0 ) {
         xerbla('SDRVVX', -INFO );
         return;
      }

      // If nothing to do check on NIUNIT

      if (NSIZES == 0 || NTYPES == 0) GO TO 160;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = SLAMCH( 'Precision' );
      ULPINV = ONE / ULP;
      RTULP = sqrt( ULP );
      RTULPI = ONE / RTULP;

      // Loop over sizes, types

      NERRS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 150
         N = NN( JSIZE );
         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 140
            if( !DOTYPE( JTYPE ) ) GO TO 140;

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

            slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 70
                  A[JCOL][JCOL] = ANORM;
               } // 70

            } else if ( ITYPE == 3 ) {

               // Jordan Block

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL][JCOL] = ANORM;
                  if (JCOL > 1) A( JCOL, JCOL-1 ) = ONE;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // Symmetric, eigenvalues specified

               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

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
               slatme(N, 'S', ISEED, WORK, IMODE, COND, ONE, ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random eigenvalues

               slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // General, random eigenvalues

               slatmr(N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );
               if ( N >= 4 ) {
                  slaset('Full', 2, N, ZERO, ZERO, A, LDA );
                  slaset('Full', N-3, 1, ZERO, ZERO, A( 3, 1 ), LDA );
                  slaset('Full', N-3, 2, ZERO, ZERO, A( 3, N-1 ), LDA );
                  slaset('Full', 1, N, ZERO, ZERO, A( N, 1 ), LDA );
               }

            } else if ( ITYPE == 10 ) {

               // Triangular, random eigenvalues

               slatmr(N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

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

            for (IWK = 1; IWK <= 3; IWK++) { // 130
               if ( IWK == 1 ) {
                  NNWORK = 3*N;
               } else if ( IWK == 2 ) {
                  NNWORK = 6*N + N**2;
               } else {
                  NNWORK = 6*N + 2*N**2;
               }
               NNWORK = max( NNWORK, 1 );

               // Test for all balancing options

               for (IBAL = 1; IBAL <= 4; IBAL++) { // 120
                  BALANC = BAL( IBAL );

                  // Perform tests

                  sget23( false , BALANC, JTYPE, THRESH, IOLDSD, NOUNIT, N, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, NNWORK, IWORK, INFO );

                  // Check for RESULT(j) > THRESH

                  NTEST = 0;
                  NFAIL = 0;
                  for (J = 1; J <= 9; J++) { // 100
                     if( RESULT( J ) >= ZERO ) NTEST = NTEST + 1;
                     IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1;
                  } // 100

                  if (NFAIL > 0) NTESTF = NTESTF + 1;
                  if ( NTESTF == 1 ) {
                     WRITE( NOUNIT, FMT = 9999 )PATH;
                     WRITE( NOUNIT, FMT = 9998 );
                     WRITE( NOUNIT, FMT = 9997 );
                     WRITE( NOUNIT, FMT = 9996 );
                     WRITE( NOUNIT, FMT = 9995 )THRESH;
                     NTESTF = 2;
                  }

                  for (J = 1; J <= 9; J++) { // 110
                     if ( RESULT( J ) >= THRESH ) {
                        WRITE( NOUNIT, FMT = 9994 )BALANC, N, IWK, IOLDSD, JTYPE, J, RESULT( J );
                     }
                  } // 110

                  NERRS = NERRS + NFAIL;
                  NTESTT = NTESTT + NTEST;

               } // 120
            } // 130
         } // 140
      } // 150

      } // 160

      // Read in data from file to check accuracy of condition estimation.
      // Assume input eigenvalues are sorted lexicographically (increasing
      // by real part, then decreasing by imaginary part)

      JTYPE = 0;
      } // 170
      READ( NIUNIT, FMT = *, END = 220 )N;

      // Read input data until N=0

      if (N == 0) GO TO 220;
      JTYPE = JTYPE + 1;
      ISEED[1] = JTYPE;
      for (I = 1; I <= N; I++) { // 180
         READ( NIUNIT, FMT = * )( A( I, J ), J = 1, N );
      } // 180
      for (I = 1; I <= N; I++) { // 190
         READ( NIUNIT, FMT = * )WR1( I ), WI1( I ), RCDEIN( I ), RCDVIN( I );
      } // 190
      sget23( true , 'N', 22, THRESH, ISEED, NOUNIT, N, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, 6*N+2*N**2, IWORK, INFO );

      // Check for RESULT(j) > THRESH

      NTEST = 0;
      NFAIL = 0;
      for (J = 1; J <= 11; J++) { // 200
         if( RESULT( J ) >= ZERO ) NTEST = NTEST + 1;
         IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1;
      } // 200

      if (NFAIL > 0) NTESTF = NTESTF + 1;
      if ( NTESTF == 1 ) {
         WRITE( NOUNIT, FMT = 9999 )PATH;
         WRITE( NOUNIT, FMT = 9998 );
         WRITE( NOUNIT, FMT = 9997 );
         WRITE( NOUNIT, FMT = 9996 );
         WRITE( NOUNIT, FMT = 9995 )THRESH;
         NTESTF = 2;
      }

      for (J = 1; J <= 11; J++) { // 210
         if ( RESULT( J ) >= THRESH ) {
            WRITE( NOUNIT, FMT = 9993 )N, JTYPE, J, RESULT( J );
         }
      } // 210

      NERRS = NERRS + NFAIL;
      NTESTT = NTESTT + NTEST;
      GO TO 170;
      } // 220

      // Summary

      slasum(PATH, NOUNIT, NERRS, NTESTT );

 9999 FORMAT('\n ${.a3} -- Real Eigenvalue-Eigenvector Decomposition Expert Driver\n Matrix types (see SDRVVX for details): ' );

 9998 FORMAT('\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.' );
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ' );
 9996 FORMAT( ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.    22=Matrix read from input file\n');
 9995 FORMAT( ' Tests performed with test threshold =${.f8_2}\n\n 1 = | A VR - VR W | / ( n |A| ulp ) \n 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) \n 3 = | |VR(i)| - 1 | / ulp \n 4 = | |VL(i)| - 1 | / ulp \n 5 = 0 if W same no matter if VR or VL computed, 1/ulp otherwise\n 6 = 0 if VR same no matter what else computed,  1/ulp otherwise\n 7 = 0 if VL same no matter what else computed,  1/ulp otherwise\n 8 = 0 if RCONDV same no matter what else computed,  1/ulp otherwise\n 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else computed,  1/ulp otherwise\n 10 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),\n 11 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),' );
 9994 FORMAT( ' BALANC=''${.a1}'',N=${.i4},IWK=${.i1}, seed=${i4(4, ',')}', ' type ${.i2}, test(${.i2})=${.g10_3}');
 9993 FORMAT( ' N=${.i5}, input example =${.i3},  test(${.i2})=${.g10_3}');
 9992 FORMAT( ' SDRVVX: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );

      }
