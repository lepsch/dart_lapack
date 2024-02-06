      void cdrvev(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, H, W, W1, VL, LDVL, VR, LDVR, LRE, LDLRE, RESULT, WORK, NWORK, RWORK, IWORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDLRE, LDVL, LDVR, NOUNIT, NSIZES, NTYPES, NWORK;
      double               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double               RESULT( 7 ), RWORK( * );
      Complex            A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), W1( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               TWO;
      const              TWO = 2.0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             PATH;
      int                IINFO, IMODE, ITYPE, IWK, J, JCOL, JJ, JSIZE, JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NTEST, NTESTF, NTESTT;
      double               ANORM, COND, CONDS, OVFL, RTULP, RTULPI, TNRM, ULP, ULPINV, UNFL, VMX, VRMX, VTST;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double               RES( 2 );
      Complex            DUM( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SCNRM2, SLAMCH;
      // EXTERNAL SCNRM2, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEEV, CGET22, CLACPY, CLATME, CLATMR, CLATMS, CLASET, SLASUM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, MAX, MIN, REAL, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 ];
      const KCONDS = [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0,];
      // ..
      // .. Executable Statements ..

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'EV';

      // Check for errors

      NTESTT = 0;
      NTESTF = 0;
      INFO = 0;

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
      } else if ( LDVL < 1 || LDVL < NMAX ) {
         INFO = -14;
      } else if ( LDVR < 1 || LDVR < NMAX ) {
         INFO = -16;
      } else if ( LDLRE < 1 || LDLRE < NMAX ) {
         INFO = -28;
      } else if ( 5*NMAX+2*NMAX**2 > NWORK ) {
         INFO = -21;
      }

      if ( INFO != 0 ) {
         xerbla('CDRVEV', -INFO );
         return;
      }

      // Quick return if nothing to do

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = SLAMCH( 'Precision' );
      ULPINV = ONE / ULP;
      RTULP = sqrt( ULP );
      RTULPI = ONE / RTULP;

      // Loop over sizes, types

      NERRS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 270
         N = NN( JSIZE );
         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

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

            claset('Full', LDA, N, CZERO, CZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 70
                  A[JCOL][JCOL] = CMPLX( ANORM );
               } // 70

            } else if ( ITYPE == 3 ) {

               // Jordan Block

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL][JCOL] = CMPLX( ANORM );
                  if (JCOL > 1) A( JCOL, JCOL-1 ) = CONE;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // Hermitian, eigenvalues specified

               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 6 ) {

               // General, eigenvalues specified

               if ( KCONDS( JTYPE ) == 1 ) {
                  CONDS = ONE;
               } else if ( KCONDS( JTYPE ) == 2 ) {
                  CONDS = RTULPI;
               } else {
                  CONDS = ZERO;
               }

               clatme(N, 'D', ISEED, WORK, IMODE, COND, CONE, 'T', 'T', 'T', RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               clatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random eigenvalues

               clatmr(N, N, 'D', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // General, random eigenvalues

               clatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );
               if ( N >= 4 ) {
                  claset('Full', 2, N, CZERO, CZERO, A, LDA );
                  claset('Full', N-3, 1, CZERO, CZERO, A( 3, 1 ), LDA );
                  claset('Full', N-3, 2, CZERO, CZERO, A( 3, N-1 ), LDA );
                  claset('Full', 1, N, CZERO, CZERO, A( N, 1 ), LDA );
               }

            } else if ( ITYPE == 10 ) {

               // Triangular, random eigenvalues

               clatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else {

               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9993 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 90

            // Test for minimal and generous workspace

            for (IWK = 1; IWK <= 2; IWK++) { // 250
               if ( IWK == 1 ) {
                  NNWORK = 2*N;
               } else {
                  NNWORK = 5*N + 2*N**2;
               }
               NNWORK = max( NNWORK, 1 );

               // Initialize RESULT

               for (J = 1; J <= 7; J++) { // 100
                  RESULT[J] = -ONE;
               } // 100

               // Compute eigenvalues and eigenvectors, and test them

               clacpy('F', N, N, A, LDA, H, LDA );
               cgeev('V', 'V', N, H, LDA, W, VL, LDVL, VR, LDVR, WORK, NNWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  RESULT[1] = ULPINV;
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV1', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  GO TO 220;
               }

               // Do Test (1)

               cget22('N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES );
               RESULT[1] = RES( 1 );

               // Do Test (2)

               cget22('C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES );
               RESULT[2] = RES( 1 );

               // Do Test (3)

               for (J = 1; J <= N; J++) { // 120
                  TNRM = SCNRM2( N, VR( 1, J ), 1 );
                  RESULT[3] = max( RESULT( 3 ), min( ULPINV, ( TNRM-ONE ).abs() / ULP ) );
                  VMX = ZERO;
                  VRMX = ZERO;
                  for (JJ = 1; JJ <= N; JJ++) { // 110
                     VTST = ( VR( JJ, J ) ).abs();
                     if (VTST > VMX) VMX = VTST;
                     IF( AIMAG( VR( JJ, J ) ) == ZERO && ABS( REAL( VR( JJ, J ) ) ) > VRMX ) VRMX = ABS( double( VR( JJ, J ) ) );
                  } // 110
                  if (VRMX / VMX < ONE-TWO*ULP) RESULT( 3 ) = ULPINV;
               } // 120

               // Do Test (4)

               for (J = 1; J <= N; J++) { // 140
                  TNRM = SCNRM2( N, VL( 1, J ), 1 );
                  RESULT[4] = max( RESULT( 4 ), min( ULPINV, ( TNRM-ONE ).abs() / ULP ) );
                  VMX = ZERO;
                  VRMX = ZERO;
                  for (JJ = 1; JJ <= N; JJ++) { // 130
                     VTST = ( VL( JJ, J ) ).abs();
                     if (VTST > VMX) VMX = VTST;
                     IF( AIMAG( VL( JJ, J ) ) == ZERO && ABS( REAL( VL( JJ, J ) ) ) > VRMX ) VRMX = ABS( double( VL( JJ, J ) ) );
                  } // 130
                  if (VRMX / VMX < ONE-TWO*ULP) RESULT( 4 ) = ULPINV;
               } // 140

               // Compute eigenvalues only, and test them

               clacpy('F', N, N, A, LDA, H, LDA );
               cgeev('N', 'N', N, H, LDA, W1, DUM, 1, DUM, 1, WORK, NNWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  RESULT[1] = ULPINV;
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV2', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  GO TO 220;
               }

               // Do Test (5)

               for (J = 1; J <= N; J++) { // 150
                  if( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV;
               } // 150

               // Compute eigenvalues and right eigenvectors, and test them

               clacpy('F', N, N, A, LDA, H, LDA );
               cgeev('N', 'V', N, H, LDA, W1, DUM, 1, LRE, LDLRE, WORK, NNWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  RESULT[1] = ULPINV;
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV3', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  GO TO 220;
               }

               // Do Test (5) again

               for (J = 1; J <= N; J++) { // 160
                  if( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV;
               } // 160

               // Do Test (6)

               for (J = 1; J <= N; J++) { // 180
                  for (JJ = 1; JJ <= N; JJ++) { // 170
                     if( VR( J, JJ ) != LRE( J, JJ ) ) RESULT( 6 ) = ULPINV;
                  } // 170
               } // 180

               // Compute eigenvalues and left eigenvectors, and test them

               clacpy('F', N, N, A, LDA, H, LDA );
               cgeev('V', 'N', N, H, LDA, W1, LRE, LDLRE, DUM, 1, WORK, NNWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  RESULT[1] = ULPINV;
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV4', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  GO TO 220;
               }

               // Do Test (5) again

               for (J = 1; J <= N; J++) { // 190
                  if( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV;
               } // 190

               // Do Test (7)

               for (J = 1; J <= N; J++) { // 210
                  for (JJ = 1; JJ <= N; JJ++) { // 200
                     if( VL( J, JJ ) != LRE( J, JJ ) ) RESULT( 7 ) = ULPINV;
                  } // 200
               } // 210

               // End of Loop -- Check for RESULT(j) > THRESH

               } // 220

               NTEST = 0;
               NFAIL = 0;
               for (J = 1; J <= 7; J++) { // 230
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
                  NTESTF = 2;
               }

               for (J = 1; J <= 7; J++) { // 240
                  if ( RESULT( J ) >= THRESH ) {
                     WRITE( NOUNIT, FMT = 9994 )N, IWK, IOLDSD, JTYPE, J, RESULT( J );
                  }
               } // 240

               NERRS = NERRS + NFAIL;
               NTESTT = NTESTT + NTEST;

            } // 250
         } // 260
      } // 270

      // Summary

      slasum(PATH, NOUNIT, NERRS, NTESTT );

 9999 FORMAT( / 1X, A3, ' -- Complex Eigenvalue-Eigenvector ', 'Decomposition Driver', / ' Matrix types (see CDRVEV for details): ' );

 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ', '           ', '  5=Diagonal: geometr. spaced entries.', / '  2=Identity matrix.                    ', '  6=Diagona', 'l: clustered entries.', / '  3=Transposed Jordan block.  ', '          ', '  7=Diagonal: large, evenly spaced.', / '  ', '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s', 'mall, evenly spaced.' );
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev', 'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e', 'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ', ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond', 'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp', 'lex ', A6, / ' 12=Well-cond., random complex ', A6, '   ', ' 17=Ill-cond., large rand. complx ', A4, / ' 13=Ill-condi', 'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.', ' complx ', A4 );
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ', 'with small random entries.', / ' 20=Matrix with large ran', 'dom entries.   ', / );
 9995 FORMAT( ' Tests performed with test threshold =', F8.2, / / ' 1 = | A VR - VR W | / ( n |A| ulp ) ', / ' 2 = | conj-trans(A) VL - VL conj-trans(W) | /', ' ( n |A| ulp ) ', / ' 3 = | |VR(i)| - 1 | / ulp ', / ' 4 = | |VL(i)| - 1 | / ulp ', / ' 5 = 0 if W same no matter if VR or VL computed,', ' 1/ulp otherwise', / ' 6 = 0 if VR same no matter if VL computed,', '  1/ulp otherwise', / ' 7 = 0 if VL same no matter if VR computed,', '  1/ulp otherwise', / );
 9994 FORMAT( ' N=', I5, ', IWK=', I2, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 );
 9993 FORMAT( ' CDRVEV: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

      return;
      }
