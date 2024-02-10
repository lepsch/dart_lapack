      void zchkhs(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, H, T1, T2, U, LDU, Z, UZ, W1, W3, EVECTL, EVECTR, EVECTY, EVECTX, UU, TAU, WORK, NWORK, RWORK, IWORK, SELECT, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK;
      double             THRESH;
      bool               DOTYPE( * ), SELECT( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             RESULT( 16 ), RWORK( * );
      Complex         A( LDA, * ), EVECTL( LDU, * ), EVECTR( LDU, * ), EVECTX( LDU, * ), EVECTY( LDU, * ), H( LDA, * ), T1( LDA, * ), T2( LDA, * ), TAU( * ), U( LDU, * ), UU( LDU, * ), UZ( LDU, * ), W1( * ), W3( * ), WORK( * ), Z( LDU, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      bool               BADNN, MATCH;
      int                I, IHI, IINFO, ILO, IMODE, IN, ITYPE, J, JCOL, JJ, JSIZE, JTYPE, K, MTYPES, N, N1, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ANINV, ANORM, COND, CONDS, OVFL, RTOVFL, RTULP, RTULPI, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL;
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double             DUMMA( 4 );
      Complex         CDUMMA( 4 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAFTS, DLASUM, XERBLA, ZCOPY, ZGEHRD, ZGEMM, ZGET10, ZGET22, ZHSEIN, ZHSEQR, ZHST01, ZLACPY, ZLASET, ZLATME, ZLATMR, ZLATMS, ZTREVC, ZTREVC3, ZUNGHR, ZUNMHR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 ];
      const KCONDS = [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0,];

      // Check for errors

      NTESTT = 0;
      INFO = 0;

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
      } else if ( LDA <= 1 || LDA < NMAX ) {
         INFO = -9;
      } else if ( LDU <= 1 || LDU < NMAX ) {
         INFO = -14;
      } else if ( 4*NMAX*NMAX+2 > NWORK ) {
         INFO = -26;
      }

      if ( INFO != 0 ) {
         xerbla('ZCHKHS', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      // More important constants

      UNFL = dlamch( 'Safe minimum' );
      OVFL = dlamch( 'Overflow' );
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );
      RTULP = sqrt( ULP );
      RTULPI = ONE / RTULP;

      // Loop over sizes, types

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 260
         N = NN( JSIZE );
         if (N == 0) GO TO 260;
         N1 = max( 1, N );
         ANINV = ONE / N1.toDouble();

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 250
            if( !DOTYPE( JTYPE ) ) GO TO 250;
            NMATS = NMATS + 1;
            NTEST = 0;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            // Initialize RESULT

            for (J = 1; J <= 14; J++) { // 30
               RESULT[J] = ZERO;
            } // 30

            // Compute "A"

            // Control parameters:

            // KMAGN  KCONDS  KMODE        KTYPE
        // =1  O(1)   1       clustered 1  zero
        // =2  large  large   clustered 2  identity
        // =3  small          exponential  Jordan
        // =4                 arithmetic   diagonal, (w/ eigenvalues)
        // =5                 random log   hermitian, w/ eigenvalues
        // =6                 random       general, w/ eigenvalues
        // =7                              random diagonal
        // =8                              random hermitian
        // =9                              random general
        // =10                             random triangular

            if (MTYPES > MAXTYP) GO TO 100;

            ITYPE = KTYPE( JTYPE );
            IMODE = KMODE( JTYPE );

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE );

            } // 40
            ANORM = ONE;
            GO TO 70;

            } // 50
            ANORM = ( RTOVFL*ULP )*ANINV;
            GO TO 70;

            } // 60
            ANORM = RTUNFL*N*ULPINV;
            GO TO 70;

            } // 70

            zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices

            if ( ITYPE == 1 ) {

               // Zero

               IINFO = 0;
            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL][JCOL] = ANORM;
               } // 80

            } else if ( ITYPE == 3 ) {

               // Jordan Block

               for (JCOL = 1; JCOL <= N; JCOL++) { // 90
                  A[JCOL][JCOL] = ANORM;
                  if (JCOL > 1) A( JCOL, JCOL-1 ) = ONE;
               } // 90

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               zlatmr(N, N, 'D', ISEED, 'N', WORK, IMODE, COND, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 5 ) {

               // Hermitian, eigenvalues specified

               zlatms(N, N, 'D', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 6 ) {

               // General, eigenvalues specified

               if ( KCONDS( JTYPE ) == 1 ) {
                  CONDS = ONE;
               } else if ( KCONDS( JTYPE ) == 2 ) {
                  CONDS = RTULPI;
               } else {
                  CONDS = ZERO;
               }

               zlatme(N, 'D', ISEED, WORK, IMODE, COND, CONE, 'T', 'T', 'T', RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Hermitian, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // General, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 10 ) {

               // Triangular, random eigenvalues

               zlatmr(N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else {

               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 100

            // Call ZGEHRD to compute H and U, do tests.

            zlacpy(' ', N, N, A, LDA, H, LDA );
            NTEST = 1;

            ILO = 1;
            IHI = N;

            zgehrd(N, ILO, IHI, H, LDA, WORK, WORK( N+1 ), NWORK-N, IINFO );

            if ( IINFO != 0 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'ZGEHRD', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            for (J = 1; J <= N - 1; J++) { // 120
               UU[J+1][J] = CZERO;
               for (I = J + 2; I <= N; I++) { // 110
                  U[I][J] = H( I, J );
                  UU[I][J] = H( I, J );
                  H[I][J] = CZERO;
               } // 110
            } // 120
            zcopy(N-1, WORK, 1, TAU, 1 );
            zunghr(N, ILO, IHI, U, LDU, WORK, WORK( N+1 ), NWORK-N, IINFO );
            NTEST = 2;

            zhst01(N, ILO, IHI, A, LDA, H, LDA, U, LDU, WORK, NWORK, RWORK, RESULT( 1 ) );

            // Call ZHSEQR to compute T1, T2 and Z, do tests.

            // Eigenvalues only (W3)

            zlacpy(' ', N, N, H, LDA, T2, LDA );
            NTEST = 3;
            RESULT[3] = ULPINV;

            zhseqr('E', 'N', N, ILO, IHI, T2, LDA, W3, UZ, LDU, WORK, NWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHSEQR(E)', IINFO, N, JTYPE, IOLDSD;
               if ( IINFO <= N+2 ) {
                  INFO = ( IINFO ).abs();
                  GO TO 240;
               }
            }

            // Eigenvalues (W1) and Full Schur Form (T2)

            zlacpy(' ', N, N, H, LDA, T2, LDA );

            zhseqr('S', 'N', N, ILO, IHI, T2, LDA, W1, UZ, LDU, WORK, NWORK, IINFO );
            if ( IINFO != 0 && IINFO <= N+2 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHSEQR(S)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            // Eigenvalues (W1), Schur Form (T1), and Schur Vectors (UZ)

            zlacpy(' ', N, N, H, LDA, T1, LDA );
            zlacpy(' ', N, N, U, LDU, UZ, LDU );

            zhseqr('S', 'V', N, ILO, IHI, T1, LDA, W1, UZ, LDU, WORK, NWORK, IINFO );
            if ( IINFO != 0 && IINFO <= N+2 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHSEQR(V)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            // Compute Z = U' UZ

            zgemm('C', 'N', N, N, N, CONE, U, LDU, UZ, LDU, CZERO, Z, LDU );
            NTEST = 8;

            // Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
            //      and 4: | I - Z Z' | / ( n ulp )

            zhst01(N, ILO, IHI, H, LDA, T1, LDA, Z, LDU, WORK, NWORK, RWORK, RESULT( 3 ) );

            // Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
            //      and 6: | I - UZ (UZ)' | / ( n ulp )

            zhst01(N, ILO, IHI, A, LDA, T1, LDA, UZ, LDU, WORK, NWORK, RWORK, RESULT( 5 ) );

            // Do Test 7: | T2 - T1 | / ( |T| n ulp )

            zget10(N, N, T2, LDA, T1, LDA, WORK, RWORK, RESULT( 7 ) );

            // Do Test 8: | W3 - W1 | / ( max(|W1|,|W3|) ulp )

            TEMP1 = ZERO;
            TEMP2 = ZERO;
            for (J = 1; J <= N; J++) { // 130
               TEMP1 = max( TEMP1, ( W1( J ) ).abs(), ( W3( J ) ).abs() );
               TEMP2 = max( TEMP2, ABS( W1( J )-W3( J ) ) );
            } // 130

            RESULT[8] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

            // Compute the Left and Right Eigenvectors of T

            // Compute the Right eigenvector Matrix:

            NTEST = 9;
            RESULT[9] = ULPINV;

            // Select every other eigenvector

            for (J = 1; J <= N; J++) { // 140
               SELECT[J] = false;
            } // 140
            for (J = 1; J <= N; J += 2) { // 150
               SELECT[J] = true;
            } // 150
            ztrevc('Right', 'All', SELECT, N, T1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(R,A)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            // Test 9:  | TR - RW | / ( |T| |R| ulp )

            zget22('N', 'N', 'N', N, T1, LDA, EVECTR, LDU, W1, WORK, RWORK, DUMMA( 1 ) );
            RESULT[9] = DUMMA( 1 );
            if ( DUMMA( 2 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZTREVC', DUMMA( 2 ), N, JTYPE, IOLDSD;
            }

            // Compute selected right eigenvectors and confirm that
            // they agree with previous right eigenvectors

            ztrevc('Right', 'Some', SELECT, N, T1, LDA, CDUMMA, LDU, EVECTL, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(R,S)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            K = 1;
            MATCH = true;
            for (J = 1; J <= N; J++) { // 170
               if ( SELECT( J ) ) {
                  for (JJ = 1; JJ <= N; JJ++) { // 160
                     if ( EVECTR( JJ, J ) != EVECTL( JJ, K ) ) {
                        MATCH = false;
                        GO TO 180;
                     }
                  } // 160
                  K = K + 1;
               }
            } // 170
            } // 180
            if ( !MATCH) WRITE( NOUNIT, FMT = 9997 )'Right', 'ZTREVC', N, JTYPE, IOLDSD;

            // Compute the Left eigenvector Matrix:

            NTEST = 10;
            RESULT[10] = ULPINV;
            ztrevc('Left', 'All', SELECT, N, T1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(L,A)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            // Test 10:  | LT - WL | / ( |T| |L| ulp )

            zget22('C', 'N', 'C', N, T1, LDA, EVECTL, LDU, W1, WORK, RWORK, DUMMA( 3 ) );
            RESULT[10] = DUMMA( 3 );
            if ( DUMMA( 4 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZTREVC', DUMMA( 4 ), N, JTYPE, IOLDSD;
            }

            // Compute selected left eigenvectors and confirm that
            // they agree with previous left eigenvectors

            ztrevc('Left', 'Some', SELECT, N, T1, LDA, EVECTR, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(L,S)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 240;
            }

            K = 1;
            MATCH = true;
            for (J = 1; J <= N; J++) { // 200
               if ( SELECT( J ) ) {
                  for (JJ = 1; JJ <= N; JJ++) { // 190
                     if ( EVECTL( JJ, J ) != EVECTR( JJ, K ) ) {
                        MATCH = false;
                        GO TO 210;
                     }
                  } // 190
                  K = K + 1;
               }
            } // 200
            } // 210
            if ( !MATCH) WRITE( NOUNIT, FMT = 9997 )'Left', 'ZTREVC', N, JTYPE, IOLDSD;

            // Call ZHSEIN for Right eigenvectors of H, do test 11

            NTEST = 11;
            RESULT[11] = ULPINV;
            for (J = 1; J <= N; J++) { // 220
               SELECT[J] = true;
            } // 220

            zhsein('Right', 'Qr', 'Ninitv', SELECT, N, H, LDA, W3, CDUMMA, LDU, EVECTX, LDU, N1, IN, WORK, RWORK, IWORK, IWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHSEIN(R)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if (IINFO < 0) GO TO 240;
            } else {

               // Test 11:  | HX - XW | / ( |H| |X| ulp )

                         // (from inverse iteration)

               CALL ZGET22( 'N', 'N', 'N', N, H, LDA, EVECTX, LDU, W3, WORK, RWORK, DUMMA( 1 ) )                IF( DUMMA( 1 ) < ULPINV ) RESULT( 11 ) = DUMMA( 1 )*ANINV;
               if ( DUMMA( 2 ) > THRESH ) {
                  WRITE( NOUNIT, FMT = 9998 )'Right', 'ZHSEIN', DUMMA( 2 ), N, JTYPE, IOLDSD;
               }
            }

            // Call ZHSEIN for Left eigenvectors of H, do test 12

            NTEST = 12;
            RESULT[12] = ULPINV;
            for (J = 1; J <= N; J++) { // 230
               SELECT[J] = true;
            } // 230

            zhsein('Left', 'Qr', 'Ninitv', SELECT, N, H, LDA, W3, EVECTY, LDU, CDUMMA, LDU, N1, IN, WORK, RWORK, IWORK, IWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHSEIN(L)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if (IINFO < 0) GO TO 240;
            } else {

               // Test 12:  | YH - WY | / ( |H| |Y| ulp )

                         // (from inverse iteration)

               CALL ZGET22( 'C', 'N', 'C', N, H, LDA, EVECTY, LDU, W3, WORK, RWORK, DUMMA( 3 ) )                IF( DUMMA( 3 ) < ULPINV ) RESULT( 12 ) = DUMMA( 3 )*ANINV;
               if ( DUMMA( 4 ) > THRESH ) {
                  WRITE( NOUNIT, FMT = 9998 )'Left', 'ZHSEIN', DUMMA( 4 ), N, JTYPE, IOLDSD;
               }
            }

            // Call ZUNMHR for Right eigenvectors of A, do test 13

            NTEST = 13;
            RESULT[13] = ULPINV;

            zunmhr('Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTX, LDU, WORK, NWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZUNMHR(L)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if (IINFO < 0) GO TO 240;
            } else {

               // Test 13:  | AX - XW | / ( |A| |X| ulp )

                         // (from inverse iteration)

               CALL ZGET22( 'N', 'N', 'N', N, A, LDA, EVECTX, LDU, W3, WORK, RWORK, DUMMA( 1 ) )                IF( DUMMA( 1 ) < ULPINV ) RESULT( 13 ) = DUMMA( 1 )*ANINV;
            }

            // Call ZUNMHR for Left eigenvectors of A, do test 14

            NTEST = 14;
            RESULT[14] = ULPINV;

            zunmhr('Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTY, LDU, WORK, NWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZUNMHR(L)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if (IINFO < 0) GO TO 240;
            } else {

               // Test 14:  | YA - WY | / ( |A| |Y| ulp )

                         // (from inverse iteration)

               CALL ZGET22( 'C', 'N', 'C', N, A, LDA, EVECTY, LDU, W3, WORK, RWORK, DUMMA( 3 ) )                IF( DUMMA( 3 ) < ULPINV ) RESULT( 14 ) = DUMMA( 3 )*ANINV;
            }

            // Compute Left and Right Eigenvectors of A

            // Compute a Right eigenvector matrix:

            NTEST = 15;
            RESULT[15] = ULPINV;

            zlacpy(' ', N, N, UZ, LDU, EVECTR, LDU );

            ztrevc3('Right', 'Back', SELECT, N, T1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, NWORK, RWORK, N, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC3(R,B)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 250;
            }

            // Test 15:  | AR - RW | / ( |A| |R| ulp )

                      // (from Schur decomposition)

            zget22('N', 'N', 'N', N, A, LDA, EVECTR, LDU, W1, WORK, RWORK, DUMMA( 1 ) );
            RESULT[15] = DUMMA( 1 );
            if ( DUMMA( 2 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZTREVC3', DUMMA( 2 ), N, JTYPE, IOLDSD;
            }

            // Compute a Left eigenvector matrix:

            NTEST = 16;
            RESULT[16] = ULPINV;

            zlacpy(' ', N, N, UZ, LDU, EVECTL, LDU );

            ztrevc3('Left', 'Back', SELECT, N, T1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, NWORK, RWORK, N, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC3(L,B)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 250;
            }

            // Test 16:  | LA - WL | / ( |A| |L| ulp )

                      // (from Schur decomposition)

            zget22('Conj', 'N', 'Conj', N, A, LDA, EVECTL, LDU, W1, WORK, RWORK, DUMMA( 3 ) );
            RESULT[16] = DUMMA( 3 );
            if ( DUMMA( 4 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZTREVC3', DUMMA( 4 ), N, JTYPE, IOLDSD;
            }

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 240

            NTESTT = NTESTT + NTEST;
            dlafts('ZHS', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );

         } // 250
      } // 260

      // Summary

      dlasum('ZHS', NOUNIT, NERRS, NTESTT );

      return;

 9999 FORMAT( ' ZCHKHS: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );
 9998 FORMAT( ' ZCHKHS: ${} Eigenvectors from ${} incorrectly normalized.\n Bits of error=${.g10_3},${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );
 9997 FORMAT( ' ZCHKHS: Selected ${} Eigenvectors from ${} do not match other eigenvectors ${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );
      }
