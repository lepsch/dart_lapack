      void cdrgev3(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, QE, LDQE, ALPHA, BETA, ALPHA1, BETA1, WORK, LWORK, RWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDQ, LDQE, LWORK, NOUNIT, NSIZES, NTYPES;
      double               THRESH;
      bool               DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      double               RESULT( * ), RWORK( * );
      Complex            A( LDA, * ), ALPHA( * ), ALPHA1( * ), B( LDA, * ), BETA( * ), BETA1( * ), Q( LDQ, * ), QE( LDQE, * ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      bool               BADNN;
      int                I, IADD, IERR, IN, J, JC, JR, JSIZE, JTYPE, MAXWRK, MINWRK, MTYPES, N, N1, NB, NERRS, NMATS, NMAX, NTESTT;
      double               SAFMAX, SAFMIN, ULP, ULPINV;
      Complex            CTEMP;
      bool               LASIGN( MAXTYP ), LBSIGN( MAXTYP );
      int                IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      double               RMAGN( 0: 3 );
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SLAMCH;
      //- COMPLEX            CLARND;
      // EXTERNAL ILAENV, SLAMCH, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, CGET52, CGGEV3, CLACPY, CLARFG, CLASET, CLATM4, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, REAL, SIGN
      // ..
      // .. Data statements ..
      const KCLASS = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,];
      const KZ1 = [ 0, 1, 2, 1, 3, 3 ];
      const KZ2 = [ 0, 0, 1, 2, 1, 1 ];
      const KADD = [ 0, 0, 0, 0, 3, 2 ];
      const KATYPE = [ 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0 ];
      const KBTYPE = [ 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0 ];
      const KAZERO = [ 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1 ];
      const KBZERO = [ 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1 ];
      const KAMAGN = [ 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1 ];
      const KBMAGN = [ 1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1 ];
      const KTRIAN = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,];
      const LASIGN = [ false, false, false, false, false, false, true , false , true, true, false, false, true, true, true, false , true , false, false, false, true, true, true, true, true, false ];
      const LBSIGN = [ false, false, false, false, false, false, false, true , false, false, true, true, false, false, true , false , true , false, false, false, false, false, false, false, false, false,];

      // Check for errors

      INFO = 0;

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
      } // 10

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
      } else if ( LDQ <= 1 || LDQ < NMAX ) {
         INFO = -14;
      } else if ( LDQE <= 1 || LDQE < NMAX ) {
         INFO = -17;
      }

      // Compute workspace
      //  (Note: Comments in the code beginning "Workspace:" describe the
      //   minimal amount of workspace needed at that point in the code,
      //   as well as the preferred amount for good performance.
      //   NB refers to the optimal block size for the immediately
      //   following subroutine, as returned by ILAENV.

      MINWRK = 1;
      if ( INFO == 0 && LWORK >= 1 ) {
         MINWRK = NMAX*( NMAX+1 );
         NB = max( 1, ilaenv( 1, 'CGEQRF', ' ', NMAX, NMAX, -1, -1 ), ilaenv( 1, 'CUNMQR', 'LC', NMAX, NMAX, NMAX, -1 ), ilaenv( 1, 'CUNGQR', ' ', NMAX, NMAX, NMAX, -1 ) );
         MAXWRK = max( 2*NMAX, NMAX*( NB+1 ), NMAX*( NMAX+1 ) );
         WORK[1] = MAXWRK;
      }

      if (LWORK < MINWRK) INFO = -23;

      if ( INFO != 0 ) {
         xerbla('CDRGEV3', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      ULP = SLAMCH( 'Precision' );
      SAFMIN = SLAMCH( 'Safe minimum' );
      SAFMIN = SAFMIN / ULP;
      SAFMAX = ONE / SAFMIN;
      ULPINV = ONE / ULP;

      // The values RMAGN(2:3) depend on N, see below.

      RMAGN[0] = ZERO;
      RMAGN[1] = ONE;

      // Loop over sizes, types

      NTESTT = 0;
      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 220
         N = NN( JSIZE );
         N1 = max( 1, N );
         RMAGN[2] = SAFMAX*ULP / REAL( N1 );
         RMAGN[3] = SAFMIN*ULPINV*N1;

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 210
            if( !DOTYPE( JTYPE ) ) GO TO 210;
            NMATS = NMATS + 1;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            // Generate test matrices A and B

            // Description of control parameters:

            // KCLASS: =1 means w/o rotation, =2 means w/ rotation,
            //         =3 means random.
            // KATYPE: the "type" to be passed to CLATM4 for computing A.
            // KAZERO: the pattern of zeros on the diagonal for A:
            //         =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
            //         =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
            //         =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
            //         non-zero entries.)
            // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
            //         =2: large, =3: small.
            // LASIGN: true if the diagonal elements of A are to be
            //         multiplied by a random magnitude 1 number.
            // KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN: used to implement KAMAGN and KBMAGN.

            if (MTYPES > MAXTYP) GO TO 100;
            IERR = 0;
            if ( KCLASS( JTYPE ) < 3 ) {

               // Generate A (w/o rotation)

               if ( ( KATYPE( JTYPE ) ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) claset( 'Full', N, N, CZERO, CZERO, A, LDA );
               } else {
                  IN = N;
               }
               clatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) );
               if (IADD > 0 && IADD <= N) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) );

               // Generate B (w/o rotation)

               if ( ( KBTYPE( JTYPE ) ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) claset( 'Full', N, N, CZERO, CZERO, B, LDA );
               } else {
                  IN = N;
               }
               clatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) );
               if (IADD != 0 && IADD <= N) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) );

               if ( KCLASS( JTYPE ) == 2 && N > 0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  for (JC = 1; JC <= N - 1; JC++) { // 40
                     for (JR = JC; JR <= N; JR++) { // 30
                        Q[JR][JC] = CLARND( 3, ISEED );
                        Z[JR][JC] = CLARND( 3, ISEED );
                     } // 30
                     clarfg(N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK[2*N+JC] = sign( ONE, double( Q( JC, JC ) ) );
                     Q[JC][JC] = CONE;
                     clarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK[3*N+JC] = sign( ONE, double( Z( JC, JC ) ) );
                     Z[JC][JC] = CONE;
                  } // 40
                  CTEMP = CLARND( 3, ISEED );
                  Q[N][N] = CONE;
                  WORK[N] = CZERO;
                  WORK[3*N] = CTEMP / ( CTEMP ).abs();
                  CTEMP = CLARND( 3, ISEED );
                  Z[N][N] = CONE;
                  WORK[2*N] = CZERO;
                  WORK[4*N] = CTEMP / ( CTEMP ).abs();

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 60
                     for (JR = 1; JR <= N; JR++) { // 50
                        A[JR][JC] = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* B( JR, JC );
                     } // 50
                  } // 60
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
                  CALL CUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
                  CALL CUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 80
                  for (JR = 1; JR <= N; JR++) { // 70
                     A[JR][JC] = RMAGN( KAMAGN( JTYPE ) )* CLARND( 4, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* CLARND( 4, ISEED );
                  } // 70
               } // 80
            }

            } // 90

            if ( IERR != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IERR, N, JTYPE, IOLDSD;
               INFO = ( IERR ).abs();
               return;
            }

            } // 100

            for (I = 1; I <= 7; I++) { // 110
               RESULT[I] = -ONE;
            } // 110

            // Call XLAENV to set the parameters used in CLAQZ0

            xlaenv(12, 10 );
            xlaenv(13, 12 );
            xlaenv(14, 13 );
            xlaenv(15, 2 );
            xlaenv(17, 10 );

            // Call CGGEV3 to compute eigenvalues and eigenvectors.

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('V', 'V', N, S, LDA, T, LDA, ALPHA, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'CGGEV31', IERR, N, JTYPE, IOLDSD;
               INFO = ( IERR ).abs();
               GO TO 190;
            }

            // Do the tests (1) and (2)

            cget52( true , N, A, LDA, B, LDA, Q, LDQ, ALPHA, BETA, WORK, RWORK, RESULT( 1 ) );
            if ( RESULT( 2 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'CGGEV31', RESULT( 2 ), N, JTYPE, IOLDSD;
            }

            // Do the tests (3) and (4)

            cget52( false , N, A, LDA, B, LDA, Z, LDQ, ALPHA, BETA, WORK, RWORK, RESULT( 3 ) );
            if ( RESULT( 4 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'CGGEV31', RESULT( 4 ), N, JTYPE, IOLDSD;
            }

            // Do test (5)

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('N', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'CGGEV32', IERR, N, JTYPE, IOLDSD;
               INFO = ( IERR ).abs();
               GO TO 190;
            }

            for (J = 1; J <= N; J++) { // 120
               if( ALPHA( J ) != ALPHA1( J ) || BETA( J ) != BETA1( J ) ) RESULT( 5 ) = ULPINV;
            } // 120

            // Do the test (6): Compute eigenvalues and left eigenvectors,
            // and test them

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('V', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, QE, LDQE, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'CGGEV33', IERR, N, JTYPE, IOLDSD;
               INFO = ( IERR ).abs();
               GO TO 190;
            }


            for (J = 1; J <= N; J++) { // 130
               if ( ALPHA( J ) != ALPHA1( J ) || BETA( J ) != BETA1( J ) ) {
                  RESULT[6] = ULPINV;
               }
            } // 130

            for (J = 1; J <= N; J++) { // 150
               for (JC = 1; JC <= N; JC++) { // 140
                  if ( Q( J, JC ) != QE( J, JC ) ) {
                     RESULT[6] = ULPINV;
                  }
               } // 140
            } // 150

            // DO the test (7): Compute eigenvalues and right eigenvectors,
            // and test them

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('N', 'V', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, QE, LDQE, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'CGGEV34', IERR, N, JTYPE, IOLDSD;
               INFO = ( IERR ).abs();
               GO TO 190;
            }

            for (J = 1; J <= N; J++) { // 160
               if( ALPHA( J ) != ALPHA1( J ) || BETA( J ) != BETA1( J ) )RESULT( 7 ) = ULPINV;
            } // 160

            for (J = 1; J <= N; J++) { // 180
               for (JC = 1; JC <= N; JC++) { // 170
                  if( Z( J, JC ) != QE( J, JC ) ) RESULT( 7 ) = ULPINV;
               } // 170
            } // 180

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 190

            NTESTT = NTESTT + 7;

            // Print out tests which fail.

            for (JR = 1; JR <= 7; JR++) { // 200
               if ( RESULT( JR ) >= THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9997 )'CGV';

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 );
                     WRITE( NOUNIT, FMT = 9995 );
                     WRITE( NOUNIT, FMT = 9994 )'Orthogonal';

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 );

                  }
                  NERRS = NERRS + 1;
                  if ( RESULT( JR ) < 10000.0 ) {
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  } else {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  }
               }
            } // 200

         } // 210
      } // 220

      // Summary

      alasvm('CGV3', NOUNIT, NERRS, NTESTT, 0 );

      WORK[1] = MAXWRK;

      return;

 9999 FORMAT( ' CDRGEV3: ${} returned INFO=${.i6}.\n${' ' * 3}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );

 9998 FORMAT( ' CDRGEV3: ${} Eigenvectors from ${} incorrectly normalized.\n Bits of error=${.g10_3},${' ' * 3}N=${.i4}, JTYPE=${.i3}, ISEED=(${i4(3, ',')}', I5, ')' );

 9997 FORMAT('\n ${.a3} -- Complex Generalized eigenvalue problem driver' );

 9996 FORMAT( ' Matrix types (see CDRGEV3 for details): ' );

 9995 FORMAT( ' Special Matrices:${' ' * 23}(J''=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  6=(diag(J'',I), diag(I,J''))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)' );
 9994 FORMAT( ' Matrices Rotated by Random ${} Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.' );

 9993 FORMAT('\n Tests performed:    \n 1 = max | ( b A - a B )''*l | / const.,\n 2 = | |VR(i)| - 1 | / ulp,\n 3 = max | ( b A - a B )*r | / const.\n 4 = | |VL(i)| - 1 | / ulp,\n 5 = 0 if W same no matter if r or l computed,\n 6 = 0 if l same no matter if l computed,\n 7 = 0 if r same no matter if r computed,/n ');
 9992 FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i2} is${.f8_2}');
 9991 FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i2} is', 1P, E10.3 );
      }
