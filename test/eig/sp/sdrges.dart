      void sdrges(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHAR, ALPHAI, BETA, WORK, LWORK, RESULT, BWORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES;
      double               THRESH;
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      double               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDA, * ), BETA( * ), Q( LDQ, * ), RESULT( 13 ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      bool               BADNN, ILABAD;
      String             SORT;
      int                I, I1, IADD, IERR, IINFO, IN, ISORT, J, JC, JR, JSIZE, JTYPE, KNTEIG, MAXWRK, MINWRK, MTYPES, N, N1, NB, NERRS, NMATS, NMAX, NTEST, NTESTT, RSUB, SDIM;
      double               SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV;
      int                IASIGN( MAXTYP ), IBSIGN( MAXTYP ), IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      double               RMAGN( 0: 3 );
      // ..
      // .. External Functions ..
      //- bool               SLCTES;
      //- int                ILAENV;
      //- REAL               SLAMCH, SLARND;
      // EXTERNAL SLCTES, ILAENV, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, SGET51, SGET53, SGET54, SGGES, SLACPY, SLARFG, SLASET, SLATM4, SORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SIGN
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
      const IASIGN = [ 0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0 ];
      const IBSIGN = [ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,];

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
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      MINWRK = 1;
      if ( INFO == 0 && LWORK >= 1 ) {
         MINWRK = max( 10*( NMAX+1 ), 3*NMAX*NMAX );
         NB = max( 1, ilaenv( 1, 'SGEQRF', ' ', NMAX, NMAX, -1, -1 ), ilaenv( 1, 'SORMQR', 'LT', NMAX, NMAX, NMAX, -1 ), ilaenv( 1, 'SORGQR', ' ', NMAX, NMAX, NMAX, -1 ) );
         MAXWRK = max( 10*( NMAX+1 ), 2*NMAX+NMAX*NB, 3*NMAX*NMAX );
         WORK[1] = MAXWRK;
      }

      if (LWORK < MINWRK) INFO = -20;

      if ( INFO != 0 ) {
         xerbla('SDRGES', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      SAFMIN = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
      SAFMIN = SAFMIN / ULP;
      SAFMAX = ONE / SAFMIN;
      ULPINV = ONE / ULP;

      // The values RMAGN(2:3) depend on N, see below.

      RMAGN[0] = ZERO;
      RMAGN[1] = ONE;

      // Loop over matrix sizes

      NTESTT = 0;
      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 190
         N = NN( JSIZE );
         N1 = max( 1, N );
         RMAGN[2] = SAFMAX*ULP / REAL( N1 );
         RMAGN[3] = SAFMIN*ULPINV*double( N1 );

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         // Loop over matrix types

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 180
            if( !DOTYPE( JTYPE ) ) GO TO 180;
            NMATS = NMATS + 1;
            NTEST = 0;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            // Initialize RESULT

            for (J = 1; J <= 13; J++) { // 30
               RESULT[J] = ZERO;
            } // 30

            // Generate test matrices A and B

            // Description of control parameters:

            // KCLASS: =1 means w/o rotation, =2 means w/ rotation,
                    // =3 means random.
            // KATYPE: the "type" to be passed to SLATM4 for computing A.
            // KAZERO: the pattern of zeros on the diagonal for A:
                    // =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
                    // =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
                    // =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
                    // non-zero entries.)
            // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
                    // =2: large, =3: small.
            // IASIGN: 1 if the diagonal elements of A are to be
                    // multiplied by a random magnitude 1 number, =2 if
                    // randomly chosen diagonal blocks are to be rotated
                    // to form 2x2 blocks.
            // KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN: used to implement KAMAGN and KBMAGN.

            if (MTYPES > MAXTYP) GO TO 110;
            IINFO = 0;
            if ( KCLASS( JTYPE ) < 3 ) {

               // Generate A (w/o rotation)

               if ( ( KATYPE( JTYPE ) ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) slaset( 'Full', N, N, ZERO, ZERO, A, LDA );
               } else {
                  IN = N;
               }
               slatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), IASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) );
               if (IADD > 0 && IADD <= N) A( IADD, IADD ) = ONE;

               // Generate B (w/o rotation)

               if ( ( KBTYPE( JTYPE ) ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) slaset( 'Full', N, N, ZERO, ZERO, B, LDA );
               } else {
                  IN = N;
               }
               slatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), IBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) );
               if (IADD != 0 && IADD <= N) B( IADD, IADD ) = ONE;

               if ( KCLASS( JTYPE ) == 2 && N > 0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  for (JC = 1; JC <= N - 1; JC++) { // 50
                     for (JR = JC; JR <= N; JR++) { // 40
                        Q[JR][JC] = SLARND( 3, ISEED );
                        Z[JR][JC] = SLARND( 3, ISEED );
                     } // 40
                     slarfg(N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK[2*N+JC] = sign( ONE, Q( JC, JC ) );
                     Q[JC][JC] = ONE;
                     slarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK[3*N+JC] = sign( ONE, Z( JC, JC ) );
                     Z[JC][JC] = ONE;
                  } // 50
                  Q[N][N] = ONE;
                  WORK[N] = ZERO;
                  WORK[3*N] = sign( ONE, SLARND( 2, ISEED ) );
                  Z[N][N] = ONE;
                  WORK[2*N] = ZERO;
                  WORK[4*N] = sign( ONE, SLARND( 2, ISEED ) );

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 70
                     for (JR = 1; JR <= N; JR++) { // 60
                        A[JR][JC] = WORK( 2*N+JR )*WORK( 3*N+JC )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* B( JR, JC );
                     } // 60
                  } // 70
                  CALL SORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL SORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL SORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL SORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 90
                  for (JR = 1; JR <= N; JR++) { // 80
                     A[JR][JC] = RMAGN( KAMAGN( JTYPE ) )* SLARND( 2, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* SLARND( 2, ISEED );
                  } // 80
               } // 90
            }

            } // 100

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 110

            for (I = 1; I <= 13; I++) { // 120
               RESULT[I] = -ONE;
            } // 120

            // Test with and without sorting of eigenvalues

            for (ISORT = 0; ISORT <= 1; ISORT++) { // 150
               if ( ISORT == 0 ) {
                  SORT = 'N';
                  RSUB = 0;
               } else {
                  SORT = 'S';
                  RSUB = 5;
               }

               // Call SGGES to compute H, T, Q, Z, alpha, and beta.

               slacpy('Full', N, N, A, LDA, S, LDA );
               slacpy('Full', N, N, B, LDA, T, LDA );
               NTEST = 1 + RSUB + ISORT;
               RESULT[1+RSUB+ISORT] = ULPINV;
               sgges('V', 'V', SORT, SLCTES, N, S, LDA, T, LDA, SDIM, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, BWORK, IINFO );
               if ( IINFO != 0 && IINFO != N+2 ) {
                  RESULT[1+RSUB+ISORT] = ULPINV;
                  WRITE( NOUNIT, FMT = 9999 )'SGGES', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  GO TO 160;
               }

               NTEST = 4 + RSUB;

               // Do tests 1--4 (or tests 7--9 when reordering )

               if ( ISORT == 0 ) {
                  sget51(1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 1 ) );
                  sget51(1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 2 ) );
               } else {
                  sget54(N, A, LDA, B, LDA, S, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 7 ) );
               }
               sget51(3, N, A, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK, RESULT( 3+RSUB ) );
               sget51(3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK, RESULT( 4+RSUB ) );

               // Do test 5 and 6 (or Tests 10 and 11 when reordering):
               // check Schur form of A and compare eigenvalues with
               // diagonals.

               NTEST = 6 + RSUB;
               TEMP1 = ZERO;

               for (J = 1; J <= N; J++) { // 130
                  ILABAD = false;
                  if ( ALPHAI( J ) == ZERO ) {
                     TEMP2 = ( ABS( ALPHAR( J )-S( J, J ) ) / max( SAFMIN, ( ALPHAR( J ) ).abs(), ( S( J, J ) ).abs() )+ABS( BETA( J )-T( J, J ) ) / max( SAFMIN, ( BETA( J ) ).abs(), ( T( J, J ) ).abs() ) ) / ULP;

                     if ( J < N ) {
                        if ( S( J+1, J ) != ZERO ) {
                           ILABAD = true;
                           RESULT[5+RSUB] = ULPINV;
                        }
                     }
                     if ( J > 1 ) {
                        if ( S( J, J-1 ) != ZERO ) {
                           ILABAD = true;
                           RESULT[5+RSUB] = ULPINV;
                        }
                     }

                  } else {
                     if ( ALPHAI( J ) > ZERO ) {
                        I1 = J;
                     } else {
                        I1 = J - 1;
                     }
                     if ( I1 <= 0 || I1 >= N ) {
                        ILABAD = true;
                     } else if ( I1 < N-1 ) {
                        if ( S( I1+2, I1+1 ) != ZERO ) {
                           ILABAD = true;
                           RESULT[5+RSUB] = ULPINV;
                        }
                     } else if ( I1 > 1 ) {
                        if ( S( I1, I1-1 ) != ZERO ) {
                           ILABAD = true;
                           RESULT[5+RSUB] = ULPINV;
                        }
                     }
                     if ( !ILABAD ) {
                        sget53(S( I1, I1 ), LDA, T( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IERR );
                        if ( IERR >= 3 ) {
                           WRITE( NOUNIT, FMT = 9998 )IERR, J, N, JTYPE, IOLDSD;
                           INFO = ( IERR ).abs();
                        }
                     } else {
                        TEMP2 = ULPINV;
                     }

                  }
                  TEMP1 = max( TEMP1, TEMP2 );
                  if ( ILABAD ) {
                     WRITE( NOUNIT, FMT = 9997 )J, N, JTYPE, IOLDSD;
                  }
               } // 130
               RESULT[6+RSUB] = TEMP1;

               if ( ISORT >= 1 ) {

                  // Do test 12

                  NTEST = 12;
                  RESULT[12] = ZERO;
                  KNTEIG = 0;
                  for (I = 1; I <= N; I++) { // 140
                     if ( SLCTES( ALPHAR( I ), ALPHAI( I ), BETA( I ) ) || SLCTES( ALPHAR( I ), -ALPHAI( I ), BETA( I ) ) ) {
                        KNTEIG = KNTEIG + 1;
                     }
                     if ( I < N ) {
                        if ( ( SLCTES( ALPHAR( I+1 ), ALPHAI( I+1 ), BETA( I+1 ) ) || SLCTES( ALPHAR( I+1 ), -ALPHAI( I+1 ), BETA( I+1 ) ) ) && ( !( SLCTES( ALPHAR( I ), ALPHAI( I ), BETA( I ) ) || SLCTES( ALPHAR( I ), -ALPHAI( I ), BETA( I ) ) ) ) && IINFO != N+2 ) {
                           RESULT[12] = ULPINV;
                        }
                     }
                  } // 140
                  if ( SDIM != KNTEIG ) {
                     RESULT[12] = ULPINV;
                  }
               }

            } // 150

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 160

            NTESTT = NTESTT + NTEST;

            // Print out tests which fail.

            for (JR = 1; JR <= NTEST; JR++) { // 170
               if ( RESULT( JR ) >= THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9996 )'SGS';

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9995 );
                     WRITE( NOUNIT, FMT = 9994 );
                     WRITE( NOUNIT, FMT = 9993 )'Orthogonal';

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', J = 1, 8 );

                  }
                  NERRS = NERRS + 1;
                  if ( RESULT( JR ) < 10000.0 ) {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  } else {
                     WRITE( NOUNIT, FMT = 9990 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  }
               }
            } // 170

         } // 180
      } // 190

      // Summary

      alasvm('SGS', NOUNIT, NERRS, NTESTT, 0 );

      WORK[1] = MAXWRK;

      return;

 9999 FORMAT( ' SDRGES: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${i4(4, ',')}', I5, ')' );

 9998 FORMAT( ' SDRGES: SGET53 returned INFO=${.i1} for eigenvalue ${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${i4(4, ',')}', I5, ')' );

 9997 FORMAT( ' SDRGES: S not in Schur form at eigenvalue ${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${i5(3, ',')}', I5, ')' );

 9996 FORMAT('\n ${.a3} -- Real Generalized Schur form driver' );

 9995 FORMAT( ' Matrix types (see SDRGES for details): ' );

 9994 FORMAT( ' Special Matrices:${' ' * 23}(J''=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  6=(diag(J'',I), diag(I,J''))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)' );
 9993 FORMAT( ' Matrices Rotated by Random ${} Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.' );

 9992 FORMAT('\n Tests performed:  (S is Schur, T is triangular, Q and Z are ${},\n${' ' * 19}l and r are the appropriate left and right\n${' ' * 19}eigenvectors, resp., a is alpha, b is beta, and', / 19X, A, ' means ${}.)\n Without ordering: \n  1 = | A - Q S Z${} | / ( |A| n ulp )      2 = | B - Q T Z${} | / ( |B| n ulp )\n  3 = | I - QQ${} | / ( n ulp )             4 = | I - ZZ${} | / ( n ulp )\n  5 = A is in Schur form S\n  6 = difference between (alpha,beta) and diagonals of (S,T)\n With ordering: \n  7 = | (A,B) - Q (S,T) Z${} | / ( |(A,B)| n ulp )  \n  8 = | I - QQ${} | / ( n ulp )            9 = | I - ZZ${} | / ( n ulp )\n 10 = A is in Schur form S\n 11 = difference between (alpha,beta) and diagonals of (S,T)\n 12 = SDIM is the correct number of selected eigenvalues\n');
 9991 FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i2} is' F8.2 );
 9990 FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i2} is', 1P, E10.3 );
      }
