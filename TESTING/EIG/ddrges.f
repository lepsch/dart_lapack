      SUBROUTINE DDRGES( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHAR, ALPHAI, BETA, WORK, LWORK, RESULT, BWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDA, * ), BETA( * ), Q( LDQ, * ), RESULT( 13 ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, ILABAD;
      String             SORT;
      int                I, I1, IADD, IERR, IINFO, IN, ISORT, J, JC, JR, JSIZE, JTYPE, KNTEIG, MAXWRK, MINWRK, MTYPES, N, N1, NB, NERRS, NMATS, NMAX, NTEST, NTESTT, RSUB, SDIM;
      double             SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV;
      // ..
      // .. Local Arrays ..
      int                IASIGN( MAXTYP ), IBSIGN( MAXTYP ), IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      double             RMAGN( 0: 3 );
      // ..
      // .. External Functions ..
      bool               DLCTES;
      int                ILAENV;
      double             DLAMCH, DLARND;
      // EXTERNAL DLCTES, ILAENV, DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, DGET51, DGET53, DGET54, DGGES, DLACPY, DLARFG, DLASET, DLATM4, DORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SIGN
      // ..
      // .. Data statements ..
      DATA               KCLASS / 15*1, 10*2, 1*3 /
      DATA               KZ1 / 0, 1, 2, 1, 3, 3 /
      DATA               KZ2 / 0, 0, 1, 2, 1, 1 /
      DATA               KADD / 0, 0, 0, 0, 3, 2 /
      DATA               KATYPE / 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4*4, 0 /
      DATA               KBTYPE / 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8*8, 0 /
      DATA               KAZERO / 6*1, 2, 1, 2*2, 2*1, 2*2, 3, 1, 3, 4*5, 4*3, 1 /
      DATA               KBZERO / 6*1, 1, 2, 2*1, 2*2, 2*1, 4, 1, 4, 4*6, 4*4, 1 /
      DATA               KAMAGN / 8*1, 2, 3, 2, 3, 2, 3, 7*1, 2, 3, 3, 2, 1 /
      DATA               KBMAGN / 8*1, 3, 2, 3, 2, 2, 3, 7*1, 3, 2, 3, 2, 1 /
      DATA               KTRIAN / 16*0, 10*1 /
      DATA               IASIGN / 6*0, 2, 0, 2*2, 2*0, 3*2, 0, 2, 3*0, 5*2, 0 /
      DATA               IBSIGN / 7*0, 2, 2*0, 2*2, 2*0, 2, 0, 2, 9*0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0

      BADNN = false;
      NMAX = 1
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ) < 0 ) BADNN = true;
      } // 10

      if ( NSIZES < 0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES < 0 ) {
         INFO = -3
      } else if ( THRESH < ZERO ) {
         INFO = -6
      } else if ( LDA.LE.1 || LDA < NMAX ) {
         INFO = -9
      } else if ( LDQ.LE.1 || LDQ < NMAX ) {
         INFO = -14
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      MINWRK = 1
      if ( INFO == 0 && LWORK.GE.1 ) {
         MINWRK = MAX( 10*( NMAX+1 ), 3*NMAX*NMAX )
         NB = MAX( 1, ILAENV( 1, 'DGEQRF', ' ', NMAX, NMAX, -1, -1 ), ILAENV( 1, 'DORMQR', 'LT', NMAX, NMAX, NMAX, -1 ), ILAENV( 1, 'DORGQR', ' ', NMAX, NMAX, NMAX, -1 ) )
         MAXWRK = MAX( 10*( NMAX+1 ), 2*NMAX+NMAX*NB, 3*NMAX*NMAX )
         WORK( 1 ) = MAXWRK
      }

      if (LWORK < MINWRK) INFO = -20;

      if ( INFO != 0 ) {
         xerbla('DDRGES', -INFO );
         RETURN
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) RETURN;

      SAFMIN = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SAFMIN = SAFMIN / ULP
      SAFMAX = ONE / SAFMIN
      ULPINV = ONE / ULP

      // The values RMAGN(2:3) depend on N, see below.

      RMAGN( 0 ) = ZERO
      RMAGN( 1 ) = ONE

      // Loop over matrix sizes

      NTESTT = 0
      NERRS = 0
      NMATS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 190
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / DBLE( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*DBLE( N1 )

         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         // Loop over matrix types

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 180
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 180
            NMATS = NMATS + 1
            NTEST = 0

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD( J ) = ISEED( J )
            } // 20

            // Initialize RESULT

            for (J = 1; J <= 13; J++) { // 30
               RESULT( J ) = ZERO
            } // 30

            // Generate test matrices A and B

            // Description of control parameters:

            // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
                    // =3 means random.
            // KATYPE: the "type" to be passed to DLATM4 for computing A.
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

            if (MTYPES.GT.MAXTYP) GO TO 110;
            IINFO = 0
            if ( KCLASS( JTYPE ) < 3 ) {

               // Generate A (w/o rotation)

               if ( ABS( KATYPE( JTYPE ) ) == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  if (IN != N) CALL DLASET( 'Full', N, N, ZERO, ZERO, A, LDA );
               } else {
                  IN = N
               }
               dlatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), IASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) )
               if (IADD.GT.0 && IADD.LE.N) A( IADD, IADD ) = ONE;

               // Generate B (w/o rotation)

               if ( ABS( KBTYPE( JTYPE ) ) == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  if (IN != N) CALL DLASET( 'Full', N, N, ZERO, ZERO, B, LDA );
               } else {
                  IN = N
               }
               dlatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), IBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) )
               if (IADD != 0 && IADD.LE.N) B( IADD, IADD ) = ONE;

               if ( KCLASS( JTYPE ) == 2 && N.GT.0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  for (JC = 1; JC <= N - 1; JC++) { // 50
                     for (JR = JC; JR <= N; JR++) { // 40
                        Q( JR, JC ) = DLARND( 3, ISEED )
                        Z( JR, JC ) = DLARND( 3, ISEED )
                     } // 40
                     dlarfg(N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK( 2*N+JC ) = SIGN( ONE, Q( JC, JC ) )
                     Q( JC, JC ) = ONE
                     dlarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK( 3*N+JC ) = SIGN( ONE, Z( JC, JC ) )
                     Z( JC, JC ) = ONE
                  } // 50
                  Q( N, N ) = ONE
                  WORK( N ) = ZERO
                  WORK( 3*N ) = SIGN( ONE, DLARND( 2, ISEED ) )
                  Z( N, N ) = ONE
                  WORK( 2*N ) = ZERO
                  WORK( 4*N ) = SIGN( ONE, DLARND( 2, ISEED ) )

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 70
                     for (JR = 1; JR <= N; JR++) { // 60
                        A( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* B( JR, JC )
                     } // 60
                  } // 70
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 90
                  for (JR = 1; JR <= N; JR++) { // 80
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* DLARND( 2, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* DLARND( 2, ISEED )
                  } // 80
               } // 90
            }

            } // 100

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

            } // 110

            for (I = 1; I <= 13; I++) { // 120
               RESULT( I ) = -ONE
            } // 120

            // Test with and without sorting of eigenvalues

            for (ISORT = 0; ISORT <= 1; ISORT++) { // 150
               if ( ISORT == 0 ) {
                  SORT = 'N'
                  RSUB = 0
               } else {
                  SORT = 'S'
                  RSUB = 5
               }

               // Call DGGES to compute H, T, Q, Z, alpha, and beta.

               dlacpy('Full', N, N, A, LDA, S, LDA );
               dlacpy('Full', N, N, B, LDA, T, LDA );
               NTEST = 1 + RSUB + ISORT
               RESULT( 1+RSUB+ISORT ) = ULPINV
               dgges('V', 'V', SORT, DLCTES, N, S, LDA, T, LDA, SDIM, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, BWORK, IINFO );
               if ( IINFO != 0 && IINFO != N+2 ) {
                  RESULT( 1+RSUB+ISORT ) = ULPINV
                  WRITE( NOUNIT, FMT = 9999 )'DGGES', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 160
               }

               NTEST = 4 + RSUB

               // Do tests 1--4 (or tests 7--9 when reordering )

               if ( ISORT == 0 ) {
                  dget51(1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 1 ) );
                  dget51(1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 2 ) );
               } else {
                  dget54(N, A, LDA, B, LDA, S, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 7 ) );
               }
               dget51(3, N, A, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK, RESULT( 3+RSUB ) );
               dget51(3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK, RESULT( 4+RSUB ) );

               // Do test 5 and 6 (or Tests 10 and 11 when reordering):
               // check Schur form of A and compare eigenvalues with
               // diagonals.

               NTEST = 6 + RSUB
               TEMP1 = ZERO

               for (J = 1; J <= N; J++) { // 130
                  ILABAD = false;
                  if ( ALPHAI( J ) == ZERO ) {
                     TEMP2 = ( ABS( ALPHAR( J )-S( J, J ) ) / MAX( SAFMIN, ABS( ALPHAR( J ) ), ABS( S( J, J ) ) )+ABS( BETA( J )-T( J, J ) ) / MAX( SAFMIN, ABS( BETA( J ) ), ABS( T( J, J ) ) ) ) / ULP

                     if ( J < N ) {
                        if ( S( J+1, J ) != ZERO ) {
                           ILABAD = true;
                           RESULT( 5+RSUB ) = ULPINV
                        }
                     }
                     if ( J.GT.1 ) {
                        if ( S( J, J-1 ) != ZERO ) {
                           ILABAD = true;
                           RESULT( 5+RSUB ) = ULPINV
                        }
                     }

                  } else {
                     if ( ALPHAI( J ).GT.ZERO ) {
                        I1 = J
                     } else {
                        I1 = J - 1
                     }
                     if ( I1.LE.0 || I1.GE.N ) {
                        ILABAD = true;
                     } else if ( I1 < N-1 ) {
                        if ( S( I1+2, I1+1 ) != ZERO ) {
                           ILABAD = true;
                           RESULT( 5+RSUB ) = ULPINV
                        }
                     } else if ( I1.GT.1 ) {
                        if ( S( I1, I1-1 ) != ZERO ) {
                           ILABAD = true;
                           RESULT( 5+RSUB ) = ULPINV
                        }
                     }
                     if ( .NOT.ILABAD ) {
                        dget53(S( I1, I1 ), LDA, T( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IERR );
                        if ( IERR.GE.3 ) {
                           WRITE( NOUNIT, FMT = 9998 )IERR, J, N, JTYPE, IOLDSD
                           INFO = ABS( IERR )
                        }
                     } else {
                        TEMP2 = ULPINV
                     }

                  }
                  TEMP1 = MAX( TEMP1, TEMP2 )
                  if ( ILABAD ) {
                     WRITE( NOUNIT, FMT = 9997 )J, N, JTYPE, IOLDSD
                  }
               } // 130
               RESULT( 6+RSUB ) = TEMP1

               if ( ISORT.GE.1 ) {

                  // Do test 12

                  NTEST = 12
                  RESULT( 12 ) = ZERO
                  KNTEIG = 0
                  for (I = 1; I <= N; I++) { // 140
                     if ( DLCTES( ALPHAR( I ), ALPHAI( I ), BETA( I ) ) || DLCTES( ALPHAR( I ), -ALPHAI( I ), BETA( I ) ) ) {
                        KNTEIG = KNTEIG + 1
                     }
                     if ( I < N ) {
                        if ( ( DLCTES( ALPHAR( I+1 ), ALPHAI( I+1 ), BETA( I+1 ) ) || DLCTES( ALPHAR( I+1 ), -ALPHAI( I+1 ), BETA( I+1 ) ) ) && ( .NOT.( DLCTES( ALPHAR( I ), ALPHAI( I ), BETA( I ) ) || DLCTES( ALPHAR( I ), -ALPHAI( I ), BETA( I ) ) ) ) && IINFO != N+2 ) {
                           RESULT( 12 ) = ULPINV
                        }
                     }
                  } // 140
                  if ( SDIM != KNTEIG ) {
                     RESULT( 12 ) = ULPINV
                  }
               }

            } // 150

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 160

            NTESTT = NTESTT + NTEST

            // Print out tests which fail.

            for (JR = 1; JR <= NTEST; JR++) { // 170
               if ( RESULT( JR ).GE.THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9996 )'DGS'

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )
                     WRITE( NOUNIT, FMT = 9993 )'Orthogonal'

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', J = 1, 8 )

                  }
                  NERRS = NERRS + 1
                  if ( RESULT( JR ) < 10000.0D0 ) {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  } else {
                     WRITE( NOUNIT, FMT = 9990 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  }
               }
            } // 170

         } // 180
      } // 190

      // Summary

      alasvm('DGS', NOUNIT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK

      RETURN

 9999 FORMAT( ' DDRGES: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' )

 9998 FORMAT( ' DDRGES: DGET53 returned INFO=', I1, ' for eigenvalue ', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' )

 9997 FORMAT( ' DDRGES: S not in Schur form at eigenvalue ', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9996 FORMAT( / 1X, A3, ' -- Real Generalized Schur form driver' )

 9995 FORMAT( ' Matrix types (see DDRGES for details): ' )

 9994 FORMAT( ' Special Matrices:', 23X, '(J''=transposed Jordan block)', / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ', '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ', 'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I', ')  11=(large*I, small*D)  13=(large*D, large*I)', / '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ', ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9993 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:', / '  16=Transposed Jordan Blocks             19=geometric ', 'alpha, beta=0,1', / '  17=arithm. alpha&beta             ', '      20=arithmetic alpha, beta=0,1', / '  18=clustered ', 'alpha, beta=0,1            21=random alpha, beta=0,1', / ' Large & Small Matrices:', / '  22=(large, small)   ', '23=(small,large)    24=(small,small)    25=(large,large)', / '  26=random O(1) matrices.' )

 9992 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ', 'Q and Z are ', A, ',', / 19X, 'l and r are the appropriate left and right', / 19X, 'eigenvectors, resp., a is alpha, b is beta, and', / 19X, A, ' means ', A, '.)', / ' Without ordering: ', / '  1 = | A - Q S Z', A, ' | / ( |A| n ulp )      2 = | B - Q T Z', A, ' | / ( |B| n ulp )', / '  3 = | I - QQ', A, ' | / ( n ulp )             4 = | I - ZZ', A, ' | / ( n ulp )', / '  5 = A is in Schur form S', / '  6 = difference between (alpha,beta)', ' and diagonals of (S,T)', / ' With ordering: ', / '  7 = | (A,B) - Q (S,T) Z', A, ' | / ( |(A,B)| n ulp )  ', / '  8 = | I - QQ', A, ' | / ( n ulp )            9 = | I - ZZ', A, ' | / ( n ulp )', / ' 10 = A is in Schur form S', / ' 11 = difference between (alpha,beta) and diagonals', ' of (S,T)', / ' 12 = SDIM is the correct number of ', 'selected eigenvalues', / )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9990 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 1P, D10.3 )

      // End of DDRGES

      }
