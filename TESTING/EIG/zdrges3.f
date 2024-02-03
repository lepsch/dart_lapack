      SUBROUTINE ZDRGES3( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHA, BETA, WORK, LWORK, RWORK, RESULT, BWORK, INFO );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      double             RESULT( 13 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDA, * ), BETA( * ), Q( LDQ, * ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, ILABAD;
      String             SORT;
      int                I, IADD, IINFO, IN, ISORT, J, JC, JR, JSIZE, JTYPE, KNTEIG, MAXWRK, MINWRK, MTYPES, N, N1, NB, NERRS, NMATS, NMAX, NTEST, NTESTT, RSUB, SDIM;
      double             SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV;
      COMPLEX*16         CTEMP, X;
      // ..
      // .. Local Arrays ..
      bool               LASIGN( MAXTYP ), LBSIGN( MAXTYP );
      int                IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      double             RMAGN( 0: 3 );
      // ..
      // .. External Functions ..
      bool               ZLCTES;
      int                ILAENV;
      double             DLAMCH;
      COMPLEX*16         ZLARND;
      // EXTERNAL ZLCTES, ILAENV, DLAMCH, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, XERBLA, ZGET51, ZGET54, ZGGES3, ZLACPY, ZLARFG, ZLASET, ZLATM4, ZUNM2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SIGN
      // ..
      // .. Statement Functions ..
      double             ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) );
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
      // ..
      // .. Executable Statements ..

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
         MINWRK = 3*NMAX*NMAX;
         NB = max( 1, ILAENV( 1, 'ZGEQRF', ' ', NMAX, NMAX, -1, -1 ), ILAENV( 1, 'ZUNMQR', 'LC', NMAX, NMAX, NMAX, -1 ), ILAENV( 1, 'ZUNGQR', ' ', NMAX, NMAX, NMAX, -1 ) );
         MAXWRK = max( NMAX+NMAX*NB, 3*NMAX*NMAX );
         WORK( 1 ) = MAXWRK;
      }

      if (LWORK < MINWRK) INFO = -19;

      if ( INFO != 0 ) {
         xerbla('ZDRGES3', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      ULP = DLAMCH( 'Precision' );
      SAFMIN = DLAMCH( 'Safe minimum' );
      SAFMIN = SAFMIN / ULP;
      SAFMAX = ONE / SAFMIN;
      ULPINV = ONE / ULP;

      // The values RMAGN(2:3) depend on N, see below.

      RMAGN( 0 ) = ZERO;
      RMAGN( 1 ) = ONE;

      // Loop over matrix sizes

      NTESTT = 0;
      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 190
         N = NN( JSIZE );
         N1 = max( 1, N );
         RMAGN( 2 ) = SAFMAX*ULP / DBLE( N1 );
         RMAGN( 3 ) = SAFMIN*ULPINV*DBLE( N1 );

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
               IOLDSD( J ) = ISEED( J );
            } // 20

            // Initialize RESULT

            for (J = 1; J <= 13; J++) { // 30
               RESULT( J ) = ZERO;
            } // 30

            // Generate test matrices A and B

            // Description of control parameters:

            // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
                    // =3 means random.
            // KATYPE: the "type" to be passed to ZLATM4 for computing A.
            // KAZERO: the pattern of zeros on the diagonal for A:
                    // =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
                    // =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
                    // =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
                    // non-zero entries.)
            // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
                    // =2: large, =3: small.
            // LASIGN: true if the diagonal elements of A are to be
                    // multiplied by a random magnitude 1 number.
            // KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN: used to implement KAMAGN and KBMAGN.

            if (MTYPES > MAXTYP) GO TO 110;
            IINFO = 0;
            if ( KCLASS( JTYPE ) < 3 ) {

               // Generate A (w/o rotation)

               if ( ABS( KATYPE( JTYPE ) ) == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) CALL ZLASET( 'Full', N, N, CZERO, CZERO, A, LDA );
               } else {
                  IN = N;
               }
               zlatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) );
               if (IADD > 0 && IADD <= N) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) );

               // Generate B (w/o rotation)

               if ( ABS( KBTYPE( JTYPE ) ) == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) CALL ZLASET( 'Full', N, N, CZERO, CZERO, B, LDA );
               } else {
                  IN = N;
               }
               zlatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) );
               if (IADD != 0 && IADD <= N) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) );

               if ( KCLASS( JTYPE ) == 2 && N > 0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  for (JC = 1; JC <= N - 1; JC++) { // 50
                     for (JR = JC; JR <= N; JR++) { // 40
                        Q( JR, JC ) = ZLARND( 3, ISEED );
                        Z( JR, JC ) = ZLARND( 3, ISEED );
                     } // 40
                     zlarfg(N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK( 2*N+JC ) = SIGN( ONE, DBLE( Q( JC, JC ) ) );
                     Q( JC, JC ) = CONE;
                     zlarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK( 3*N+JC ) = SIGN( ONE, DBLE( Z( JC, JC ) ) );
                     Z( JC, JC ) = CONE;
                  } // 50
                  CTEMP = ZLARND( 3, ISEED );
                  Q( N, N ) = CONE;
                  WORK( N ) = CZERO;
                  WORK( 3*N ) = CTEMP / ABS( CTEMP );
                  CTEMP = ZLARND( 3, ISEED );
                  Z( N, N ) = CONE;
                  WORK( 2*N ) = CZERO;
                  WORK( 4*N ) = CTEMP / ABS( CTEMP );

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 70
                     for (JR = 1; JR <= N; JR++) { // 60
                        A( JR, JC ) = WORK( 2*N+JR )* DCONJG( WORK( 3*N+JC ) )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )* DCONJG( WORK( 3*N+JC ) )* B( JR, JC );
                     } // 60
                  } // 70
                  CALL ZUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL ZUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL ZUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL ZUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 90
                  for (JR = 1; JR <= N; JR++) { // 80
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* ZLARND( 4, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* ZLARND( 4, ISEED );
                  } // 80
               } // 90
            }

            } // 100

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ABS( IINFO );
               return;
            }

            } // 110

            for (I = 1; I <= 13; I++) { // 120
               RESULT( I ) = -ONE;
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

               // Call XLAENV to set the parameters used in ZLAQZ0

               xlaenv(12, 10 );
               xlaenv(13, 12 );
               xlaenv(14, 13 );
               xlaenv(15, 2 );
               xlaenv(17, 10 );

               // Call ZGGES3 to compute H, T, Q, Z, alpha, and beta.

               zlacpy('Full', N, N, A, LDA, S, LDA );
               zlacpy('Full', N, N, B, LDA, T, LDA );
               NTEST = 1 + RSUB + ISORT;
               RESULT( 1+RSUB+ISORT ) = ULPINV;
               zgges3('V', 'V', SORT, ZLCTES, N, S, LDA, T, LDA, SDIM, ALPHA, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, BWORK, IINFO );
               if ( IINFO != 0 && IINFO != N+2 ) {
                  RESULT( 1+RSUB+ISORT ) = ULPINV;
                  WRITE( NOUNIT, FMT = 9999 )'ZGGES3', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  GO TO 160;
               }

               NTEST = 4 + RSUB;

               // Do tests 1--4 (or tests 7--9 when reordering )

               if ( ISORT == 0 ) {
                  zget51(1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK, RWORK, RESULT( 1 ) );
                  zget51(1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RWORK, RESULT( 2 ) );
               } else {
                  zget54(N, A, LDA, B, LDA, S, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 2+RSUB ) );
               }

               zget51(3, N, B, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK, RWORK, RESULT( 3+RSUB ) );
               zget51(3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK, RWORK, RESULT( 4+RSUB ) );

               // Do test 5 and 6 (or Tests 10 and 11 when reordering):
               // check Schur form of A and compare eigenvalues with
               // diagonals.

               NTEST = 6 + RSUB;
               TEMP1 = ZERO;

               for (J = 1; J <= N; J++) { // 130
                  ILABAD = false;
                  TEMP2 = ( ABS1( ALPHA( J )-S( J, J ) ) / max( SAFMIN, ABS1( ALPHA( J ) ), ABS1( S( J, J ) ) )+ABS1( BETA( J )-T( J, J ) ) / max( SAFMIN, ABS1( BETA( J ) ), ABS1( T( J, J ) ) ) ) / ULP;

                  if ( J < N ) {
                     if ( S( J+1, J ) != ZERO ) {
                        ILABAD = true;
                        RESULT( 5+RSUB ) = ULPINV;
                     }
                  }
                  if ( J > 1 ) {
                     if ( S( J, J-1 ) != ZERO ) {
                        ILABAD = true;
                        RESULT( 5+RSUB ) = ULPINV;
                     }
                  }
                  TEMP1 = max( TEMP1, TEMP2 );
                  if ( ILABAD ) {
                     WRITE( NOUNIT, FMT = 9998 )J, N, JTYPE, IOLDSD;
                  }
               } // 130
               RESULT( 6+RSUB ) = TEMP1;

               if ( ISORT >= 1 ) {

                  // Do test 12

                  NTEST = 12;
                  RESULT( 12 ) = ZERO;
                  KNTEIG = 0;
                  for (I = 1; I <= N; I++) { // 140
                     if( ZLCTES( ALPHA( I ), BETA( I ) ) ) KNTEIG = KNTEIG + 1;
                  } // 140
                  if (SDIM != KNTEIG) RESULT( 13 ) = ULPINV;
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
                     WRITE( NOUNIT, FMT = 9997 )'ZGS';

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 );
                     WRITE( NOUNIT, FMT = 9995 );
                     WRITE( NOUNIT, FMT = 9994 )'Unitary';

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 )'unitary', '''', 'transpose', ( '''', J = 1, 8 );

                  }
                  NERRS = NERRS + 1;
                  if ( RESULT( JR ) < 10000.0 ) {
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  } else {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  }
               }
            } // 170

         } // 180
      } // 190

      // Summary

      alasvm('ZGS', NOUNIT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK;

      return;

 9999 FORMAT( ' ZDRGES3: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' );

 9998 FORMAT( ' ZDRGES3: S not in Schur form at eigenvalue ', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

 9997 FORMAT( / 1X, A3, ' -- Complex Generalized Schur from problem ', 'driver' );

 9996 FORMAT( ' Matrix types (see ZDRGES3 for details): ' );

 9995 FORMAT( ' Special Matrices:', 23X, '(J''=transposed Jordan block)', / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ', '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ', 'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I', ')  11=(large*I, small*D)  13=(large*D, large*I)', / '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ', ' 14=(small*D, small*I)', / '  15=(D, reversed D)' );
 9994 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:', / '  16=Transposed Jordan Blocks             19=geometric ', 'alpha, beta=0,1', / '  17=arithm. alpha&beta             ', '      20=arithmetic alpha, beta=0,1', / '  18=clustered ', 'alpha, beta=0,1            21=random alpha, beta=0,1', / ' Large & Small Matrices:', / '  22=(large, small)   ', '23=(small,large)    24=(small,small)    25=(large,large)', / '  26=random O(1) matrices.' );

 9993 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ', 'Q and Z are ', A, ',', / 19X, 'l and r are the appropriate left and right', / 19X, 'eigenvectors, resp., a is alpha, b is beta, and', / 19X, A, ' means ', A, '.)', / ' Without ordering: ', / '  1 = | A - Q S Z', A, ' | / ( |A| n ulp )      2 = | B - Q T Z', A, ' | / ( |B| n ulp )', / '  3 = | I - QQ', A, ' | / ( n ulp )             4 = | I - ZZ', A, ' | / ( n ulp )', / '  5 = A is in Schur form S', / '  6 = difference between (alpha,beta)', ' and diagonals of (S,T)', / ' With ordering: ', / '  7 = | (A,B) - Q (S,T) Z', A, ' | / ( |(A,B)| n ulp )', / '  8 = | I - QQ', A, ' | / ( n ulp )             9 = | I - ZZ', A, ' | / ( n ulp )', / ' 10 = A is in Schur form S', / ' 11 = difference between (alpha,beta) and diagonals', ' of (S,T)', / ' 12 = SDIM is the correct number of ', 'selected eigenvalues', / );
 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 );
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 1P, D10.3 );

      // End of ZDRGES3

      }
