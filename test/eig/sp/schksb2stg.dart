      void schksb2stg(NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, SD, SE, D1, D2, D3, U, LDU, WORK, LWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, NWDTHS;
      double               THRESH;
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), NN( * );
      double               A( LDA, * ), RESULT( * ), SD( * ), SE( * ), D1( * ), D2( * ), D3( * ), U( LDU, * ), WORK( * );
      // ..

      double               ZERO, ONE, TWO, TEN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0 ;
      double               HALF;
      const              HALF = ONE / TWO ;
      int                MAXTYP;
      const              MAXTYP = 15 ;
      bool               BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KMAX, LH, LW, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double               ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, TEMP4, ULP, ULPINV, UNFL;
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLASET, SLASUM, SLATMR, SLATMS, SSBT21, SSBTRD, XERBLA, SSYTRD_SB2ST, SSTEQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8,];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 ];

      // Check for errors

      NTESTT = 0;
      INFO = 0;

      // Important constants

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
      } // 10

      BADNNB = false;
      KMAX = 0;
      for (J = 1; J <= NSIZES; J++) { // 20
         KMAX = max( KMAX, KK( J ) );
         if( KK( J ) < 0 ) BADNNB = true;
      } // 20
      KMAX = min( NMAX-1, KMAX );

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADNN ) {
         INFO = -2;
      } else if ( NWDTHS < 0 ) {
         INFO = -3;
      } else if ( BADNNB ) {
         INFO = -4;
      } else if ( NTYPES < 0 ) {
         INFO = -5;
      } else if ( LDA < KMAX+1 ) {
         INFO = -11;
      } else if ( LDU < NMAX ) {
         INFO = -15;
      } else if ( ( max( LDA, NMAX )+1 )*NMAX > LWORK ) {
         INFO = -17;
      }

      if ( INFO != 0 ) {
         xerbla('SCHKSB2STG', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      // Loop over sizes, types

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 190
         N = NN( JSIZE );
         ANINV = ONE / REAL( max( 1, N ) );

         for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) { // 180
            K = KK( JWIDTH );
            if (K > N) GO TO 180;
            K = max( 0, min( N-1, K ) );

            if ( NSIZES != 1 ) {
               MTYPES = min( MAXTYP, NTYPES );
            } else {
               MTYPES = min( MAXTYP+1, NTYPES );
            }

            for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 170
               if( !DOTYPE( JTYPE ) ) GO TO 170;
               NMATS = NMATS + 1;
               NTEST = 0;

               for (J = 1; J <= 4; J++) { // 30
                  IOLDSD[J] = ISEED( J );
               } // 30

               // Compute "A".
               // Store as "Upper"; later, we will copy to other format.

               // Control parameters:

                   // KMAGN  KMODE        KTYPE
               // =1  O(1)   clustered 1  zero
               // =2  large  clustered 2  identity
               // =3  small  exponential  (none)
               // =4         arithmetic   diagonal, (w/ eigenvalues)
               // =5         random log   symmetric, w/ eigenvalues
               // =6         random       (none)
               // =7                      random diagonal
               // =8                      random symmetric
               // =9                      positive definite
               // =10                     diagonally dominant tridiagonal

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

               slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
               IINFO = 0;
               if ( JTYPE <= 15 ) {
                  COND = ULPINV;
               } else {
                  COND = ULPINV*ANINV / TEN;
               }

               // Special Matrices -- Identity & Jordan block

                  // Zero

               if ( ITYPE == 1 ) {
                  IINFO = 0;

               } else if ( ITYPE == 2 ) {

                  // Identity

                  for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                     A[K+1][JCOL] = ANORM;
                  } // 80

               } else if ( ITYPE == 4 ) {

                  // Diagonal Matrix, [Eigen]values Specified

                  slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q', A( K+1, 1 ), LDA, WORK( N+1 ), IINFO );

               } else if ( ITYPE == 5 ) {

                  // Symmetric, eigenvalues specified

                  slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO );

               } else if ( ITYPE == 7 ) {

                  // Diagonal, random eigenvalues

                  slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'Q', A( K+1, 1 ), LDA, IDUMMA, IINFO );

               } else if ( ITYPE == 8 ) {

                  // Symmetric, random eigenvalues

                  slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, K, K, ZERO, ANORM, 'Q', A, LDA, IDUMMA, IINFO );

               } else if ( ITYPE == 9 ) {

                  // Positive definite, eigenvalues specified.

                  slatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO );

               } else if ( ITYPE == 10 ) {

                  // Positive definite tridiagonal, eigenvalues specified.

                  if (N > 1) K = max( 1, K );
                  slatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q', A( K, 1 ), LDA, WORK( N+1 ), IINFO );
                  for (I = 2; I <= N; I++) { // 90
                     TEMP1 = ( A( K, I ) ).abs() / sqrt( ABS( A( K+1, I-1 )*A( K+1, I ) ) );
                     if ( TEMP1 > HALF ) {
                        A[K][I] = HALF*sqrt( ABS( A( K+1, I-1 )*A( K+1, I ) ) );
                     }
                  } // 90

               } else {

                  IINFO = 1;
               }

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               } // 100

               // Call SSBTRD to compute S and U from upper triangle.

               slacpy(' ', K+1, N, A, LDA, WORK, LDA );

               NTEST = 1;
               ssbtrd('V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBTRD(U)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[1] = ULPINV;
                     GO TO 150;
                  }
               }

               // Do tests 1 and 2

               ssbt21('Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 1 ) );

               // Before converting A into lower for SSBTRD, run SSYTRD_SB2ST
               // otherwise matrix A will be converted to lower and then need
               // to be converted back to upper in order to run the upper case
               // ofSSYTRD_SB2ST

               // Compute D1 the eigenvalues resulting from the tridiagonal
               // form using the SSBTRD and used as reference to compare
               // with the SSYTRD_SB2ST routine

               // Compute D1 from the SSBTRD and used as reference for the
               // SSYTRD_SB2ST

               scopy(N, SD, 1, D1, 1 );
               if (N > 0) scopy( N-1, SE, 1, WORK, 1 );

               ssteqr('N', N, D1, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[5] = ULPINV;
                     GO TO 150;
                  }
               }

               // SSYTRD_SB2ST Upper case is used to compute D2.
               // Note to set SD and SE to zero to be sure not reusing
               // the one from above. Compare it with D1 computed
               // using the SSBTRD.

               slaset('Full', N, 1, ZERO, ZERO, SD, N );
               slaset('Full', N, 1, ZERO, ZERO, SE, N );
               slacpy(' ', K+1, N, A, LDA, U, LDU );
               LH = max(1, 4*N);
               LW = LWORK - LH;
               ssytrd_sb2st('N', 'N', "U", N, K, U, LDU, SD, SE, WORK, LH, WORK( LH+1 ), LW, IINFO );

               // Compute D2 from the SSYTRD_SB2ST Upper case

               scopy(N, SD, 1, D2, 1 );
               if (N > 0) scopy( N-1, SE, 1, WORK, 1 );

               ssteqr('N', N, D2, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[5] = ULPINV;
                     GO TO 150;
                  }
               }

               // Convert A from Upper-Triangle-Only storage to
               // Lower-Triangle-Only storage.

               for (JC = 1; JC <= N; JC++) { // 120
                  for (JR = 0; JR <= min( K, N-JC ); JR++) { // 110
                     A[JR+1][JC] = A( K+1-JR, JC+JR );
                  } // 110
               } // 120
               for (JC = N + 1 - K; JC <= N; JC++) { // 140
                  for (JR = min( K, N-JC ) + 1; JR <= K; JR++) { // 130
                     A[JR+1][JC] = ZERO;
                  } // 130
               } // 140

               // Call SSBTRD to compute S and U from lower triangle

               slacpy(' ', K+1, N, A, LDA, WORK, LDA );

               NTEST = 3;
               ssbtrd('V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBTRD(L)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[3] = ULPINV;
                     GO TO 150;
                  }
               }
               NTEST = 4;

               // Do tests 3 and 4

               ssbt21('Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 3 ) );

               // SSYTRD_SB2ST Lower case is used to compute D3.
               // Note to set SD and SE to zero to be sure not reusing
               // the one from above. Compare it with D1 computed
               // using the SSBTRD.

               slaset('Full', N, 1, ZERO, ZERO, SD, N );
               slaset('Full', N, 1, ZERO, ZERO, SE, N );
               slacpy(' ', K+1, N, A, LDA, U, LDU );
               LH = max(1, 4*N);
               LW = LWORK - LH;
               ssytrd_sb2st('N', 'N', "L", N, K, U, LDU, SD, SE, WORK, LH, WORK( LH+1 ), LW, IINFO );

               // Compute D3 from the 2-stage Upper case

               scopy(N, SD, 1, D3, 1 );
               if (N > 0) scopy( N-1, SE, 1, WORK, 1 );

               ssteqr('N', N, D3, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[6] = ULPINV;
                     GO TO 150;
                  }
               }


               // Do Tests 3 and 4 which are similar to 11 and 12 but with the
               // D1 computed using the standard 1-stage reduction as reference

               NTEST = 6;
               TEMP1 = ZERO;
               TEMP2 = ZERO;
               TEMP3 = ZERO;
               TEMP4 = ZERO;

               for (J = 1; J <= N; J++) { // 151
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ).abs() );
                  TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
                  TEMP3 = max( TEMP3, ( D1( J ) ).abs(), ( D3( J ) ).abs() );
                  TEMP4 = max( TEMP4, ABS( D1( J )-D3( J ) ) );
               } // 151

               RESULT[5] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );
               RESULT[6] = TEMP4 / max( UNFL, ULP*max( TEMP3, TEMP4 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               } // 150
               NTESTT = NTESTT + NTEST;

               // Print out tests which fail.

               for (JR = 1; JR <= NTEST; JR++) { // 160
                  if ( RESULT( JR ) >= THRESH ) {

                     // If this is the first test to fail,
                     // print a header to the data file.

                     if ( NERRS == 0 ) {
                        WRITE( NOUNIT, FMT = 9998 )'SSB';
                        WRITE( NOUNIT, FMT = 9997 );
                        WRITE( NOUNIT, FMT = 9996 );
                        WRITE( NOUNIT, FMT = 9995 )'Symmetric';
                        WRITE( NOUNIT, FMT = 9994 )'orthogonal', '''', 'transpose', ( '''', J = 1, 6 );
                     }
                     NERRS = NERRS + 1;
                     WRITE( NOUNIT, FMT = 9993 )N, K, IOLDSD, JTYPE, JR, RESULT( JR );
                  }
               } // 160

            } // 170
         } // 180
      } // 190

      // Summary

      slasum('SSB', NOUNIT, NERRS, NTESTT );
      return;

 9999 FORMAT( ' SCHKSB2STG: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${i5(3, ',')}', I5, ')' );

 9998 FORMAT('\n ${.a3} -- Real Symmetric Banded Tridiagonal Reduction Routines' );
 9997 FORMAT( ' Matrix types (see SCHKSB2STG for details): ' );

 9996 FORMAT('\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.' );
 9995 FORMAT( ' Dense ${} Banded Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.' );

 9994 FORMAT('\n Tests performed:   (S is Tridiag,  U is ${},', / 20X, A, ' means ${}.\n UPLO=''U'':\n  1= | A - U S U${.a1} | / ( |A| n ulp )       2= | I - U U${.a1} | / ( n ulp )\n UPLO=''L'':\n  3= | A - U S U${.a1} | / ( |A| n ulp )       4= | I - U U${.a1} | / ( n ulp )' / ' Eig check:\n  5= | D1 - D2 | / ( |D1| ulp )           6= | D1 - D3 | / ( |D1| ulp )          ' );
 9993 FORMAT( ' N=${.i5}, K=${.i4}, seed=${i4(4, ',')}', ' type ${.i2}, test(${.i2})=${.g10_3}');
      }
