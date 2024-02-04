      void dchksb(NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, SD, SE, U, LDU, WORK, LWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, NWDTHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), NN( * );
      double             A( LDA, * ), RESULT( * ), SD( * ), SE( * ), U( LDU, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, TEN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0 ;
      double             HALF;
      const              HALF = ONE / TWO ;
      int                MAXTYP;
      const              MAXTYP = 15 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KMAX, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACPY, DLASET, DLASUM, DLATMR, DLATMS, DSBT21, DSBTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8,];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 ];
      // ..
      // .. Executable Statements ..

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
         xerbla('DCHKSB', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      // Loop over sizes, types

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 190
         N = NN( JSIZE );
         ANINV = ONE / (max( 1, N )).toDouble();

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

               dlaset('Full', LDA, N, ZERO, ZERO, A, LDA );
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
                     A[K+1, JCOL] = ANORM;
                  } // 80

               } else if ( ITYPE == 4 ) {

                  // Diagonal Matrix, [Eigen]values Specified

                  dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q', A( K+1, 1 ), LDA, WORK( N+1 ), IINFO );

               } else if ( ITYPE == 5 ) {

                  // Symmetric, eigenvalues specified

                  dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO );

               } else if ( ITYPE == 7 ) {

                  // Diagonal, random eigenvalues

                  dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'Q', A( K+1, 1 ), LDA, IDUMMA, IINFO );

               } else if ( ITYPE == 8 ) {

                  // Symmetric, random eigenvalues

                  dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, K, K, ZERO, ANORM, 'Q', A, LDA, IDUMMA, IINFO );

               } else if ( ITYPE == 9 ) {

                  // Positive definite, eigenvalues specified.

                  dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO );

               } else if ( ITYPE == 10 ) {

                  // Positive definite tridiagonal, eigenvalues specified.

                  if (N > 1) K = max( 1, K );
                  dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q', A( K, 1 ), LDA, WORK( N+1 ), IINFO );
                  for (I = 2; I <= N; I++) { // 90
                     TEMP1 = ( A( K, I ) ).abs() / sqrt( ABS( A( K+1, I-1 )*A( K+1, I ) ) );
                     if ( TEMP1 > HALF ) {
                        A[K, I] = HALF*sqrt( ABS( A( K+1, I-1 )*A( K+1, I ) ) );
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

               // Call DSBTRD to compute S and U from upper triangle.

               dlacpy(' ', K+1, N, A, LDA, WORK, LDA );

               NTEST = 1;
               dsbtrd('V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSBTRD(U)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[1] = ULPINV;
                     GO TO 150;
                  }
               }

               // Do tests 1 and 2

               dsbt21('Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 1 ) );

               // Convert A from Upper-Triangle-Only storage to
               // Lower-Triangle-Only storage.

               for (JC = 1; JC <= N; JC++) { // 120
                  for (JR = 0; JR <= min( K, N-JC ); JR++) { // 110
                     A[JR+1, JC] = A( K+1-JR, JC+JR );
                  } // 110
               } // 120
               for (JC = N + 1 - K; JC <= N; JC++) { // 140
                  for (JR = min( K, N-JC ) + 1; JR <= K; JR++) { // 130
                     A[JR+1, JC] = ZERO;
                  } // 130
               } // 140

               // Call DSBTRD to compute S and U from lower triangle

               dlacpy(' ', K+1, N, A, LDA, WORK, LDA );

               NTEST = 3;
               dsbtrd('V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSBTRD(L)', IINFO, N, JTYPE, IOLDSD;
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

               dsbt21('Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 3 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               } // 150
               NTESTT = NTESTT + NTEST;

               // Print out tests which fail.

               for (JR = 1; JR <= NTEST; JR++) { // 160
                  if ( RESULT( JR ) >= THRESH ) {

                     // If this is the first test to fail,
                     // print a header to the data file.

                     if ( NERRS == 0 ) {
                        WRITE( NOUNIT, FMT = 9998 )'DSB';
                        WRITE( NOUNIT, FMT = 9997 );
                        WRITE( NOUNIT, FMT = 9996 );
                        WRITE( NOUNIT, FMT = 9995 )'Symmetric';
                        WRITE( NOUNIT, FMT = 9994 )'orthogonal', '''', 'transpose', ( '''', J = 1, 4 );
                     }
                     NERRS = NERRS + 1;
                     WRITE( NOUNIT, FMT = 9993 )N, K, IOLDSD, JTYPE, JR, RESULT( JR );
                  }
               } // 160

            } // 170
         } // 180
      } // 190

      // Summary

      dlasum('DSB', NOUNIT, NERRS, NTESTT );
      return;

 9999 FORMAT( ' DCHKSB: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

 9998 FORMAT( / 1X, A3, ' -- Real Symmetric Banded Tridiagonal Reduction Routines' );
 9997 FORMAT( ' Matrix types (see DCHKSB for details): ' );

 9996 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.                        ', '  5=Diagonal: clustered entries.', / '  2=Identity matrix.                    ', '  6=Diagonal: large, evenly spaced.', / '  3=Diagonal: evenly spaced entries.    ', '  7=Diagonal: small, evenly spaced.', / '  4=Diagonal: geometr. spaced entries.' );
 9995 FORMAT( ' Dense ', A, ' Banded Matrices:', / '  8=Evenly spaced eigenvals.            ', ' 12=Small, evenly spaced eigenvals.', / '  9=Geometrically spaced eigenvals.     ', ' 13=Matrix with random O(1) entries.', / ' 10=Clustered eigenvalues.              ', ' 14=Matrix with large random entries.', / ' 11=Large, evenly spaced eigenvals.     ', ' 15=Matrix with small random entries.' );

 9994 FORMAT( / ' Tests performed:   (S is Tridiag,  U is ', A, ',', / 20X, A, ' means ', A, '.', / ' UPLO=''U'':', / '  1= | A - U S U', A1, ' | / ( |A| n ulp )     ', '  2= | I - U U', A1, ' | / ( n ulp )', / ' UPLO=''L'':', / '  3= | A - U S U', A1, ' | / ( |A| n ulp )     ', '  4= | I - U U', A1, ' | / ( n ulp )' );
 9993 FORMAT( ' N=', I5, ', K=', I4, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 );
      }