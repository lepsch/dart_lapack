      import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dchksb(final int NSIZES,
          final Array<int> NN, final int NWDTHS,
          final Array<int> KK, final int NTYPES,
          final Array<bool> DOTYPE,
          final Array<int> ISEED, final int THRESH, final int NOUNIT,
          final Matrix<double> A, final int LDA,
          final Array<double> SD,
          final Array<double> SE,
          final Matrix<double> U, final int LDU,
          final Array<double> WORK, final int LWORK,
          final Array<double> RESULT, final Box<int> INFO, ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      // int                INFO.value, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, NWDTHS;
      // double             THRESH;
      // bool               DOTYPE[ * ];
      // int                ISEED[ 4 ], KK[ * ], NN[ * ];
      // double             A[ LDA][ * ], RESULT[ * ], SD[ * ], SE[ * ], U[ LDU][ * ], WORK[ * ];
      // ..

      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0 ;
      const              HALF = ONE / TWO ;
      const              MAXTYP = 15 ;
      bool               BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KMAX, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, ULP, ULPINV, UNFL;
      final                IDUMMA=Array<int>( 1 ), IOLDSD=Array<int>( 4 );
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
      final KTYPE = Array.fromList([ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, //
      ]);
      final KMAGN = Array.fromList([ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 //
      ]);
      final KMODE = Array.fromList([ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 //
      ]);

      // Check for errors

      NTESTT = 0;
      INFO.value = 0;

      // Important constants

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN[ J ] );
         if( NN[ J ] < 0 ) BADNN = true;
      } // 10

      BADNNB = false;
      KMAX = 0;
      for (J = 1; J <= NSIZES; J++) { // 20
         KMAX = max( KMAX, KK[ J ] );
         if( KK[ J ] < 0 ) BADNNB = true;
      } // 20
      KMAX = min( NMAX-1, KMAX );

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO.value = -1;
      } else if ( BADNN ) {
         INFO.value = -2;
      } else if ( NWDTHS < 0 ) {
         INFO.value = -3;
      } else if ( BADNNB ) {
         INFO.value = -4;
      } else if ( NTYPES < 0 ) {
         INFO.value = -5;
      } else if ( LDA < KMAX+1 ) {
         INFO.value = -11;
      } else if ( LDU < NMAX ) {
         INFO.value = -15;
      } else if ( ( max( LDA, NMAX )+1 )*NMAX > LWORK ) {
         INFO.value = -17;
      }

      if ( INFO.value != 0 ) {
         xerbla('DCHKSB', -INFO.value );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

      // More Important constants

      UNFL = dlamch( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      // Loop over sizes, types

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 190
         N = NN[ JSIZE ];
         ANINV = ONE / (max( 1, N )).toDouble();

         for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) { // 180
            K = KK[ JWIDTH ];
            if (K > N) continue;
            K = max( 0, min( N-1, K ) );

            if ( NSIZES != 1 ) {
               MTYPES = min( MAXTYP, NTYPES );
            } else {
               MTYPES = min( MAXTYP+1, NTYPES );
            }

            for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 170
               if( !DOTYPE[ JTYPE ] ) continue;
               NMATS = NMATS + 1;
               NTEST = 0;

               for (J = 1; J <= 4; J++) { // 30
                  IOLDSD[J] = ISEED[ J ];
               } // 30

               // Compute "A".
               // Store as "Upper"; later, we will copy to other format.
               //
               // Control parameters:
               //
               //     KMAGN  KMODE        KTYPE
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

               if (MTYPES > MAXTYP) GOTO100;

               ITYPE = KTYPE[ JTYPE ];
               IMODE = KMODE[ JTYPE ];

               // Compute norm

              //  GOTO( 40, 50, 60 )KMAGN[ JTYPE ];
              switch(KMAGN[ JTYPE ]){

               case 1:
               ANORM = ONE;
               break;

               case 2:
               ANORM = ( RTOVFL*ULP )*ANINV;
               break;

               case 3:
               ANORM = RTUNFL*N*ULPINV;
               break;

               }

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
                     A[K+1][JCOL] = ANORM;
                  } // 80

               } else if ( ITYPE == 4 ) {

                  // Diagonal Matrix, [Eigen]values Specified

                  dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q', A[ K+1][ 1 ], LDA, WORK[ N+1 ], IINFO );

               } else if ( ITYPE == 5 ) {

                  // Symmetric, eigenvalues specified

                  dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK[ N+1 ], IINFO );

               } else if ( ITYPE == 7 ) {

                  // Diagonal, random eigenvalues

                  dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK[ N+1 ], 1, ONE, WORK[ 2*N+1 ], 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'Q', A[ K+1][ 1 ], LDA, IDUMMA, IINFO );

               } else if ( ITYPE == 8 ) {

                  // Symmetric, random eigenvalues

                  dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK[ N+1 ], 1, ONE, WORK[ 2*N+1 ], 1, ONE, 'N', IDUMMA, K, K, ZERO, ANORM, 'Q', A, LDA, IDUMMA, IINFO );

               } else if ( ITYPE == 9 ) {

                  // Positive definite, eigenvalues specified.

                  dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK[ N+1 ], IINFO );

               } else if ( ITYPE == 10 ) {

                  // Positive definite tridiagonal, eigenvalues specified.

                  if (N > 1) K = max( 1, K );
                  dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q', A[ K][ 1 ], LDA, WORK[ N+1 ], IINFO );
                  for (I = 2; I <= N; I++) { // 90
                     TEMP1 = ( A[ K][ I ] ).abs() / sqrt( ABS( A[ K+1][ I-1 ]*A[ K+1][ I ] ) );
                     if ( TEMP1 > HALF ) {
                        A[K][I] = HALF*sqrt( ABS( A[ K+1][ I-1 ]*A[ K+1][ I ] ) );
                     }
                  } // 90

               } else {

                  IINFO = 1;
               }

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
                  INFO.value = ( IINFO ).abs();
                  return;
               }

               } // 100

               // Call DSBTRD to compute S and U from upper triangle.

               dlacpy(' ', K+1, N, A, LDA, WORK, LDA );

               NTEST = 1;
               dsbtrd('V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, WORK[ LDA*N+1 ], IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSBTRD(U)', IINFO, N, JTYPE, IOLDSD;
                  INFO.value = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[1] = ULPINV;
                     GOTO150;
                  }
               }

               // Do tests 1 and 2

               dsbt21('Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT[ 1 ] );

               // Convert A from Upper-Triangle-Only storage to
               // Lower-Triangle-Only storage.

               for (JC = 1; JC <= N; JC++) { // 120
                  for (JR = 0; JR <= min( K, N-JC ); JR++) { // 110
                     A[JR+1, JC] = A[ K+1-JR][ JC+JR ];
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
               dsbtrd('V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, WORK[ LDA*N+1 ], IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSBTRD(L)', IINFO, N, JTYPE, IOLDSD;
                  INFO.value = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[3] = ULPINV;
                     GOTO150;
                  }
               }
               NTEST = 4;

               // Do tests 3 and 4

               dsbt21('Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT[ 3 ] );

               // End of Loop -- Check for RESULT[j] > THRESH

               } // 150
               NTESTT = NTESTT + NTEST;

               // Print out tests which fail.

               for (JR = 1; JR <= NTEST; JR++) { // 160
                  if ( RESULT[ JR ] >= THRESH ) {

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
                     WRITE( NOUNIT, FMT = 9993 )N, K, IOLDSD, JTYPE, JR, RESULT[ JR ];
                  }
               } // 160

            } // 170
         } // 180
      } // 190

      // Summary

      dlasum('DSB', NOUNIT, NERRS, NTESTT );
      return;

 9999 FORMAT( ' DCHKSB: ${} returned INFO.value=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${i5(3, ',')}', I5, ')' );

 9998 FORMAT( / 1X, '${.a3} -- Real Symmetric Banded Tridiagonal Reduction Routines' );
 9997 FORMAT( ' Matrix types (see DCHKSB for details): ' );

 9996 FORMAT('\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.' );
 9995 FORMAT( ' Dense ${} Banded Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.' );

 9994 FORMAT('\n Tests performed:   (S is Tridiag,  U is ${},', / 20X, A, ' means ${}.\n UPLO=''U'':\n  1= | A - U S U${.a1} | / ( |A| n ulp )       2= | I - U U${.a1} | / ( n ulp )\n UPLO=''L'':\n  3= | A - U S U${.a1} | / ( |A| n ulp )       4= | I - U U${.a1} | / ( n ulp )' );
 9993 FORMAT( ' N=${.i5}, K=${.i4}, seed=${i4(4, ',')}', ' type ${.i2}, test(${.i2})=${.g10_3}');
      }
