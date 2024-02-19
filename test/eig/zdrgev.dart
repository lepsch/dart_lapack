import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zggev.dart';
import 'package:lapack/src/zunm2r.dart';

import '../matgen/zlarnd.dart';
import 'alasvm.dart';
import 'zget52.dart';
import 'zlatm4.dart';

void zdrgev(final int NSIZES, final Array<int> NN_, final int NTYPES, final Array<bool> DOTYPE_, final Array<int> ISEED_, final double THRESH,
final int NOUNIT, final Matrix<Complex> A_, final int LDA, final Matrix<Complex> B_, final Matrix<Complex> S_, final Matrix<Complex> T_, final Matrix<Complex> Q_, final int LDQ,
final Matrix<Complex> Z_, final Matrix<Complex> QE_, final int LDQE, final Array<Complex> ALPHA_, final Array<Complex> BETA_, final Array<Complex> ALPHA1_, final Array<Complex> BETA1_,
final Array<Complex> WORK_, final int LWORK, final Array<double> RWORK_, final Array<double> RESULT_, final Box<int> INFO,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.dim();
  final ISEED = ISEED_.dim(4);
  final NN = NN_.dim();
  final A = A_.dim(LDA);
  final B = B_.dim(LDA);
  final S = S_.dim(LDA);
  final T = T_.dim(LDA);
  final Q = Q_.dim(LDQ);
  final QE = QE_.dim(LDQE);
  final Z = Z_.dim(LDQ);
final ALPHA=ALPHA_.dim();
final BETA=BETA_.dim();
final ALPHA1=ALPHA1_.dim();
final BETA1=BETA1_.dim();
  final RESULT = RESULT_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

      const              ZERO = 0.0, ONE = 1.0 ;
      const              MAXTYP = 26 ;
      bool               BADNN;
      int                I, IADD, IERR, IN, J, JC, JR, JSIZE, JTYPE, MAXWRK, MINWRK, MTYPES, N=0, N1, NB, NERRS, NMATS, NMAX, NTESTT;
      double             SAFMAX, SAFMIN, ULP, ULPINV;
      Complex         CTEMP;
      final                IOLDSD=Array<int>( 4 );
      final             RMAGN=Array<double>( 4)(1,offset: 1);
      final KCLASS = Array.fromList([
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, //
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,]);
      final KZ1 = Array.fromList([ 0, 1, 2, 1, 3, 3 ]);
      final KZ2 = Array.fromList([ 0, 0, 1, 2, 1, 1 ]);
      final KADD = Array.fromList([ 0, 0, 0, 0, 3, 2 ]);
      final KATYPE = Array.fromList([
        0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, //
        4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0 ]);
      final KBTYPE = Array.fromList([
        0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, //
        1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0 ]);
      final KAZERO = Array.fromList([
        1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, //
        2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1 ]);
      final KBZERO = Array.fromList([
        1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, //
        1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1 ]);
      final KAMAGN = Array.fromList([
        1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, //
        3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1 ]);
      final KBMAGN = Array.fromList([
        1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, //
        3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1 ]);
      final KTRIAN = Array.fromList([
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
        0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,]);
      final LASIGN = Array.fromList([
        false, false, false, false, false, false, true , false , true, true, //
        false, false, true, true, true, false , true , false, false, false, //
        true, true, true, true, true, false ]);
      final LBSIGN = Array.fromList([
        false, false, false, false, false, false, false, true , false, false, //
        true, true, false, false, true , false , true , false, false, false, //
        false, false, false, false, false, false,]);

      // Check for errors

      INFO.value = 0;

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN[J] );
         if( NN[J] < 0 ) BADNN = true;
      } // 10

      if ( NSIZES < 0 ) {
         INFO.value = -1;
      } else if ( BADNN ) {
         INFO.value = -2;
      } else if ( NTYPES < 0 ) {
         INFO.value = -3;
      } else if ( THRESH < ZERO ) {
         INFO.value = -6;
      } else if ( LDA <= 1 || LDA < NMAX ) {
         INFO.value = -9;
      } else if ( LDQ <= 1 || LDQ < NMAX ) {
         INFO.value = -14;
      } else if ( LDQE <= 1 || LDQE < NMAX ) {
         INFO.value = -17;
      }

      // Compute workspace
      //  (Note: Comments in the code beginning "Workspace:" describe the
      //   minimal amount of workspace needed at that point in the code,
      //   as well as the preferred amount for good performance.
      //   NB refers to the optimal block size for the immediately
      //   following subroutine, as returned by ILAENV.

      MINWRK = 1;
      if ( INFO.value == 0 && LWORK >= 1 ) {
         MINWRK = NMAX*( NMAX+1 );
         NB = max( max(1, ilaenv( 1, 'ZGEQRF', ' ', NMAX, NMAX, -1, -1 )), max(ilaenv( 1, 'ZUNMQR', 'LC', NMAX, NMAX, NMAX, -1 ), ilaenv( 1, 'ZUNGQR', ' ', NMAX, NMAX, NMAX, -1 )) );
         MAXWRK = max( 2*NMAX, max(NMAX*( NB+1 ), NMAX*( NMAX+1 )) );
         WORK[1] = MAXWRK.toComplex();
      }

      if (LWORK < MINWRK) INFO.value = -23;

      if ( INFO.value != 0 ) {
         xerbla('ZDRGEV', -INFO.value );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      ULP = dlamch( 'Precision' );
      SAFMIN = dlamch( 'Safe minimum' );
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
         N = NN[ JSIZE ];
         N1 = max( 1, N );
         RMAGN[2] = SAFMAX*ULP / N1.toDouble();
         RMAGN[3] = SAFMIN*ULPINV*N1;

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 210
            if( !DOTYPE[ JTYPE ] ) continue;
            NMATS = NMATS + 1;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED[J];
            } // 20

            // Generate test matrices A and B

            // Description of control parameters:

            // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
            //         =3 means random.
            // KATYPE: the "type" to be passed to ZLATM4 for computing A.
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
            if ( KCLASS[ JTYPE ] < 3 ) {

               // Generate A (w/o rotation)

               if ( ( KATYPE[ JTYPE ] ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) ~/ 2 ) + 1;
                  if (IN != N) zlaset( 'Full', N, N, Complex.zero, Complex.zero, A, LDA );
               } else {
                  IN = N;
               }
               zlatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) );
               if (IADD > 0 && IADD <= N) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) );

               // Generate B (w/o rotation)

               if ( ( KBTYPE[ JTYPE ] ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) ~/ 2 ) + 1;
                  if (IN != N) zlaset( 'Full', N, N, Complex.zero, Complex.zero, B, LDA );
               } else {
                  IN = N;
               }
               zlatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) );
               if (IADD != 0 && IADD <= N) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) );

               if ( KCLASS[ JTYPE ] == 2 && N > 0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  for (JC = 1; JC <= N - 1; JC++) { // 40
                     for (JR = JC; JR <= N; JR++) { // 30
                        Q[JR][JC] = zlarnd( 3, ISEED );
                        Z[JR][JC] = zlarnd( 3, ISEED );
                     } // 30
                     zlarfg(N+1-JC, Q[JC][ JC] , Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK[2*N+JC] = sign( ONE, (Q[JC][ JC] ).toDouble() );
                     Q[JC][JC] = Complex.one;
                     zlarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK[3*N+JC] = sign( ONE, (Z[JC][ JC] ).toDouble() );
                     Z[JC][JC] = Complex.one;
                  } // 40
                  CTEMP = zlarnd( 3, ISEED );
                  Q[N][N] = Complex.one;
                  WORK[N] = Complex.zero;
                  WORK[3*N] = CTEMP / ( CTEMP ).abs();
                  CTEMP = zlarnd( 3, ISEED );
                  Z[N][N] = Complex.one;
                  WORK[2*N] = Complex.zero;
                  WORK[4*N] = CTEMP / ( CTEMP ).abs();

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 60
                     for (JR = 1; JR <= N; JR++) { // 50
                        A[JR][JC] = WORK( 2*N+JR )* dconjg( WORK( 3*N+JC ) )* A( JR, JC )                         ;
                        B( JR, JC ) = WORK( 2*N+JR )* dconjg( WORK( 3*N+JC ) )* B( JR, JC );
                     } // 50
                  zunm2r( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IERR )                   ;
                  } // 60
                  if( IERR != 0 ) GO TO 90;
                  zunm2r( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IERR )                   ;
                  if( IERR != 0 ) GO TO 90;
                  zunm2r( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IERR )                   ;
                  if( IERR != 0 ) GO TO 90;
                  zunm2r( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IERR )                   ;
                  if( IERR != 0 ) GO TO 90;
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 80
                  for (JR = 1; JR <= N; JR++) { // 70
                     A[JR][JC] = RMAGN( KAMAGN[ JTYPE ] )* zlarnd( 4, ISEED )                      ;
                     B[JR][ JC]  = RMAGN( KBMAGN[ JTYPE ] )* zlarnd( 4, ISEED );
                  } // 70
               } // 80
            }

            // } // 90

            if ( IERR != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IERR, N, JTYPE, IOLDSD;
               INFO.value = ( IERR ).abs();
               return;
            }

            // } // 100

            for (I = 1; I <= 7; I++) { // 110
               RESULT[I] = -ONE;
            } // 110

            // Call ZGGEV to compute eigenvalues and eigenvectors.

            zlacpy(' ', N, N, A, LDA, S, LDA );
            zlacpy(' ', N, N, B, LDA, T, LDA );
            zggev('V', 'V', N, S, LDA, T, LDA, ALPHA, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'ZGGEV1', IERR, N, JTYPE, IOLDSD;
               INFO.value = ( IERR ).abs();
               GO TO 190;
            }

            // Do the tests (1) and (2)

            zget52( true , N, A, LDA, B, LDA, Q, LDQ, ALPHA, BETA, WORK, RWORK, RESULT( 1 ) );
            if ( RESULT( 2 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZGGEV1', RESULT( 2 ), N, JTYPE, IOLDSD;
            }

            // Do the tests (3) and (4)

            zget52( false , N, A, LDA, B, LDA, Z, LDQ, ALPHA, BETA, WORK, RWORK, RESULT( 3 ) );
            if ( RESULT( 4 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZGGEV1', RESULT( 4 ), N, JTYPE, IOLDSD;
            }

            // Do test (5)

            zlacpy(' ', N, N, A, LDA, S, LDA );
            zlacpy(' ', N, N, B, LDA, T, LDA );
            zggev('N', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'ZGGEV2', IERR, N, JTYPE, IOLDSD;
               INFO.value = ( IERR ).abs();
               GO TO 190;
            }

            for (J = 1; J <= N; J++) { // 120
               if( ALPHA[J] != ALPHA1[J] || BETA[J] != BETA1[J] )RESULT( 5 ) = ULPINV;
            } // 120

            // Do test (6): Compute eigenvalues and left eigenvectors,
            // and test them

            zlacpy(' ', N, N, A, LDA, S, LDA );
            zlacpy(' ', N, N, B, LDA, T, LDA );
            zggev('V', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, QE, LDQE, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'ZGGEV3', IERR, N, JTYPE, IOLDSD;
               INFO.value = ( IERR ).abs();
               GO TO 190;
            }

            for (J = 1; J <= N; J++) { // 130
               if( ALPHA[J] != ALPHA1[J] || BETA[J] != BETA1[J] )RESULT( 6 ) = ULPINV;
            } // 130

            for (J = 1; J <= N; J++) { // 150
               for (JC = 1; JC <= N; JC++) { // 140
                  if( Q( J, JC ) != QE( J, JC ) ) RESULT( 6 ) = ULPINV;
               } // 140
            } // 150

            // Do test (7): Compute eigenvalues and right eigenvectors,
            // and test them

            zlacpy(' ', N, N, A, LDA, S, LDA );
            zlacpy(' ', N, N, B, LDA, T, LDA );
            zggev('N', 'V', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, QE, LDQE, WORK, LWORK, RWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT[1] = ULPINV;
               WRITE( NOUNIT, FMT = 9999 )'ZGGEV4', IERR, N, JTYPE, IOLDSD;
               INFO.value = ( IERR ).abs();
               GO TO 190;
            }

            for (J = 1; J <= N; J++) { // 160
               if( ALPHA[J] != ALPHA1[J] || BETA[J] != BETA1[J] )RESULT[7] = ULPINV;
            } // 160

            for (J = 1; J <= N; J++) { // 180
               for (JC = 1; JC <= N; JC++) { // 170
                  if( Z( J, JC ) != QE( J, JC ) ) RESULT[7] = ULPINV;
               } // 170
            } // 180

            // End of Loop -- Check for RESULT(j) > THRESH

            // } // 190

            NTESTT = NTESTT + 7;

            // Print out tests which fail.

            for (JR = 1; JR <= 7; JR++) { // 200
               if ( RESULT[ JR ] >= THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9997 )'ZGV';

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

      alasvm('ZGV', NOUNIT, NERRS, NTESTT, 0 );

      WORK[1] = MAXWRK;

      return;

 9999 FORMAT( ' ZDRGEV: ${} returned INFO.value=${.i6}.\n${' ' * 3}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );

 9998 FORMAT( ' ZDRGEV: ${} Eigenvectors from ${} incorrectly normalized.\n Bits of error=${.g10_3},${' ' * 3}N=${.i4}, JTYPE=${.i3}, ISEED=(${i4(3, ',')}', I5, ')' );

 9997 FORMAT('\n ${.a3} -- Complex Generalized eigenvalue problem driver' );

 9996 FORMAT( ' Matrix types (see ZDRGEV for details): ' );

 9995 FORMAT( ' Special Matrices:${' ' * 23}(J''=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  6=(diag(J'',I), diag(I,J''))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)' );
 9994 FORMAT( ' Matrices Rotated by Random ${} Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.' );

 9993 FORMAT('\n Tests performed:    \n 1 = max | ( b A - a B )''*l | / const.,\n 2 = | |VR(i)| - 1 | / ulp,\n 3 = max | ( b A - a B )*r | / const.\n 4 = | |VL(i)| - 1 | / ulp,\n 5 = 0 if W same no matter if r or l computed,\n 6 = 0 if l same no matter if l computed,\n 7 = 0 if r same no matter if r computed,/n ');
 9992 FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i2} is${.f8_2}');
 9991 FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i2} is${( * 10).d10_3}');
      }
