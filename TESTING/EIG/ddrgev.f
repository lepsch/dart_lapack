      SUBROUTINE DDRGEV( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, QE, LDQE, ALPHAR, ALPHAI, BETA, ALPHR1, ALPHI1, BETA1, WORK, LWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDQ, LDQE, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), ALPHI1( * ), ALPHR1( * ), B( LDA, * ), BETA( * ), BETA1( * ), Q( LDQ, * ), QE( LDQE, * ), RESULT( * ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      int                I, IADD, IERR, IN, J, JC, JR, JSIZE, JTYPE, MAXWRK, MINWRK, MTYPES, N, N1, NERRS, NMATS, NMAX, NTESTT;
      double             SAFMAX, SAFMIN, ULP, ULPINV;
      // ..
      // .. Local Arrays ..
      int                IASIGN( MAXTYP ), IBSIGN( MAXTYP ), IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      double             RMAGN( 0: 3 );
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, DLARND;
      // EXTERNAL ILAENV, DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, DGET52, DGGEV, DLACPY, DLARFG, DLASET, DLATM4, DORM2R, XERBLA
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
         IF( NN( J ).LT.0 ) BADNN = true;
      } // 10

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( THRESH.LT.ZERO ) {
         INFO = -6
      } else if ( LDA.LE.1 || LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDQ.LE.1 || LDQ.LT.NMAX ) {
         INFO = -14
      } else if ( LDQE.LE.1 || LDQE.LT.NMAX ) {
         INFO = -17
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      MINWRK = 1
      if ( INFO == 0 && LWORK.GE.1 ) {
         MINWRK = MAX( 1, 8*NMAX, NMAX*( NMAX+1 ) )
         MAXWRK = 7*NMAX + NMAX*ILAENV( 1, 'DGEQRF', ' ', NMAX, 1, NMAX, 0 )
         MAXWRK = MAX( MAXWRK, NMAX*( NMAX+1 ) )
         WORK( 1 ) = MAXWRK
      }

      if (LWORK.LT.MINWRK) INFO = -25;

      if ( INFO != 0 ) {
         xerbla('DDRGEV', -INFO );
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

      // Loop over sizes, types

      NTESTT = 0
      NERRS = 0
      NMATS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 220
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / DBLE( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*N1

         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 210
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 210
            NMATS = NMATS + 1

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD( J ) = ISEED( J )
            } // 20

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

            if (MTYPES.GT.MAXTYP) GO TO 100;
            IERR = 0
            if ( KCLASS( JTYPE ).LT.3 ) {

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

                  for (JC = 1; JC <= N - 1; JC++) { // 40
                     for (JR = JC; JR <= N; JR++) { // 30
                        Q( JR, JC ) = DLARND( 3, ISEED )
                        Z( JR, JC ) = DLARND( 3, ISEED )
                     } // 30
                     dlarfg(N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK( 2*N+JC ) = SIGN( ONE, Q( JC, JC ) )
                     Q( JC, JC ) = ONE
                     dlarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK( 3*N+JC ) = SIGN( ONE, Z( JC, JC ) )
                     Z( JC, JC ) = ONE
                  } // 40
                  Q( N, N ) = ONE
                  WORK( N ) = ZERO
                  WORK( 3*N ) = SIGN( ONE, DLARND( 2, ISEED ) )
                  Z( N, N ) = ONE
                  WORK( 2*N ) = ZERO
                  WORK( 4*N ) = SIGN( ONE, DLARND( 2, ISEED ) )

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 60
                     for (JR = 1; JR <= N; JR++) { // 50
                        A( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* B( JR, JC )
                     } // 50
                  } // 60
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
                  CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90;
                  CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR != 0 ) GO TO 90
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 80
                  for (JR = 1; JR <= N; JR++) { // 70
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* DLARND( 2, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* DLARND( 2, ISEED )
                  } // 70
               } // 80
            }

            } // 90

            if ( IERR != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               RETURN
            }

            } // 100

            for (I = 1; I <= 7; I++) { // 110
               RESULT( I ) = -ONE
            } // 110

            // Call DGGEV to compute eigenvalues and eigenvectors.

            dlacpy(' ', N, N, A, LDA, S, LDA );
            dlacpy(' ', N, N, B, LDA, T, LDA );
            dggev('V', 'V', N, S, LDA, T, LDA, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV1', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            // Do the tests (1) and (2)

            dget52( true , N, A, LDA, B, LDA, Q, LDQ, ALPHAR, ALPHAI, BETA, WORK, RESULT( 1 ) );
            if ( RESULT( 2 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'DGGEV1', RESULT( 2 ), N, JTYPE, IOLDSD
            }

            // Do the tests (3) and (4)

            dget52( false , N, A, LDA, B, LDA, Z, LDQ, ALPHAR, ALPHAI, BETA, WORK, RESULT( 3 ) );
            if ( RESULT( 4 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'DGGEV1', RESULT( 4 ), N, JTYPE, IOLDSD
            }

            // Do the test (5)

            dlacpy(' ', N, N, A, LDA, S, LDA );
            dlacpy(' ', N, N, B, LDA, T, LDA );
            dggev('N', 'N', N, S, LDA, T, LDA, ALPHR1, ALPHI1, BETA1, Q, LDQ, Z, LDQ, WORK, LWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV2', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            for (J = 1; J <= N; J++) { // 120
               IF( ALPHAR( J ) != ALPHR1( J ) || ALPHAI( J ) != ALPHI1( J ) || BETA( J ) != BETA1( J ) )RESULT( 5 ) = ULPINV
            } // 120

            // Do the test (6): Compute eigenvalues and left eigenvectors,
            // and test them

            dlacpy(' ', N, N, A, LDA, S, LDA );
            dlacpy(' ', N, N, B, LDA, T, LDA );
            dggev('V', 'N', N, S, LDA, T, LDA, ALPHR1, ALPHI1, BETA1, QE, LDQE, Z, LDQ, WORK, LWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV3', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            for (J = 1; J <= N; J++) { // 130
               IF( ALPHAR( J ) != ALPHR1( J ) || ALPHAI( J ) != ALPHI1( J ) || BETA( J ) != BETA1( J ) )RESULT( 6 ) = ULPINV
            } // 130

            for (J = 1; J <= N; J++) { // 150
               for (JC = 1; JC <= N; JC++) { // 140
                  IF( Q( J, JC ) != QE( J, JC ) ) RESULT( 6 ) = ULPINV
               } // 140
            } // 150

            // DO the test (7): Compute eigenvalues and right eigenvectors,
            // and test them

            dlacpy(' ', N, N, A, LDA, S, LDA );
            dlacpy(' ', N, N, B, LDA, T, LDA );
            dggev('N', 'V', N, S, LDA, T, LDA, ALPHR1, ALPHI1, BETA1, Q, LDQ, QE, LDQE, WORK, LWORK, IERR );
            if ( IERR != 0 && IERR != N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV4', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            for (J = 1; J <= N; J++) { // 160
               IF( ALPHAR( J ) != ALPHR1( J ) || ALPHAI( J ) != ALPHI1( J ) || BETA( J ) != BETA1( J ) )RESULT( 7 ) = ULPINV
            } // 160

            for (J = 1; J <= N; J++) { // 180
               for (JC = 1; JC <= N; JC++) { // 170
                  IF( Z( J, JC ) != QE( J, JC ) ) RESULT( 7 ) = ULPINV
               } // 170
            } // 180

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 190

            NTESTT = NTESTT + 7

            // Print out tests which fail.

            for (JR = 1; JR <= 7; JR++) { // 200
               if ( RESULT( JR ).GE.THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9997 )'DGV'

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )'Orthogonal'

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 )

                  }
                  NERRS = NERRS + 1
                  if ( RESULT( JR ).LT.10000.0D0 ) {
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  } else {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  }
               }
            } // 200

         } // 210
      } // 220

      // Summary

      alasvm('DGV', NOUNIT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK

      RETURN

 9999 FORMAT( ' DDRGEV: ', A, ' returned INFO=', I6, '.', / 3X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' )

 9998 FORMAT( ' DDRGEV: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 3X, 'N=', I4, ', JTYPE=', I3, ', ISEED=(', 4( I4, ',' ), I5, ')' )

 9997 FORMAT( / 1X, A3, ' -- Real Generalized eigenvalue problem driver' )

 9996 FORMAT( ' Matrix types (see DDRGEV for details): ' )

 9995 FORMAT( ' Special Matrices:', 23X, '(J''=transposed Jordan block)', / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ', '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ', 'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I', ')  11=(large*I, small*D)  13=(large*D, large*I)', / '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ', ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9994 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:', / '  16=Transposed Jordan Blocks             19=geometric ', 'alpha, beta=0,1', / '  17=arithm. alpha&beta             ', '      20=arithmetic alpha, beta=0,1', / '  18=clustered ', 'alpha, beta=0,1            21=random alpha, beta=0,1', / ' Large & Small Matrices:', / '  22=(large, small)   ', '23=(small,large)    24=(small,small)    25=(large,large)', / '  26=random O(1) matrices.' )

 9993 FORMAT( / ' Tests performed:    ', / ' 1 = max | ( b A - a B )''*l | / const.,', / ' 2 = | |VR(i)| - 1 | / ulp,', / ' 3 = max | ( b A - a B )*r | / const.', / ' 4 = | |VL(i)| - 1 | / ulp,', / ' 5 = 0 if W same no matter if r or l computed,', / ' 6 = 0 if l same no matter if l computed,', / ' 7 = 0 if r same no matter if r computed,', / 1X )
 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 1P, D10.3 )

      // End of DDRGEV

      }
