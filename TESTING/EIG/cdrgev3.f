      SUBROUTINE CDRGEV3( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, QE, LDQE, ALPHA, BETA, ALPHA1, BETA1, WORK, LWORK, RWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDQ, LDQE, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), ALPHA1( * ), B( LDA, * ), BETA( * ), BETA1( * ), Q( LDQ, * ), QE( LDQE, * ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      int                I, IADD, IERR, IN, J, JC, JR, JSIZE, JTYPE, MAXWRK, MINWRK, MTYPES, N, N1, NB, NERRS, NMATS, NMAX, NTESTT;
      REAL               SAFMAX, SAFMIN, ULP, ULPINV
      COMPLEX            CTEMP
      // ..
      // .. Local Arrays ..
      bool               LASIGN( MAXTYP ), LBSIGN( MAXTYP );
      int                IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      REAL               RMAGN( 0: 3 )
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SLAMCH
      COMPLEX            CLARND
      // EXTERNAL ILAENV, SLAMCH, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, CGET52, CGGEV3, CLACPY, CLARFG, CLASET, CLATM4, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, REAL, SIGN
      // ..
      // .. Data statements ..
      DATA               KCLASS / 15*1, 10*2, 1*3 /
      DATA               KZ1 / 0, 1, 2, 1, 3, 3 /
      DATA               KZ2 / 0, 0, 1, 2, 1, 1 /
      DATA               KADD / 0, 0, 0, 0, 3, 2 /
      DATA               KATYPE / 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4*4, 0 /       DATA               KBTYPE / 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8*8, 0 /       DATA               KAZERO / 6*1, 2, 1, 2*2, 2*1, 2*2, 3, 1, 3, 4*5, 4*3, 1 /       DATA               KBZERO / 6*1, 1, 2, 2*1, 2*2, 2*1, 4, 1, 4, 4*6, 4*4, 1 /       DATA               KAMAGN / 8*1, 2, 3, 2, 3, 2, 3, 7*1, 2, 3, 3, 2, 1 /       DATA               KBMAGN / 8*1, 3, 2, 3, 2, 2, 3, 7*1, 3, 2, 3, 2, 1 /
      DATA               KTRIAN / 16*0, 10*1 /
      DATA               LASIGN / 6*.FALSE., .TRUE., .FALSE., 2*.TRUE., 2*.FALSE., 3*.TRUE., .FALSE., .TRUE., 3*.FALSE., 5*.TRUE., .FALSE. /       DATA               LBSIGN / 7*.FALSE., .TRUE., 2*.FALSE., 2*.TRUE., 2*.FALSE., .TRUE., .FALSE., .TRUE., 9*.FALSE. /
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0

      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( THRESH.LT.ZERO ) {
         INFO = -6
      } else if ( LDA.LE.1 .OR. LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDQ.LE.1 .OR. LDQ.LT.NMAX ) {
         INFO = -14
      } else if ( LDQE.LE.1 .OR. LDQE.LT.NMAX ) {
         INFO = -17
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      MINWRK = 1
      if ( INFO.EQ.0 .AND. LWORK.GE.1 ) {
         MINWRK = NMAX*( NMAX+1 )
         NB = MAX( 1, ILAENV( 1, 'CGEQRF', ' ', NMAX, NMAX, -1, -1 ), ILAENV( 1, 'CUNMQR', 'LC', NMAX, NMAX, NMAX, -1 ), ILAENV( 1, 'CUNGQR', ' ', NMAX, NMAX, NMAX, -1 ) )
         MAXWRK = MAX( 2*NMAX, NMAX*( NB+1 ), NMAX*( NMAX+1 ) )
         WORK( 1 ) = MAXWRK
      }

      IF( LWORK.LT.MINWRK ) INFO = -23

      if ( INFO.NE.0 ) {
         xerbla('CDRGEV3', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      ULP = SLAMCH( 'Precision' )
      SAFMIN = SLAMCH( 'Safe minimum' )
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

      DO 220 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / REAL( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*N1

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 210 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 210
            NMATS = NMATS + 1

            // Save ISEED in case of an error.

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Generate test matrices A and B

            // Description of control parameters:

            // KCLASS: =1 means w/o rotation, =2 means w/ rotation,
                    // =3 means random.
            // KATYPE: the "type" to be passed to CLATM4 for computing A.
            // KAZERO: the pattern of zeros on the diagonal for A:
                    // =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
                    // =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
                    // =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
                    // non-zero entries.)
            // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
                    // =2: large, =3: small.
            // LASIGN: .TRUE. if the diagonal elements of A are to be
                    // multiplied by a random magnitude 1 number.
            // KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN: used to implement KAMAGN and KBMAGN.

            IF( MTYPES.GT.MAXTYP ) GO TO 100
            IERR = 0
            if ( KCLASS( JTYPE ).LT.3 ) {

               // Generate A (w/o rotation)

               if ( ABS( KATYPE( JTYPE ) ).EQ.3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL CLASET( 'Full', N, N, CZERO, CZERO, A, LDA )
               } else {
                  IN = N
               }
               clatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) )
               IF( IADD.GT.0 .AND. IADD.LE.N ) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) )

               // Generate B (w/o rotation)

               if ( ABS( KBTYPE( JTYPE ) ).EQ.3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL CLASET( 'Full', N, N, CZERO, CZERO, B, LDA )
               } else {
                  IN = N
               }
               clatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) )
               IF( IADD.NE.0 .AND. IADD.LE.N ) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) )

               if ( KCLASS( JTYPE ).EQ.2 .AND. N.GT.0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  DO 40 JC = 1, N - 1
                     DO 30 JR = JC, N
                        Q( JR, JC ) = CLARND( 3, ISEED )
                        Z( JR, JC ) = CLARND( 3, ISEED )
   30                CONTINUE
                     clarfg(N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) );
                     WORK( 2*N+JC ) = SIGN( ONE, REAL( Q( JC, JC ) ) )
                     Q( JC, JC ) = CONE
                     clarfg(N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK( 3*N+JC ) = SIGN( ONE, REAL( Z( JC, JC ) ) )
                     Z( JC, JC ) = CONE
   40             CONTINUE
                  CTEMP = CLARND( 3, ISEED )
                  Q( N, N ) = CONE
                  WORK( N ) = CZERO
                  WORK( 3*N ) = CTEMP / ABS( CTEMP )
                  CTEMP = CLARND( 3, ISEED )
                  Z( N, N ) = CONE
                  WORK( 2*N ) = CZERO
                  WORK( 4*N ) = CTEMP / ABS( CTEMP )

                  // Apply the diagonal matrices

                  DO 60 JC = 1, N
                     DO 50 JR = 1, N
                        A( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* B( JR, JC )
   50                CONTINUE
   60             CONTINUE
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR.NE.0 ) GO TO 90                   CALL CUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR.NE.0 ) GO TO 90                   CALL CUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR.NE.0 ) GO TO 90                   CALL CUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IERR )                   IF( IERR.NE.0 ) GO TO 90
               }
            } else {

               // Random matrices

               DO 80 JC = 1, N
                  DO 70 JR = 1, N
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* CLARND( 4, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* CLARND( 4, ISEED )
   70             CONTINUE
   80          CONTINUE
            }

   90       CONTINUE

            if ( IERR.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               RETURN
            }

  100       CONTINUE

            DO 110 I = 1, 7
               RESULT( I ) = -ONE
  110       CONTINUE

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
            if ( IERR.NE.0 .AND. IERR.NE.N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'CGGEV31', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            // Do the tests (1) and (2)

            cget52(.TRUE., N, A, LDA, B, LDA, Q, LDQ, ALPHA, BETA, WORK, RWORK, RESULT( 1 ) );
            if ( RESULT( 2 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'CGGEV31', RESULT( 2 ), N, JTYPE, IOLDSD
            }

            // Do the tests (3) and (4)

            cget52(.FALSE., N, A, LDA, B, LDA, Z, LDQ, ALPHA, BETA, WORK, RWORK, RESULT( 3 ) );
            if ( RESULT( 4 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'CGGEV31', RESULT( 4 ), N, JTYPE, IOLDSD
            }

            // Do test (5)

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('N', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR.NE.0 .AND. IERR.NE.N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'CGGEV32', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            DO 120 J = 1, N
               IF( ALPHA( J ).NE.ALPHA1( J ) .OR. BETA( J ).NE. BETA1( J ) ) RESULT( 5 ) = ULPINV
  120       CONTINUE

            // Do the test (6): Compute eigenvalues and left eigenvectors,
            // and test them

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('V', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, QE, LDQE, Z, LDQ, WORK, LWORK, RWORK, IERR );
            if ( IERR.NE.0 .AND. IERR.NE.N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'CGGEV33', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }


            DO 130 J = 1, N
               if ( ALPHA( J ).NE.ALPHA1( J ) .OR. BETA( J ).NE.BETA1( J ) ) {
                  RESULT( 6 ) = ULPINV
               ENDIF
  130       CONTINUE

            DO 150 J = 1, N
               DO 140 JC = 1, N
                  if ( Q( J, JC ).NE.QE( J, JC ) ) {
                     RESULT( 6 ) = ULPINV
                  }
  140          CONTINUE
  150       CONTINUE

            // DO the test (7): Compute eigenvalues and right eigenvectors,
            // and test them

            clacpy(' ', N, N, A, LDA, S, LDA );
            clacpy(' ', N, N, B, LDA, T, LDA );
            cggev3('N', 'V', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, QE, LDQE, WORK, LWORK, RWORK, IERR );
            if ( IERR.NE.0 .AND. IERR.NE.N+1 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'CGGEV34', IERR, N, JTYPE, IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            }

            DO 160 J = 1, N
               IF( ALPHA( J ).NE.ALPHA1( J ) .OR. BETA( J ).NE. BETA1( J ) )RESULT( 7 ) = ULPINV
  160       CONTINUE

            DO 180 J = 1, N
               DO 170 JC = 1, N
                  IF( Z( J, JC ).NE.QE( J, JC ) ) RESULT( 7 ) = ULPINV
  170          CONTINUE
  180       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

  190       CONTINUE

            NTESTT = NTESTT + 7

            // Print out tests which fail.

            DO 200 JR = 1, 7
               if ( RESULT( JR ).GE.THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS.EQ.0 ) {
                     WRITE( NOUNIT, FMT = 9997 )'CGV'

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )'Orthogonal'

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 )

                  }
                  NERRS = NERRS + 1
                  if ( RESULT( JR ).LT.10000.0 ) {
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  } else {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  }
               }
  200       CONTINUE

  210    CONTINUE
  220 CONTINUE

      // Summary

      alasvm('CGV3', NOUNIT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK

      RETURN

 9999 FORMAT( ' CDRGEV3: ', A, ' returned INFO=', I6, '.', / 3X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9998 FORMAT( ' CDRGEV3: ', A, ' Eigenvectors from ', A, ' incorrectly normalized.', / ' Bits of error=', 0P, G10.3, ',', 3X, 'N=', I4, ', JTYPE=', I3, ', ISEED=(', 3( I4, ',' ), I5, ')' )

 9997 FORMAT( / 1X, A3, ' -- Complex Generalized eigenvalue problem ', 'driver' )

 9996 FORMAT( ' Matrix types (see CDRGEV3 for details): ' )

 9995 FORMAT( ' Special Matrices:', 23X, '(J''=transposed Jordan block)', / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ', '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ', 'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I', ')  11=(large*I, small*D)  13=(large*D, large*I)', / '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ', ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9994 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:', / '  16=Transposed Jordan Blocks             19=geometric ', 'alpha, beta=0,1', / '  17=arithm. alpha&beta             ', '      20=arithmetic alpha, beta=0,1', / '  18=clustered ', 'alpha, beta=0,1            21=random alpha, beta=0,1', / ' Large & Small Matrices:', / '  22=(large, small)   ', '23=(small,large)    24=(small,small)    25=(large,large)', / '  26=random O(1) matrices.' )

 9993 FORMAT( / ' Tests performed:    ', / ' 1 = max | ( b A - a B )''*l | / const.,', / ' 2 = | |VR(i)| - 1 | / ulp,', / ' 3 = max | ( b A - a B )*r | / const.', / ' 4 = | |VL(i)| - 1 | / ulp,', / ' 5 = 0 if W same no matter if r or l computed,', / ' 6 = 0 if l same no matter if l computed,', / ' 7 = 0 if r same no matter if r computed,', / 1X )
 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 1P, E10.3 )

      // End of CDRGEV3

      }
