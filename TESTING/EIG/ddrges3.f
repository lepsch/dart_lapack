      SUBROUTINE DDRGES3( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHAR, ALPHAI, BETA, WORK, LWORK, RESULT, BWORK, INFO )

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
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 26 )
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
      // EXTERNAL ALASVM, DGET51, DGET53, DGET54, DGGES3, DLACPY, DLARFG, DLASET, DLATM4, DORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SIGN
      // ..
      // .. Data statements ..
      DATA               KCLASS / 15*1, 10*2, 1*3 /
      DATA               KZ1 / 0, 1, 2, 1, 3, 3 /
      DATA               KZ2 / 0, 0, 1, 2, 1, 1 /
      DATA               KADD / 0, 0, 0, 0, 3, 2 /
      DATA               KATYPE / 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4*4, 0 /       DATA               KBTYPE / 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8*8, 0 /       DATA               KAZERO / 6*1, 2, 1, 2*2, 2*1, 2*2, 3, 1, 3, 4*5, 4*3, 1 /       DATA               KBZERO / 6*1, 1, 2, 2*1, 2*2, 2*1, 4, 1, 4, 4*6, 4*4, 1 /       DATA               KAMAGN / 8*1, 2, 3, 2, 3, 2, 3, 7*1, 2, 3, 3, 2, 1 /       DATA               KBMAGN / 8*1, 3, 2, 3, 2, 2, 3, 7*1, 3, 2, 3, 2, 1 /
      DATA               KTRIAN / 16*0, 10*1 /
      DATA               IASIGN / 6*0, 2, 0, 2*2, 2*0, 3*2, 0, 2, 3*0, 5*2, 0 /
      DATA               IBSIGN / 7*0, 2, 2*0, 2*2, 2*0, 2, 0, 2, 9*0 /
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

      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LE.1 .OR. LDA.LT.NMAX ) THEN
         INFO = -9
      ELSE IF( LDQ.LE.1 .OR. LDQ.LT.NMAX ) THEN
         INFO = -14
      END IF

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
         MINWRK = MAX( 10*( NMAX+1 ), 3*NMAX*NMAX )
         NB = MAX( 1, ILAENV( 1, 'DGEQRF', ' ', NMAX, NMAX, -1, -1 ), ILAENV( 1, 'DORMQR', 'LT', NMAX, NMAX, NMAX, -1 ), ILAENV( 1, 'DORGQR', ' ', NMAX, NMAX, NMAX, -1 ) )
         MAXWRK = MAX( 10*( NMAX+1 ), 2*NMAX+NMAX*NB, 3*NMAX*NMAX )
         WORK( 1 ) = MAXWRK
      END IF

      IF( LWORK.LT.MINWRK ) INFO = -20

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DDRGES3', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

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

      DO 190 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / DBLE( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*DBLE( N1 )

         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF

         // Loop over matrix types

         DO 180 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 180
            NMATS = NMATS + 1
            NTEST = 0

            // Save ISEED in case of an error.

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Initialize RESULT

            DO 30 J = 1, 13
               RESULT( J ) = ZERO
   30       CONTINUE

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
                   t // o form 2x2 blocks.
            // KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN: used to implement KAMAGN and KBMAGN.

            IF( MTYPES.GT.MAXTYP ) GO TO 110
            IINFO = 0
            IF( KCLASS( JTYPE ).LT.3 ) THEN

               // Generate A (w/o rotation)

               IF( ABS( KATYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL DLASET( 'Full', N, N, ZERO, ZERO, A, LDA )
               ELSE
                  IN = N
               END IF
               CALL DLATM4( KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), IASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA )
               IADD = KADD( KAZERO( JTYPE ) )
               IF( IADD.GT.0 .AND. IADD.LE.N ) A( IADD, IADD ) = ONE

               // Generate B (w/o rotation)

               IF( ABS( KBTYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL DLASET( 'Full', N, N, ZERO, ZERO, B, LDA )
               ELSE
                  IN = N
               END IF
               CALL DLATM4( KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), IBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA )
               IADD = KADD( KBZERO( JTYPE ) )
               IF( IADD.NE.0 .AND. IADD.LE.N ) B( IADD, IADD ) = ONE

               IF( KCLASS( JTYPE ).EQ.2 .AND. N.GT.0 ) THEN

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  DO 50 JC = 1, N - 1
                     DO 40 JR = JC, N
                        Q( JR, JC ) = DLARND( 3, ISEED )
                        Z( JR, JC ) = DLARND( 3, ISEED )
   40                CONTINUE
                     CALL DLARFG( N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, Q( JC, JC ) )
                     Q( JC, JC ) = ONE
                     CALL DLARFG( N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, Z( JC, JC ) )
                     Z( JC, JC ) = ONE
   50             CONTINUE
                  Q( N, N ) = ONE
                  WORK( N ) = ZERO
                  WORK( 3*N ) = SIGN( ONE, DLARND( 2, ISEED ) )
                  Z( N, N ) = ONE
                  WORK( 2*N ) = ZERO
                  WORK( 4*N ) = SIGN( ONE, DLARND( 2, ISEED ) )

                  // Apply the diagonal matrices

                  DO 70 JC = 1, N
                     DO 60 JR = 1, N
                        A( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )* B( JR, JC )
   60                CONTINUE
   70             CONTINUE
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100
               END IF
            ELSE

               // Random matrices

               DO 90 JC = 1, N
                  DO 80 JR = 1, N
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* DLARND( 2, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* DLARND( 2, ISEED )
   80             CONTINUE
   90          CONTINUE
            END IF

  100       CONTINUE

            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF

  110       CONTINUE

            DO 120 I = 1, 13
               RESULT( I ) = -ONE
  120       CONTINUE

            // Test with and without sorting of eigenvalues

            DO 150 ISORT = 0, 1
               IF( ISORT.EQ.0 ) THEN
                  SORT = 'N'
                  RSUB = 0
               ELSE
                  SORT = 'S'
                  RSUB = 5
               END IF

               // Call XLAENV to set the parameters used in DLAQZ0

               CALL XLAENV( 12, 10 )
               CALL XLAENV( 13, 12 )
               CALL XLAENV( 14, 13 )
               CALL XLAENV( 15, 2 )
               CALL XLAENV( 17, 10 )

               // Call DGGES3 to compute H, T, Q, Z, alpha, and beta.

               CALL DLACPY( 'Full', N, N, A, LDA, S, LDA )
               CALL DLACPY( 'Full', N, N, B, LDA, T, LDA )
               NTEST = 1 + RSUB + ISORT
               RESULT( 1+RSUB+ISORT ) = ULPINV
               CALL DGGES3( 'V', 'V', SORT, DLCTES, N, S, LDA, T, LDA, SDIM, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, BWORK, IINFO )
               IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
                  RESULT( 1+RSUB+ISORT ) = ULPINV
                  WRITE( NOUNIT, FMT = 9999 )'DGGES3', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 160
               END IF

               NTEST = 4 + RSUB

               // Do tests 1--4 (or tests 7--9 when reordering )

               IF( ISORT.EQ.0 ) THEN
                  CALL DGET51( 1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 1 ) )                   CALL DGET51( 1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 2 ) )
               ELSE
                  CALL DGET54( N, A, LDA, B, LDA, S, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 7 ) )
               END IF
               CALL DGET51( 3, N, A, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK, RESULT( 3+RSUB ) )                CALL DGET51( 3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK, RESULT( 4+RSUB ) )

               // Do test 5 and 6 (or Tests 10 and 11 when reordering):
               // check Schur form of A and compare eigenvalues with
               // diagonals.

               NTEST = 6 + RSUB
               TEMP1 = ZERO

               DO 130 J = 1, N
                  ILABAD = .FALSE.
                  IF( ALPHAI( J ).EQ.ZERO ) THEN
                     TEMP2 = ( ABS( ALPHAR( J )-S( J, J ) ) / MAX( SAFMIN, ABS( ALPHAR( J ) ), ABS( S( J, J ) ) )+ABS( BETA( J )-T( J, J ) ) / MAX( SAFMIN, ABS( BETA( J ) ), ABS( T( J, J ) ) ) ) / ULP

                     IF( J.LT.N ) THEN
                        IF( S( J+1, J ).NE.ZERO ) THEN
                           ILABAD = .TRUE.
                           RESULT( 5+RSUB ) = ULPINV
                        END IF
                     END IF
                     IF( J.GT.1 ) THEN
                        IF( S( J, J-1 ).NE.ZERO ) THEN
                           ILABAD = .TRUE.
                           RESULT( 5+RSUB ) = ULPINV
                        END IF
                     END IF

                  ELSE
                     IF( ALPHAI( J ).GT.ZERO ) THEN
                        I1 = J
                     ELSE
                        I1 = J - 1
                     END IF
                     IF( I1.LE.0 .OR. I1.GE.N ) THEN
                        ILABAD = .TRUE.
                     ELSE IF( I1.LT.N-1 ) THEN
                        IF( S( I1+2, I1+1 ).NE.ZERO ) THEN
                           ILABAD = .TRUE.
                           RESULT( 5+RSUB ) = ULPINV
                        END IF
                     ELSE IF( I1.GT.1 ) THEN
                        IF( S( I1, I1-1 ).NE.ZERO ) THEN
                           ILABAD = .TRUE.
                           RESULT( 5+RSUB ) = ULPINV
                        END IF
                     END IF
                     IF( .NOT.ILABAD ) THEN
                        CALL DGET53( S( I1, I1 ), LDA, T( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IERR )
                        IF( IERR.GE.3 ) THEN
                           WRITE( NOUNIT, FMT = 9998 )IERR, J, N, JTYPE, IOLDSD
                           INFO = ABS( IERR )
                        END IF
                     ELSE
                        TEMP2 = ULPINV
                     END IF

                  END IF
                  TEMP1 = MAX( TEMP1, TEMP2 )
                  IF( ILABAD ) THEN
                     WRITE( NOUNIT, FMT = 9997 )J, N, JTYPE, IOLDSD
                  END IF
  130          CONTINUE
               RESULT( 6+RSUB ) = TEMP1

               IF( ISORT.GE.1 ) THEN

                  // Do test 12

                  NTEST = 12
                  RESULT( 12 ) = ZERO
                  KNTEIG = 0
                  DO 140 I = 1, N
                     IF( DLCTES( ALPHAR( I ), ALPHAI( I ), BETA( I ) ) .OR. DLCTES( ALPHAR( I ), -ALPHAI( I ), BETA( I ) ) ) THEN
                        KNTEIG = KNTEIG + 1
                     END IF
                     IF( I.LT.N ) THEN
                        IF( ( DLCTES( ALPHAR( I+1 ), ALPHAI( I+1 ), BETA( I+1 ) ) .OR. DLCTES( ALPHAR( I+1 ), -ALPHAI( I+1 ), BETA( I+1 ) ) ) .AND. ( .NOT.( DLCTES( ALPHAR( I ), ALPHAI( I ), BETA( I ) ) .OR. DLCTES( ALPHAR( I ), -ALPHAI( I ), BETA( I ) ) ) ) .AND. IINFO.NE.N+2 ) THEN
                           RESULT( 12 ) = ULPINV
                        END IF
                     END IF
  140             CONTINUE
                  IF( SDIM.NE.KNTEIG ) THEN
                     RESULT( 12 ) = ULPINV
                  END IF
               END IF

  150       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

  160       CONTINUE

            NTESTT = NTESTT + NTEST

            // Print out tests which fail.

            DO 170 JR = 1, NTEST
               IF( RESULT( JR ).GE.THRESH ) THEN

                  // If this is the first test to fail,
                  // print a header to the data file.

                  IF( NERRS.EQ.0 ) THEN
                     WRITE( NOUNIT, FMT = 9996 )'DGS'

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )
                     WRITE( NOUNIT, FMT = 9993 )'Orthogonal'

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', J = 1, 8 )

                  END IF
                  NERRS = NERRS + 1
                  IF( RESULT( JR ).LT.10000.0D0 ) THEN
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  ELSE
                     WRITE( NOUNIT, FMT = 9990 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  END IF
               END IF
  170       CONTINUE

  180    CONTINUE
  190 CONTINUE

      // Summary

      CALL ALASVM( 'DGS', NOUNIT, NERRS, NTESTT, 0 )

      WORK( 1 ) = MAXWRK

      RETURN

 9999 FORMAT( ' DDRGES3: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' )

 9998 FORMAT( ' DDRGES3: DGET53 returned INFO=', I1, ' for eigenvalue ',
     $      I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(',
     $      4( I4, ',' ), I5, ')' )

 9997 FORMAT( ' DDRGES3: S not in Schur form at eigenvalue ', I6, '.',
     $      / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )

 9996 FORMAT( / 1X, A3, ' -- Real Generalized Schur form driver' )

 9995 FORMAT( ' Matrix types (see DDRGES3 for details): ' )

 9994 FORMAT( ' Special Matrices:', 23X,
     $      '(J''=transposed Jordan block)',
     $      / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',
     $      '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ',
     $      'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I',
     $      ')  11=(large*I, small*D)  13=(large*D, large*I)', /
     $      '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ',
     $      ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9993 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:',
     $      / '  16=Transposed Jordan Blocks             19=geometric ',
     $      'alpha, beta=0,1', / '  17=arithm. alpha&beta             ',
     $      '      20=arithmetic alpha, beta=0,1', / '  18=clustered ',
     $      'alpha, beta=0,1            21=random alpha, beta=0,1',
     $      / ' Large & Small Matrices:', / '  22=(large, small)   ',
     $      '23=(small,large)    24=(small,small)    25=(large,large)',
     $      / '  26=random O(1) matrices.' )

 9992 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ',
     $      'Q and Z are ', A, ',', / 19X,
     $      'l and r are the appropriate left and right', / 19X,
     $      'eigenvectors, resp., a is alpha, b is beta, and', / 19X, A,
     $      ' means ', A, '.)', / ' Without ordering: ',
     $      / '  1 = | A - Q S Z', A,
     $      ' | / ( |A| n ulp )      2 = | B - Q T Z', A,
     $      ' | / ( |B| n ulp )', / '  3 = | I - QQ', A,
     $      ' | / ( n ulp )             4 = | I - ZZ', A,
     $      ' | / ( n ulp )', / '  5 = A is in Schur form S',
     $      / '  6 = difference between (alpha,beta)',
     $      ' and diagonals of (S,T)', / ' With ordering: ',
     $      / '  7 = | (A,B) - Q (S,T) Z', A,
     $      ' | / ( |(A,B)| n ulp )  ', / '  8 = | I - QQ', A,
     $      ' | / ( n ulp )            9 = | I - ZZ', A,
     $      ' | / ( n ulp )', / ' 10 = A is in Schur form S',
     $      / ' 11 = difference between (alpha,beta) and diagonals',
     $      ' of (S,T)', / ' 12 = SDIM is the correct number of ',
     $      'selected eigenvalues', / )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9990 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I2, ' is', 1P, D10.3 )

      // End of DDRGES3

      END
