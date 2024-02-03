      SUBROUTINE CDRGES( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHA, BETA, WORK, LWORK, RWORK, RESULT, BWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               BWORK( * ), DOTYPE( * );
      int                ISEED( 4 ), NN( * );
      REAL               RESULT( 13 ), RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDA, * ), BETA( * ), Q( LDQ, * ), S( LDA, * ), T( LDA, * ), WORK( * ), Z( LDQ, * )
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
      bool               BADNN, ILABAD;
      String             SORT;
      int                I, IADD, IINFO, IN, ISORT, J, JC, JR, JSIZE, JTYPE, KNTEIG, MAXWRK, MINWRK, MTYPES, N, N1, NB, NERRS, NMATS, NMAX, NTEST, NTESTT, RSUB, SDIM;
      REAL               SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV
      COMPLEX            CTEMP, X
      // ..
      // .. Local Arrays ..
      bool               LASIGN( MAXTYP ), LBSIGN( MAXTYP );
      int                IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      REAL               RMAGN( 0: 3 )
      // ..
      // .. External Functions ..
      bool               CLCTES;
      int                ILAENV;
      REAL               SLAMCH
      COMPLEX            CLARND
      // EXTERNAL CLCTES, ILAENV, SLAMCH, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, CGET51, CGET54, CGGES, CLACPY, CLARFG, CLASET, CLATM4, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, MIN, REAL, SIGN
      // ..
      // .. Statement Functions ..
      REAL               ABS1
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( REAL( X ) ) + ABS( AIMAG( X ) )
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
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      MINWRK = 1
      if ( INFO.EQ.0 .AND. LWORK.GE.1 ) {
         MINWRK = 3*NMAX*NMAX
         NB = MAX( 1, ILAENV( 1, 'CGEQRF', ' ', NMAX, NMAX, -1, -1 ), ILAENV( 1, 'CUNMQR', 'LC', NMAX, NMAX, NMAX, -1 ), ILAENV( 1, 'CUNGQR', ' ', NMAX, NMAX, NMAX, -1 ) )
         MAXWRK = MAX( NMAX+NMAX*NB, 3*NMAX*NMAX )
         WORK( 1 ) = MAXWRK
      }

      IF( LWORK.LT.MINWRK ) INFO = -19

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CDRGES', -INFO )
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

      // Loop over matrix sizes

      NTESTT = 0
      NERRS = 0
      NMATS = 0

      DO 190 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / REAL( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*REAL( N1 )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

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

            IF( MTYPES.GT.MAXTYP ) GO TO 110
            IINFO = 0
            if ( KCLASS( JTYPE ).LT.3 ) {

               // Generate A (w/o rotation)

               if ( ABS( KATYPE( JTYPE ) ).EQ.3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL CLASET( 'Full', N, N, CZERO, CZERO, A, LDA )
               } else {
                  IN = N
               }
               CALL CLATM4( KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2, ISEED, A, LDA )
               IADD = KADD( KAZERO( JTYPE ) )
               IF( IADD.GT.0 .AND. IADD.LE.N ) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) )

               // Generate B (w/o rotation)

               if ( ABS( KBTYPE( JTYPE ) ).EQ.3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL CLASET( 'Full', N, N, CZERO, CZERO, B, LDA )
               } else {
                  IN = N
               }
               CALL CLATM4( KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2, ISEED, B, LDA )
               IADD = KADD( KBZERO( JTYPE ) )
               IF( IADD.NE.0 .AND. IADD.LE.N ) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) )

               if ( KCLASS( JTYPE ).EQ.2 .AND. N.GT.0 ) {

                  // Include rotations

                  // Generate Q, Z as Householder transformations times
                  // a diagonal matrix.

                  DO 50 JC = 1, N - 1
                     DO 40 JR = JC, N
                        Q( JR, JC ) = CLARND( 3, ISEED )
                        Z( JR, JC ) = CLARND( 3, ISEED )
   40                CONTINUE
                     CALL CLARFG( N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1, WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, REAL( Q( JC, JC ) ) )
                     Q( JC, JC ) = CONE
                     CALL CLARFG( N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1, WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, REAL( Z( JC, JC ) ) )
                     Z( JC, JC ) = CONE
   50             CONTINUE
                  CTEMP = CLARND( 3, ISEED )
                  Q( N, N ) = CONE
                  WORK( N ) = CZERO
                  WORK( 3*N ) = CTEMP / ABS( CTEMP )
                  CTEMP = CLARND( 3, ISEED )
                  Z( N, N ) = CONE
                  WORK( 2*N ) = CZERO
                  WORK( 4*N ) = CTEMP / ABS( CTEMP )

                  // Apply the diagonal matrices

                  DO 70 JC = 1, N
                     DO 60 JR = 1, N
                        A( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* B( JR, JC )
   60                CONTINUE
   70             CONTINUE
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL CUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL CUNM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL CUNM2R( 'R', 'C', N, N, N-1, Z, LDQ, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100
               }
            } else {

               // Random matrices

               DO 90 JC = 1, N
                  DO 80 JR = 1, N
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* CLARND( 4, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* CLARND( 4, ISEED )
   80             CONTINUE
   90          CONTINUE
            }

  100       CONTINUE

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

  110       CONTINUE

            DO 120 I = 1, 13
               RESULT( I ) = -ONE
  120       CONTINUE

            // Test with and without sorting of eigenvalues

            DO 150 ISORT = 0, 1
               if ( ISORT.EQ.0 ) {
                  SORT = 'N'
                  RSUB = 0
               } else {
                  SORT = 'S'
                  RSUB = 5
               }

               // Call CGGES to compute H, T, Q, Z, alpha, and beta.

               CALL CLACPY( 'Full', N, N, A, LDA, S, LDA )
               CALL CLACPY( 'Full', N, N, B, LDA, T, LDA )
               NTEST = 1 + RSUB + ISORT
               RESULT( 1+RSUB+ISORT ) = ULPINV
               CALL CGGES( 'V', 'V', SORT, CLCTES, N, S, LDA, T, LDA, SDIM, ALPHA, BETA, Q, LDQ, Z, LDQ, WORK, LWORK, RWORK, BWORK, IINFO )
               if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
                  RESULT( 1+RSUB+ISORT ) = ULPINV
                  WRITE( NOUNIT, FMT = 9999 )'CGGES', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 160
               }

               NTEST = 4 + RSUB

               // Do tests 1--4 (or tests 7--9 when reordering )

               if ( ISORT.EQ.0 ) {
                  CALL CGET51( 1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK, RWORK, RESULT( 1 ) )                   CALL CGET51( 1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RWORK, RESULT( 2 ) )
               } else {
                  CALL CGET54( N, A, LDA, B, LDA, S, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT( 2+RSUB ) )
               }

               CALL CGET51( 3, N, B, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK, RWORK, RESULT( 3+RSUB ) )                CALL CGET51( 3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK, RWORK, RESULT( 4+RSUB ) )

               // Do test 5 and 6 (or Tests 10 and 11 when reordering):
               // check Schur form of A and compare eigenvalues with
               // diagonals.

               NTEST = 6 + RSUB
               TEMP1 = ZERO

               DO 130 J = 1, N
                  ILABAD = .FALSE.
                  TEMP2 = ( ABS1( ALPHA( J )-S( J, J ) ) / MAX( SAFMIN, ABS1( ALPHA( J ) ), ABS1( S( J, J ) ) )+ABS1( BETA( J )-T( J, J ) ) / MAX( SAFMIN, ABS1( BETA( J ) ), ABS1( T( J, J ) ) ) ) / ULP

                  if ( J.LT.N ) {
                     if ( S( J+1, J ).NE.ZERO ) {
                        ILABAD = .TRUE.
                        RESULT( 5+RSUB ) = ULPINV
                     }
                  }
                  if ( J.GT.1 ) {
                     if ( S( J, J-1 ).NE.ZERO ) {
                        ILABAD = .TRUE.
                        RESULT( 5+RSUB ) = ULPINV
                     }
                  }
                  TEMP1 = MAX( TEMP1, TEMP2 )
                  if ( ILABAD ) {
                     WRITE( NOUNIT, FMT = 9998 )J, N, JTYPE, IOLDSD
                  }
  130          CONTINUE
               RESULT( 6+RSUB ) = TEMP1

               if ( ISORT.GE.1 ) {

                  // Do test 12

                  NTEST = 12
                  RESULT( 12 ) = ZERO
                  KNTEIG = 0
                  DO 140 I = 1, N
                     IF( CLCTES( ALPHA( I ), BETA( I ) ) ) KNTEIG = KNTEIG + 1
  140             CONTINUE
                  IF( SDIM.NE.KNTEIG ) RESULT( 13 ) = ULPINV
               }

  150       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

  160       CONTINUE

            NTESTT = NTESTT + NTEST

            // Print out tests which fail.

            DO 170 JR = 1, NTEST
               if ( RESULT( JR ).GE.THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS.EQ.0 ) {
                     WRITE( NOUNIT, FMT = 9997 )'CGS'

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )'Unitary'

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 )'unitary', '''', 'transpose', ( '''', J = 1, 8 )

                  }
                  NERRS = NERRS + 1
                  if ( RESULT( JR ).LT.10000.0 ) {
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  } else {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  }
               }
  170       CONTINUE

  180    CONTINUE
  190 CONTINUE

      // Summary

      CALL ALASVM( 'CGS', NOUNIT, NERRS, NTESTT, 0 )

      WORK( 1 ) = MAXWRK

      RETURN

 9999 FORMAT( ' CDRGES: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' )

 9998 FORMAT( ' CDRGES: S not in Schur form at eigenvalue ', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9997 FORMAT( / 1X, A3, ' -- Complex Generalized Schur from problem ', 'driver' )

 9996 FORMAT( ' Matrix types (see CDRGES for details): ' )

 9995 FORMAT( ' Special Matrices:', 23X, '(J''=transposed Jordan block)', / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ', '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ', 'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I', ')  11=(large*I, small*D)  13=(large*D, large*I)', / '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ', ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9994 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:', / '  16=Transposed Jordan Blocks             19=geometric ', 'alpha, beta=0,1', / '  17=arithm. alpha&beta             ', '      20=arithmetic alpha, beta=0,1', / '  18=clustered ', 'alpha, beta=0,1            21=random alpha, beta=0,1', / ' Large & Small Matrices:', / '  22=(large, small)   ', '23=(small,large)    24=(small,small)    25=(large,large)', / '  26=random O(1) matrices.' )

 9993 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ', 'Q and Z are ', A, ',', / 19X, 'l and r are the appropriate left and right', / 19X, 'eigenvectors, resp., a is alpha, b is beta, and', / 19X, A, ' means ', A, '.)', / ' Without ordering: ', / '  1 = | A - Q S Z', A, ' | / ( |A| n ulp )      2 = | B - Q T Z', A, ' | / ( |B| n ulp )', / '  3 = | I - QQ', A, ' | / ( n ulp )             4 = | I - ZZ', A, ' | / ( n ulp )', / '  5 = A is in Schur form S', / '  6 = difference between (alpha,beta)', ' and diagonals of (S,T)', / ' With ordering: ', / '  7 = | (A,B) - Q (S,T) Z', A, ' | / ( |(A,B)| n ulp )', / '  8 = | I - QQ', A, ' | / ( n ulp )             9 = | I - ZZ', A, ' | / ( n ulp )', / ' 10 = A is in Schur form S', / ' 11 = difference between (alpha,beta) and diagonals', ' of (S,T)', / ' 12 = SDIM is the correct number of ', 'selected eigenvalues', / )
 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 1P, E10.3 )

      // End of CDRGES

      }
