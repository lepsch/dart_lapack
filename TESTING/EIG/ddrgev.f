      SUBROUTINE DDRGEV( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, QE, LDQE,
     $                   ALPHAR, ALPHAI, BETA, ALPHR1, ALPHI1, BETA1,
     $                   WORK, LWORK, RESULT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDQ, LDQE, LWORK, NOUNIT, NSIZES,
     $                   NTYPES
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), NN( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   ALPHI1( * ), ALPHR1( * ), B( LDA, * ),
     $                   BETA( * ), BETA1( * ), Q( LDQ, * ),
     $                   QE( LDQE, * ), RESULT( * ), S( LDA, * ),
     $                   T( LDA, * ), WORK( * ), Z( LDQ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 26 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN
      INTEGER            I, IADD, IERR, IN, J, JC, JR, JSIZE, JTYPE,
     $                   MAXWRK, MINWRK, MTYPES, N, N1, NERRS, NMATS,
     $                   NMAX, NTESTT
      DOUBLE PRECISION   SAFMAX, SAFMIN, ULP, ULPINV
*     ..
*     .. Local Arrays ..
      INTEGER            IASIGN( MAXTYP ), IBSIGN( MAXTYP ),
     $                   IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ),
     $                   KATYPE( MAXTYP ), KAZERO( MAXTYP ),
     $                   KBMAGN( MAXTYP ), KBTYPE( MAXTYP ),
     $                   KBZERO( MAXTYP ), KCLASS( MAXTYP ),
     $                   KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 )
      DOUBLE PRECISION   RMAGN( 0: 3 )
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLARND
      EXTERNAL           ILAENV, DLAMCH, DLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, DGET52, DGGEV, DLACPY, DLARFG,
     $                   DLASET, DLATM4, DORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN
*     ..
*     .. Data statements ..
      DATA               KCLASS / 15*1, 10*2, 1*3 /
      DATA               KZ1 / 0, 1, 2, 1, 3, 3 /
      DATA               KZ2 / 0, 0, 1, 2, 1, 1 /
      DATA               KADD / 0, 0, 0, 0, 3, 2 /
      DATA               KATYPE / 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4,
     $                   4, 4, 2, 4, 5, 8, 7, 9, 4*4, 0 /
      DATA               KBTYPE / 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4,
     $                   1, 1, -4, 2, -4, 8*8, 0 /
      DATA               KAZERO / 6*1, 2, 1, 2*2, 2*1, 2*2, 3, 1, 3,
     $                   4*5, 4*3, 1 /
      DATA               KBZERO / 6*1, 1, 2, 2*1, 2*2, 2*1, 4, 1, 4,
     $                   4*6, 4*4, 1 /
      DATA               KAMAGN / 8*1, 2, 3, 2, 3, 2, 3, 7*1, 2, 3, 3,
     $                   2, 1 /
      DATA               KBMAGN / 8*1, 3, 2, 3, 2, 2, 3, 7*1, 3, 2, 3,
     $                   2, 1 /
      DATA               KTRIAN / 16*0, 10*1 /
      DATA               IASIGN / 6*0, 2, 0, 2*2, 2*0, 3*2, 0, 2, 3*0,
     $                   5*2, 0 /
      DATA               IBSIGN / 7*0, 2, 2*0, 2*2, 2*0, 2, 0, 2, 9*0 /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      INFO = 0
*
      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 )
     $      BADNN = .TRUE.
   10 CONTINUE
*
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
      ELSE IF( LDQE.LE.1 .OR. LDQE.LT.NMAX ) THEN
         INFO = -17
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.
*
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
         MINWRK = MAX( 1, 8*NMAX, NMAX*( NMAX+1 ) )
         MAXWRK = 7*NMAX + NMAX*ILAENV( 1, 'DGEQRF', ' ', NMAX, 1, NMAX,
     $            0 )
         MAXWRK = MAX( MAXWRK, NMAX*( NMAX+1 ) )
         WORK( 1 ) = MAXWRK
      END IF
*
      IF( LWORK.LT.MINWRK )
     $   INFO = -25
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DDRGEV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SAFMIN = SAFMIN / ULP
      SAFMAX = ONE / SAFMIN
      ULPINV = ONE / ULP
*
*     The values RMAGN(2:3) depend on N, see below.
*
      RMAGN( 0 ) = ZERO
      RMAGN( 1 ) = ONE
*
*     Loop over sizes, types
*
      NTESTT = 0
      NERRS = 0
      NMATS = 0
*
      DO 220 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / DBLE( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*N1
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 210 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 210
            NMATS = NMATS + 1
*
*           Save ISEED in case of an error.
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
*           Generate test matrices A and B
*
*           Description of control parameters:
*
*           KZLASS: =1 means w/o rotation, =2 means w/ rotation,
*                   =3 means random.
*           KATYPE: the "type" to be passed to DLATM4 for computing A.
*           KAZERO: the pattern of zeros on the diagonal for A:
*                   =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
*                   =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
*                   =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
*                   non-zero entries.)
*           KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
*                   =2: large, =3: small.
*           IASIGN: 1 if the diagonal elements of A are to be
*                   multiplied by a random magnitude 1 number, =2 if
*                   randomly chosen diagonal blocks are to be rotated
*                   to form 2x2 blocks.
*           KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
*           KTRIAN: =0: don't fill in the upper triangle, =1: do.
*           KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
*           RMAGN: used to implement KAMAGN and KBMAGN.
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 100
            IERR = 0
            IF( KCLASS( JTYPE ).LT.3 ) THEN
*
*              Generate A (w/o rotation)
*
               IF( ABS( KATYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N )
     $               CALL DLASET( 'Full', N, N, ZERO, ZERO, A, LDA )
               ELSE
                  IN = N
               END IF
               CALL DLATM4( KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ),
     $                      KZ2( KAZERO( JTYPE ) ), IASIGN( JTYPE ),
     $                      RMAGN( KAMAGN( JTYPE ) ), ULP,
     $                      RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2,
     $                      ISEED, A, LDA )
               IADD = KADD( KAZERO( JTYPE ) )
               IF( IADD.GT.0 .AND. IADD.LE.N )
     $            A( IADD, IADD ) = ONE
*
*              Generate B (w/o rotation)
*
               IF( ABS( KBTYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N )
     $               CALL DLASET( 'Full', N, N, ZERO, ZERO, B, LDA )
               ELSE
                  IN = N
               END IF
               CALL DLATM4( KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ),
     $                      KZ2( KBZERO( JTYPE ) ), IBSIGN( JTYPE ),
     $                      RMAGN( KBMAGN( JTYPE ) ), ONE,
     $                      RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2,
     $                      ISEED, B, LDA )
               IADD = KADD( KBZERO( JTYPE ) )
               IF( IADD.NE.0 .AND. IADD.LE.N )
     $            B( IADD, IADD ) = ONE
*
               IF( KCLASS( JTYPE ).EQ.2 .AND. N.GT.0 ) THEN
*
*                 Include rotations
*
*                 Generate Q, Z as Householder transformations times
*                 a diagonal matrix.
*
                  DO 40 JC = 1, N - 1
                     DO 30 JR = JC, N
                        Q( JR, JC ) = DLARND( 3, ISEED )
                        Z( JR, JC ) = DLARND( 3, ISEED )
   30                CONTINUE
                     CALL DLARFG( N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1,
     $                            WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, Q( JC, JC ) )
                     Q( JC, JC ) = ONE
                     CALL DLARFG( N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1,
     $                            WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, Z( JC, JC ) )
                     Z( JC, JC ) = ONE
   40             CONTINUE
                  Q( N, N ) = ONE
                  WORK( N ) = ZERO
                  WORK( 3*N ) = SIGN( ONE, DLARND( 2, ISEED ) )
                  Z( N, N ) = ONE
                  WORK( 2*N ) = ZERO
                  WORK( 4*N ) = SIGN( ONE, DLARND( 2, ISEED ) )
*
*                 Apply the diagonal matrices
*
                  DO 60 JC = 1, N
                     DO 50 JR = 1, N
                        A( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )*
     $                                A( JR, JC )
                        B( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )*
     $                                B( JR, JC )
   50                CONTINUE
   60             CONTINUE
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A,
     $                         LDA, WORK( 2*N+1 ), IERR )
                  IF( IERR.NE.0 )
     $               GO TO 90
                  CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ),
     $                         A, LDA, WORK( 2*N+1 ), IERR )
                  IF( IERR.NE.0 )
     $               GO TO 90
                  CALL DORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B,
     $                         LDA, WORK( 2*N+1 ), IERR )
                  IF( IERR.NE.0 )
     $               GO TO 90
                  CALL DORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ),
     $                         B, LDA, WORK( 2*N+1 ), IERR )
                  IF( IERR.NE.0 )
     $               GO TO 90
               END IF
            ELSE
*
*              Random matrices
*
               DO 80 JC = 1, N
                  DO 70 JR = 1, N
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )*
     $                             DLARND( 2, ISEED )
                     B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )*
     $                             DLARND( 2, ISEED )
   70             CONTINUE
   80          CONTINUE
            END IF
*
   90       CONTINUE
*
            IF( IERR.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IERR, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IERR )
               RETURN
            END IF
*
  100       CONTINUE
*
            DO 110 I = 1, 7
               RESULT( I ) = -ONE
  110       CONTINUE
*
*           Call DGGEV to compute eigenvalues and eigenvectors.
*
            CALL DLACPY( ' ', N, N, A, LDA, S, LDA )
            CALL DLACPY( ' ', N, N, B, LDA, T, LDA )
            CALL DGGEV( 'V', 'V', N, S, LDA, T, LDA, ALPHAR, ALPHAI,
     $                  BETA, Q, LDQ, Z, LDQ, WORK, LWORK, IERR )
            IF( IERR.NE.0 .AND. IERR.NE.N+1 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV1', IERR, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            END IF
*
*           Do the tests (1) and (2)
*
            CALL DGET52( .TRUE., N, A, LDA, B, LDA, Q, LDQ, ALPHAR,
     $                   ALPHAI, BETA, WORK, RESULT( 1 ) )
            IF( RESULT( 2 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'DGGEV1',
     $            RESULT( 2 ), N, JTYPE, IOLDSD
            END IF
*
*           Do the tests (3) and (4)
*
            CALL DGET52( .FALSE., N, A, LDA, B, LDA, Z, LDQ, ALPHAR,
     $                   ALPHAI, BETA, WORK, RESULT( 3 ) )
            IF( RESULT( 4 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'DGGEV1',
     $            RESULT( 4 ), N, JTYPE, IOLDSD
            END IF
*
*           Do the test (5)
*
            CALL DLACPY( ' ', N, N, A, LDA, S, LDA )
            CALL DLACPY( ' ', N, N, B, LDA, T, LDA )
            CALL DGGEV( 'N', 'N', N, S, LDA, T, LDA, ALPHR1, ALPHI1,
     $                  BETA1, Q, LDQ, Z, LDQ, WORK, LWORK, IERR )
            IF( IERR.NE.0 .AND. IERR.NE.N+1 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV2', IERR, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            END IF
*
            DO 120 J = 1, N
               IF( ALPHAR( J ).NE.ALPHR1( J ) .OR. ALPHAI( J ).NE.
     $             ALPHI1( J ) .OR. BETA( J ).NE.BETA1( J ) )RESULT( 5 )
     $              = ULPINV
  120       CONTINUE
*
*           Do the test (6): Compute eigenvalues and left eigenvectors,
*           and test them
*
            CALL DLACPY( ' ', N, N, A, LDA, S, LDA )
            CALL DLACPY( ' ', N, N, B, LDA, T, LDA )
            CALL DGGEV( 'V', 'N', N, S, LDA, T, LDA, ALPHR1, ALPHI1,
     $                  BETA1, QE, LDQE, Z, LDQ, WORK, LWORK, IERR )
            IF( IERR.NE.0 .AND. IERR.NE.N+1 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV3', IERR, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            END IF
*
            DO 130 J = 1, N
               IF( ALPHAR( J ).NE.ALPHR1( J ) .OR. ALPHAI( J ).NE.
     $             ALPHI1( J ) .OR. BETA( J ).NE.BETA1( J ) )RESULT( 6 )
     $              = ULPINV
  130       CONTINUE
*
            DO 150 J = 1, N
               DO 140 JC = 1, N
                  IF( Q( J, JC ).NE.QE( J, JC ) )
     $               RESULT( 6 ) = ULPINV
  140          CONTINUE
  150       CONTINUE
*
*           DO the test (7): Compute eigenvalues and right eigenvectors,
*           and test them
*
            CALL DLACPY( ' ', N, N, A, LDA, S, LDA )
            CALL DLACPY( ' ', N, N, B, LDA, T, LDA )
            CALL DGGEV( 'N', 'V', N, S, LDA, T, LDA, ALPHR1, ALPHI1,
     $                  BETA1, Q, LDQ, QE, LDQE, WORK, LWORK, IERR )
            IF( IERR.NE.0 .AND. IERR.NE.N+1 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGGEV4', IERR, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IERR )
               GO TO 190
            END IF
*
            DO 160 J = 1, N
               IF( ALPHAR( J ).NE.ALPHR1( J ) .OR. ALPHAI( J ).NE.
     $             ALPHI1( J ) .OR. BETA( J ).NE.BETA1( J ) )RESULT( 7 )
     $              = ULPINV
  160       CONTINUE
*
            DO 180 J = 1, N
               DO 170 JC = 1, N
                  IF( Z( J, JC ).NE.QE( J, JC ) )
     $               RESULT( 7 ) = ULPINV
  170          CONTINUE
  180       CONTINUE
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
  190       CONTINUE
*
            NTESTT = NTESTT + 7
*
*           Print out tests which fail.
*
            DO 200 JR = 1, 7
               IF( RESULT( JR ).GE.THRESH ) THEN
*
*                 If this is the first test to fail,
*                 print a header to the data file.
*
                  IF( NERRS.EQ.0 ) THEN
                     WRITE( NOUNIT, FMT = 9997 )'DGV'
*
*                    Matrix types
*
                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )'Orthogonal'
*
*                    Tests performed
*
                     WRITE( NOUNIT, FMT = 9993 )
*
                  END IF
                  NERRS = NERRS + 1
                  IF( RESULT( JR ).LT.10000.0D0 ) THEN
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR,
     $                  RESULT( JR )
                  ELSE
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR,
     $                  RESULT( JR )
                  END IF
               END IF
  200       CONTINUE
*
  210    CONTINUE
  220 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'DGV', NOUNIT, NERRS, NTESTT, 0 )
*
      WORK( 1 ) = MAXWRK
*
      RETURN
*
 9999 FORMAT( ' DDRGEV: ', A, ' returned INFO=', I6, '.', / 3X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 4( I4, ',' ), I5, ')' )
*
 9998 FORMAT( ' DDRGEV: ', A, ' Eigenvectors from ', A, ' incorrectly ',
     $      'normalized.', / ' Bits of error=', 0P, G10.3, ',', 3X,
     $      'N=', I4, ', JTYPE=', I3, ', ISEED=(', 4( I4, ',' ), I5,
     $      ')' )
*
 9997 FORMAT( / 1X, A3, ' -- Real Generalized eigenvalue problem driver'
     $       )
*
 9996 FORMAT( ' Matrix types (see DDRGEV for details): ' )
*
 9995 FORMAT( ' Special Matrices:', 23X,
     $      '(J''=transposed Jordan block)',
     $      / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',
     $      '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ',
     $      'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I',
     $      ')  11=(large*I, small*D)  13=(large*D, large*I)', /
     $      '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ',
     $      ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9994 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:',
     $      / '  16=Transposed Jordan Blocks             19=geometric ',
     $      'alpha, beta=0,1', / '  17=arithm. alpha&beta             ',
     $      '      20=arithmetic alpha, beta=0,1', / '  18=clustered ',
     $      'alpha, beta=0,1            21=random alpha, beta=0,1',
     $      / ' Large & Small Matrices:', / '  22=(large, small)   ',
     $      '23=(small,large)    24=(small,small)    25=(large,large)',
     $      / '  26=random O(1) matrices.' )
*
 9993 FORMAT( / ' Tests performed:    ',
     $      / ' 1 = max | ( b A - a B )''*l | / const.,',
     $      / ' 2 = | |VR(i)| - 1 | / ulp,',
     $      / ' 3 = max | ( b A - a B )*r | / const.',
     $      / ' 4 = | |VL(i)| - 1 | / ulp,',
     $      / ' 5 = 0 if W same no matter if r or l computed,',
     $      / ' 6 = 0 if l same no matter if l computed,',
     $      / ' 7 = 0 if r same no matter if r computed,', / 1X )
 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I2, ' is', 1P, D10.3 )
*
*     End of DDRGEV
*
      END