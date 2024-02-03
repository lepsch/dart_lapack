      SUBROUTINE CCHKGG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, TSTDIF, THRSHN, NOUNIT, A, LDA, B, H, T, S1, S2, P1, P2, U, LDU, V, Q, Z, ALPHA1, BETA1, ALPHA3, BETA3, EVECTL, EVECTR, WORK, LWORK, RWORK, LLWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTDIF;
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH, THRSHN
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * ), LLWORK( * );
      int                ISEED( 4 ), NN( * );
      REAL               RESULT( 15 ), RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA1( * ), ALPHA3( * ), B( LDA, * ), BETA1( * ), BETA3( * ), EVECTL( LDU, * ), EVECTR( LDU, * ), H( LDA, * ), P1( LDA, * ), P2( LDA, * ), Q( LDU, * ), S1( LDA, * ), S2( LDA, * ), T( LDA, * ), U( LDU, * ), V( LDU, * ), WORK( * ), Z( LDU, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 26 )
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      int                I1, IADD, IINFO, IN, J, JC, JR, JSIZE, JTYPE, LWKOPT, MTYPES, N, N1, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ANORM, BNORM, SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV
      COMPLEX            CTEMP
      // ..
      // .. Local Arrays ..
      bool               LASIGN( MAXTYP ), LBSIGN( MAXTYP );
      int                IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      REAL               DUMMA( 4 ), RMAGN( 0: 3 )
      COMPLEX            CDUMMA( 4 )
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      COMPLEX            CLARND
      // EXTERNAL CLANGE, SLAMCH, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQR2, CGET51, CGET52, CGGHRD, CHGEQZ, CLACPY, CLARFG, CLASET, CLATM4, CTGEVC, CUNM2R, SLASUM, XERBLA
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

      LWKOPT = MAX( 2*NMAX*NMAX, 4*NMAX, 1 )

      // Check for errors

      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LE.1 .OR. LDA.LT.NMAX ) THEN
         INFO = -10
      ELSE IF( LDU.LE.1 .OR. LDU.LT.NMAX ) THEN
         INFO = -19
      ELSE IF( LWKOPT.GT.LWORK ) THEN
         INFO = -30
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CCHKGG', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      SAFMIN = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
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

      DO 240 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / REAL( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*N1

         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF

         DO 230 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 230
            NMATS = NMATS + 1
            NTEST = 0

            // Save ISEED in case of an error.

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Initialize RESULT

            DO 30 J = 1, 15
               RESULT( J ) = ZERO
   30       CONTINUE

            // Compute A and B

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
            // RMAGN:  used to implement KAMAGN and KBMAGN.

            IF( MTYPES.GT.MAXTYP ) GO TO 110
            IINFO = 0
            IF( KCLASS( JTYPE ).LT.3 ) THEN

               // Generate A (w/o rotation)

               IF( ABS( KATYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL CLASET( 'Full', N, N, CZERO, CZERO, A, LDA )
               ELSE
                  IN = N
               END IF
               CALL CLATM4( KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 4, ISEED, A, LDA )
               IADD = KADD( KAZERO( JTYPE ) )
               IF( IADD.GT.0 .AND. IADD.LE.N ) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) )

               // Generate B (w/o rotation)

               IF( ABS( KBTYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N ) CALL CLASET( 'Full', N, N, CZERO, CZERO, B, LDA )
               ELSE
                  IN = N
               END IF
               CALL CLATM4( KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 4, ISEED, B, LDA )
               IADD = KADD( KBZERO( JTYPE ) )
               IF( IADD.NE.0 ) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) )

               IF( KCLASS( JTYPE ).EQ.2 .AND. N.GT.0 ) THEN

                  // Include rotations

                  // Generate U, V as Householder transformations times a
                  // diagonal matrix.  (Note that CLARFG makes U(j,j) and
                  // V(j,j) real.)

                  DO 50 JC = 1, N - 1
                     DO 40 JR = JC, N
                        U( JR, JC ) = CLARND( 3, ISEED )
                        V( JR, JC ) = CLARND( 3, ISEED )
   40                CONTINUE
                     CALL CLARFG( N+1-JC, U( JC, JC ), U( JC+1, JC ), 1, WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, REAL( U( JC, JC ) ) )
                     U( JC, JC ) = CONE
                     CALL CLARFG( N+1-JC, V( JC, JC ), V( JC+1, JC ), 1, WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, REAL( V( JC, JC ) ) )
                     V( JC, JC ) = CONE
   50             CONTINUE
                  CTEMP = CLARND( 3, ISEED )
                  U( N, N ) = CONE
                  WORK( N ) = CZERO
                  WORK( 3*N ) = CTEMP / ABS( CTEMP )
                  CTEMP = CLARND( 3, ISEED )
                  V( N, N ) = CONE
                  WORK( 2*N ) = CZERO
                  WORK( 4*N ) = CTEMP / ABS( CTEMP )

                  // Apply the diagonal matrices

                  DO 70 JC = 1, N
                     DO 60 JR = 1, N
                        A( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )* CONJG( WORK( 3*N+JC ) )* B( JR, JC )
   60                CONTINUE
   70             CONTINUE
                  CALL CUNM2R( 'L', 'N', N, N, N-1, U, LDU, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL CUNM2R( 'R', 'C', N, N, N-1, V, LDU, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL CUNM2R( 'L', 'N', N, N, N-1, U, LDU, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100                   CALL CUNM2R( 'R', 'C', N, N, N-1, V, LDU, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO.NE.0 ) GO TO 100
               END IF
            ELSE

               // Random matrices

               DO 90 JC = 1, N
                  DO 80 JR = 1, N
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )* CLARND( 4, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* CLARND( 4, ISEED )
   80             CONTINUE
   90          CONTINUE
            END IF

            ANORM = CLANGE( '1', N, N, A, LDA, RWORK )
            BNORM = CLANGE( '1', N, N, B, LDA, RWORK )

  100       CONTINUE

            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF

  110       CONTINUE

            // Call CGEQR2, CUNM2R, and CGGHRD to compute H, T, U, and V

            CALL CLACPY( ' ', N, N, A, LDA, H, LDA )
            CALL CLACPY( ' ', N, N, B, LDA, T, LDA )
            NTEST = 1
            RESULT( 1 ) = ULPINV

            CALL CGEQR2( N, N, T, LDA, WORK, WORK( N+1 ), IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CGEQR2', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CUNM2R( 'L', 'C', N, N, N, T, LDA, WORK, H, LDA, WORK( N+1 ), IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CUNM2R', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CLASET( 'Full', N, N, CZERO, CONE, U, LDU )
            CALL CUNM2R( 'R', 'N', N, N, N, T, LDA, WORK, U, LDU, WORK( N+1 ), IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CUNM2R', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CGGHRD( 'V', 'I', N, 1, N, H, LDA, T, LDA, U, LDU, V, LDU, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CGGHRD', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF
            NTEST = 4

            // Do tests 1--4

            CALL CGET51( 1, N, A, LDA, H, LDA, U, LDU, V, LDU, WORK, RWORK, RESULT( 1 ) )             CALL CGET51( 1, N, B, LDA, T, LDA, U, LDU, V, LDU, WORK, RWORK, RESULT( 2 ) )             CALL CGET51( 3, N, B, LDA, T, LDA, U, LDU, U, LDU, WORK, RWORK, RESULT( 3 ) )             CALL CGET51( 3, N, B, LDA, T, LDA, V, LDU, V, LDU, WORK, RWORK, RESULT( 4 ) )

            // Call CHGEQZ to compute S1, P1, S2, P2, Q, and Z, do tests.

            // Compute T1 and UZ

            // Eigenvalues only

            CALL CLACPY( ' ', N, N, H, LDA, S2, LDA )
            CALL CLACPY( ' ', N, N, T, LDA, P2, LDA )
            NTEST = 5
            RESULT( 5 ) = ULPINV

            CALL CHGEQZ( 'E', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHA3, BETA3, Q, LDU, Z, LDU, WORK, LWORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CHGEQZ(E)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            // Eigenvalues and Full Schur Form

            CALL CLACPY( ' ', N, N, H, LDA, S2, LDA )
            CALL CLACPY( ' ', N, N, T, LDA, P2, LDA )

            CALL CHGEQZ( 'S', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHA1, BETA1, Q, LDU, Z, LDU, WORK, LWORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CHGEQZ(S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            // Eigenvalues, Schur Form, and Schur Vectors

            CALL CLACPY( ' ', N, N, H, LDA, S1, LDA )
            CALL CLACPY( ' ', N, N, T, LDA, P1, LDA )

            CALL CHGEQZ( 'S', 'I', 'I', N, 1, N, S1, LDA, P1, LDA, ALPHA1, BETA1, Q, LDU, Z, LDU, WORK, LWORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CHGEQZ(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            NTEST = 8

            // Do Tests 5--8

            CALL CGET51( 1, N, H, LDA, S1, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 5 ) )             CALL CGET51( 1, N, T, LDA, P1, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 6 ) )             CALL CGET51( 3, N, T, LDA, P1, LDA, Q, LDU, Q, LDU, WORK, RWORK, RESULT( 7 ) )             CALL CGET51( 3, N, T, LDA, P1, LDA, Z, LDU, Z, LDU, WORK, RWORK, RESULT( 8 ) )

            // Compute the Left and Right Eigenvectors of (S1,P1)

            // 9: Compute the left eigenvector Matrix without
               // back transforming:

            NTEST = 9
            RESULT( 9 ) = ULPINV

            // To test "SELECT" option, compute half of the eigenvectors
            // in one call, and half in another

            I1 = N / 2
            DO 120 J = 1, I1
               LLWORK( J ) = .TRUE.
  120       CONTINUE
            DO 130 J = I1 + 1, N
               LLWORK( J ) = .FALSE.
  130       CONTINUE

            CALL CTGEVC( 'L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CTGEVC(L,S1)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            I1 = IN
            DO 140 J = 1, I1
               LLWORK( J ) = .FALSE.
  140       CONTINUE
            DO 150 J = I1 + 1, N
               LLWORK( J ) = .TRUE.
  150       CONTINUE

            CALL CTGEVC( 'L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL( 1, I1+1 ), LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CTGEVC(L,S2)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CGET52( .TRUE., N, S1, LDA, P1, LDA, EVECTL, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) )
            RESULT( 9 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRSHN ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'CTGEVC(HOWMNY=S)', DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF

            // 10: Compute the left eigenvector Matrix with
                // back transforming:

            NTEST = 10
            RESULT( 10 ) = ULPINV
            CALL CLACPY( 'F', N, N, Q, LDU, EVECTL, LDU )
            CALL CTGEVC( 'L', 'B', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CTGEVC(L,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CGET52( .TRUE., N, H, LDA, T, LDA, EVECTL, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) )
            RESULT( 10 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRSHN ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'CTGEVC(HOWMNY=B)', DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF

            // 11: Compute the right eigenvector Matrix without
                // back transforming:

            NTEST = 11
            RESULT( 11 ) = ULPINV

            // To test "SELECT" option, compute half of the eigenvectors
            // in one call, and half in another

            I1 = N / 2
            DO 160 J = 1, I1
               LLWORK( J ) = .TRUE.
  160       CONTINUE
            DO 170 J = I1 + 1, N
               LLWORK( J ) = .FALSE.
  170       CONTINUE

            CALL CTGEVC( 'R', 'S', LLWORK, N, S1, LDA, P1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CTGEVC(R,S1)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            I1 = IN
            DO 180 J = 1, I1
               LLWORK( J ) = .FALSE.
  180       CONTINUE
            DO 190 J = I1 + 1, N
               LLWORK( J ) = .TRUE.
  190       CONTINUE

            CALL CTGEVC( 'R', 'S', LLWORK, N, S1, LDA, P1, LDA, CDUMMA, LDU, EVECTR( 1, I1+1 ), LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CTGEVC(R,S2)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CGET52( .FALSE., N, S1, LDA, P1, LDA, EVECTR, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) )
            RESULT( 11 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'CTGEVC(HOWMNY=S)', DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF

            // 12: Compute the right eigenvector Matrix with
                // back transforming:

            NTEST = 12
            RESULT( 12 ) = ULPINV
            CALL CLACPY( 'F', N, N, Z, LDU, EVECTR, LDU )
            CALL CTGEVC( 'R', 'B', LLWORK, N, S1, LDA, P1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'CTGEVC(R,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 210
            END IF

            CALL CGET52( .FALSE., N, H, LDA, T, LDA, EVECTR, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) )
            RESULT( 12 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'CTGEVC(HOWMNY=B)', DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF

            // Tests 13--15 are done only on request

            IF( TSTDIF ) THEN

               // Do Tests 13--14

               CALL CGET51( 2, N, S1, LDA, S2, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 13 ) )                CALL CGET51( 2, N, P1, LDA, P2, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 14 ) )

               // Do Test 15

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 200 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( ALPHA1( J )-ALPHA3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( BETA1( J )-BETA3( J ) ) )
  200          CONTINUE

               TEMP1 = TEMP1 / MAX( SAFMIN, ULP*MAX( TEMP1, ANORM ) )
               TEMP2 = TEMP2 / MAX( SAFMIN, ULP*MAX( TEMP2, BNORM ) )
               RESULT( 15 ) = MAX( TEMP1, TEMP2 )
               NTEST = 15
            ELSE
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               RESULT( 15 ) = ZERO
               NTEST = 12
            END IF

            // End of Loop -- Check for RESULT(j) > THRESH

  210       CONTINUE

            NTESTT = NTESTT + NTEST

            // Print out tests which fail.

            DO 220 JR = 1, NTEST
               IF( RESULT( JR ).GE.THRESH ) THEN

                  // If this is the first test to fail,
                  // print a header to the data file.

                  IF( NERRS.EQ.0 ) THEN
                     WRITE( NOUNIT, FMT = 9997 )'CGG'

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )
                     WRITE( NOUNIT, FMT = 9994 )'Unitary'

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 )'unitary', '*', 'conjugate transpose', ( '*', J = 1, 10 )

                  END IF
                  NERRS = NERRS + 1
                  IF( RESULT( JR ).LT.10000.0 ) THEN
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  ELSE
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  END IF
               END IF
  220       CONTINUE

  230    CONTINUE
  240 CONTINUE

      // Summary

      CALL SLASUM( 'CGG', NOUNIT, NERRS, NTESTT )
      RETURN

 9999 FORMAT( ' CCHKGG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9998 FORMAT( ' CCHKGG: ', A, ' Eigenvectors from ', A, ' incorrectly ',
     $      'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X,
     $      'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5,
     $      ')' )

 9997 FORMAT( 1X, A3, ' -- Complex Generalized eigenvalue problem' )

 9996 FORMAT( ' Matrix types (see CCHKGG for details): ' )

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

 9993 FORMAT( / ' Tests performed:   (H is Hessenberg, S is Schur, B, ',
     $      'T, P are triangular,', / 20X, 'U, V, Q, and Z are ', A,
     $      ', l and r are the', / 20X,
     $      'appropriate left and right eigenvectors, resp., a is',
     $      / 20X, 'alpha, b is beta, and ', A, ' means ', A, '.)',
     $      / ' 1 = | A - U H V', A,
     $      ' | / ( |A| n ulp )      2 = | B - U T V', A,
     $      ' | / ( |B| n ulp )', / ' 3 = | I - UU', A,
     $      ' | / ( n ulp )             4 = | I - VV', A,
     $      ' | / ( n ulp )', / ' 5 = | H - Q S Z', A,
     $      ' | / ( |H| n ulp )', 6X, '6 = | T - Q P Z', A,
     $      ' | / ( |T| n ulp )', / ' 7 = | I - QQ', A,
     $      ' | / ( n ulp )             8 = | I - ZZ', A,
     $      ' | / ( n ulp )', / ' 9 = max | ( b S - a P )', A,
     $      ' l | / const.  10 = max | ( b H - a T )', A,
     $      ' l | / const.', /
     $      ' 11= max | ( b S - a P ) r | / const.   12 = max | ( b H',
     $      ' - a T ) r | / const.', / 1X )

 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 )
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I2, ' is', 1P, E10.3 )

      // End of CCHKGG

      END
