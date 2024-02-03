      SUBROUTINE ZCHKHS( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, H, T1, T2, U, LDU, Z, UZ, W1, W3, EVECTL, EVECTR, EVECTY, EVECTX, UU, TAU, WORK, NWORK, RWORK, IWORK, SELECT, RESULT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      bool               DOTYPE( * ), SELECT( * );
      int                ISEED( 4 ), IWORK( * ), NN( * )
      DOUBLE PRECISION   RESULT( 16 ), RWORK( * )
      COMPLEX*16         A( LDA, * ), EVECTL( LDU, * ), EVECTR( LDU, * ), EVECTX( LDU, * ), EVECTY( LDU, * ), H( LDA, * ), T1( LDA, * ), T2( LDA, * ), TAU( * ), U( LDU, * ), UU( LDU, * ), UZ( LDU, * ), W1( * ), W3( * ), WORK( * ), Z( LDU, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      int                MAXTYP
      PARAMETER          ( MAXTYP = 21 )
*     ..
*     .. Local Scalars ..
      bool               BADNN, MATCH;
      int                I, IHI, IINFO, ILO, IMODE, IN, ITYPE, J, JCOL, JJ, JSIZE, JTYPE, K, MTYPES, N, N1, NERRS, NMATS, NMAX, NTEST, NTESTT
      DOUBLE PRECISION   ANINV, ANORM, COND, CONDS, OVFL, RTOVFL, RTULP, RTULPI, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL
*     ..
*     .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP )
      DOUBLE PRECISION   DUMMA( 4 )
      COMPLEX*16         CDUMMA( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAFTS, DLASUM, XERBLA, ZCOPY, ZGEHRD, ZGEMM, ZGET10, ZGET22, ZHSEIN, ZHSEQR, ZHST01, ZLACPY, ZLASET, ZLATME, ZLATMR, ZLATMS, ZTREVC, ZTREVC3, ZUNGHR, ZUNMHR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
      DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
      DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      NTESTT = 0
      INFO = 0
*
      BADNN = .FALSE.
      NMAX = 0
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE
*
*     Check for errors
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
      ELSE IF( LDU.LE.1 .OR. LDU.LT.NMAX ) THEN
         INFO = -14
      ELSE IF( 4*NMAX*NMAX+2.GT.NWORK ) THEN
         INFO = -26
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZCHKHS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN
*
*     More important constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
      RTULP = SQRT( ULP )
      RTULPI = ONE / RTULP
*
*     Loop over sizes, types
*
      NERRS = 0
      NMATS = 0
*
      DO 260 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( N.EQ.0 ) GO TO 260
         N1 = MAX( 1, N )
         ANINV = ONE / DBLE( N1 )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 250 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 250
            NMATS = NMATS + 1
            NTEST = 0
*
*           Save ISEED in case of an error.
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
*           Initialize RESULT
*
            DO 30 J = 1, 14
               RESULT( J ) = ZERO
   30       CONTINUE
*
*           Compute "A"
*
*           Control parameters:
*
*           KMAGN  KCONDS  KMODE        KTYPE
*       =1  O(1)   1       clustered 1  zero
*       =2  large  large   clustered 2  identity
*       =3  small          exponential  Jordan
*       =4                 arithmetic   diagonal, (w/ eigenvalues)
*       =5                 random log   hermitian, w/ eigenvalues
*       =6                 random       general, w/ eigenvalues
*       =7                              random diagonal
*       =8                              random hermitian
*       =9                              random general
*       =10                             random triangular
*
            IF( MTYPES.GT.MAXTYP ) GO TO 100
*
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
*
*           Compute norm
*
            GO TO ( 40, 50, 60 )KMAGN( JTYPE )
*
   40       CONTINUE
            ANORM = ONE
            GO TO 70
*
   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70
*
   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70
*
   70       CONTINUE
*
            CALL ZLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
            IINFO = 0
            COND = ULPINV
*
*           Special Matrices
*
            IF( ITYPE.EQ.1 ) THEN
*
*              Zero
*
               IINFO = 0
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE
*
            ELSE IF( ITYPE.EQ.3 ) THEN
*
*              Jordan Block
*
               DO 90 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
                  IF( JCOL.GT.1 ) A( JCOL, JCOL-1 ) = ONE
   90          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              Diagonal Matrix, [Eigen]values Specified
*
               CALL ZLATMR( N, N, 'D', ISEED, 'N', WORK, IMODE, COND, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Hermitian, eigenvalues specified
*
               CALL ZLATMS( N, N, 'D', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO )
*
            ELSE IF( ITYPE.EQ.6 ) THEN
*
*              General, eigenvalues specified
*
               IF( KCONDS( JTYPE ).EQ.1 ) THEN
                  CONDS = ONE
               ELSE IF( KCONDS( JTYPE ).EQ.2 ) THEN
                  CONDS = RTULPI
               ELSE
                  CONDS = ZERO
               END IF
*
               CALL ZLATME( N, 'D', ISEED, WORK, IMODE, COND, CONE, 'T', 'T', 'T', RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK( N+1 ), IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              Diagonal, random eigenvalues
*
               CALL ZLATMR( N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              Hermitian, random eigenvalues
*
               CALL ZLATMR( N, N, 'D', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              General, random eigenvalues
*
               CALL ZLATMR( N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.10 ) THEN
*
*              Triangular, random eigenvalues
*
               CALL ZLATMR( N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE
*
               IINFO = 1
            END IF
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
  100       CONTINUE
*
*           Call ZGEHRD to compute H and U, do tests.
*
            CALL ZLACPY( ' ', N, N, A, LDA, H, LDA )
            NTEST = 1
*
            ILO = 1
            IHI = N
*
            CALL ZGEHRD( N, ILO, IHI, H, LDA, WORK, WORK( N+1 ), NWORK-N, IINFO )
*
            IF( IINFO.NE.0 ) THEN
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'ZGEHRD', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
            DO 120 J = 1, N - 1
               UU( J+1, J ) = CZERO
               DO 110 I = J + 2, N
                  U( I, J ) = H( I, J )
                  UU( I, J ) = H( I, J )
                  H( I, J ) = CZERO
  110          CONTINUE
  120       CONTINUE
            CALL ZCOPY( N-1, WORK, 1, TAU, 1 )
            CALL ZUNGHR( N, ILO, IHI, U, LDU, WORK, WORK( N+1 ), NWORK-N, IINFO )
            NTEST = 2
*
            CALL ZHST01( N, ILO, IHI, A, LDA, H, LDA, U, LDU, WORK, NWORK, RWORK, RESULT( 1 ) )
*
*           Call ZHSEQR to compute T1, T2 and Z, do tests.
*
*           Eigenvalues only (W3)
*
            CALL ZLACPY( ' ', N, N, H, LDA, T2, LDA )
            NTEST = 3
            RESULT( 3 ) = ULPINV
*
            CALL ZHSEQR( 'E', 'N', N, ILO, IHI, T2, LDA, W3, UZ, LDU, WORK, NWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHSEQR(E)', IINFO, N, JTYPE, IOLDSD
               IF( IINFO.LE.N+2 ) THEN
                  INFO = ABS( IINFO )
                  GO TO 240
               END IF
            END IF
*
*           Eigenvalues (W1) and Full Schur Form (T2)
*
            CALL ZLACPY( ' ', N, N, H, LDA, T2, LDA )
*
            CALL ZHSEQR( 'S', 'N', N, ILO, IHI, T2, LDA, W1, UZ, LDU, WORK, NWORK, IINFO )
            IF( IINFO.NE.0 .AND. IINFO.LE.N+2 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHSEQR(S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
*           Eigenvalues (W1), Schur Form (T1), and Schur Vectors (UZ)
*
            CALL ZLACPY( ' ', N, N, H, LDA, T1, LDA )
            CALL ZLACPY( ' ', N, N, U, LDU, UZ, LDU )
*
            CALL ZHSEQR( 'S', 'V', N, ILO, IHI, T1, LDA, W1, UZ, LDU, WORK, NWORK, IINFO )
            IF( IINFO.NE.0 .AND. IINFO.LE.N+2 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHSEQR(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
*           Compute Z = U' UZ
*
            CALL ZGEMM( 'C', 'N', N, N, N, CONE, U, LDU, UZ, LDU, CZERO, Z, LDU )
            NTEST = 8
*
*           Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
*                and 4: | I - Z Z' | / ( n ulp )
*
            CALL ZHST01( N, ILO, IHI, H, LDA, T1, LDA, Z, LDU, WORK, NWORK, RWORK, RESULT( 3 ) )
*
*           Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
*                and 6: | I - UZ (UZ)' | / ( n ulp )
*
            CALL ZHST01( N, ILO, IHI, A, LDA, T1, LDA, UZ, LDU, WORK, NWORK, RWORK, RESULT( 5 ) )
*
*           Do Test 7: | T2 - T1 | / ( |T| n ulp )
*
            CALL ZGET10( N, N, T2, LDA, T1, LDA, WORK, RWORK, RESULT( 7 ) )
*
*           Do Test 8: | W3 - W1 | / ( max(|W1|,|W3|) ulp )
*
            TEMP1 = ZERO
            TEMP2 = ZERO
            DO 130 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( W1( J ) ), ABS( W3( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( W1( J )-W3( J ) ) )
  130       CONTINUE
*
            RESULT( 8 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
*           Compute the Left and Right Eigenvectors of T
*
*           Compute the Right eigenvector Matrix:
*
            NTEST = 9
            RESULT( 9 ) = ULPINV
*
*           Select every other eigenvector
*
            DO 140 J = 1, N
               SELECT( J ) = .FALSE.
  140       CONTINUE
            DO 150 J = 1, N, 2
               SELECT( J ) = .TRUE.
  150       CONTINUE
            CALL ZTREVC( 'Right', 'All', SELECT, N, T1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(R,A)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
*           Test 9:  | TR - RW | / ( |T| |R| ulp )
*
            CALL ZGET22( 'N', 'N', 'N', N, T1, LDA, EVECTR, LDU, W1, WORK, RWORK, DUMMA( 1 ) )
            RESULT( 9 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZTREVC', DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF
*
*           Compute selected right eigenvectors and confirm that
*           they agree with previous right eigenvectors
*
            CALL ZTREVC( 'Right', 'Some', SELECT, N, T1, LDA, CDUMMA, LDU, EVECTL, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(R,S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
            K = 1
            MATCH = .TRUE.
            DO 170 J = 1, N
               IF( SELECT( J ) ) THEN
                  DO 160 JJ = 1, N
                     IF( EVECTR( JJ, J ).NE.EVECTL( JJ, K ) ) THEN
                        MATCH = .FALSE.
                        GO TO 180
                     END IF
  160             CONTINUE
                  K = K + 1
               END IF
  170       CONTINUE
  180       CONTINUE
            IF( .NOT.MATCH ) WRITE( NOUNIT, FMT = 9997 )'Right', 'ZTREVC', N, JTYPE, IOLDSD
*
*           Compute the Left eigenvector Matrix:
*
            NTEST = 10
            RESULT( 10 ) = ULPINV
            CALL ZTREVC( 'Left', 'All', SELECT, N, T1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(L,A)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
*           Test 10:  | LT - WL | / ( |T| |L| ulp )
*
            CALL ZGET22( 'C', 'N', 'C', N, T1, LDA, EVECTL, LDU, W1, WORK, RWORK, DUMMA( 3 ) )
            RESULT( 10 ) = DUMMA( 3 )
            IF( DUMMA( 4 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZTREVC', DUMMA( 4 ), N, JTYPE, IOLDSD
            END IF
*
*           Compute selected left eigenvectors and confirm that
*           they agree with previous left eigenvectors
*
            CALL ZTREVC( 'Left', 'Some', SELECT, N, T1, LDA, EVECTR, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC(L,S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 240
            END IF
*
            K = 1
            MATCH = .TRUE.
            DO 200 J = 1, N
               IF( SELECT( J ) ) THEN
                  DO 190 JJ = 1, N
                     IF( EVECTL( JJ, J ).NE.EVECTR( JJ, K ) ) THEN
                        MATCH = .FALSE.
                        GO TO 210
                     END IF
  190             CONTINUE
                  K = K + 1
               END IF
  200       CONTINUE
  210       CONTINUE
            IF( .NOT.MATCH ) WRITE( NOUNIT, FMT = 9997 )'Left', 'ZTREVC', N, JTYPE, IOLDSD
*
*           Call ZHSEIN for Right eigenvectors of H, do test 11
*
            NTEST = 11
            RESULT( 11 ) = ULPINV
            DO 220 J = 1, N
               SELECT( J ) = .TRUE.
  220       CONTINUE
*
            CALL ZHSEIN( 'Right', 'Qr', 'Ninitv', SELECT, N, H, LDA, W3, CDUMMA, LDU, EVECTX, LDU, N1, IN, WORK, RWORK, IWORK, IWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHSEIN(R)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 240
            ELSE
*
*              Test 11:  | HX - XW | / ( |H| |X| ulp )
*
*                        (from inverse iteration)
*
               CALL ZGET22( 'N', 'N', 'N', N, H, LDA, EVECTX, LDU, W3, WORK, RWORK, DUMMA( 1 ) )                IF( DUMMA( 1 ).LT.ULPINV ) RESULT( 11 ) = DUMMA( 1 )*ANINV
               IF( DUMMA( 2 ).GT.THRESH ) THEN
                  WRITE( NOUNIT, FMT = 9998 )'Right', 'ZHSEIN', DUMMA( 2 ), N, JTYPE, IOLDSD
               END IF
            END IF
*
*           Call ZHSEIN for Left eigenvectors of H, do test 12
*
            NTEST = 12
            RESULT( 12 ) = ULPINV
            DO 230 J = 1, N
               SELECT( J ) = .TRUE.
  230       CONTINUE
*
            CALL ZHSEIN( 'Left', 'Qr', 'Ninitv', SELECT, N, H, LDA, W3, EVECTY, LDU, CDUMMA, LDU, N1, IN, WORK, RWORK, IWORK, IWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZHSEIN(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 240
            ELSE
*
*              Test 12:  | YH - WY | / ( |H| |Y| ulp )
*
*                        (from inverse iteration)
*
               CALL ZGET22( 'C', 'N', 'C', N, H, LDA, EVECTY, LDU, W3, WORK, RWORK, DUMMA( 3 ) )                IF( DUMMA( 3 ).LT.ULPINV ) RESULT( 12 ) = DUMMA( 3 )*ANINV
               IF( DUMMA( 4 ).GT.THRESH ) THEN
                  WRITE( NOUNIT, FMT = 9998 )'Left', 'ZHSEIN', DUMMA( 4 ), N, JTYPE, IOLDSD
               END IF
            END IF
*
*           Call ZUNMHR for Right eigenvectors of A, do test 13
*
            NTEST = 13
            RESULT( 13 ) = ULPINV
*
            CALL ZUNMHR( 'Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTX, LDU, WORK, NWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZUNMHR(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 240
            ELSE
*
*              Test 13:  | AX - XW | / ( |A| |X| ulp )
*
*                        (from inverse iteration)
*
               CALL ZGET22( 'N', 'N', 'N', N, A, LDA, EVECTX, LDU, W3, WORK, RWORK, DUMMA( 1 ) )                IF( DUMMA( 1 ).LT.ULPINV ) RESULT( 13 ) = DUMMA( 1 )*ANINV
            END IF
*
*           Call ZUNMHR for Left eigenvectors of A, do test 14
*
            NTEST = 14
            RESULT( 14 ) = ULPINV
*
            CALL ZUNMHR( 'Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTY, LDU, WORK, NWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZUNMHR(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 240
            ELSE
*
*              Test 14:  | YA - WY | / ( |A| |Y| ulp )
*
*                        (from inverse iteration)
*
               CALL ZGET22( 'C', 'N', 'C', N, A, LDA, EVECTY, LDU, W3, WORK, RWORK, DUMMA( 3 ) )                IF( DUMMA( 3 ).LT.ULPINV ) RESULT( 14 ) = DUMMA( 3 )*ANINV
            END IF
*
*           Compute Left and Right Eigenvectors of A
*
*           Compute a Right eigenvector matrix:
*
            NTEST = 15
            RESULT( 15 ) = ULPINV
*
            CALL ZLACPY( ' ', N, N, UZ, LDU, EVECTR, LDU )
*
            CALL ZTREVC3( 'Right', 'Back', SELECT, N, T1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, NWORK, RWORK, N, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC3(R,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            END IF
*
*           Test 15:  | AR - RW | / ( |A| |R| ulp )
*
*                     (from Schur decomposition)
*
            CALL ZGET22( 'N', 'N', 'N', N, A, LDA, EVECTR, LDU, W1, WORK, RWORK, DUMMA( 1 ) )
            RESULT( 15 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZTREVC3', DUMMA( 2 ), N, JTYPE, IOLDSD
            END IF
*
*           Compute a Left eigenvector matrix:
*
            NTEST = 16
            RESULT( 16 ) = ULPINV
*
            CALL ZLACPY( ' ', N, N, UZ, LDU, EVECTL, LDU )
*
            CALL ZTREVC3( 'Left', 'Back', SELECT, N, T1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, NWORK, RWORK, N, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'ZTREVC3(L,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            END IF
*
*           Test 16:  | LA - WL | / ( |A| |L| ulp )
*
*                     (from Schur decomposition)
*
            CALL ZGET22( 'Conj', 'N', 'Conj', N, A, LDA, EVECTL, LDU, W1, WORK, RWORK, DUMMA( 3 ) )
            RESULT( 16 ) = DUMMA( 3 )
            IF( DUMMA( 4 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZTREVC3', DUMMA( 4 ), N, JTYPE, IOLDSD
            END IF
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
  240       CONTINUE
*
            NTESTT = NTESTT + NTEST
            CALL DLAFTS( 'ZHS', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS )
*
  250    CONTINUE
  260 CONTINUE
*
*     Summary
*
      CALL DLASUM( 'ZHS', NOUNIT, NERRS, NTESTT )
*
      RETURN
*
 9999 FORMAT( ' ZCHKHS: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' ZCHKHS: ', A, ' Eigenvectors from ', A, ' incorrectly ',
     $      'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X,
     $      'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5,
     $      ')' )
 9997 FORMAT( ' ZCHKHS: Selected ', A, ' Eigenvectors from ', A,
     $      ' do not match other eigenvectors ', 9X, 'N=', I6,
     $      ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
*     End of ZCHKHS
*
      END
