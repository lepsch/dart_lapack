      SUBROUTINE CDRVEV( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, H, W, W1, VL, LDVL, VR, LDVR, LRE, LDLRE, RESULT, WORK, NWORK, RWORK, IWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDLRE, LDVL, LDVR, NOUNIT, NSIZES, NTYPES, NWORK;
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               RESULT( 7 ), RWORK( * )
      COMPLEX            A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), W1( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 21 )
*     ..
*     .. Local Scalars ..
      bool               BADNN;
      String             PATH;
      int                IINFO, IMODE, ITYPE, IWK, J, JCOL, JJ, JSIZE, JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NTEST, NTESTF, NTESTT;
      REAL               ANORM, COND, CONDS, OVFL, RTULP, RTULPI, TNRM, ULP, ULPINV, UNFL, VMX, VRMX, VTST
*     ..
*     .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      REAL               RES( 2 )
      COMPLEX            DUM( 1 )
*     ..
*     .. External Functions ..
      REAL               SCNRM2, SLAMCH
      // EXTERNAL SCNRM2, SLAMCH
*     ..
*     .. External Subroutines ..
      // EXTERNAL CGEEV, CGET22, CLACPY, CLATME, CLATMR, CLATMS, CLASET, SLASUM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, MAX, MIN, REAL, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
      DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
      DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
*     ..
*     .. Executable Statements ..
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'EV'
*
*     Check for errors
*
      NTESTT = 0
      NTESTF = 0
      INFO = 0
*
*     Important constants
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
      ELSE IF( NOUNIT.LE.0 ) THEN
         INFO = -7
      ELSE IF( LDA.LT.1 .OR. LDA.LT.NMAX ) THEN
         INFO = -9
      ELSE IF( LDVL.LT.1 .OR. LDVL.LT.NMAX ) THEN
         INFO = -14
      ELSE IF( LDVR.LT.1 .OR. LDVR.LT.NMAX ) THEN
         INFO = -16
      ELSE IF( LDLRE.LT.1 .OR. LDLRE.LT.NMAX ) THEN
         INFO = -28
      ELSE IF( 5*NMAX+2*NMAX**2.GT.NWORK ) THEN
         INFO = -21
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CDRVEV', -INFO )
         RETURN
      END IF
*
*     Quick return if nothing to do
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN
*
*     More Important constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      RTULP = SQRT( ULP )
      RTULPI = ONE / RTULP
*
*     Loop over sizes, types
*
      NERRS = 0
*
      DO 270 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 260 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 260
*
*           Save ISEED in case of an error.
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
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
*       =5                 random log   symmetric, w/ eigenvalues
*       =6                 random       general, w/ eigenvalues
*       =7                              random diagonal
*       =8                              random symmetric
*       =9                              random general
*       =10                             random triangular
*
            IF( MTYPES.GT.MAXTYP ) GO TO 90
*
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
*
*           Compute norm
*
            GO TO ( 30, 40, 50 )KMAGN( JTYPE )
*
   30       CONTINUE
            ANORM = ONE
            GO TO 60
*
   40       CONTINUE
            ANORM = OVFL*ULP
            GO TO 60
*
   50       CONTINUE
            ANORM = UNFL*ULPINV
            GO TO 60
*
   60       CONTINUE
*
            CALL CLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
            IINFO = 0
            COND = ULPINV
*
*           Special Matrices -- Identity & Jordan block
*
*              Zero
*
            IF( ITYPE.EQ.1 ) THEN
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 70 JCOL = 1, N
                  A( JCOL, JCOL ) = CMPLX( ANORM )
   70          CONTINUE
*
            ELSE IF( ITYPE.EQ.3 ) THEN
*
*              Jordan Block
*
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = CMPLX( ANORM )
                  IF( JCOL.GT.1 ) A( JCOL, JCOL-1 ) = CONE
   80          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              Diagonal Matrix, [Eigen]values Specified
*
               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Hermitian, eigenvalues specified
*
               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO )
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
               CALL CLATME( N, 'D', ISEED, WORK, IMODE, COND, CONE, 'T', 'T', 'T', RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              Diagonal, random eigenvalues
*
               CALL CLATMR( N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              Symmetric, random eigenvalues
*
               CALL CLATMR( N, N, 'D', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              General, random eigenvalues
*
               CALL CLATMR( N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
               IF( N.GE.4 ) THEN
                  CALL CLASET( 'Full', 2, N, CZERO, CZERO, A, LDA )
                  CALL CLASET( 'Full', N-3, 1, CZERO, CZERO, A( 3, 1 ), LDA )                   CALL CLASET( 'Full', N-3, 2, CZERO, CZERO, A( 3, N-1 ), LDA )                   CALL CLASET( 'Full', 1, N, CZERO, CZERO, A( N, 1 ), LDA )
               END IF
*
            ELSE IF( ITYPE.EQ.10 ) THEN
*
*              Triangular, random eigenvalues
*
               CALL CLATMR( N, N, 'D', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE
*
               IINFO = 1
            END IF
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9993 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
   90       CONTINUE
*
*           Test for minimal and generous workspace
*
            DO 250 IWK = 1, 2
               IF( IWK.EQ.1 ) THEN
                  NNWORK = 2*N
               ELSE
                  NNWORK = 5*N + 2*N**2
               END IF
               NNWORK = MAX( NNWORK, 1 )
*
*              Initialize RESULT
*
               DO 100 J = 1, 7
                  RESULT( J ) = -ONE
  100          CONTINUE
*
*              Compute eigenvalues and eigenvectors, and test them
*
               CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
               CALL CGEEV( 'V', 'V', N, H, LDA, W, VL, LDVL, VR, LDVR, WORK, NNWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  RESULT( 1 ) = ULPINV
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV1', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 220
               END IF
*
*              Do Test (1)
*
               CALL CGET22( 'N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES )
               RESULT( 1 ) = RES( 1 )
*
*              Do Test (2)
*
               CALL CGET22( 'C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES )
               RESULT( 2 ) = RES( 1 )
*
*              Do Test (3)
*
               DO 120 J = 1, N
                  TNRM = SCNRM2( N, VR( 1, J ), 1 )
                  RESULT( 3 ) = MAX( RESULT( 3 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
                  VMX = ZERO
                  VRMX = ZERO
                  DO 110 JJ = 1, N
                     VTST = ABS( VR( JJ, J ) )
                     IF( VTST.GT.VMX ) VMX = VTST                      IF( AIMAG( VR( JJ, J ) ).EQ.ZERO .AND. ABS( REAL( VR( JJ, J ) ) ).GT.VRMX ) VRMX = ABS( REAL( VR( JJ, J ) ) )
  110             CONTINUE
                  IF( VRMX / VMX.LT.ONE-TWO*ULP ) RESULT( 3 ) = ULPINV
  120          CONTINUE
*
*              Do Test (4)
*
               DO 140 J = 1, N
                  TNRM = SCNRM2( N, VL( 1, J ), 1 )
                  RESULT( 4 ) = MAX( RESULT( 4 ), MIN( ULPINV, ABS( TNRM-ONE ) / ULP ) )
                  VMX = ZERO
                  VRMX = ZERO
                  DO 130 JJ = 1, N
                     VTST = ABS( VL( JJ, J ) )
                     IF( VTST.GT.VMX ) VMX = VTST                      IF( AIMAG( VL( JJ, J ) ).EQ.ZERO .AND. ABS( REAL( VL( JJ, J ) ) ).GT.VRMX ) VRMX = ABS( REAL( VL( JJ, J ) ) )
  130             CONTINUE
                  IF( VRMX / VMX.LT.ONE-TWO*ULP ) RESULT( 4 ) = ULPINV
  140          CONTINUE
*
*              Compute eigenvalues only, and test them
*
               CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
               CALL CGEEV( 'N', 'N', N, H, LDA, W1, DUM, 1, DUM, 1, WORK, NNWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  RESULT( 1 ) = ULPINV
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV2', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 220
               END IF
*
*              Do Test (5)
*
               DO 150 J = 1, N
                  IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
  150          CONTINUE
*
*              Compute eigenvalues and right eigenvectors, and test them
*
               CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
               CALL CGEEV( 'N', 'V', N, H, LDA, W1, DUM, 1, LRE, LDLRE, WORK, NNWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  RESULT( 1 ) = ULPINV
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV3', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 220
               END IF
*
*              Do Test (5) again
*
               DO 160 J = 1, N
                  IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
  160          CONTINUE
*
*              Do Test (6)
*
               DO 180 J = 1, N
                  DO 170 JJ = 1, N
                     IF( VR( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 6 ) = ULPINV
  170             CONTINUE
  180          CONTINUE
*
*              Compute eigenvalues and left eigenvectors, and test them
*
               CALL CLACPY( 'F', N, N, A, LDA, H, LDA )
               CALL CGEEV( 'V', 'N', N, H, LDA, W1, LRE, LDLRE, DUM, 1, WORK, NNWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  RESULT( 1 ) = ULPINV
                  WRITE( NOUNIT, FMT = 9993 )'CGEEV4', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 220
               END IF
*
*              Do Test (5) again
*
               DO 190 J = 1, N
                  IF( W( J ).NE.W1( J ) ) RESULT( 5 ) = ULPINV
  190          CONTINUE
*
*              Do Test (7)
*
               DO 210 J = 1, N
                  DO 200 JJ = 1, N
                     IF( VL( J, JJ ).NE.LRE( J, JJ ) ) RESULT( 7 ) = ULPINV
  200             CONTINUE
  210          CONTINUE
*
*              End of Loop -- Check for RESULT(j) > THRESH
*
  220          CONTINUE
*
               NTEST = 0
               NFAIL = 0
               DO 230 J = 1, 7
                  IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1                   IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  230          CONTINUE
*
               IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
               IF( NTESTF.EQ.1 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )PATH
                  WRITE( NOUNIT, FMT = 9998 )
                  WRITE( NOUNIT, FMT = 9997 )
                  WRITE( NOUNIT, FMT = 9996 )
                  WRITE( NOUNIT, FMT = 9995 )THRESH
                  NTESTF = 2
               END IF
*
               DO 240 J = 1, 7
                  IF( RESULT( J ).GE.THRESH ) THEN
                     WRITE( NOUNIT, FMT = 9994 )N, IWK, IOLDSD, JTYPE, J, RESULT( J )
                  END IF
  240          CONTINUE
*
               NERRS = NERRS + NFAIL
               NTESTT = NTESTT + NTEST
*
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
*
*     Summary
*
      CALL SLASUM( PATH, NOUNIT, NERRS, NTESTT )
*
 9999 FORMAT( / 1X, A3, ' -- Complex Eigenvalue-Eigenvector ',
     $      'Decomposition Driver', /
     $      ' Matrix types (see CDRVEV for details): ' )
*
 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ',
     $      '           ', '  5=Diagonal: geometr. spaced entries.',
     $      / '  2=Identity matrix.                    ', '  6=Diagona',
     $      'l: clustered entries.', / '  3=Transposed Jordan block.  ',
     $      '          ', '  7=Diagonal: large, evenly spaced.', / '  ',
     $      '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s',
     $      'mall, evenly spaced.' )
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev',
     $      'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e',
     $      'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ',
     $      ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond',
     $      'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp',
     $      'lex ', A6, / ' 12=Well-cond., random complex ', A6, '   ',
     $      ' 17=Ill-cond., large rand. complx ', A4, / ' 13=Ill-condi',
     $      'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.',
     $      ' complx ', A4 )
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ',
     $      'with small random entries.', / ' 20=Matrix with large ran',
     $      'dom entries.   ', / )
 9995 FORMAT( ' Tests performed with test threshold =', F8.2,
     $      / / ' 1 = | A VR - VR W | / ( n |A| ulp ) ',
     $      / ' 2 = | conj-trans(A) VL - VL conj-trans(W) | /',
     $      ' ( n |A| ulp ) ', / ' 3 = | |VR(i)| - 1 | / ulp ',
     $      / ' 4 = | |VL(i)| - 1 | / ulp ',
     $      / ' 5 = 0 if W same no matter if VR or VL computed,',
     $      ' 1/ulp otherwise', /
     $      ' 6 = 0 if VR same no matter if VL computed,',
     $      '  1/ulp otherwise', /
     $      ' 7 = 0 if VL same no matter if VR computed,',
     $      '  1/ulp otherwise', / )
 9994 FORMAT( ' N=', I5, ', IWK=', I2, ', seed=', 4( I4, ',' ),
     $      ' type ', I2, ', test(', I2, ')=', G10.3 )
 9993 FORMAT( ' CDRVEV: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of CDRVEV
*
      END
