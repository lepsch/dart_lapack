      SUBROUTINE CDRVST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, D1, D2, D3, WA1, WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LRWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               D1( * ), D2( * ), D3( * ), RESULT( * ), RWORK( * ), WA1( * ), WA2( * ), WA3( * )       COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ), V( LDU, * ), WORK( * ), Z( LDU, * )
      // ..
*
*  =====================================================================
*
*
      // .. Parameters ..
      REAL               ZERO, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0, TEN = 10.0E+0 )
      REAL               HALF
      PARAMETER          ( HALF = ONE / TWO )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 18 )
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IDIAG, IHBW, IINFO, IL, IMODE, INDWRK, INDX, IROW, ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL, JSIZE, JTYPE, KD, LGN, LIWEDC, LRWEDC, LWEDC, M, M2, M3, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLARND, SSXT1
      // EXTERNAL SLAMCH, SLARND, SSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, SLAFTS, XERBLA, CHBEV, CHBEVD, CHBEVX, CHEEV, CHEEVD, CHEEVR, CHEEVX, CHET21, CHET22, CHPEV, CHPEVD, CHPEVX, CLACPY, CLASET, CHEEVD_2STAGE, CHEEVR_2STAGE, CHEEVX_2STAGE, CHEEV_2STAGE, CHBEV_2STAGE, CHBEVD_2STAGE, CHBEVX_2STAGE, CLATMR, CLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 3*9 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4 /
      // ..
      // .. Executable Statements ..
*
      // 1)      Check for errors
*
      NTESTT = 0
      INFO = 0
*
      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE
*
      // Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.NMAX ) THEN
         INFO = -9
      ELSE IF( LDU.LT.NMAX ) THEN
         INFO = -16
      ELSE IF( 2*MAX( 2, NMAX )**2.GT.LWORK ) THEN
         INFO = -22
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CDRVST2STG', -INFO )
         RETURN
      END IF
*
      // Quick return if nothing to do
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN
*
      // More Important constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
      // Loop over sizes, types
*
      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
         ISEED3( I ) = ISEED( I )
   20 CONTINUE
*
      NERRS = 0
      NMATS = 0
*
      DO 1220 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( N.GT.0 ) THEN
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1
            LWEDC = MAX( 2*N+N*N, 2*N*N )
            LRWEDC = 1 + 4*N + 2*N*LGN + 3*N**2
            LIWEDC = 3 + 5*N
         ELSE
            LWEDC = 2
            LRWEDC = 8
            LIWEDC = 8
         END IF
         ANINV = ONE / REAL( MAX( 1, N ) )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 1210 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 1210
            NMATS = NMATS + 1
            NTEST = 0
*
            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE
*
            // 2)      Compute "A"
*
                    // Control parameters:
*
                // KMAGN  KMODE        KTYPE
            // =1  O(1)   clustered 1  zero
            // =2  large  clustered 2  identity
            // =3  small  exponential  (none)
            // =4         arithmetic   diagonal, (w/ eigenvalues)
            // =5         random log   Hermitian, w/ eigenvalues
            // =6         random       (none)
            // =7                      random diagonal
            // =8                      random Hermitian
            // =9                      band Hermitian, w/ eigenvalues
*
            IF( MTYPES.GT.MAXTYP ) GO TO 110
*
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
*
            // Compute norm
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
            CALL CLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
            IINFO = 0
            COND = ULPINV
*
            // Special Matrices -- Identity & Jordan block
*
                    // Zero
*
            IF( ITYPE.EQ.1 ) THEN
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
               // Identity
*
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
               // Diagonal Matrix, [Eigen]values Specified
*
               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
               // Hermitian, eigenvalues specified
*
               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
               // Diagonal, random eigenvalues
*
               CALL CLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
               // Hermitian, random eigenvalues
*
               CALL CLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
               // Hermitian banded, eigenvalues specified
*
               IHBW = INT( ( N-1 )*SLARND( 1, ISEED3 ) )
               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, IHBW, IHBW, 'Z', U, LDU, WORK, IINFO )
*
               // Store as dense matrix for most routines.
*
               CALL CLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
               DO 100 IDIAG = -IHBW, IHBW
                  IROW = IHBW - IDIAG + 1
                  J1 = MAX( 1, IDIAG+1 )
                  J2 = MIN( N, N+IDIAG )
                  DO 90 J = J1, J2
                     I = J - IDIAG
                     A( I, J ) = U( IROW, J )
   90             CONTINUE
  100          CONTINUE
            ELSE
               IINFO = 1
            END IF
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
  110       CONTINUE
*
            ABSTOL = UNFL + UNFL
            IF( N.LE.1 ) THEN
               IL = 1
               IU = N
            ELSE
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IF( IL.GT.IU ) THEN
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               END IF
            END IF
*
            // Perform tests storing upper or lower triangular
            // part of matrix.
*
            DO 1200 IUPLO = 0, 1
               IF( IUPLO.EQ.0 ) THEN
                  UPLO = 'L'
               ELSE
                  UPLO = 'U'
               END IF
*
               // Call CHEEVD and CHEEVX.
*
               CALL CLACPY( ' ', N, N, A, LDA, V, LDU )
*
               NTEST = NTEST + 1
               CALL CHEEVD( 'V', UPLO, N, A, LDU, D1, WORK, LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 130
                  END IF
               END IF
*
               // Do tests 1 and 2.
*
               CALL CHET21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 2
               CALL CHEEVD_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, LWORK, RWORK, LRWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVD_2STAGE(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 130
                  END IF
               END IF
*
               // Do test 3.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 120 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  120          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
  130          CONTINUE
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 1
*
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
                  IF( IL.NE.1 ) THEN
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
                  IF( IU.NE.N ) THEN
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
               ELSE
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               END IF
*
               CALL CHEEVX( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
               // Do tests 4 and 5.
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL CHET21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL CHEEVX_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVX_2STAGE(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
               // Do test 6.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 140 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  140          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
  150          CONTINUE
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 1
*
               CALL CHEEVX( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 160
                  END IF
               END IF
*
               // Do tests 7 and 8.
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               CALL CHEEVX_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVX_2STAGE(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 160
                  END IF
               END IF
*
               // Do test 9.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  160          CONTINUE
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 1
*
               CALL CHEEVX( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 170
                  END IF
               END IF
*
               // Do tests 10 and 11.
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               CALL CHEEVX_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVX_2STAGE(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 170
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 170
               END IF
*
               // Do test 12.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  170          CONTINUE
*
               // Call CHPEVD and CHPEVX.
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 190 J = 1, N
                     DO 180 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  180                CONTINUE
  190             CONTINUE
               ELSE
                  INDX = 1
                  DO 210 J = 1, N
                     DO 200 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  200                CONTINUE
  210             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               INDWRK = N*( N+1 ) / 2 + 1
               CALL CHPEVD( 'V', UPLO, N, WORK, D1, Z, LDU, WORK( INDWRK ), LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 270
                  END IF
               END IF
*
               // Do tests 13 and 14.
*
               CALL CHET21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 230 J = 1, N
                     DO 220 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  220                CONTINUE
  230             CONTINUE
               ELSE
                  INDX = 1
                  DO 250 J = 1, N
                     DO 240 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  240                CONTINUE
  250             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               INDWRK = N*( N+1 ) / 2 + 1
               CALL CHPEVD( 'N', UPLO, N, WORK, D3, Z, LDU, WORK( INDWRK ), LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 270
                  END IF
               END IF
*
               // Do test 15.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 260 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  260          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
               // Load array WORK with the upper or lower triangular part
               // of the matrix in packed form.
*
  270          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 290 J = 1, N
                     DO 280 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  280                CONTINUE
  290             CONTINUE
               ELSE
                  INDX = 1
                  DO 310 J = 1, N
                     DO 300 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  300                CONTINUE
  310             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
                  IF( IL.NE.1 ) THEN
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
                  IF( IU.NE.N ) THEN
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  ELSE IF( N.GT.0 ) THEN
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  END IF
               ELSE
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               END IF
*
               CALL CHPEVX( 'V', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 370
                  END IF
               END IF
*
               // Do tests 16 and 17.
*
               CALL CHET21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 330 J = 1, N
                     DO 320 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  320                CONTINUE
  330             CONTINUE
               ELSE
                  INDX = 1
                  DO 350 J = 1, N
                     DO 340 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  340                CONTINUE
  350             CONTINUE
               END IF
*
               CALL CHPEVX( 'N', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 370
                  END IF
               END IF
*
               // Do test 18.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 360 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  360          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
  370          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 390 J = 1, N
                     DO 380 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  380                CONTINUE
  390             CONTINUE
               ELSE
                  INDX = 1
                  DO 410 J = 1, N
                     DO 400 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  400                CONTINUE
  410             CONTINUE
               END IF
*
               CALL CHPEVX( 'V', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 460
                  END IF
               END IF
*
               // Do tests 19 and 20.
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 430 J = 1, N
                     DO 420 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  420                CONTINUE
  430             CONTINUE
               ELSE
                  INDX = 1
                  DO 450 J = 1, N
                     DO 440 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  440                CONTINUE
  450             CONTINUE
               END IF
*
               CALL CHPEVX( 'N', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 460
                  END IF
               END IF
*
               // Do test 21.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  460          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 480 J = 1, N
                     DO 470 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  470                CONTINUE
  480             CONTINUE
               ELSE
                  INDX = 1
                  DO 500 J = 1, N
                     DO 490 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  490                CONTINUE
  500             CONTINUE
               END IF
*
               CALL CHPEVX( 'V', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 550
                  END IF
               END IF
*
               // Do tests 22 and 23.
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 520 J = 1, N
                     DO 510 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  510                CONTINUE
  520             CONTINUE
               ELSE
                  INDX = 1
                  DO 540 J = 1, N
                     DO 530 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  530                CONTINUE
  540             CONTINUE
               END IF
*
               CALL CHPEVX( 'N', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 550
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 550
               END IF
*
               // Do test 24.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  550          CONTINUE
*
               // Call CHBEVD and CHBEVX.
*
               IF( JTYPE.LE.7 ) THEN
                  KD = 0
               ELSE IF( JTYPE.GE.8 .AND. JTYPE.LE.15 ) THEN
                  KD = MAX( N-1, 0 )
               ELSE
                  KD = IHBW
               END IF
*
               // Load array V with the upper or lower triangular part
               // of the matrix in band form.
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 570 J = 1, N
                     DO 560 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  560                CONTINUE
  570             CONTINUE
               ELSE
                  DO 590 J = 1, N
                     DO 580 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  580                CONTINUE
  590             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               CALL CHBEVD( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 )'CHBEVD(V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 650
                  END IF
               END IF
*
               // Do tests 25 and 26.
*
               CALL CHET21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 610 J = 1, N
                     DO 600 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  600                CONTINUE
  610             CONTINUE
               ELSE
                  DO 630 J = 1, N
                     DO 620 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  620                CONTINUE
  630             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               CALL CHBEVD_2STAGE( 'N', UPLO, N, KD, V, LDU, D3,  Z, LDU, WORK, LWORK, RWORK, LRWEDC, IWORK, LIWEDC, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 ) 'CHBEVD_2STAGE(N,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 650
                  END IF
               END IF
*
               // Do test 27.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 640 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  640          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
               // Load array V with the upper or lower triangular part
               // of the matrix in band form.
*
  650          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  DO 670 J = 1, N
                     DO 660 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  660                CONTINUE
  670             CONTINUE
               ELSE
                  DO 690 J = 1, N
                     DO 680 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  680                CONTINUE
  690             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               CALL CHBEVX( 'V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHBEVX(V,A,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 750
                  END IF
               END IF
*
               // Do tests 28 and 29.
*
               CALL CHET21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 710 J = 1, N
                     DO 700 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  700                CONTINUE
  710             CONTINUE
               ELSE
                  DO 730 J = 1, N
                     DO 720 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  720                CONTINUE
  730             CONTINUE
               END IF
*
               CALL CHBEVX_2STAGE( 'N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 ) 'CHBEVX_2STAGE(N,A,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 750
                  END IF
               END IF
*
               // Do test 30.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 740 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  740          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
               // Load array V with the upper or lower triangular part
               // of the matrix in band form.
*
  750          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  DO 770 J = 1, N
                     DO 760 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  760                CONTINUE
  770             CONTINUE
               ELSE
                  DO 790 J = 1, N
                     DO 780 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  780                CONTINUE
  790             CONTINUE
               END IF
*
               CALL CHBEVX( 'V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 )'CHBEVX(V,I,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 840
                  END IF
               END IF
*
               // Do tests 31 and 32.
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 810 J = 1, N
                     DO 800 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  800                CONTINUE
  810             CONTINUE
               ELSE
                  DO 830 J = 1, N
                     DO 820 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  820                CONTINUE
  830             CONTINUE
               END IF
               CALL CHBEVX_2STAGE( 'N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 ) 'CHBEVX_2STAGE(N,I,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 840
                  END IF
               END IF
*
               // Do test 33.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
               // Load array V with the upper or lower triangular part
               // of the matrix in band form.
*
  840          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  DO 860 J = 1, N
                     DO 850 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  850                CONTINUE
  860             CONTINUE
               ELSE
                  DO 880 J = 1, N
                     DO 870 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  870                CONTINUE
  880             CONTINUE
               END IF
               CALL CHBEVX( 'V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 )'CHBEVX(V,V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 930
                  END IF
               END IF
*
               // Do tests 34 and 35.
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 900 J = 1, N
                     DO 890 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  890                CONTINUE
  900             CONTINUE
               ELSE
                  DO 920 J = 1, N
                     DO 910 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  910                CONTINUE
  920             CONTINUE
               END IF
               CALL CHBEVX_2STAGE( 'N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 ) 'CHBEVX_2STAGE(N,V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 930
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 930
               END IF
*
               // Do test 36.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  930          CONTINUE
*
               // Call CHEEV
*
               CALL CLACPY( ' ', N, N, A, LDA, V, LDU )
*
               NTEST = NTEST + 1
               CALL CHEEV( 'V', UPLO, N, A, LDU, D1, WORK, LWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 950
                  END IF
               END IF
*
               // Do tests 37 and 38
*
               CALL CHET21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 2
               CALL CHEEV_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, LWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEV_2STAGE(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 950
                  END IF
               END IF
*
               // Do test 39
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 940 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  940          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
  950          CONTINUE
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               // Call CHPEV
*
               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 970 J = 1, N
                     DO 960 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  960                CONTINUE
  970             CONTINUE
               ELSE
                  INDX = 1
                  DO 990 J = 1, N
                     DO 980 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  980                CONTINUE
  990             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               INDWRK = N*( N+1 ) / 2 + 1
               CALL CHPEV( 'V', UPLO, N, WORK, D1, Z, LDU, WORK( INDWRK ), RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1050
                  END IF
               END IF
*
               // Do tests 40 and 41.
*
               CALL CHET21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 1010 J = 1, N
                     DO 1000 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1000                CONTINUE
 1010             CONTINUE
               ELSE
                  INDX = 1
                  DO 1030 J = 1, N
                     DO 1020 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1020                CONTINUE
 1030             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               INDWRK = N*( N+1 ) / 2 + 1
               CALL CHPEV( 'N', UPLO, N, WORK, D3, Z, LDU, WORK( INDWRK ), RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHPEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1050
                  END IF
               END IF
*
               // Do test 42
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1040 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1040          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
 1050          CONTINUE
*
               // Call CHBEV
*
               IF( JTYPE.LE.7 ) THEN
                  KD = 0
               ELSE IF( JTYPE.GE.8 .AND. JTYPE.LE.15 ) THEN
                  KD = MAX( N-1, 0 )
               ELSE
                  KD = IHBW
               END IF
*
               // Load array V with the upper or lower triangular part
               // of the matrix in band form.
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1070 J = 1, N
                     DO 1060 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1060                CONTINUE
 1070             CONTINUE
               ELSE
                  DO 1090 J = 1, N
                     DO 1080 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1080                CONTINUE
 1090             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               CALL CHBEV( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 )'CHBEV(V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1140
                  END IF
               END IF
*
               // Do tests 43 and 44.
*
               CALL CHET21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1110 J = 1, N
                     DO 1100 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1100                CONTINUE
 1110             CONTINUE
               ELSE
                  DO 1130 J = 1, N
                     DO 1120 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1120                CONTINUE
 1130             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               CALL CHBEV_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9998 ) 'CHBEV_2STAGE(N,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1140
                  END IF
               END IF
*
 1140          CONTINUE
*
               // Do test 45.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1150 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1150          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
               CALL CLACPY( ' ', N, N, A, LDA, V, LDU )
               NTEST = NTEST + 1
               CALL CHEEVR( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVR(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1170
                  END IF
               END IF
*
               // Do tests 45 and 46 (or ... )
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL CHET21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL CHEEVR_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVR_2STAGE(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1170
                  END IF
               END IF
*
               // Do test 47 (or ... )
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1160 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
 1160          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
*
 1170          CONTINUE
*
               NTEST = NTEST + 1
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL CHEEVR( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVR(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1180
                  END IF
               END IF
*
               // Do tests 48 and 49 (or +??)
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL CHEEVR_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVR_2STAGE(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1180
                  END IF
               END IF
*
               // Do test 50 (or +??)
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
 1180          CONTINUE
*
               NTEST = NTEST + 1
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL CHEEVR( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'CHEEVR(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1190
                  END IF
               END IF
*
               // Do tests 51 and 52 (or +??)
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL CHET22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL CHEEVR_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 ) 'CHEEVR_2STAGE(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1190
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 1190
               END IF
*
               // Do test 52 (or +??)
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
               CALL CLACPY( ' ', N, N, V, LDU, A, LDA )
*
*
*
*
               // Load array V with the upper or lower triangular part
               // of the matrix in band form.
*
 1190          CONTINUE
*
 1200       CONTINUE
*
            // End of Loop -- Check for RESULT(j) > THRESH
*
            NTESTT = NTESTT + NTEST
            CALL SLAFTS( 'CST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS )
*
 1210    CONTINUE
 1220 CONTINUE
*
      // Summary
*
      CALL ALASVM( 'CST', NOUNIT, NERRS, NTESTT, 0 )
*
 9999 FORMAT( ' CDRVST2STG: ', A, ' returned INFO=', I6, / 9X, 'N=', I6,
     $      ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' CDRVST2STG: ', A, ' returned INFO=', I6, / 9X, 'N=', I6,
     $      ', KD=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5,
     $      ')' )
*
      RETURN
*
      // End of CDRVST2STG
*
      END
