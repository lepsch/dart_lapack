      SUBROUTINE DCHKBD( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS,
     $                   ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX,
     $                   Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK,
     $                   IWORK, NOUT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS,
     $                   NSIZES, NTYPES
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * )
      DOUBLE PRECISION   A( LDA, * ), BD( * ), BE( * ), PT( LDPT, * ),
     $                   Q( LDQ, * ), S1( * ), S2( * ), U( LDPT, * ),
     $                   VT( LDPT, * ), WORK( * ), X( LDX, * ),
     $                   Y( LDX, * ), Z( LDX, * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 0.5D0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 16 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADMM, BADNN, BIDIAG
      CHARACTER          UPLO
      CHARACTER*3        PATH
      INTEGER            I, IINFO, IL, IMODE, ITEMP, ITYPE, IU, IWBD,
     $                   IWBE, IWBS, IWBZ, IWWORK, J, JCOL, JSIZE,
     $                   JTYPE, LOG2UI, M, MINWRK, MMAX, MNMAX, MNMIN,
     $                   MNMIN2, MQ, MTYPES, N, NFAIL, NMAX,
     $                   NS1, NS2, NTEST
      DOUBLE PRECISION   ABSTOL, AMNINV, ANORM, COND, OVFL, RTOVFL,
     $                   RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL,
     $                   VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 1 ), IOLDSD( 4 ), ISEED2( 4 ),
     $                   KMAGN( MAXTYP ), KMODE( MAXTYP ),
     $                   KTYPE( MAXTYP )
      DOUBLE PRECISION   DUM( 1 ), DUMMA( 1 ), RESULT( 40 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLARND, DSXT1
      EXTERNAL           DLAMCH, DLARND, DSXT1
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASUM, DBDSDC, DBDSQR, DBDSVDX, DBDT01,
     $                   DBDT02, DBDT03, DBDT04, DCOPY, DGEBRD,
     $                   DGEMM, DLACPY, DLAHD2, DLASET, DLATMR,
     $                   DLATMS, DORGBR, DORT01, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, EXP, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA            KTYPE / 1, 2, 5*4, 5*6, 3*9, 10 /
      DATA            KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3, 0 /
      DATA            KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0,
     $                0, 0, 0 /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      INFO = 0
*
      BADMM = .FALSE.
      BADNN = .FALSE.
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      MINWRK = 1
      DO 10 J = 1, NSIZES
         MMAX = MAX( MMAX, MVAL( J ) )
         IF( MVAL( J ).LT.0 )
     $      BADMM = .TRUE.
         NMAX = MAX( NMAX, NVAL( J ) )
         IF( NVAL( J ).LT.0 )
     $      BADNN = .TRUE.
         MNMAX = MAX( MNMAX, MIN( MVAL( J ), NVAL( J ) ) )
         MINWRK = MAX( MINWRK, 3*( MVAL( J )+NVAL( J ) ),
     $            MVAL( J )*( MVAL( J )+MAX( MVAL( J ), NVAL( J ),
     $            NRHS )+1 )+NVAL( J )*MIN( NVAL( J ), MVAL( J ) ) )
   10 CONTINUE
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADMM ) THEN
         INFO = -2
      ELSE IF( BADNN ) THEN
         INFO = -3
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MMAX ) THEN
         INFO = -11
      ELSE IF( LDX.LT.MMAX ) THEN
         INFO = -17
      ELSE IF( LDQ.LT.MMAX ) THEN
         INFO = -21
      ELSE IF( LDPT.LT.MNMAX ) THEN
         INFO = -23
      ELSE IF( MINWRK.GT.LWORK ) THEN
         INFO = -27
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DCHKBD', -INFO )
         RETURN
      END IF
*
*     Initialize constants
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'BD'
      NFAIL = 0
      NTEST = 0
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) )
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
      INFOT = 0
      ABSTOL = 2*UNFL
*
*     Loop over sizes, types
*
      DO 300 JSIZE = 1, NSIZES
         M = MVAL( JSIZE )
         N = NVAL( JSIZE )
         MNMIN = MIN( M, N )
         AMNINV = ONE / MAX( M, N, 1 )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 290 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 290
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
            DO 30 J = 1, 34
               RESULT( J ) = -ONE
   30       CONTINUE
*
            UPLO = ' '
*
*           Compute "A"
*
*           Control parameters:
*
*           KMAGN  KMODE        KTYPE
*       =1  O(1)   clustered 1  zero
*       =2  large  clustered 2  identity
*       =3  small  exponential  (none)
*       =4         arithmetic   diagonal, (w/ eigenvalues)
*       =5         random       symmetric, w/ eigenvalues
*       =6                      nonsymmetric, w/ singular values
*       =7                      random diagonal
*       =8                      random symmetric
*       =9                      random nonsymmetric
*       =10                     random bidiagonal (log. distrib.)
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 100
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
            ANORM = ( RTOVFL*ULP )*AMNINV
            GO TO 70
*
   60       CONTINUE
            ANORM = RTUNFL*MAX( M, N )*ULPINV
            GO TO 70
*
   70       CONTINUE
*
            CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            COND = ULPINV
*
            BIDIAG = .FALSE.
            IF( ITYPE.EQ.1 ) THEN
*
*              Zero matrix
*
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 80 JCOL = 1, MNMIN
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              Diagonal Matrix, [Eigen]values Specified
*
               CALL DLATMS( MNMIN, MNMIN, 'S', ISEED, 'N', WORK, IMODE,
     $                      COND, ANORM, 0, 0, 'N', A, LDA,
     $                      WORK( MNMIN+1 ), IINFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Symmetric, eigenvalues specified
*
               CALL DLATMS( MNMIN, MNMIN, 'S', ISEED, 'S', WORK, IMODE,
     $                      COND, ANORM, M, N, 'N', A, LDA,
     $                      WORK( MNMIN+1 ), IINFO )
*
            ELSE IF( ITYPE.EQ.6 ) THEN
*
*              Nonsymmetric, singular values specified
*
               CALL DLATMS( M, N, 'S', ISEED, 'N', WORK, IMODE, COND,
     $                      ANORM, M, N, 'N', A, LDA, WORK( MNMIN+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              Diagonal, random entries
*
               CALL DLATMR( MNMIN, MNMIN, 'S', ISEED, 'N', WORK, 6, ONE,
     $                      ONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE,
     $                      WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, 0, 0,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              Symmetric, random entries
*
               CALL DLATMR( MNMIN, MNMIN, 'S', ISEED, 'S', WORK, 6, ONE,
     $                      ONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE,
     $                      WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              Nonsymmetric, random entries
*
               CALL DLATMR( M, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( MNMIN+1 ), 1, ONE,
     $                      WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.10 ) THEN
*
*              Bidiagonal, random entries
*
               TEMP1 = -TWO*LOG( ULP )
               DO 90 J = 1, MNMIN
                  BD( J ) = EXP( TEMP1*DLARND( 2, ISEED ) )
                  IF( J.LT.MNMIN )
     $               BE( J ) = EXP( TEMP1*DLARND( 2, ISEED ) )
   90          CONTINUE
*
               IINFO = 0
               BIDIAG = .TRUE.
               IF( M.GE.N ) THEN
                  UPLO = 'U'
               ELSE
                  UPLO = 'L'
               END IF
            ELSE
               IINFO = 1
            END IF
*
            IF( IINFO.EQ.0 ) THEN
*
*              Generate Right-Hand Side
*
               IF( BIDIAG ) THEN
                  CALL DLATMR( MNMIN, NRHS, 'S', ISEED, 'N', WORK, 6,
     $                         ONE, ONE, 'T', 'N', WORK( MNMIN+1 ), 1,
     $                         ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N',
     $                         IWORK, MNMIN, NRHS, ZERO, ONE, 'NO', Y,
     $                         LDX, IWORK, IINFO )
               ELSE
                  CALL DLATMR( M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE,
     $                         ONE, 'T', 'N', WORK( M+1 ), 1, ONE,
     $                         WORK( 2*M+1 ), 1, ONE, 'N', IWORK, M,
     $                         NRHS, ZERO, ONE, 'NO', X, LDX, IWORK,
     $                         IINFO )
               END IF
            END IF
*
*           Error Exit
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'Generator', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
  100       CONTINUE
*
*           Call DGEBRD and DORGBR to compute B, Q, and P, do tests.
*
            IF( .NOT.BIDIAG ) THEN
*
*              Compute transformations to reduce A to bidiagonal form:
*              B := Q' * A * P.
*
               CALL DLACPY( ' ', M, N, A, LDA, Q, LDQ )
               CALL DGEBRD( M, N, Q, LDQ, BD, BE, WORK, WORK( MNMIN+1 ),
     $                      WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )
*
*              Check error code from DGEBRD.
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9998 )'DGEBRD', IINFO, M, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
               CALL DLACPY( ' ', M, N, Q, LDQ, PT, LDPT )
               IF( M.GE.N ) THEN
                  UPLO = 'U'
               ELSE
                  UPLO = 'L'
               END IF
*
*              Generate Q
*
               MQ = M
               IF( NRHS.LE.0 )
     $            MQ = MNMIN
               CALL DORGBR( 'Q', M, MQ, N, Q, LDQ, WORK,
     $                      WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )
*
*              Check error code from DORGBR.
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9998 )'DORGBR(Q)', IINFO, M, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
*              Generate P'
*
               CALL DORGBR( 'P', MNMIN, N, M, PT, LDPT, WORK( MNMIN+1 ),
     $                      WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )
*
*              Check error code from DORGBR.
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9998 )'DORGBR(P)', IINFO, M, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
*              Apply Q' to an M by NRHS matrix X:  Y := Q' * X.
*
               CALL DGEMM( 'Transpose', 'No transpose', M, NRHS, M, ONE,
     $                     Q, LDQ, X, LDX, ZERO, Y, LDX )
*
*              Test 1:  Check the decomposition A := Q * B * PT
*                   2:  Check the orthogonality of Q
*                   3:  Check the orthogonality of PT
*
               CALL DBDT01( M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT,
     $                      WORK, RESULT( 1 ) )
               CALL DORT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK,
     $                      RESULT( 2 ) )
               CALL DORT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK,
     $                      RESULT( 3 ) )
            END IF
*
*           Use DBDSQR to form the SVD of the bidiagonal matrix B:
*           B := U * S1 * VT, and compute Z = U' * Y.
*
            CALL DCOPY( MNMIN, BD, 1, S1, 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK, 1 )
            CALL DLACPY( ' ', M, NRHS, Y, LDX, Z, LDX )
            CALL DLASET( 'Full', MNMIN, MNMIN, ZERO, ONE, U, LDPT )
            CALL DLASET( 'Full', MNMIN, MNMIN, ZERO, ONE, VT, LDPT )
*
            CALL DBDSQR( UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, WORK, VT,
     $                   LDPT, U, LDPT, Z, LDX, WORK( MNMIN+1 ), IINFO )
*
*           Check error code from DBDSQR.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSQR(vects)', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 4 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Use DBDSQR to compute only the singular values of the
*           bidiagonal matrix B;  U, VT, and Z should not be modified.
*
            CALL DCOPY( MNMIN, BD, 1, S2, 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK, 1 )
*
            CALL DBDSQR( UPLO, MNMIN, 0, 0, 0, S2, WORK, VT, LDPT, U,
     $                   LDPT, Z, LDX, WORK( MNMIN+1 ), IINFO )
*
*           Check error code from DBDSQR.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSQR(values)', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 9 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Test 4:  Check the decomposition B := U * S1 * VT
*                5:  Check the computation Z := U' * Y
*                6:  Check the orthogonality of U
*                7:  Check the orthogonality of VT
*
            CALL DBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT,
     $                   WORK, RESULT( 4 ) )
            CALL DBDT02( MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK,
     $                   RESULT( 5 ) )
            CALL DORT01( 'Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK,
     $                   RESULT( 6 ) )
            CALL DORT01( 'Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK,
     $                   RESULT( 7 ) )
*
*           Test 8:  Check that the singular values are sorted in
*                    non-increasing order and are non-negative
*
            RESULT( 8 ) = ZERO
            DO 110 I = 1, MNMIN - 1
               IF( S1( I ).LT.S1( I+1 ) )
     $            RESULT( 8 ) = ULPINV
               IF( S1( I ).LT.ZERO )
     $            RESULT( 8 ) = ULPINV
  110       CONTINUE
            IF( MNMIN.GE.1 ) THEN
               IF( S1( MNMIN ).LT.ZERO )
     $            RESULT( 8 ) = ULPINV
            END IF
*
*           Test 9:  Compare DBDSQR with and without singular vectors
*
            TEMP2 = ZERO
*
            DO 120 J = 1, MNMIN
               TEMP1 = ABS( S1( J )-S2( J ) ) /
     $                 MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ),
     $                 ULP*MAX( ABS( S1( J ) ), ABS( S2( J ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  120       CONTINUE
*
            RESULT( 9 ) = TEMP2
*
*           Test 10:  Sturm sequence test of singular values
*                     Go up by factors of two until it succeeds
*
            TEMP1 = THRESH*( HALF-ULP )
*
            DO 130 J = 0, LOG2UI
*               CALL DSVDCH( MNMIN, BD, BE, S1, TEMP1, IINFO )
               IF( IINFO.EQ.0 )
     $            GO TO 140
               TEMP1 = TEMP1*TWO
  130       CONTINUE
*
  140       CONTINUE
            RESULT( 10 ) = TEMP1
*
*           Use DBDSQR to form the decomposition A := (QU) S (VT PT)
*           from the bidiagonal form A := Q B PT.
*
            IF( .NOT.BIDIAG ) THEN
               CALL DCOPY( MNMIN, BD, 1, S2, 1 )
               IF( MNMIN.GT.0 )
     $            CALL DCOPY( MNMIN-1, BE, 1, WORK, 1 )
*
               CALL DBDSQR( UPLO, MNMIN, N, M, NRHS, S2, WORK, PT, LDPT,
     $                      Q, LDQ, Y, LDX, WORK( MNMIN+1 ), IINFO )
*
*              Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
*                   12:  Check the computation Z := U' * Q' * X
*                   13:  Check the orthogonality of Q*U
*                   14:  Check the orthogonality of VT*PT
*
               CALL DBDT01( M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT,
     $                      LDPT, WORK, RESULT( 11 ) )
               CALL DBDT02( M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK,
     $                      RESULT( 12 ) )
               CALL DORT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK,
     $                      RESULT( 13 ) )
               CALL DORT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK,
     $                      RESULT( 14 ) )
            END IF
*
*           Use DBDSDC to form the SVD of the bidiagonal matrix B:
*           B := U * S1 * VT
*
            CALL DCOPY( MNMIN, BD, 1, S1, 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK, 1 )
            CALL DLASET( 'Full', MNMIN, MNMIN, ZERO, ONE, U, LDPT )
            CALL DLASET( 'Full', MNMIN, MNMIN, ZERO, ONE, VT, LDPT )
*
            CALL DBDSDC( UPLO, 'I', MNMIN, S1, WORK, U, LDPT, VT, LDPT,
     $                   DUM, IDUM, WORK( MNMIN+1 ), IWORK, IINFO )
*
*           Check error code from DBDSDC.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSDC(vects)', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 15 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Use DBDSDC to compute only the singular values of the
*           bidiagonal matrix B;  U and VT should not be modified.
*
            CALL DCOPY( MNMIN, BD, 1, S2, 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK, 1 )
*
            CALL DBDSDC( UPLO, 'N', MNMIN, S2, WORK, DUM, 1, DUM, 1,
     $                   DUM, IDUM, WORK( MNMIN+1 ), IWORK, IINFO )
*
*           Check error code from DBDSDC.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSDC(values)', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 18 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Test 15:  Check the decomposition B := U * S1 * VT
*                16:  Check the orthogonality of U
*                17:  Check the orthogonality of VT
*
            CALL DBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT,
     $                   WORK, RESULT( 15 ) )
            CALL DORT01( 'Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK,
     $                   RESULT( 16 ) )
            CALL DORT01( 'Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK,
     $                   RESULT( 17 ) )
*
*           Test 18:  Check that the singular values are sorted in
*                     non-increasing order and are non-negative
*
            RESULT( 18 ) = ZERO
            DO 150 I = 1, MNMIN - 1
               IF( S1( I ).LT.S1( I+1 ) )
     $            RESULT( 18 ) = ULPINV
               IF( S1( I ).LT.ZERO )
     $            RESULT( 18 ) = ULPINV
  150       CONTINUE
            IF( MNMIN.GE.1 ) THEN
               IF( S1( MNMIN ).LT.ZERO )
     $            RESULT( 18 ) = ULPINV
            END IF
*
*           Test 19:  Compare DBDSQR with and without singular vectors
*
            TEMP2 = ZERO
*
            DO 160 J = 1, MNMIN
               TEMP1 = ABS( S1( J )-S2( J ) ) /
     $                 MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ),
     $                 ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  160       CONTINUE
*
            RESULT( 19 ) = TEMP2
*
*
*           Use DBDSVDX to compute the SVD of the bidiagonal matrix B:
*           B := U * S1 * VT
*
            IF( JTYPE.EQ.10 .OR. JTYPE.EQ.16 ) THEN
*              =================================
*              Matrix types temporarily disabled
*              =================================
               RESULT( 20:34 ) = ZERO
               GO TO 270
            END IF
*
            IWBS = 1
            IWBD = IWBS + MNMIN
            IWBE = IWBD + MNMIN
            IWBZ = IWBE + MNMIN
            IWWORK = IWBZ + 2*MNMIN*(MNMIN+1)
            MNMIN2 = MAX( 1,MNMIN*2 )
*
            CALL DCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
*
            CALL DBDSVDX( UPLO, 'V', 'A', MNMIN, WORK( IWBD ),
     $                    WORK( IWBE ), ZERO, ZERO, 0, 0, NS1, S1,
     $                    WORK( IWBZ ), MNMIN2, WORK( IWWORK ),
     $                    IWORK, IINFO)
*
*           Check error code from DBDSVDX.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSVDX(vects,A)', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 20 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
            J = IWBZ
            DO 170 I = 1, NS1
               CALL DCOPY( MNMIN, WORK( J ), 1, U( 1,I ), 1 )
               J = J + MNMIN
               CALL DCOPY( MNMIN, WORK( J ), 1, VT( I,1 ), LDPT )
               J = J + MNMIN
  170       CONTINUE
*
*           Use DBDSVDX to compute only the singular values of the
*           bidiagonal matrix B;  U and VT should not be modified.
*
            IF( JTYPE.EQ.9 ) THEN
*              =================================
*              Matrix types temporarily disabled
*              =================================
               RESULT( 24 ) = ZERO
               GO TO 270
            END IF
*
            CALL DCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
*
            CALL DBDSVDX( UPLO, 'N', 'A', MNMIN, WORK( IWBD ),
     $                    WORK( IWBE ), ZERO, ZERO, 0, 0, NS2, S2,
     $                    WORK( IWBZ ), MNMIN2, WORK( IWWORK ),
     $                    IWORK, IINFO )
*
*           Check error code from DBDSVDX.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSVDX(values,A)', IINFO,
     $            M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 24 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Save S1 for tests 30-34.
*
            CALL DCOPY( MNMIN, S1, 1, WORK( IWBS ), 1 )
*
*           Test 20:  Check the decomposition B := U * S1 * VT
*                21:  Check the orthogonality of U
*                22:  Check the orthogonality of VT
*                23:  Check that the singular values are sorted in
*                     non-increasing order and are non-negative
*                24:  Compare DBDSVDX with and without singular vectors
*
            CALL DBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT,
     $                   LDPT, WORK( IWBS+MNMIN ), RESULT( 20 ) )
            CALL DORT01( 'Columns', MNMIN, MNMIN, U, LDPT,
     $                   WORK( IWBS+MNMIN ), LWORK-MNMIN,
     $                   RESULT( 21 ) )
            CALL DORT01( 'Rows', MNMIN, MNMIN, VT, LDPT,
     $                   WORK( IWBS+MNMIN ), LWORK-MNMIN,
     $                   RESULT( 22) )
*
            RESULT( 23 ) = ZERO
            DO 180 I = 1, MNMIN - 1
               IF( S1( I ).LT.S1( I+1 ) )
     $            RESULT( 23 ) = ULPINV
               IF( S1( I ).LT.ZERO )
     $            RESULT( 23 ) = ULPINV
  180       CONTINUE
            IF( MNMIN.GE.1 ) THEN
               IF( S1( MNMIN ).LT.ZERO )
     $            RESULT( 23 ) = ULPINV
            END IF
*
            TEMP2 = ZERO
            DO 190 J = 1, MNMIN
               TEMP1 = ABS( S1( J )-S2( J ) ) /
     $                 MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ),
     $                 ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  190       CONTINUE
            RESULT( 24 ) = TEMP2
            ANORM = S1( 1 )
*
*           Use DBDSVDX with RANGE='I': choose random values for IL and
*           IU, and ask for the IL-th through IU-th singular values
*           and corresponding vectors.
*
            DO 200 I = 1, 4
               ISEED2( I ) = ISEED( I )
  200       CONTINUE
            IF( MNMIN.LE.1 ) THEN
               IL = 1
               IU = MNMIN
            ELSE
               IL = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
               IF( IU.LT.IL ) THEN
                  ITEMP = IU
                  IU = IL
                  IL = ITEMP
               END IF
            END IF
*
            CALL DCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
*
            CALL DBDSVDX( UPLO, 'V', 'I', MNMIN, WORK( IWBD ),
     $                    WORK( IWBE ), ZERO, ZERO, IL, IU, NS1, S1,
     $                    WORK( IWBZ ), MNMIN2, WORK( IWWORK ),
     $                    IWORK, IINFO)
*
*           Check error code from DBDSVDX.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSVDX(vects,I)', IINFO,
     $            M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 25 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
            J = IWBZ
            DO 210 I = 1, NS1
               CALL DCOPY( MNMIN, WORK( J ), 1, U( 1,I ), 1 )
               J = J + MNMIN
               CALL DCOPY( MNMIN, WORK( J ), 1, VT( I,1 ), LDPT )
               J = J + MNMIN
  210       CONTINUE
*
*           Use DBDSVDX to compute only the singular values of the
*           bidiagonal matrix B;  U and VT should not be modified.
*
            CALL DCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
*
            CALL DBDSVDX( UPLO, 'N', 'I', MNMIN, WORK( IWBD ),
     $                    WORK( IWBE ), ZERO, ZERO, IL, IU, NS2, S2,
     $                    WORK( IWBZ ), MNMIN2, WORK( IWWORK ),
     $                    IWORK, IINFO )
*
*           Check error code from DBDSVDX.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSVDX(values,I)', IINFO,
     $            M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 29 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Test 25:  Check S1 - U' * B * VT'
*                26:  Check the orthogonality of U
*                27:  Check the orthogonality of VT
*                28:  Check that the singular values are sorted in
*                     non-increasing order and are non-negative
*                29:  Compare DBDSVDX with and without singular vectors
*
            CALL DBDT04( UPLO, MNMIN, BD, BE, S1, NS1, U,
     $                   LDPT, VT, LDPT, WORK( IWBS+MNMIN ),
     $                   RESULT( 25 ) )
            CALL DORT01( 'Columns', MNMIN, NS1, U, LDPT,
     $                   WORK( IWBS+MNMIN ), LWORK-MNMIN,
     $                   RESULT( 26 ) )
            CALL DORT01( 'Rows', NS1, MNMIN, VT, LDPT,
     $                   WORK( IWBS+MNMIN ), LWORK-MNMIN,
     $                   RESULT( 27 ) )
*
            RESULT( 28 ) = ZERO
            DO 220 I = 1, NS1 - 1
               IF( S1( I ).LT.S1( I+1 ) )
     $            RESULT( 28 ) = ULPINV
               IF( S1( I ).LT.ZERO )
     $            RESULT( 28 ) = ULPINV
  220       CONTINUE
            IF( NS1.GE.1 ) THEN
               IF( S1( NS1 ).LT.ZERO )
     $            RESULT( 28 ) = ULPINV
            END IF
*
            TEMP2 = ZERO
            DO 230 J = 1, NS1
               TEMP1 = ABS( S1( J )-S2( J ) ) /
     $                 MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ),
     $                 ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  230       CONTINUE
            RESULT( 29 ) = TEMP2
*
*           Use DBDSVDX with RANGE='V': determine the values VL and VU
*           of the IL-th and IU-th singular values and ask for all
*           singular values in this range.
*
            CALL DCOPY( MNMIN, WORK( IWBS ), 1, S1, 1 )
*
            IF( MNMIN.GT.0 ) THEN
               IF( IL.NE.1 ) THEN
                  VU = S1( IL ) + MAX( HALF*ABS( S1( IL )-S1( IL-1 ) ),
     $                 ULP*ANORM, TWO*RTUNFL )
               ELSE
                  VU = S1( 1 ) + MAX( HALF*ABS( S1( MNMIN )-S1( 1 ) ),
     $                 ULP*ANORM, TWO*RTUNFL )
               END IF
               IF( IU.NE.NS1 ) THEN
                  VL = S1( IU ) - MAX( ULP*ANORM, TWO*RTUNFL,
     $                 HALF*ABS( S1( IU+1 )-S1( IU ) ) )
               ELSE
                  VL = S1( NS1 ) - MAX( ULP*ANORM, TWO*RTUNFL,
     $                 HALF*ABS( S1( MNMIN )-S1( 1 ) ) )
               END IF
               VL = MAX( VL,ZERO )
               VU = MAX( VU,ZERO )
               IF( VL.GE.VU ) VU = MAX( VU*2, VU+VL+HALF )
            ELSE
               VL = ZERO
               VU = ONE
            END IF
*
            CALL DCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
*
            CALL DBDSVDX( UPLO, 'V', 'V', MNMIN, WORK( IWBD ),
     $                    WORK( IWBE ), VL, VU, 0, 0, NS1, S1,
     $                    WORK( IWBZ ), MNMIN2, WORK( IWWORK ),
     $                    IWORK, IINFO )
*
*           Check error code from DBDSVDX.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSVDX(vects,V)', IINFO,
     $            M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 30 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
            J = IWBZ
            DO 240 I = 1, NS1
               CALL DCOPY( MNMIN, WORK( J ), 1, U( 1,I ), 1 )
               J = J + MNMIN
               CALL DCOPY( MNMIN, WORK( J ), 1, VT( I,1 ), LDPT )
               J = J + MNMIN
  240       CONTINUE
*
*           Use DBDSVDX to compute only the singular values of the
*           bidiagonal matrix B;  U and VT should not be modified.
*
            CALL DCOPY( MNMIN, BD, 1, WORK( IWBD ), 1 )
            IF( MNMIN.GT.0 )
     $         CALL DCOPY( MNMIN-1, BE, 1, WORK( IWBE ), 1 )
*
            CALL DBDSVDX( UPLO, 'N', 'V', MNMIN, WORK( IWBD ),
     $                    WORK( IWBE ), VL, VU, 0, 0, NS2, S2,
     $                    WORK( IWBZ ), MNMIN2, WORK( IWWORK ),
     $                    IWORK, IINFO )
*
*           Check error code from DBDSVDX.
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'DBDSVDX(values,V)', IINFO,
     $            M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               ELSE
                  RESULT( 34 ) = ULPINV
                  GO TO 270
               END IF
            END IF
*
*           Test 30:  Check S1 - U' * B * VT'
*                31:  Check the orthogonality of U
*                32:  Check the orthogonality of VT
*                33:  Check that the singular values are sorted in
*                     non-increasing order and are non-negative
*                34:  Compare DBDSVDX with and without singular vectors
*
            CALL DBDT04( UPLO, MNMIN, BD, BE, S1, NS1, U,
     $                   LDPT, VT, LDPT, WORK( IWBS+MNMIN ),
     $                   RESULT( 30 ) )
            CALL DORT01( 'Columns', MNMIN, NS1, U, LDPT,
     $                   WORK( IWBS+MNMIN ), LWORK-MNMIN,
     $                   RESULT( 31 ) )
            CALL DORT01( 'Rows', NS1, MNMIN, VT, LDPT,
     $                   WORK( IWBS+MNMIN ), LWORK-MNMIN,
     $                   RESULT( 32 ) )
*
            RESULT( 33 ) = ZERO
            DO 250 I = 1, NS1 - 1
               IF( S1( I ).LT.S1( I+1 ) )
     $            RESULT( 28 ) = ULPINV
               IF( S1( I ).LT.ZERO )
     $            RESULT( 28 ) = ULPINV
  250       CONTINUE
            IF( NS1.GE.1 ) THEN
               IF( S1( NS1 ).LT.ZERO )
     $            RESULT( 28 ) = ULPINV
            END IF
*
            TEMP2 = ZERO
            DO 260 J = 1, NS1
               TEMP1 = ABS( S1( J )-S2( J ) ) /
     $                 MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ),
     $                 ULP*MAX( ABS( S1( 1 ) ), ABS( S2( 1 ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  260       CONTINUE
            RESULT( 34 ) = TEMP2
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
  270       CONTINUE
*
            DO 280 J = 1, 34
               IF( RESULT( J ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 )
     $               CALL DLAHD2( NOUT, PATH )
                  WRITE( NOUT, FMT = 9999 )M, N, JTYPE, IOLDSD, J,
     $               RESULT( J )
                  NFAIL = NFAIL + 1
               END IF
  280       CONTINUE
            IF( .NOT.BIDIAG ) THEN
               NTEST = NTEST + 34
            ELSE
               NTEST = NTEST + 30
            END IF
*
  290    CONTINUE
  300 CONTINUE
*
*     Summary
*
      CALL ALASUM( PATH, NOUT, NFAIL, NTEST, 0 )
*
      RETURN
*
*     End of DCHKBD
*
 9999 FORMAT( ' M=', I5, ', N=', I5, ', type ', I2, ', seed=',
     $      4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9998 FORMAT( ' DCHKBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )
*
      END