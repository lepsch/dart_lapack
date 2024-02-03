      SUBROUTINE ZCHKBD( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS, ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX, Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK, RWORK, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), MVAL( * ), NVAL( * );
      double             BD( * ), BE( * ), RWORK( * ), S1( * ), S2( * );
      COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ), U( LDPT, * ), VT( LDPT, * ), WORK( * ), X( LDX, * ), Y( LDX, * ), Z( LDX, * )
      // ..

* ======================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 16 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN, BIDIAG;
      String             UPLO;
      String             PATH;
      int                I, IINFO, IMODE, ITYPE, J, JCOL, JSIZE, JTYPE, LOG2UI, M, MINWRK, MMAX, MNMAX, MNMIN, MQ, MTYPES, N, NFAIL, NMAX, NTEST;
      double             AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      int                IOLDSD( 4 ), IWORK( 1 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double             DUMMA( 1 ), RESULT( 14 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLARND;
      // EXTERNAL DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASUM, DCOPY, DLAHD2, DSVDCH, XERBLA, ZBDSQR, ZBDT01, ZBDT02, ZBDT03, ZGEBRD, ZGEMM, ZLACPY, ZLASET, ZLATMR, ZLATMS, ZUNGBR, ZUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, EXP, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA            KTYPE / 1, 2, 5*4, 5*6, 3*9, 10 /
      DATA            KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3, 0 /
      DATA            KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0

      BADMM = .FALSE.
      BADNN = .FALSE.
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      MINWRK = 1
      DO 10 J = 1, NSIZES
         MMAX = MAX( MMAX, MVAL( J ) )
         IF( MVAL( J ).LT.0 ) BADMM = .TRUE.
         NMAX = MAX( NMAX, NVAL( J ) )
         IF( NVAL( J ).LT.0 ) BADNN = .TRUE.
         MNMAX = MAX( MNMAX, MIN( MVAL( J ), NVAL( J ) ) )
         MINWRK = MAX( MINWRK, 3*( MVAL( J )+NVAL( J ) ), MVAL( J )*( MVAL( J )+MAX( MVAL( J ), NVAL( J ), NRHS )+1 )+NVAL( J )*MIN( NVAL( J ), MVAL( J ) ) )
   10 CONTINUE

      // Check for errors

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

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZCHKBD', -INFO )
         RETURN
      END IF

      // Initialize constants

      PATH( 1: 1 ) = 'Zomplex precision'
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

      // Loop over sizes, types

      DO 180 JSIZE = 1, NSIZES
         M = MVAL( JSIZE )
         N = NVAL( JSIZE )
         MNMIN = MIN( M, N )
         AMNINV = ONE / MAX( M, N, 1 )

         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF

         DO 170 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 170

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            DO 30 J = 1, 14
               RESULT( J ) = -ONE
   30       CONTINUE

            UPLO = ' '

            // Compute "A"

            // Control parameters:

            // KMAGN  KMODE        KTYPE
        // =1  O(1)   clustered 1  zero
        // =2  large  clustered 2  identity
        // =3  small  exponential  (none)
        // =4         arithmetic   diagonal, (w/ eigenvalues)
        // =5         random       symmetric, w/ eigenvalues
        // =6                      nonsymmetric, w/ singular values
        // =7                      random diagonal
        // =8                      random symmetric
        // =9                      random nonsymmetric
        // =10                     random bidiagonal (log. distrib.)

            IF( MTYPES.GT.MAXTYP ) GO TO 100

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40       CONTINUE
            ANORM = ONE
            GO TO 70

   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*AMNINV
            GO TO 70

   60       CONTINUE
            ANORM = RTUNFL*MAX( M, N )*ULPINV
            GO TO 70

   70       CONTINUE

            CALL ZLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
            IINFO = 0
            COND = ULPINV

            BIDIAG = .FALSE.
            IF( ITYPE.EQ.1 ) THEN

               // Zero matrix

               IINFO = 0

            ELSE IF( ITYPE.EQ.2 ) THEN

               // Identity

               DO 80 JCOL = 1, MNMIN
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            ELSE IF( ITYPE.EQ.4 ) THEN

               // Diagonal Matrix, [Eigen]values Specified

               CALL ZLATMS( MNMIN, MNMIN, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO )

            ELSE IF( ITYPE.EQ.5 ) THEN

               // Symmetric, eigenvalues specified

               CALL ZLATMS( MNMIN, MNMIN, 'S', ISEED, 'S', RWORK, IMODE, COND, ANORM, M, N, 'N', A, LDA, WORK, IINFO )

            ELSE IF( ITYPE.EQ.6 ) THEN

               // Nonsymmetric, singular values specified

               CALL ZLATMS( M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, M, N, 'N', A, LDA, WORK, IINFO )

            ELSE IF( ITYPE.EQ.7 ) THEN

               // Diagonal, random entries

               CALL ZLATMR( MNMIN, MNMIN, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            ELSE IF( ITYPE.EQ.8 ) THEN

               // Symmetric, random entries

               CALL ZLATMR( MNMIN, MNMIN, 'S', ISEED, 'S', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            ELSE IF( ITYPE.EQ.9 ) THEN

               // Nonsymmetric, random entries

               CALL ZLATMR( M, N, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            ELSE IF( ITYPE.EQ.10 ) THEN

               // Bidiagonal, random entries

               TEMP1 = -TWO*LOG( ULP )
               DO 90 J = 1, MNMIN
                  BD( J ) = EXP( TEMP1*DLARND( 2, ISEED ) )
                  IF( J.LT.MNMIN ) BE( J ) = EXP( TEMP1*DLARND( 2, ISEED ) )
   90          CONTINUE

               IINFO = 0
               BIDIAG = .TRUE.
               IF( M.GE.N ) THEN
                  UPLO = 'U'
               } else {
                  UPLO = 'L'
               END IF
            } else {
               IINFO = 1
            END IF

            IF( IINFO.EQ.0 ) THEN

               // Generate Right-Hand Side

               IF( BIDIAG ) THEN
                  CALL ZLATMR( MNMIN, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, MNMIN, NRHS, ZERO, ONE, 'NO', Y, LDX, IWORK, IINFO )
               } else {
                  CALL ZLATMR( M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( M+1 ), 1, ONE, WORK( 2*M+1 ), 1, ONE, 'N', IWORK, M, NRHS, ZERO, ONE, 'NO', X, LDX, IWORK, IINFO )
               END IF
            END IF

            // Error Exit

            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'Generator', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF

  100       CONTINUE

            // Call ZGEBRD and ZUNGBR to compute B, Q, and P, do tests.

            IF( .NOT.BIDIAG ) THEN

               // Compute transformations to reduce A to bidiagonal form:
               // B := Q' * A * P.

               CALL ZLACPY( ' ', M, N, A, LDA, Q, LDQ )
               CALL ZGEBRD( M, N, Q, LDQ, BD, BE, WORK, WORK( MNMIN+1 ), WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )

               // Check error code from ZGEBRD.

               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9998 )'ZGEBRD', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               CALL ZLACPY( ' ', M, N, Q, LDQ, PT, LDPT )
               IF( M.GE.N ) THEN
                  UPLO = 'U'
               } else {
                  UPLO = 'L'
               END IF

               // Generate Q

               MQ = M
               IF( NRHS.LE.0 ) MQ = MNMIN                CALL ZUNGBR( 'Q', M, MQ, N, Q, LDQ, WORK, WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )

               // Check error code from ZUNGBR.

               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9998 )'ZUNGBR(Q)', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               // Generate P'

               CALL ZUNGBR( 'P', MNMIN, N, M, PT, LDPT, WORK( MNMIN+1 ), WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )

               // Check error code from ZUNGBR.

               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9998 )'ZUNGBR(P)', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               // Apply Q' to an M by NRHS matrix X:  Y := Q' * X.

               CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, NRHS, M, CONE, Q, LDQ, X, LDX, CZERO, Y, LDX )

               // Test 1:  Check the decomposition A := Q * B * PT
                    // 2:  Check the orthogonality of Q
                    // 3:  Check the orthogonality of PT

               CALL ZBDT01( M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT, WORK, RWORK, RESULT( 1 ) )                CALL ZUNT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) )                CALL ZUNT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RWORK, RESULT( 3 ) )
            END IF

            // Use ZBDSQR to form the SVD of the bidiagonal matrix B:
            // B := U * S1 * VT, and compute Z = U' * Y.

            CALL DCOPY( MNMIN, BD, 1, S1, 1 )
            IF( MNMIN.GT.0 ) CALL DCOPY( MNMIN-1, BE, 1, RWORK, 1 )
            CALL ZLACPY( ' ', M, NRHS, Y, LDX, Z, LDX )
            CALL ZLASET( 'Full', MNMIN, MNMIN, CZERO, CONE, U, LDPT )
            CALL ZLASET( 'Full', MNMIN, MNMIN, CZERO, CONE, VT, LDPT )

            CALL ZBDSQR( UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, RWORK, VT, LDPT, U, LDPT, Z, LDX, RWORK( MNMIN+1 ), IINFO )

            // Check error code from ZBDSQR.

            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'ZBDSQR(vects)', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               } else {
                  RESULT( 4 ) = ULPINV
                  GO TO 150
               END IF
            END IF

            // Use ZBDSQR to compute only the singular values of the
            // bidiagonal matrix B;  U, VT, and Z should not be modified.

            CALL DCOPY( MNMIN, BD, 1, S2, 1 )
            IF( MNMIN.GT.0 ) CALL DCOPY( MNMIN-1, BE, 1, RWORK, 1 )

            CALL ZBDSQR( UPLO, MNMIN, 0, 0, 0, S2, RWORK, VT, LDPT, U, LDPT, Z, LDX, RWORK( MNMIN+1 ), IINFO )

            // Check error code from ZBDSQR.

            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9998 )'ZBDSQR(values)', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) THEN
                  RETURN
               } else {
                  RESULT( 9 ) = ULPINV
                  GO TO 150
               END IF
            END IF

            // Test 4:  Check the decomposition B := U * S1 * VT
                 // 5:  Check the computation Z := U' * Y
                 // 6:  Check the orthogonality of U
                 // 7:  Check the orthogonality of VT

            CALL ZBDT03( UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK, RESULT( 4 ) )             CALL ZBDT02( MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK, RWORK, RESULT( 5 ) )             CALL ZUNT01( 'Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RWORK, RESULT( 6 ) )             CALL ZUNT01( 'Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RWORK, RESULT( 7 ) )

            // Test 8:  Check that the singular values are sorted in
                     // non-increasing order and are non-negative

            RESULT( 8 ) = ZERO
            DO 110 I = 1, MNMIN - 1
               IF( S1( I ).LT.S1( I+1 ) ) RESULT( 8 ) = ULPINV                IF( S1( I ).LT.ZERO ) RESULT( 8 ) = ULPINV
  110       CONTINUE
            IF( MNMIN.GE.1 ) THEN
               IF( S1( MNMIN ).LT.ZERO ) RESULT( 8 ) = ULPINV
            END IF

            // Test 9:  Compare ZBDSQR with and without singular vectors

            TEMP2 = ZERO

            DO 120 J = 1, MNMIN
               TEMP1 = ABS( S1( J )-S2( J ) ) / MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ), ULP*MAX( ABS( S1( J ) ), ABS( S2( J ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  120       CONTINUE

            RESULT( 9 ) = TEMP2

            // Test 10:  Sturm sequence test of singular values
                      // Go up by factors of two until it succeeds

            TEMP1 = THRESH*( HALF-ULP )

            DO 130 J = 0, LOG2UI
               CALL DSVDCH( MNMIN, BD, BE, S1, TEMP1, IINFO )
               IF( IINFO.EQ.0 ) GO TO 140
               TEMP1 = TEMP1*TWO
  130       CONTINUE

  140       CONTINUE
            RESULT( 10 ) = TEMP1

            // Use ZBDSQR to form the decomposition A := (QU) S (VT PT)
            // from the bidiagonal form A := Q B PT.

            IF( .NOT.BIDIAG ) THEN
               CALL DCOPY( MNMIN, BD, 1, S2, 1 )
               IF( MNMIN.GT.0 ) CALL DCOPY( MNMIN-1, BE, 1, RWORK, 1 )

               CALL ZBDSQR( UPLO, MNMIN, N, M, NRHS, S2, RWORK, PT, LDPT, Q, LDQ, Y, LDX, RWORK( MNMIN+1 ), IINFO )

               // Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
                    // 12:  Check the computation Z := U' * Q' * X
                    // 13:  Check the orthogonality of Q*U
                    // 14:  Check the orthogonality of VT*PT

               CALL ZBDT01( M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT, LDPT, WORK, RWORK, RESULT( 11 ) )                CALL ZBDT02( M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK, RWORK, RESULT( 12 ) )                CALL ZUNT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK, RWORK, RESULT( 13 ) )                CALL ZUNT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RWORK, RESULT( 14 ) )
            END IF

            // End of Loop -- Check for RESULT(j) > THRESH

  150       CONTINUE
            DO 160 J = 1, 14
               IF( RESULT( J ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 ) CALL DLAHD2( NOUT, PATH )                   WRITE( NOUT, FMT = 9999 )M, N, JTYPE, IOLDSD, J, RESULT( J )
                  NFAIL = NFAIL + 1
               END IF
  160       CONTINUE
            IF( .NOT.BIDIAG ) THEN
               NTEST = NTEST + 14
            } else {
               NTEST = NTEST + 5
            END IF

  170    CONTINUE
  180 CONTINUE

      // Summary

      CALL ALASUM( PATH, NOUT, NFAIL, NTEST, 0 )

      RETURN

      // End of ZCHKBD

 9999 FORMAT( ' M=', I5, ', N=', I5, ', type ', I2, ', seed=',
     $      4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9998 FORMAT( ' ZCHKBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )

      }
