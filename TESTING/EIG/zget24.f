      SUBROUTINE ZGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, H, HT, W, WT, WTMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, ISRT, RESULT, WORK, LWORK, RWORK, BWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               COMP;
      int                INFO, ISRT, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT;
      double             RCDEIN, RCDVIN, THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                ISEED( 4 ), ISLCT( * );
      double             RESULT( 17 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), H( LDA, * ), HT( LDA, * ), VS( LDVS, * ), VS1( LDVS, * ), W( * ), WORK( * ), WT( * ), WTMP( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      double             EPSIN;
      PARAMETER          ( EPSIN = 5.9605D-8 )
      // ..
      // .. Local Scalars ..
      String             SORT;
      int                I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, RSUB, SDIM, SDIM1       double             ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, SMLNUM, TOL, TOLIN, ULP, ULPINV, V, VRICMP, VRIMIN, WNORM;;
      COMPLEX*16         CTMP
      // ..
      // .. Local Arrays ..
      int                IPNT( 20 );
      // ..
      // .. External Functions ..
      bool               ZSLECT;
      double             DLAMCH, ZLANGE;
      // EXTERNAL ZSLECT, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGEESX, ZGEMM, ZLACPY, ZUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Executable Statements ..
*
      // Check for errors
*
      INFO = 0
      IF( THRESH.LT.ZERO ) THEN
         INFO = -3
      ELSE IF( NOUNIT.LE.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.1 .OR. LDA.LT.N ) THEN
         INFO = -8
      ELSE IF( LDVS.LT.1 .OR. LDVS.LT.N ) THEN
         INFO = -15
      ELSE IF( LWORK.LT.2*N ) THEN
         INFO = -24
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGET24', -INFO )
         RETURN
      END IF
*
      // Quick return if nothing to do
*
      DO 10 I = 1, 17
         RESULT( I ) = -ONE
   10 CONTINUE
*
      IF( N.EQ.0 ) RETURN
*
      // Important constants
*
      SMLNUM = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
*
      // Perform tests (1)-(13)
*
      SELOPT = 0
      DO 90 ISORT = 0, 1
         IF( ISORT.EQ.0 ) THEN
            SORT = 'N'
            RSUB = 0
         ELSE
            SORT = 'S'
            RSUB = 6
         END IF
*
         // Compute Schur form and Schur vectors, and test them
*
         CALL ZLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL ZGEESX( 'V', SORT, ZSLECT, 'N', N, H, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 1+RSUB ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX1', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX1', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            RETURN
         END IF
         IF( ISORT.EQ.0 ) THEN
            CALL ZCOPY( N, W, 1, WTMP, 1 )
         END IF
*
         // Do Test (1) or Test (7)
*
         RESULT( 1+RSUB ) = ZERO
         DO 30 J = 1, N - 1
            DO 20 I = J + 1, N
               IF( H( I, J ).NE.CZERO ) RESULT( 1+RSUB ) = ULPINV
   20       CONTINUE
   30    CONTINUE
*
         // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)
*
         // Copy A to VS1, used as workspace
*
         CALL ZLACPY( ' ', N, N, A, LDA, VS1, LDVS )
*
         // Compute Q*H and store in HT.
*
         CALL ZGEMM( 'No transpose', 'No transpose', N, N, N, CONE, VS, LDVS, H, LDA, CZERO, HT, LDA )
*
         // Compute A - Q*H*Q'
*
         CALL ZGEMM( 'No transpose', 'Conjugate transpose', N, N, N, -CONE, HT, LDA, VS, LDVS, CONE, VS1, LDVS )
*
         ANORM = MAX( ZLANGE( '1', N, N, A, LDA, RWORK ), SMLNUM )
         WNORM = ZLANGE( '1', N, N, VS1, LDVS, RWORK )
*
         IF( ANORM.GT.WNORM ) THEN
            RESULT( 2+RSUB ) = ( WNORM / ANORM ) / ( N*ULP )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESULT( 2+RSUB ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            ELSE
               RESULT( 2+RSUB ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
            END IF
         END IF
*
         // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )
*
         CALL ZUNT01( 'Columns', N, N, VS, LDVS, WORK, LWORK, RWORK, RESULT( 3+RSUB ) )
*
         // Do Test (4) or Test (10)
*
         RESULT( 4+RSUB ) = ZERO
         DO 40 I = 1, N
            IF( H( I, I ).NE.W( I ) ) RESULT( 4+RSUB ) = ULPINV
   40    CONTINUE
*
         // Do Test (5) or Test (11)
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'N', SORT, ZSLECT, 'N', N, HT, LDA, SDIM, WT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 5+RSUB ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX2', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX2', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         RESULT( 5+RSUB ) = ZERO
         DO 60 J = 1, N
            DO 50 I = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV
   50       CONTINUE
   60    CONTINUE
*
         // Do Test (6) or Test (12)
*
         RESULT( 6+RSUB ) = ZERO
         DO 70 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 6+RSUB ) = ULPINV
   70    CONTINUE
*
         // Do Test (13)
*
         IF( ISORT.EQ.1 ) THEN
            RESULT( 13 ) = ZERO
            KNTEIG = 0
            DO 80 I = 1, N
               IF( ZSLECT( W( I ) ) ) KNTEIG = KNTEIG + 1
               IF( I.LT.N ) THEN
                  IF( ZSLECT( W( I+1 ) ) .AND. ( .NOT.ZSLECT( W( I ) ) ) )RESULT( 13 ) = ULPINV
               END IF
   80       CONTINUE
            IF( SDIM.NE.KNTEIG ) RESULT( 13 ) = ULPINV
         END IF
*
   90 CONTINUE
*
      // If there is enough workspace, perform tests (14) and (15)
      // as well as (10) through (13)
*
      IF( LWORK.GE.( N*( N+1 ) ) / 2 ) THEN
*
         // Compute both RCONDE and RCONDV with VS
*
         SORT = 'S'
         RESULT( 14 ) = ZERO
         RESULT( 15 ) = ZERO
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'V', SORT, ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 14 ) = ULPINV
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX3', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX3', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 110 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 10 ) = ULPINV
            DO 100 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  100       CONTINUE
  110    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute both RCONDE and RCONDV without VS, and compare
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'N', SORT, ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 14 ) = ULPINV
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX4', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX4', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         // Perform tests (14) and (15)
*
         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV          IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 130 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 10 ) = ULPINV
            DO 120 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  120       CONTINUE
  130    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDE with VS, and compare
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'V', SORT, ZSLECT, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 14 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX5', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX5', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         // Perform test (14)
*
         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 150 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 10 ) = ULPINV
            DO 140 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  140       CONTINUE
  150    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDE without VS, and compare
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'N', SORT, ZSLECT, 'E', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 14 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX6', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX6', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         // Perform test (14)
*
         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 170 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 10 ) = ULPINV
            DO 160 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  160       CONTINUE
  170    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDV with VS, and compare
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'V', SORT, ZSLECT, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX7', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX7', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         // Perform test (15)
*
         IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 190 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 10 ) = ULPINV
            DO 180 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  180       CONTINUE
  190    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDV without VS, and compare
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'N', SORT, ZSLECT, 'V', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'ZGEESX8', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'ZGEESX8', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 220
         END IF
*
         // Perform test (15)
*
         IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 210 I = 1, N
            IF( W( I ).NE.WT( I ) ) RESULT( 10 ) = ULPINV
            DO 200 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  200       CONTINUE
  210    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
      END IF
*
  220 CONTINUE
*
      // If there are precomputed reciprocal condition numbers, compare
      // computed values with them.
*
      IF( COMP ) THEN
*
         // First set up SELOPT, SELDIM, SELVAL, SELWR and SELWI so that
        t // he logical function ZSLECT selects the eigenvalues specified
         // by NSLCT, ISLCT and ISRT.
*
         SELDIM = N
         SELOPT = 1
         EPS = MAX( ULP, EPSIN )
         DO 230 I = 1, N
            IPNT( I ) = I
            SELVAL( I ) = .FALSE.
            SELWR( I ) = DBLE( WTMP( I ) )
            SELWI( I ) = DIMAG( WTMP( I ) )
  230    CONTINUE
         DO 250 I = 1, N - 1
            KMIN = I
            IF( ISRT.EQ.0 ) THEN
               VRIMIN = DBLE( WTMP( I ) )
            ELSE
               VRIMIN = DIMAG( WTMP( I ) )
            END IF
            DO 240 J = I + 1, N
               IF( ISRT.EQ.0 ) THEN
                  VRICMP = DBLE( WTMP( J ) )
               ELSE
                  VRICMP = DIMAG( WTMP( J ) )
               END IF
               IF( VRICMP.LT.VRIMIN ) THEN
                  KMIN = J
                  VRIMIN = VRICMP
               END IF
  240       CONTINUE
            CTMP = WTMP( KMIN )
            WTMP( KMIN ) = WTMP( I )
            WTMP( I ) = CTMP
            ITMP = IPNT( I )
            IPNT( I ) = IPNT( KMIN )
            IPNT( KMIN ) = ITMP
  250    CONTINUE
         DO 260 I = 1, NSLCT
            SELVAL( IPNT( ISLCT( I ) ) ) = .TRUE.
  260    CONTINUE
*
         // Compute condition numbers
*
         CALL ZLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL ZGEESX( 'N', 'S', ZSLECT, 'B', N, HT, LDA, SDIM1, WT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 16 ) = ULPINV
            RESULT( 17 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'ZGEESX9', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 270
         END IF
*
         // Compare condition number for average of selected eigenvalues
        t // aking its condition number into account
*
         ANORM = ZLANGE( '1', N, N, A, LDA, RWORK )
         V = MAX( DBLE( N )*EPS*ANORM, SMLNUM )
         IF( ANORM.EQ.ZERO ) V = ONE
         IF( V.GT.RCONDV ) THEN
            TOL = ONE
         ELSE
            TOL = V / RCONDV
         END IF
         IF( V.GT.RCDVIN ) THEN
            TOLIN = ONE
         ELSE
            TOLIN = V / RCDVIN
         END IF
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         IF( EPS*( RCDEIN-TOLIN ).GT.RCONDE+TOL ) THEN
            RESULT( 16 ) = ULPINV
         ELSE IF( RCDEIN-TOLIN.GT.RCONDE+TOL ) THEN
            RESULT( 16 ) = ( RCDEIN-TOLIN ) / ( RCONDE+TOL )
         ELSE IF( RCDEIN+TOLIN.LT.EPS*( RCONDE-TOL ) ) THEN
            RESULT( 16 ) = ULPINV
         ELSE IF( RCDEIN+TOLIN.LT.RCONDE-TOL ) THEN
            RESULT( 16 ) = ( RCONDE-TOL ) / ( RCDEIN+TOLIN )
         ELSE
            RESULT( 16 ) = ONE
         END IF
*
         // Compare condition numbers for right invariant subspace
        t // aking its condition number into account
*
         IF( V.GT.RCONDV*RCONDE ) THEN
            TOL = RCONDV
         ELSE
            TOL = V / RCONDE
         END IF
         IF( V.GT.RCDVIN*RCDEIN ) THEN
            TOLIN = RCDVIN
         ELSE
            TOLIN = V / RCDEIN
         END IF
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         IF( EPS*( RCDVIN-TOLIN ).GT.RCONDV+TOL ) THEN
            RESULT( 17 ) = ULPINV
         ELSE IF( RCDVIN-TOLIN.GT.RCONDV+TOL ) THEN
            RESULT( 17 ) = ( RCDVIN-TOLIN ) / ( RCONDV+TOL )
         ELSE IF( RCDVIN+TOLIN.LT.EPS*( RCONDV-TOL ) ) THEN
            RESULT( 17 ) = ULPINV
         ELSE IF( RCDVIN+TOLIN.LT.RCONDV-TOL ) THEN
            RESULT( 17 ) = ( RCONDV-TOL ) / ( RCDVIN+TOLIN )
         ELSE
            RESULT( 17 ) = ONE
         END IF
*
  270    CONTINUE
*
      END IF
*
 9999 FORMAT( ' ZGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' ZGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
      // End of ZGET24
*
      END
