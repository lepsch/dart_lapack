      SUBROUTINE SGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, RESULT, WORK, LWORK, IWORK, BWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               COMP;
      int                INFO, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT;
      REAL               RCDEIN, RCDVIN, THRESH
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                ISEED( 4 ), ISLCT( * ), IWORK( * );
      REAL               A( LDA, * ), H( LDA, * ), HT( LDA, * ), RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ), WI( * ), WIT( * ), WITMP( * ), WORK( * ), WR( * ), WRT( * ), WRTMP( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      REAL               EPSIN
      PARAMETER          ( EPSIN = 5.9605E-8 )
      // ..
      // .. Local Scalars ..
      String             SORT;
      int                I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, LIWORK, RSUB, SDIM, SDIM1       REAL               ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, SMLNUM, TMP, TOL, TOLIN, ULP, ULPINV, V, VIMIN, VRMIN, WNORM;
      // ..
      // .. Local Arrays ..
      int                IPNT( 20 );
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      REAL               SELWI( 20 ), SELWR( 20 )
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. External Functions ..
      bool               SSLECT;
      REAL               SLAMCH, SLANGE
      // EXTERNAL SSLECT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEESX, SGEMM, SLACPY, SORT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SIGN, SQRT
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
         INFO = -18
      ELSE IF( LWORK.LT.3*N ) THEN
         INFO = -26
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGET24', -INFO )
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
      SMLNUM = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )
      ULPINV = ONE / ULP
*
      // Perform tests (1)-(13)
*
      SELOPT = 0
      LIWORK = N*N
      DO 120 ISORT = 0, 1
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
         CALL SLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL SGEESX( 'V', SORT, SSLECT, 'N', N, H, LDA, SDIM, WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 1+RSUB ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX1', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX1', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            RETURN
         END IF
         IF( ISORT.EQ.0 ) THEN
            CALL SCOPY( N, WR, 1, WRTMP, 1 )
            CALL SCOPY( N, WI, 1, WITMP, 1 )
         END IF
*
         // Do Test (1) or Test (7)
*
         RESULT( 1+RSUB ) = ZERO
         DO 30 J = 1, N - 2
            DO 20 I = J + 2, N
               IF( H( I, J ).NE.ZERO ) RESULT( 1+RSUB ) = ULPINV
   20       CONTINUE
   30    CONTINUE
         DO 40 I = 1, N - 2
            IF( H( I+1, I ).NE.ZERO .AND. H( I+2, I+1 ).NE.ZERO ) RESULT( 1+RSUB ) = ULPINV
   40    CONTINUE
         DO 50 I = 1, N - 1
            IF( H( I+1, I ).NE.ZERO ) THEN
               IF( H( I, I ).NE.H( I+1, I+1 ) .OR. H( I, I+1 ).EQ. ZERO .OR. SIGN( ONE, H( I+1, I ) ).EQ. SIGN( ONE, H( I, I+1 ) ) )RESULT( 1+RSUB ) = ULPINV
            END IF
   50    CONTINUE
*
         // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)
*
         // Copy A to VS1, used as workspace
*
         CALL SLACPY( ' ', N, N, A, LDA, VS1, LDVS )
*
         // Compute Q*H and store in HT.
*
         CALL SGEMM( 'No transpose', 'No transpose', N, N, N, ONE, VS, LDVS, H, LDA, ZERO, HT, LDA )
*
         // Compute A - Q*H*Q'
*
         CALL SGEMM( 'No transpose', 'Transpose', N, N, N, -ONE, HT, LDA, VS, LDVS, ONE, VS1, LDVS )
*
         ANORM = MAX( SLANGE( '1', N, N, A, LDA, WORK ), SMLNUM )
         WNORM = SLANGE( '1', N, N, VS1, LDVS, WORK )
*
         IF( ANORM.GT.WNORM ) THEN
            RESULT( 2+RSUB ) = ( WNORM / ANORM ) / ( N*ULP )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESULT( 2+RSUB ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            ELSE
               RESULT( 2+RSUB ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
            END IF
         END IF
*
         // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )
*
         CALL SORT01( 'Columns', N, N, VS, LDVS, WORK, LWORK, RESULT( 3+RSUB ) )
*
         // Do Test (4) or Test (10)
*
         RESULT( 4+RSUB ) = ZERO
         DO 60 I = 1, N
            IF( H( I, I ).NE.WR( I ) ) RESULT( 4+RSUB ) = ULPINV
   60    CONTINUE
         IF( N.GT.1 ) THEN
            IF( H( 2, 1 ).EQ.ZERO .AND. WI( 1 ).NE.ZERO ) RESULT( 4+RSUB ) = ULPINV             IF( H( N, N-1 ).EQ.ZERO .AND. WI( N ).NE.ZERO ) RESULT( 4+RSUB ) = ULPINV
         END IF
         DO 70 I = 1, N - 1
            IF( H( I+1, I ).NE.ZERO ) THEN
               TMP = SQRT( ABS( H( I+1, I ) ) )* SQRT( ABS( H( I, I+1 ) ) )                RESULT( 4+RSUB ) = MAX( RESULT( 4+RSUB ), ABS( WI( I )-TMP ) / MAX( ULP*TMP, SMLNUM ) )                RESULT( 4+RSUB ) = MAX( RESULT( 4+RSUB ), ABS( WI( I+1 )+TMP ) / MAX( ULP*TMP, SMLNUM ) )
            ELSE IF( I.GT.1 ) THEN
               IF( H( I+1, I ).EQ.ZERO .AND. H( I, I-1 ).EQ.ZERO .AND. WI( I ).NE.ZERO )RESULT( 4+RSUB ) = ULPINV
            END IF
   70    CONTINUE
*
         // Do Test (5) or Test (11)
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'N', SORT, SSLECT, 'N', N, HT, LDA, SDIM, WRT, WIT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 5+RSUB ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX2', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX2', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         RESULT( 5+RSUB ) = ZERO
         DO 90 J = 1, N
            DO 80 I = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV
   80       CONTINUE
   90    CONTINUE
*
         // Do Test (6) or Test (12)
*
         RESULT( 6+RSUB ) = ZERO
         DO 100 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 6+RSUB ) = ULPINV
  100    CONTINUE
*
         // Do Test (13)
*
         IF( ISORT.EQ.1 ) THEN
            RESULT( 13 ) = ZERO
            KNTEIG = 0
            DO 110 I = 1, N
               IF( SSLECT( WR( I ), WI( I ) ) .OR. SSLECT( WR( I ), -WI( I ) ) )KNTEIG = KNTEIG + 1
               IF( I.LT.N ) THEN
                  IF( ( SSLECT( WR( I+1 ), WI( I+1 ) ) .OR. SSLECT( WR( I+1 ), -WI( I+1 ) ) ) .AND. ( .NOT.( SSLECT( WR( I ), WI( I ) ) .OR. SSLECT( WR( I ), -WI( I ) ) ) ) .AND. IINFO.NE.N+2 )RESULT( 13 ) = ULPINV
               END IF
  110       CONTINUE
            IF( SDIM.NE.KNTEIG ) RESULT( 13 ) = ULPINV
         END IF
*
  120 CONTINUE
*
      // If there is enough workspace, perform tests (14) and (15)
      // as well as (10) through (13)
*
      IF( LWORK.GE.N+( N*N ) / 2 ) THEN
*
         // Compute both RCONDE and RCONDV with VS
*
         SORT = 'S'
         RESULT( 14 ) = ZERO
         RESULT( 15 ) = ZERO
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'V', SORT, SSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 14 ) = ULPINV
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX3', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX3', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 140 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 130 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  130       CONTINUE
  140    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute both RCONDE and RCONDV without VS, and compare
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'N', SORT, SSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 14 ) = ULPINV
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX4', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX4', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         // Perform tests (14) and (15)
*
         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV          IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 160 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 150 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  150       CONTINUE
  160    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDE with VS, and compare
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'V', SORT, SSLECT, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 14 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX5', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX5', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         // Perform test (14)
*
         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 180 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 170 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  170       CONTINUE
  180    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDE without VS, and compare
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'N', SORT, SSLECT, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 14 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX6', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX6', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         // Perform test (14)
*
         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 200 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 190 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  190       CONTINUE
  200    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDV with VS, and compare
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'V', SORT, SSLECT, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX7', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX7', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         // Perform test (15)
*
         IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 220 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 210 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  210       CONTINUE
  220    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
         // Compute RCONDV without VS, and compare
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'N', SORT, SSLECT, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 15 ) = ULPINV
            IF( JTYPE.NE.22 ) THEN
               WRITE( NOUNIT, FMT = 9998 )'SGEESX8', IINFO, N, JTYPE, ISEED
            ELSE
               WRITE( NOUNIT, FMT = 9999 )'SGEESX8', IINFO, N, ISEED( 1 )
            END IF
            INFO = ABS( IINFO )
            GO TO 250
         END IF
*
         // Perform test (15)
*
         IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV
*
         // Perform tests (10), (11), (12), and (13)
*
         DO 240 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 230 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  230       CONTINUE
  240    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV
*
      END IF
*
  250 CONTINUE
*
      // If there are precomputed reciprocal condition numbers, compare
      // computed values with them.
*
      IF( COMP ) THEN
*
         // First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
        t // he logical function SSLECT selects the eigenvalues specified
         // by NSLCT and ISLCT.
*
         SELDIM = N
         SELOPT = 1
         EPS = MAX( ULP, EPSIN )
         DO 260 I = 1, N
            IPNT( I ) = I
            SELVAL( I ) = .FALSE.
            SELWR( I ) = WRTMP( I )
            SELWI( I ) = WITMP( I )
  260    CONTINUE
         DO 280 I = 1, N - 1
            KMIN = I
            VRMIN = WRTMP( I )
            VIMIN = WITMP( I )
            DO 270 J = I + 1, N
               IF( WRTMP( J ).LT.VRMIN ) THEN
                  KMIN = J
                  VRMIN = WRTMP( J )
                  VIMIN = WITMP( J )
               END IF
  270       CONTINUE
            WRTMP( KMIN ) = WRTMP( I )
            WITMP( KMIN ) = WITMP( I )
            WRTMP( I ) = VRMIN
            WITMP( I ) = VIMIN
            ITMP = IPNT( I )
            IPNT( I ) = IPNT( KMIN )
            IPNT( KMIN ) = ITMP
  280    CONTINUE
         DO 290 I = 1, NSLCT
            SELVAL( IPNT( ISLCT( I ) ) ) = .TRUE.
  290    CONTINUE
*
         // Compute condition numbers
*
         CALL SLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL SGEESX( 'N', 'S', SSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         IF( IINFO.NE.0 .AND. IINFO.NE.N+2 ) THEN
            RESULT( 16 ) = ULPINV
            RESULT( 17 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'SGEESX9', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 300
         END IF
*
         // Compare condition number for average of selected eigenvalues
        t // aking its condition number into account
*
         ANORM = SLANGE( '1', N, N, A, LDA, WORK )
         V = MAX( REAL( N )*EPS*ANORM, SMLNUM )
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
  300    CONTINUE
*
      END IF
*
 9999 FORMAT( ' SGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' SGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
      // End of SGET24
*
      END
