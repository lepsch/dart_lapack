      SUBROUTINE DGET24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA, H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS, LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT, RESULT, WORK, LWORK, IWORK, BWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               COMP;
      int                INFO, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT;
      double             RCDEIN, RCDVIN, THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                ISEED( 4 ), ISLCT( * ), IWORK( * );
      double             A( LDA, * ), H( LDA, * ), HT( LDA, * ), RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ), WI( * ), WIT( * ), WITMP( * ), WORK( * ), WR( * ), WRT( * ), WRTMP( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      double             EPSIN;
      const              EPSIN = 5.9605D-8 ;
      // ..
      // .. Local Scalars ..
      String             SORT;
      int                I, IINFO, ISORT, ITMP, J, KMIN, KNTEIG, LIWORK, RSUB, SDIM, SDIM1;
      double             ANORM, EPS, RCNDE1, RCNDV1, RCONDE, RCONDV, SMLNUM, TMP, TOL, TOLIN, ULP, ULPINV, V, VIMIN, VRMIN, WNORM;
      // ..
      // .. Local Arrays ..
      int                IPNT( 20 );
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
      // .. External Functions ..
      bool               DSLECT;
      double             DLAMCH, DLANGE;
      // EXTERNAL DSLECT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEESX, DGEMM, DLACPY, DORT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0
      if ( THRESH.LT.ZERO ) {
         INFO = -3
      } else if ( NOUNIT.LE.0 ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.1 .OR. LDA.LT.N ) {
         INFO = -8
      } else if ( LDVS.LT.1 .OR. LDVS.LT.N ) {
         INFO = -18
      } else if ( LWORK.LT.3*N ) {
         INFO = -26
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DGET24', -INFO )
         RETURN
      }

      // Quick return if nothing to do

      DO 10 I = 1, 17
         RESULT( I ) = -ONE
   10 CONTINUE

      IF( N.EQ.0 ) RETURN

      // Important constants

      SMLNUM = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP

      // Perform tests (1)-(13)

      SELOPT = 0
      LIWORK = N*N
      DO 120 ISORT = 0, 1
         if ( ISORT.EQ.0 ) {
            SORT = 'N'
            RSUB = 0
         } else {
            SORT = 'S'
            RSUB = 6
         }

         // Compute Schur form and Schur vectors, and test them

         CALL DLACPY( 'F', N, N, A, LDA, H, LDA )
         CALL DGEESX( 'V', SORT, DSLECT, 'N', N, H, LDA, SDIM, WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 1+RSUB ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX1', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX1', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            RETURN
         }
         if ( ISORT.EQ.0 ) {
            CALL DCOPY( N, WR, 1, WRTMP, 1 )
            CALL DCOPY( N, WI, 1, WITMP, 1 )
         }

         // Do Test (1) or Test (7)

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
            if ( H( I+1, I ).NE.ZERO ) {
               IF( H( I, I ).NE.H( I+1, I+1 ) .OR. H( I, I+1 ).EQ. ZERO .OR. SIGN( ONE, H( I+1, I ) ).EQ. SIGN( ONE, H( I, I+1 ) ) )RESULT( 1+RSUB ) = ULPINV
            }
   50    CONTINUE

         // Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)

         // Copy A to VS1, used as workspace

         CALL DLACPY( ' ', N, N, A, LDA, VS1, LDVS )

         // Compute Q*H and store in HT.

         CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE, VS, LDVS, H, LDA, ZERO, HT, LDA )

         // Compute A - Q*H*Q'

         CALL DGEMM( 'No transpose', 'Transpose', N, N, N, -ONE, HT, LDA, VS, LDVS, ONE, VS1, LDVS )

         ANORM = MAX( DLANGE( '1', N, N, A, LDA, WORK ), SMLNUM )
         WNORM = DLANGE( '1', N, N, VS1, LDVS, WORK )

         if ( ANORM.GT.WNORM ) {
            RESULT( 2+RSUB ) = ( WNORM / ANORM ) / ( N*ULP )
         } else {
            if ( ANORM.LT.ONE ) {
               RESULT( 2+RSUB ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            } else {
               RESULT( 2+RSUB ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
            }
         }

         // Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )

         CALL DORT01( 'Columns', N, N, VS, LDVS, WORK, LWORK, RESULT( 3+RSUB ) )

         // Do Test (4) or Test (10)

         RESULT( 4+RSUB ) = ZERO
         DO 60 I = 1, N
            IF( H( I, I ).NE.WR( I ) ) RESULT( 4+RSUB ) = ULPINV
   60    CONTINUE
         if ( N.GT.1 ) {
            IF( H( 2, 1 ).EQ.ZERO .AND. WI( 1 ).NE.ZERO ) RESULT( 4+RSUB ) = ULPINV             IF( H( N, N-1 ).EQ.ZERO .AND. WI( N ).NE.ZERO ) RESULT( 4+RSUB ) = ULPINV
         }
         DO 70 I = 1, N - 1
            if ( H( I+1, I ).NE.ZERO ) {
               TMP = SQRT( ABS( H( I+1, I ) ) )* SQRT( ABS( H( I, I+1 ) ) )                RESULT( 4+RSUB ) = MAX( RESULT( 4+RSUB ), ABS( WI( I )-TMP ) / MAX( ULP*TMP, SMLNUM ) )                RESULT( 4+RSUB ) = MAX( RESULT( 4+RSUB ), ABS( WI( I+1 )+TMP ) / MAX( ULP*TMP, SMLNUM ) )
            } else if ( I.GT.1 ) {
               IF( H( I+1, I ).EQ.ZERO .AND. H( I, I-1 ).EQ.ZERO .AND. WI( I ).NE.ZERO )RESULT( 4+RSUB ) = ULPINV
            }
   70    CONTINUE

         // Do Test (5) or Test (11)

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'N', SORT, DSLECT, 'N', N, HT, LDA, SDIM, WRT, WIT, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 5+RSUB ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX2', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX2', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         RESULT( 5+RSUB ) = ZERO
         DO 90 J = 1, N
            DO 80 I = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 5+RSUB ) = ULPINV
   80       CONTINUE
   90    CONTINUE

         // Do Test (6) or Test (12)

         RESULT( 6+RSUB ) = ZERO
         DO 100 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 6+RSUB ) = ULPINV
  100    CONTINUE

         // Do Test (13)

         if ( ISORT.EQ.1 ) {
            RESULT( 13 ) = ZERO
            KNTEIG = 0
            DO 110 I = 1, N
               IF( DSLECT( WR( I ), WI( I ) ) .OR. DSLECT( WR( I ), -WI( I ) ) )KNTEIG = KNTEIG + 1
               if ( I.LT.N ) {
                  IF( ( DSLECT( WR( I+1 ), WI( I+1 ) ) .OR. DSLECT( WR( I+1 ), -WI( I+1 ) ) ) .AND. ( .NOT.( DSLECT( WR( I ), WI( I ) ) .OR. DSLECT( WR( I ), -WI( I ) ) ) ) .AND. IINFO.NE.N+2 )RESULT( 13 ) = ULPINV
               }
  110       CONTINUE
            IF( SDIM.NE.KNTEIG ) RESULT( 13 ) = ULPINV
         }

  120 CONTINUE

      // If there is enough workspace, perform tests (14) and (15)
      // as well as (10) through (13)

      if ( LWORK.GE.N+( N*N ) / 2 ) {

         // Compute both RCONDE and RCONDV with VS

         SORT = 'S'
         RESULT( 14 ) = ZERO
         RESULT( 15 ) = ZERO
         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'V', SORT, DSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 14 ) = ULPINV
            RESULT( 15 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX3', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX3', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Perform tests (10), (11), (12), and (13)

         DO 140 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 130 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  130       CONTINUE
  140    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV

         // Compute both RCONDE and RCONDV without VS, and compare

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'N', SORT, DSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 14 ) = ULPINV
            RESULT( 15 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX4', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX4', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Perform tests (14) and (15)

         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV          IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV

         // Perform tests (10), (11), (12), and (13)

         DO 160 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 150 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  150       CONTINUE
  160    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV

         // Compute RCONDE with VS, and compare

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'V', SORT, DSLECT, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 14 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX5', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX5', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Perform test (14)

         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV

         // Perform tests (10), (11), (12), and (13)

         DO 180 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 170 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  170       CONTINUE
  180    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV

         // Compute RCONDE without VS, and compare

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'N', SORT, DSLECT, 'E', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 14 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX6', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX6', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Perform test (14)

         IF( RCNDE1.NE.RCONDE ) RESULT( 14 ) = ULPINV

         // Perform tests (10), (11), (12), and (13)

         DO 200 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 190 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  190       CONTINUE
  200    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV

         // Compute RCONDV with VS, and compare

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'V', SORT, DSLECT, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 15 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX7', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX7', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Perform test (15)

         IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV

         // Perform tests (10), (11), (12), and (13)

         DO 220 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 210 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  210       CONTINUE
  220    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV

         // Compute RCONDV without VS, and compare

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'N', SORT, DSLECT, 'V', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCNDE1, RCNDV1, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 15 ) = ULPINV
            if ( JTYPE.NE.22 ) {
               WRITE( NOUNIT, FMT = 9998 )'DGEESX8', IINFO, N, JTYPE, ISEED
            } else {
               WRITE( NOUNIT, FMT = 9999 )'DGEESX8', IINFO, N, ISEED( 1 )
            }
            INFO = ABS( IINFO )
            GO TO 250
         }

         // Perform test (15)

         IF( RCNDV1.NE.RCONDV ) RESULT( 15 ) = ULPINV

         // Perform tests (10), (11), (12), and (13)

         DO 240 I = 1, N
            IF( WR( I ).NE.WRT( I ) .OR. WI( I ).NE.WIT( I ) ) RESULT( 10 ) = ULPINV
            DO 230 J = 1, N
               IF( H( I, J ).NE.HT( I, J ) ) RESULT( 11 ) = ULPINV                IF( VS( I, J ).NE.VS1( I, J ) ) RESULT( 12 ) = ULPINV
  230       CONTINUE
  240    CONTINUE
         IF( SDIM.NE.SDIM1 ) RESULT( 13 ) = ULPINV

      }

  250 CONTINUE

      // If there are precomputed reciprocal condition numbers, compare
      // computed values with them.

      if ( COMP ) {

         // First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
         // the logical function DSLECT selects the eigenvalues specified
         // by NSLCT and ISLCT.

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
               if ( WRTMP( J ).LT.VRMIN ) {
                  KMIN = J
                  VRMIN = WRTMP( J )
                  VIMIN = WITMP( J )
               }
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

         // Compute condition numbers

         CALL DLACPY( 'F', N, N, A, LDA, HT, LDA )
         CALL DGEESX( 'N', 'S', DSLECT, 'B', N, HT, LDA, SDIM1, WRT, WIT, VS1, LDVS, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, IINFO )
         if ( IINFO.NE.0 .AND. IINFO.NE.N+2 ) {
            RESULT( 16 ) = ULPINV
            RESULT( 17 ) = ULPINV
            WRITE( NOUNIT, FMT = 9999 )'DGEESX9', IINFO, N, ISEED( 1 )
            INFO = ABS( IINFO )
            GO TO 300
         }

         // Compare condition number for average of selected eigenvalues
         // taking its condition number into account

         ANORM = DLANGE( '1', N, N, A, LDA, WORK )
         V = MAX( DBLE( N )*EPS*ANORM, SMLNUM )
         IF( ANORM.EQ.ZERO ) V = ONE
         if ( V.GT.RCONDV ) {
            TOL = ONE
         } else {
            TOL = V / RCONDV
         }
         if ( V.GT.RCDVIN ) {
            TOLIN = ONE
         } else {
            TOLIN = V / RCDVIN
         }
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         if ( EPS*( RCDEIN-TOLIN ).GT.RCONDE+TOL ) {
            RESULT( 16 ) = ULPINV
         } else if ( RCDEIN-TOLIN.GT.RCONDE+TOL ) {
            RESULT( 16 ) = ( RCDEIN-TOLIN ) / ( RCONDE+TOL )
         } else if ( RCDEIN+TOLIN.LT.EPS*( RCONDE-TOL ) ) {
            RESULT( 16 ) = ULPINV
         } else if ( RCDEIN+TOLIN.LT.RCONDE-TOL ) {
            RESULT( 16 ) = ( RCONDE-TOL ) / ( RCDEIN+TOLIN )
         } else {
            RESULT( 16 ) = ONE
         }

         // Compare condition numbers for right invariant subspace
         // taking its condition number into account

         if ( V.GT.RCONDV*RCONDE ) {
            TOL = RCONDV
         } else {
            TOL = V / RCONDE
         }
         if ( V.GT.RCDVIN*RCDEIN ) {
            TOLIN = RCDVIN
         } else {
            TOLIN = V / RCDEIN
         }
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         if ( EPS*( RCDVIN-TOLIN ).GT.RCONDV+TOL ) {
            RESULT( 17 ) = ULPINV
         } else if ( RCDVIN-TOLIN.GT.RCONDV+TOL ) {
            RESULT( 17 ) = ( RCDVIN-TOLIN ) / ( RCONDV+TOL )
         } else if ( RCDVIN+TOLIN.LT.EPS*( RCONDV-TOL ) ) {
            RESULT( 17 ) = ULPINV
         } else if ( RCDVIN+TOLIN.LT.RCONDV-TOL ) {
            RESULT( 17 ) = ( RCONDV-TOL ) / ( RCDVIN+TOLIN )
         } else {
            RESULT( 17 ) = ONE
         }

  300    CONTINUE

      }

 9999 FORMAT( ' DGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', INPUT EXAMPLE NUMBER = ', I4 )
 9998 FORMAT( ' DGET24: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of DGET24

      }
