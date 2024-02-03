*
*  -- LAPACK test routine --
      SUBROUTINE CSYL01( THRESH, NFAIL, RMAX, NINFO, KNT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                KNT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                NFAIL( 3 ), NINFO( 2 );
      REAL               RMAX( 2 )
      // ..
*
*  =====================================================================
      // ..
      // .. Parameters ..
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               ONE, ZERO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      int                MAXM, MAXN, LDSWORK;
      PARAMETER          ( MAXM = 101, MAXN = 138, LDSWORK = 18 )
      // ..
      // .. Local Scalars ..
      String             TRANA, TRANB;
      int                I, INFO, IINFO, ISGN, ITRANA, ITRANB, J, KLA, KUA, KLB, KUB, M, N       REAL               ANRM, BNRM, BIGNUM, EPS, RES, RES1, SCALE, SCALE3, SMLNUM, TNRM, XNRM;
      COMPLEX            RMUL
      // ..
      // .. Local Arrays ..
      COMPLEX            DUML( MAXM ), DUMR( MAXN ), D( MAX( MAXM, MAXN ) )
      REAL               DUM( MAXN ), VM( 2 )
      int                ISEED( 4 ), IWORK( MAXM + MAXN + 2 );
      // ..
      // .. Allocatable Arrays ..
      int                AllocateStatus;
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: A, B, C, CC, X
      REAL,    DIMENSION(:,:), ALLOCATABLE :: SWORK
      // ..
      // .. External Functions ..
      bool               SISNAN;
      REAL               SLAMCH, CLANGE
      // EXTERNAL SISNAN, SLAMCH, CLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLATMR, CLACPY, CGEMM, CTRSYL, CTRSYL3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX
      // ..
      // .. Allocate memory dynamically ..
      ALLOCATE ( A( MAXM, MAXM ), STAT = AllocateStatus )
      IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
      ALLOCATE ( B( MAXN, MAXN ), STAT = AllocateStatus )
      IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
      ALLOCATE ( C( MAXM, MAXN ), STAT = AllocateStatus )
      IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
      ALLOCATE ( CC( MAXM, MAXN ), STAT = AllocateStatus )
      IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
      ALLOCATE ( X( MAXM, MAXN ), STAT = AllocateStatus )
      IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
      ALLOCATE ( SWORK( LDSWORK, 54 ), STAT = AllocateStatus )
      IF( AllocateStatus /= 0 ) STOP "*** Not enough memory ***"
      // ..
      // .. Executable Statements ..
*
      // Get machine parameters
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
*
      // Expect INFO = 0
      VM( 1 ) = ONE
      // Expect INFO = 1
      VM( 2 ) = 0.5E+0
*
      // Begin test loop
*
      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
      NFAIL( 1 ) = 0
      NFAIL( 2 ) = 0
      NFAIL( 3 ) = 0
      RMAX( 1 ) = ZERO
      RMAX( 2 ) = ZERO
      KNT = 0
      ISEED( 1 ) = 1
      ISEED( 2 ) = 1
      ISEED( 3 ) = 1
      ISEED( 4 ) = 1
      SCALE = ONE
      SCALE3 = ONE
      DO J = 1, 2
         DO ISGN = -1, 1, 2
            // Reset seed (overwritten by LATMR)
            ISEED( 1 ) = 1
            ISEED( 2 ) = 1
            ISEED( 3 ) = 1
            ISEED( 4 ) = 1
            DO M = 32, MAXM, 23
               KLA = 0
               KUA = M - 1
               CALL CLATMR( M, M, 'S', ISEED, 'N', D, 6, ONE, CONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, KLA, KUA, ZERO, ONE, 'NO', A, MAXM, IWORK, IINFO )
               DO I = 1, M
                  A( I, I ) = A( I, I ) * VM( J )
               END DO
               ANRM = CLANGE( 'M', M, M, A, MAXM, DUM )
               DO N = 51, MAXN, 29
                  KLB = 0
                  KUB = N - 1
                  CALL CLATMR( N, N, 'S', ISEED, 'N', D, 6, ONE, CONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, KLB, KUB, ZERO, ONE, 'NO', B, MAXN, IWORK, IINFO )
                  DO I = 1, N
                     B( I, I ) = B( I, I ) * VM ( J )
                  END DO
                  BNRM = CLANGE( 'M', N, N, B, MAXN, DUM )
                  TNRM = MAX( ANRM, BNRM )
                  CALL CLATMR( M, N, 'S', ISEED, 'N', D, 6, ONE, CONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, M, N, ZERO, ONE, 'NO', C, MAXM, IWORK, IINFO )
                  DO ITRANA = 1, 2
                     IF( ITRANA.EQ.1 ) TRANA = 'N'                      IF( ITRANA.EQ.2 ) TRANA = 'C'
                     DO ITRANB = 1, 2
                        IF( ITRANB.EQ.1 ) TRANB = 'N'                         IF( ITRANB.EQ.2 ) TRANB = 'C'
                        KNT = KNT + 1
*
                        CALL CLACPY( 'All', M, N, C, MAXM, X, MAXM)
                        CALL CLACPY( 'All', M, N, C, MAXM, CC, MAXM)
                        CALL CTRSYL( TRANA, TRANB, ISGN, M, N,  A, MAXM, B, MAXN, X, MAXM, SCALE, IINFO )
                        IF( IINFO.NE.0 ) NINFO( 1 ) = NINFO( 1 ) + 1
                        XNRM = CLANGE( 'M', M, N, X, MAXM, DUM )
                        RMUL = CONE
                        IF( XNRM.GT.ONE .AND. TNRM.GT.ONE ) THEN
                           IF( XNRM.GT.BIGNUM / TNRM ) THEN
                              RMUL = CONE / MAX( XNRM, TNRM )
                           END IF
                        END IF
                        CALL CGEMM( TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM, -SCALE*RMUL, CC, MAXM )                         CALL CGEMM( 'N', TRANB, M, N, N, REAL( ISGN )*RMUL, X, MAXM, B, MAXN, CONE, CC, MAXM )
                        RES1 = CLANGE( 'M', M, N, CC, MAXM, DUM )
                        RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )                         IF( RES.GT.THRESH ) NFAIL( 1 ) = NFAIL( 1 ) + 1                         IF( RES.GT.RMAX( 1 ) ) RMAX( 1 ) = RES
*
                        CALL CLACPY( 'All', M, N, C, MAXM, X, MAXM )
                        CALL CLACPY( 'All', M, N, C, MAXM, CC, MAXM )
                        CALL CTRSYL3( TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM, SCALE3, SWORK, LDSWORK, INFO)
                        IF( INFO.NE.0 ) NINFO( 2 ) = NINFO( 2 ) + 1
                        XNRM = CLANGE( 'M', M, N, X, MAXM, DUM )
                        RMUL = CONE
                        IF( XNRM.GT.ONE .AND. TNRM.GT.ONE ) THEN
                           IF( XNRM.GT.BIGNUM / TNRM ) THEN
                              RMUL = CONE / MAX( XNRM, TNRM )
                           END IF
                        END IF
                        CALL CGEMM( TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM, -SCALE3*RMUL, CC, MAXM )                         CALL CGEMM( 'N', TRANB, M, N, N, REAL( ISGN )*RMUL, X, MAXM, B, MAXN, CONE, CC, MAXM )
                        RES1 = CLANGE( 'M', M, N, CC, MAXM, DUM )
                        RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
                        // Verify that TRSYL3 only flushes if TRSYL flushes (but
                       t // here may be cases where TRSYL3 avoid flushing).
                        IF( SCALE3.EQ.ZERO .AND. SCALE.GT.ZERO .OR.  IINFO.NE.INFO ) THEN
                           NFAIL( 3 ) = NFAIL( 3 ) + 1
                        END IF
                        IF( RES.GT.THRESH .OR. SISNAN( RES ) ) NFAIL( 2 ) = NFAIL( 2 ) + 1                         IF( RES.GT.RMAX( 2 ) ) RMAX( 2 ) = RES
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
*
      DEALLOCATE (A, STAT = AllocateStatus)
      DEALLOCATE (B, STAT = AllocateStatus)
      DEALLOCATE (C, STAT = AllocateStatus)
      DEALLOCATE (CC, STAT = AllocateStatus)
      DEALLOCATE (X, STAT = AllocateStatus)
      DEALLOCATE (SWORK, STAT = AllocateStatus)
*
      RETURN
*
      // End of CSYL01
*
      END
