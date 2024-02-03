
*  -- LAPACK test routine --
      SUBROUTINE CSYL01( THRESH, NFAIL, RMAX, NINFO, KNT )
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                NFAIL( 3 ), NINFO( 2 );
      REAL               RMAX( 2 )
      // ..

*  =====================================================================
      // ..
      // .. Parameters ..
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      REAL               ONE, ZERO
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      int                MAXM, MAXN, LDSWORK;
      const              MAXM = 101, MAXN = 138, LDSWORK = 18 ;
      // ..
      // .. Local Scalars ..
      String             TRANA, TRANB;
      int                I, INFO, IINFO, ISGN, ITRANA, ITRANB, J, KLA, KUA, KLB, KUB, M, N;
      REAL               ANRM, BNRM, BIGNUM, EPS, RES, RES1, SCALE, SCALE3, SMLNUM, TNRM, XNRM;
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
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( B( MAXN, MAXN ), STAT = AllocateStatus )
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( C( MAXM, MAXN ), STAT = AllocateStatus )
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( CC( MAXM, MAXN ), STAT = AllocateStatus )
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( X( MAXM, MAXN ), STAT = AllocateStatus )
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( SWORK( LDSWORK, 54 ), STAT = AllocateStatus )
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Expect INFO = 0
      VM( 1 ) = ONE
      // Expect INFO = 1
      VM( 2 ) = 0.5E+0

      // Begin test loop

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
      for (J = 1; J <= 2; J++) {
         DO ISGN = -1, 1, 2
            // Reset seed (overwritten by LATMR)
            ISEED( 1 ) = 1
            ISEED( 2 ) = 1
            ISEED( 3 ) = 1
            ISEED( 4 ) = 1
            DO M = 32, MAXM, 23
               KLA = 0
               KUA = M - 1
               clatmr(M, M, 'S', ISEED, 'N', D, 6, ONE, CONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, KLA, KUA, ZERO, ONE, 'NO', A, MAXM, IWORK, IINFO );
               for (I = 1; I <= M; I++) {
                  A( I, I ) = A( I, I ) * VM( J )
               }
               ANRM = CLANGE( 'M', M, M, A, MAXM, DUM )
               DO N = 51, MAXN, 29
                  KLB = 0
                  KUB = N - 1
                  clatmr(N, N, 'S', ISEED, 'N', D, 6, ONE, CONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, KLB, KUB, ZERO, ONE, 'NO', B, MAXN, IWORK, IINFO );
                  for (I = 1; I <= N; I++) {
                     B( I, I ) = B( I, I ) * VM ( J )
                  }
                  BNRM = CLANGE( 'M', N, N, B, MAXN, DUM )
                  TNRM = MAX( ANRM, BNRM )
                  clatmr(M, N, 'S', ISEED, 'N', D, 6, ONE, CONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, M, N, ZERO, ONE, 'NO', C, MAXM, IWORK, IINFO );
                  for (ITRANA = 1; ITRANA <= 2; ITRANA++) {
                     if (ITRANA == 1) TRANA = 'N'                      IF( ITRANA == 2 ) TRANA = 'C';
                     for (ITRANB = 1; ITRANB <= 2; ITRANB++) {
                        if (ITRANB == 1) TRANB = 'N'                         IF( ITRANB == 2 ) TRANB = 'C';
                        KNT = KNT + 1

                        clacpy('All', M, N, C, MAXM, X, MAXM);
                        clacpy('All', M, N, C, MAXM, CC, MAXM);
                        ctrsyl(TRANA, TRANB, ISGN, M, N,  A, MAXM, B, MAXN, X, MAXM, SCALE, IINFO );
                        if (IINFO.NE.0) NINFO( 1 ) = NINFO( 1 ) + 1;
                        XNRM = CLANGE( 'M', M, N, X, MAXM, DUM )
                        RMUL = CONE
                        if ( XNRM.GT.ONE .AND. TNRM.GT.ONE ) {
                           if ( XNRM.GT.BIGNUM / TNRM ) {
                              RMUL = CONE / MAX( XNRM, TNRM )
                           }
                        }
                        cgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM, -SCALE*RMUL, CC, MAXM );
                        cgemm('N', TRANB, M, N, N, REAL( ISGN )*RMUL, X, MAXM, B, MAXN, CONE, CC, MAXM );
                        RES1 = CLANGE( 'M', M, N, CC, MAXM, DUM )
                        RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )                         IF( RES.GT.THRESH ) NFAIL( 1 ) = NFAIL( 1 ) + 1                         IF( RES.GT.RMAX( 1 ) ) RMAX( 1 ) = RES

                        clacpy('All', M, N, C, MAXM, X, MAXM );
                        clacpy('All', M, N, C, MAXM, CC, MAXM );
                        ctrsyl3(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM, SCALE3, SWORK, LDSWORK, INFO);
                        if (INFO.NE.0) NINFO( 2 ) = NINFO( 2 ) + 1;
                        XNRM = CLANGE( 'M', M, N, X, MAXM, DUM )
                        RMUL = CONE
                        if ( XNRM.GT.ONE .AND. TNRM.GT.ONE ) {
                           if ( XNRM.GT.BIGNUM / TNRM ) {
                              RMUL = CONE / MAX( XNRM, TNRM )
                           }
                        }
                        cgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM, -SCALE3*RMUL, CC, MAXM );
                        cgemm('N', TRANB, M, N, N, REAL( ISGN )*RMUL, X, MAXM, B, MAXN, CONE, CC, MAXM );
                        RES1 = CLANGE( 'M', M, N, CC, MAXM, DUM )
                        RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
                        // Verify that TRSYL3 only flushes if TRSYL flushes (but
                        // there may be cases where TRSYL3 avoid flushing).
                        if ( SCALE3 == ZERO .AND. SCALE.GT.ZERO .OR.  IINFO.NE.INFO ) {
                           NFAIL( 3 ) = NFAIL( 3 ) + 1
                        }
                        IF( RES.GT.THRESH .OR. SISNAN( RES ) ) NFAIL( 2 ) = NFAIL( 2 ) + 1                         IF( RES.GT.RMAX( 2 ) ) RMAX( 2 ) = RES
                     }
                  }
               }
            }
         }
      }

      DEALLOCATE (A, STAT = AllocateStatus)
      DEALLOCATE (B, STAT = AllocateStatus)
      DEALLOCATE (C, STAT = AllocateStatus)
      DEALLOCATE (CC, STAT = AllocateStatus)
      DEALLOCATE (X, STAT = AllocateStatus)
      DEALLOCATE (SWORK, STAT = AllocateStatus)

      RETURN

      // End of CSYL01

      }
