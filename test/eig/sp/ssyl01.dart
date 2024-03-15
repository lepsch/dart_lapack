
// -- LAPACK test routine --
      void ssyl01(final int THRESH, final int NFAIL, final int RMAX, final int NINFO, final int KNT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KNT;
      double               THRESH;
      int                NFAIL( 3 ), NINFO( 2 );
      double               RMAX( 2 );
      // ..

// =====================================================================
      // ..
      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                MAXM, MAXN, LDSWORK;
      const              MAXM = 101, MAXN = 138, LDSWORK = 18 ;
      String             TRANA, TRANB;
      int                I, INFO, IINFO, ISGN, ITRANA, ITRANB, J, KLA, KUA, KLB, KUB, LIWORK, M, N;
      double               ANRM, BNRM, BIGNUM, EPS, RES, RES1, RMUL, SCALE, SCALE3, SMLNUM, TNRM, XNRM;
      double               DUML( MAXM ), DUMR( MAXN ), D( max( MAXM, MAXN ) ), DUM( MAXN ), VM( 2 );
      int                ISEED( 4 ), IWORK( MAXM + MAXN + 2 );
      // ..
      // .. Allocatable Arrays ..
      int                AllocateStatus;
      double, DIMENSION(:,:), ALLOCATABLE :: A, B, C, CC, X, SWORK;
      // ..
      // .. External Functions ..
      //- bool               SISNAN;
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SISNAN, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLATMR, SLACPY, SGEMM, STRSYL, STRSYL3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX
      // ..
      // .. Allocate memory dynamically ..
      ALLOCATE ( A( MAXM, MAXM ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( B( MAXN, MAXN ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( C( MAXM, MAXN ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( CC( MAXM, MAXN ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( X( MAXM, MAXN ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";
      ALLOCATE ( SWORK( LDSWORK, 54 ), STAT = AllocateStatus );
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***";

      // Get machine parameters

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      VM[1] = ONE;
      VM[2] = 0.05;

      // Begin test loop

      NINFO[1] = 0;
      NINFO[2] = 0;
      NFAIL[1] = 0;
      NFAIL[2] = 0;
      NFAIL[3] = 0;
      RMAX[1] = ZERO;
      RMAX[2] = ZERO;
      KNT = 0;
      for (I = 1; I <= 4; I++) {
         ISEED[I] = 1;
      }
      SCALE = ONE;
      SCALE3 = ONE;
      LIWORK = MAXM + MAXN + 2;
      for (J = 1; J <= 2; J++) {
         for (ISGN = -1; 2 < 0 ? ISGN >= 1 : ISGN <= 1; ISGN += 2) {
            // Reset seed (overwritten by LATMR)
            for (I = 1; I <= 4; I++) {
               ISEED[I] = 1;
            }
            for (M = 32; M <= MAXM; M += 71) {
               KLA = 0;
               KUA = M - 1;
               slatmr(M, M, 'S', ISEED, 'N', D, 6, ONE, ONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, KLA, KUA, ZERO, ONE, 'NO', A, MAXM, IWORK, IINFO );
               for (I = 1; I <= M; I++) {
                  A[I][I] = A( I, I ) * VM( J );
               }
               ANRM = SLANGE( 'M', M, M, A, MAXM, DUM );
               for (N = 51; N <= MAXN; N += 47) {
                  KLB = 0;
                  KUB = N - 1;
                  slatmr(N, N, 'S', ISEED, 'N', D, 6, ONE, ONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, KLB, KUB, ZERO, ONE, 'NO', B, MAXN, IWORK, IINFO );
                  BNRM = SLANGE( 'M', N, N, B, MAXN, DUM );
                  TNRM = max( ANRM, BNRM );
                  slatmr(M, N, 'S', ISEED, 'N', D, 6, ONE, ONE, 'T', 'N', DUML, 1, ONE, DUMR, 1, ONE, 'N', IWORK, M, N, ZERO, ONE, 'NO', C, MAXM, IWORK, IINFO );
                  for (ITRANA = 1; ITRANA <= 2; ITRANA++) {
                     if ( ITRANA == 1 ) {
                        TRANA = 'N';
                     }
                     if ( ITRANA == 2 ) {
                        TRANA = 'T';
                     }
                     for (ITRANB = 1; ITRANB <= 2; ITRANB++) {
                        if ( ITRANB == 1 ) {
                           TRANB = 'N';
                        }
                        if ( ITRANB == 2 ) {
                           TRANB = 'T';
                        }
                        KNT = KNT + 1;

                        slacpy('All', M, N, C, MAXM, X, MAXM);
                        slacpy('All', M, N, C, MAXM, CC, MAXM);
                        strsyl(TRANA, TRANB, ISGN, M, N,  A, MAXM, B, MAXN, X, MAXM, SCALE, IINFO );
                        if (IINFO != 0) NINFO( 1 ) = NINFO( 1 ) + 1;
                        XNRM = SLANGE( 'M', M, N, X, MAXM, DUM );
                        RMUL = ONE;
                        if ( XNRM > ONE && TNRM > ONE ) {
                           if ( XNRM > BIGNUM / TNRM ) {
                              RMUL = ONE / max( XNRM, TNRM );
                           }
                        }
                        sgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM, -SCALE*RMUL, C, MAXM );
                        sgemm('N', TRANB, M, N, N, double( ISGN )*RMUL, X, MAXM, B, MAXN, ONE, C, MAXM );
                        RES1 = SLANGE( 'M', M, N, C, MAXM, DUM );
                        RES = RES1 / max( SMLNUM, SMLNUM*XNRM, ( ( RMUL*TNRM )*EPS )*XNRM )                         IF( RES > THRESH ) NFAIL( 1 ) = NFAIL( 1 ) + 1                         IF( RES > RMAX( 1 ) ) RMAX( 1 ) = RES;

                        slacpy('All', M, N, C, MAXM, X, MAXM );
                        slacpy('All', M, N, C, MAXM, CC, MAXM );
                        strsyl3(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM, SCALE3, IWORK, LIWORK, SWORK, LDSWORK, INFO);
                        if (INFO != 0) NINFO( 2 ) = NINFO( 2 ) + 1;
                        XNRM = SLANGE( 'M', M, N, X, MAXM, DUM );
                        RMUL = ONE;
                        if ( XNRM > ONE && TNRM > ONE ) {
                           if ( XNRM > BIGNUM / TNRM ) {
                              RMUL = ONE / max( XNRM, TNRM );
                           }
                        }
                        sgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM, -SCALE3*RMUL, CC, MAXM );
                        sgemm('N', TRANB, M, N, N, double( ISGN )*RMUL, X, MAXM, B, MAXN, ONE, CC, MAXM );
                        RES1 = SLANGE( 'M', M, N, CC, MAXM, DUM );
                        RES = RES1 / max( SMLNUM, SMLNUM*XNRM, ( ( RMUL*TNRM )*EPS )*XNRM );
                        // Verify that TRSYL3 only flushes if TRSYL flushes (but
                        // there may be cases where TRSYL3 avoid flushing).
                        if ( SCALE3 == ZERO && SCALE > ZERO || IINFO != INFO ) {
                           NFAIL[3] = NFAIL( 3 ) + 1;
                        }
                        if[RES > THRESH || SISNAN( RES ) ) NFAIL( 2] = NFAIL( 2 ) + 1;
                        IF[RES > RMAX( 2 ) ) RMAX( 2] = RES;
                     }
                  }
               }
            }
         }
      }

      DEALLOCATE (A, STAT = AllocateStatus);
      DEALLOCATE (B, STAT = AllocateStatus);
      DEALLOCATE (C, STAT = AllocateStatus);
      DEALLOCATE (CC, STAT = AllocateStatus);
      DEALLOCATE (X, STAT = AllocateStatus);
      DEALLOCATE (SWORK, STAT = AllocateStatus);

      }