      SUBROUTINE ZGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWORK, ITAU, IWORK, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SAFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD, ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVSL, 'N' ) ) {
         IJOBVL = 1
         ILVSL = .FALSE.
      } else if ( LSAME( JOBVSL, 'V' ) ) {
         IJOBVL = 2
         ILVSL = .TRUE.
      } else {
         IJOBVL = -1
         ILVSL = .FALSE.
      }

      if ( LSAME( JOBVSR, 'N' ) ) {
         IJOBVR = 1
         ILVSR = .FALSE.
      } else if ( LSAME( JOBVSR, 'V' ) ) {
         IJOBVR = 2
         ILVSR = .TRUE.
      } else {
         IJOBVR = -1
         ILVSR = .FALSE.
      }

      // Test the input arguments

      LWKMIN = MAX( 2*N, 1 )
      LWKOPT = LWKMIN
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      INFO = 0
      if ( IJOBVL.LE.0 ) {
         INFO = -1
      } else if ( IJOBVR.LE.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) {
         INFO = -11
      } else if ( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) {
         INFO = -13
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -15
      }

      if ( INFO.EQ.0 ) {
         NB1 = ILAENV( 1, 'ZGEQRF', ' ', N, N, -1, -1 )
         NB2 = ILAENV( 1, 'ZUNMQR', ' ', N, N, N, -1 )
         NB3 = ILAENV( 1, 'ZUNGQR', ' ', N, N, N, -1 )
         NB = MAX( NB1, NB2, NB3 )
         LOPT = N*( NB+1 )
         WORK( 1 ) = LOPT
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZGEGS ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SAFMIN = DLAMCH( 'S' )
      SMLNUM = N*SAFMIN / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      }

      if ( ILASCL ) {
         CALL ZLASCL( 'G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      }

      if ( ILBSCL ) {
         CALL ZLASCL( 'G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

      // Permute the matrix to make it more nearly triangular

      ILEFT = 1
      IRIGHT = N + 1
      IRWORK = IRIGHT + N
      IWORK = 1
      CALL ZGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWORK ), IINFO )
      if ( IINFO.NE.0 ) {
         INFO = N + 1
         GO TO 10
      }

      // Reduce B to triangular form, and initialize VSL and/or VSR

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) {
         INFO = N + 2
         GO TO 10
      }

      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) {
         INFO = N + 3
         GO TO 10
      }

      if ( ILVSL ) {
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL )          CALL ZUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )
         IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         if ( IINFO.NE.0 ) {
            INFO = N + 4
            GO TO 10
         }
      }

      IF( ILVSR ) CALL ZLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )

      // Reduce to generalized Hessenberg form

      CALL ZGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IINFO )
      if ( IINFO.NE.0 ) {
         INFO = N + 5
         GO TO 10
      }

      // Perform QZ algorithm, computing Schur vectors if desired

      IWORK = ITAU
      CALL ZHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWORK ), LWORK+1-IWORK, RWORK( IRWORK ), IINFO )
      IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) {
         if ( IINFO.GT.0 .AND. IINFO.LE.N ) {
            INFO = IINFO
         } else if ( IINFO.GT.N .AND. IINFO.LE.2*N ) {
            INFO = IINFO - N
         } else {
            INFO = N + 6
         }
         GO TO 10
      }

      // Apply permutation to VSL and VSR

      if ( ILVSL ) {
         CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 7
            GO TO 10
         }
      }
      if ( ILVSR ) {
         CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 8
            GO TO 10
         }
      }

      // Undo scaling

      if ( ILASCL ) {
         CALL ZLASCL( 'U', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
         CALL ZLASCL( 'G', -1, -1, ANRMTO, ANRM, N, 1, ALPHA, N, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

      if ( ILBSCL ) {
         CALL ZLASCL( 'U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
         CALL ZLASCL( 'G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

   10 CONTINUE
      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZGEGS

      }
