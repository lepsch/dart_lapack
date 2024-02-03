      SUBROUTINE CGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWORK, ITAU, IWORK, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SAFMIN, SMLNUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ, CLACPY, CLASCL, CLASET, CUNGQR, CUNMQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANGE, SLAMCH
      // EXTERNAL ILAENV, LSAME, CLANGE, SLAMCH
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
         NB1 = ILAENV( 1, 'CGEQRF', ' ', N, N, -1, -1 )
         NB2 = ILAENV( 1, 'CUNMQR', ' ', N, N, N, -1 )
         NB3 = ILAENV( 1, 'CUNGQR', ' ', N, N, N, -1 )
         NB = MAX( NB1, NB2, NB3 )
         LOPT = N*(NB+1)
         WORK( 1 ) = LOPT
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGEGS ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = SLAMCH( 'E' )*SLAMCH( 'B' )
      SAFMIN = SLAMCH( 'S' )
      SMLNUM = N*SAFMIN / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      }

      if ( ILASCL ) {
         CALL CLASCL( 'G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      }

      if ( ILBSCL ) {
         CALL CLASCL( 'G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO )
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
      CALL CGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWORK ), IINFO )
      if ( IINFO.NE.0 ) {
         INFO = N + 1
         GO TO 10
      }

      // Reduce B to triangular form, and initialize VSL and/or VSR

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL CGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) {
         INFO = N + 2
         GO TO 10
      }

      CALL CUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      if ( IINFO.NE.0 ) {
         INFO = N + 3
         GO TO 10
      }

      if ( ILVSL ) {
         CALL CLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         CALL CLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL )          CALL CUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )
         IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         if ( IINFO.NE.0 ) {
            INFO = N + 4
            GO TO 10
         }
      }

      IF( ILVSR ) CALL CLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )

      // Reduce to generalized Hessenberg form

      CALL CGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IINFO )
      if ( IINFO.NE.0 ) {
         INFO = N + 5
         GO TO 10
      }

      // Perform QZ algorithm, computing Schur vectors if desired

      IWORK = ITAU
      CALL CHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWORK ), LWORK+1-IWORK, RWORK( IRWORK ), IINFO )
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
         CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 7
            GO TO 10
         }
      }
      if ( ILVSR ) {
         CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 8
            GO TO 10
         }
      }

      // Undo scaling

      if ( ILASCL ) {
         CALL CLASCL( 'U', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
         CALL CLASCL( 'G', -1, -1, ANRMTO, ANRM, N, 1, ALPHA, N, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

      if ( ILBSCL ) {
         CALL CLASCL( 'U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
         CALL CLASCL( 'G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = N + 9
            RETURN
         }
      }

   10 CONTINUE
      WORK( 1 ) = LWKOPT

      RETURN

      // End of CGEGS

      }
