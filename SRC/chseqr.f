      SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      // ..
      // .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..

      // ==== Matrices of order NTINY or smaller must be processed by
      // .    CLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare CLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0e0, 0.0e0 ), ONE = ( 1.0e0, 0.0e0 ) ;
      REAL               RZERO
      const              RZERO = 0.0e0 ;
      // ..
      // .. Local Arrays ..
      COMPLEX            HL( NL, NL ), WORKL( NL )
      // ..
      // .. Local Scalars ..
      int                KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      int                ILAENV;
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CLAHQR, CLAQR0, CLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // ==== Decode and check the input parameters. ====

      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = CMPLX( REAL( MAX( 1, N ) ), RZERO )
      LQUERY = LWORK == -1

      INFO = 0
      if ( .NOT.LSAME( JOB, 'E' ) && .NOT.WANTT ) {
         INFO = -1
      } else if ( .NOT.LSAME( COMPZ, 'N' ) && .NOT.WANTZ ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) {
         INFO = -5
      } else if ( LDH.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDZ.LT.1 .OR. ( WANTZ && LDZ.LT.MAX( 1, N ) ) ) {
         INFO = -10
      } else if ( LWORK.LT.MAX( 1, N ) && .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO != 0 ) {

         // ==== Quick return in case of invalid argument. ====

         xerbla('CHSEQR', -INFO );
         RETURN

      } else if ( N == 0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         RETURN

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         claqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = CMPLX( MAX( REAL( WORK( 1 ) ), REAL( MAX( 1, N ) ) ), RZERO )
         RETURN

      } else {

         // ==== copy eigenvalues isolated by CGEBAL ====

         if (ILO.GT.1) CALL CCOPY( ILO-1, H, LDH+1, W, 1 )          IF( IHI.LT.N ) CALL CCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 );

         // ==== Initialize Z, if requested ====

         if (INITZ) CALL CLASET( 'A', N, N, ZERO, ONE, Z, LDZ );

         // ==== Quick return if possible ====

         if ( ILO == IHI ) {
            W( ILO ) = H( ILO, ILO )
            RETURN
         }

         // ==== CLAHQR/CLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'CHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== CLAQR0 for big matrices; CLAHQR for small ones ====

         if ( N.GT.NMIN ) {
            claqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         } else {

            // ==== Small matrix ====

            clahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, INFO );

            if ( INFO.GT.0 ) {

               // ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds
               // .    when CLAHQR fails. ====

               KBOT = INFO

               if ( N.GE.NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call CLAQR0 directly. ====

                  claqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from CLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling CLAQR0. ====

                  clacpy('A', N, N, H, LDH, HL, NL );
                  HL( N+1, N ) = ZERO
                  claset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL );
                  claqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT .OR. INFO != 0 ) CALL CLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         IF( ( WANTT .OR. INFO != 0 ) && N.GT.2 ) CALL CLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = CMPLX( MAX( REAL( MAX( 1, N ) ), REAL( WORK( 1 ) ) ), RZERO )
      }

      // ==== End of CHSEQR ====

      }
