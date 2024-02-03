      SUBROUTINE SHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      // ..
      // .. Array Arguments ..
      REAL               H( LDH, * ), WI( * ), WORK( * ), WR( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..

      // ==== Matrices of order NTINY or smaller must be processed by
      // .    SLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare SLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      REAL               ZERO, ONE
      const              ZERO = 0.0e0, ONE = 1.0e0 ;
      // ..
      // .. Local Arrays ..
      REAL               HL( NL, NL ), WORKL( NL )
      // ..
      // .. Local Scalars ..
      int                I, KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      int                ILAENV;
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLAHQR, SLAQR0, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // ==== Decode and check the input parameters. ====

      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = SROUNDUP_LWORK( MAX( 1, N ) )
      LQUERY = LWORK.EQ.-1

      INFO = 0
      if ( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) {
         INFO = -1
      } else if ( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) {
         INFO = -5
      } else if ( LDH.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) {
         INFO = -11
      } else if ( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO.NE.0 ) {

         // ==== Quick return in case of invalid argument. ====

         CALL XERBLA( 'SHSEQR', -INFO )
         RETURN

      } else if ( N.EQ.0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         RETURN

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         CALL SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = MAX( REAL( MAX( 1, N ) ), WORK( 1 ) )
         RETURN

      } else {

         // ==== copy eigenvalues isolated by SGEBAL ====

         DO 10 I = 1, ILO - 1
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   10    CONTINUE
         DO 20 I = IHI + 1, N
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   20    CONTINUE

         // ==== Initialize Z, if requested ====

         IF( INITZ ) CALL SLASET( 'A', N, N, ZERO, ONE, Z, LDZ )

         // ==== Quick return if possible ====

         if ( ILO.EQ.IHI ) {
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
         }

         // ==== SLAHQR/SLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'SHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== SLAQR0 for big matrices; SLAHQR for small ones ====

         if ( N.GT.NMIN ) {
            CALL SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         } else {

            // ==== Small matrix ====

            CALL SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, INFO )

            if ( INFO.GT.0 ) {

               // ==== A rare SLAHQR failure!  SLAQR0 sometimes succeeds
               // .    when SLAHQR fails. ====

               KBOT = INFO

               if ( N.GE.NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call SLAQR0 directly. ====

                  CALL SLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from SLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling SLAQR0. ====

                  CALL SLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL SLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL )                   CALL SLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT .OR. INFO.NE.0 ) CALL SLACPY( 'A', N, N, HL, NL, H, LDH )
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) CALL SLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = MAX( REAL( MAX( 1, N ) ), WORK( 1 ) )
      }

      // ==== End of SHSEQR ====

      }
