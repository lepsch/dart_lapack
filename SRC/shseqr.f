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
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF

      IF( INFO.NE.0 ) THEN

         // ==== Quick return in case of invalid argument. ====

         CALL XERBLA( 'SHSEQR', -INFO )
         RETURN

      ELSE IF( N.EQ.0 ) THEN

         // ==== Quick return in case N = 0; nothing to do. ====

         RETURN

      ELSE IF( LQUERY ) THEN

         // ==== Quick return in case of a workspace query ====

         CALL SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = MAX( REAL( MAX( 1, N ) ), WORK( 1 ) )
         RETURN

      ELSE

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

         IF( ILO.EQ.IHI ) THEN
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
         END IF

         // ==== SLAHQR/SLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'SHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== SLAQR0 for big matrices; SLAHQR for small ones ====

         IF( N.GT.NMIN ) THEN
            CALL SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         ELSE

            // ==== Small matrix ====

            CALL SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, INFO )

            IF( INFO.GT.0 ) THEN

               // ==== A rare SLAHQR failure!  SLAQR0 sometimes succeeds
               // .    when SLAHQR fails. ====

               KBOT = INFO

               IF( N.GE.NL ) THEN

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call SLAQR0 directly. ====

                  CALL SLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )

               ELSE

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from SLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling SLAQR0. ====

                  CALL SLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL SLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL )                   CALL SLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT .OR. INFO.NE.0 ) CALL SLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF

         // ==== Clear out the trash, if necessary. ====

         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) CALL SLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = MAX( REAL( MAX( 1, N ) ), WORK( 1 ) )
      END IF

      // ==== End of SHSEQR ====

      }
