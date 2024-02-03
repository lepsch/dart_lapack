      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      int                NTINY
      PARAMETER          ( NTINY = 15 )
*
*     ==== NL allocates some local workspace to help small matrices
*     .    through a rare DLAHQR failure.  NL > NTINY = 15 is
*     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
*     .    allows up to six simultaneous shifts and a 16-by-16
*     .    deflation window.  ====
      int                NL
      PARAMETER          ( NL = 49 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   HL( NL, NL ), WORKL( NL )
*     ..
*     .. Local Scalars ..
      int                I, KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
*     ..
*     .. External Functions ..
      int                ILAENV
      LOGICAL            LSAME
      EXTERNAL           ILAENV, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLAHQR, DLAQR0, DLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     ==== Decode and check the input parameters. ====
*
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = DBLE( MAX( 1, N ) )
      LQUERY = LWORK.EQ.-1
*
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
*
      IF( INFO.NE.0 ) THEN
*
*        ==== Quick return in case of invalid argument. ====
*
         CALL XERBLA( 'DHSEQR', -INFO )
         RETURN
*
      ELSE IF( N.EQ.0 ) THEN
*
*        ==== Quick return in case N = 0; nothing to do. ====
*
         RETURN
*
      ELSE IF( LQUERY ) THEN
*
*        ==== Quick return in case of a workspace query ====
*
         CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
         RETURN
*
      ELSE
*
*        ==== copy eigenvalues isolated by DGEBAL ====
*
         DO 10 I = 1, ILO - 1
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   10    CONTINUE
         DO 20 I = IHI + 1, N
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   20    CONTINUE
*
*        ==== Initialize Z, if requested ====
*
         IF( INITZ ) CALL DLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
*
*        ==== Quick return if possible ====
*
         IF( ILO.EQ.IHI ) THEN
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         NMIN = ILAENV( 12, 'DHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
*
         IF( N.GT.NMIN ) THEN
            CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         ELSE
*
*           ==== Small matrix ====
*
            CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, INFO )
*
            IF( INFO.GT.0 ) THEN
*
*              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
*              .    when DLAHQR fails. ====
*
               KBOT = INFO
*
               IF( N.GE.NL ) THEN
*
*                 ==== Larger matrices have enough subdiagonal scratch
*                 .    space to call DLAQR0 directly. ====
*
                  CALL DLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
*
               ELSE
*
*                 ==== Tiny matrices don't have enough subdiagonal
*                 .    scratch space to benefit from DLAQR0.  Hence,
*                 .    tiny matrices must be copied into a larger
*                 .    array before calling DLAQR0. ====
*
                  CALL DLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL DLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL )                   CALL DLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT .OR. INFO.NE.0 ) CALL DLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF
*
*        ==== Clear out the trash, if necessary. ====
*
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
*
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
*
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
      END IF
*
*     ==== End of DHSEQR ====
*
      END
