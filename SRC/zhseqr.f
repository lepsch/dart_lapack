      SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..

      // ==== Matrices of order NTINY or smaller must be processed by
      // .    ZLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      PARAMETER          ( NTINY = 15 )

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare ZLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      PARAMETER          ( NL = 49 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), ONE = ( 1.0d0, 0.0d0 ) )
      double             RZERO;
      PARAMETER          ( RZERO = 0.0d0 )
      // ..
      // .. Local Arrays ..
      COMPLEX*16         HL( NL, NL ), WORKL( NL )
      // ..
      // .. Local Scalars ..
      int                KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      int                ILAENV;
      bool               LSAME;
      // EXTERNAL ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZLAHQR, ZLAQR0, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // ==== Decode and check the input parameters. ====

      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = DCMPLX( DBLE( MAX( 1, N ) ), RZERO )
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
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF

      IF( INFO.NE.0 ) THEN

         // ==== Quick return in case of invalid argument. ====

         CALL XERBLA( 'ZHSEQR', -INFO )
         RETURN

      ELSE IF( N.EQ.0 ) THEN

         // ==== Quick return in case N = 0; nothing to do. ====

         RETURN

      ELSE IF( LQUERY ) THEN

         // ==== Quick return in case of a workspace query ====

         CALL ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = DCMPLX( MAX( DBLE( WORK( 1 ) ), DBLE( MAX( 1, N ) ) ), RZERO )
         RETURN

      ELSE

         // ==== copy eigenvalues isolated by ZGEBAL ====

         IF( ILO.GT.1 ) CALL ZCOPY( ILO-1, H, LDH+1, W, 1 )          IF( IHI.LT.N ) CALL ZCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 )

         // ==== Initialize Z, if requested ====

         IF( INITZ ) CALL ZLASET( 'A', N, N, ZERO, ONE, Z, LDZ )

         // ==== Quick return if possible ====

         IF( ILO.EQ.IHI ) THEN
            W( ILO ) = H( ILO, ILO )
            RETURN
         END IF

         // ==== ZLAHQR/ZLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'ZHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====

         IF( N.GT.NMIN ) THEN
            CALL ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
         ELSE

            // ==== Small matrix ====

            CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, INFO )

            IF( INFO.GT.0 ) THEN

               // ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
               // .    when ZLAHQR fails. ====

               KBOT = INFO

               IF( N.GE.NL ) THEN

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call ZLAQR0 directly. ====

                  CALL ZLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )

               ELSE

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from ZLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling ZLAQR0. ====

                  CALL ZLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL ZLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL )                   CALL ZLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT .OR. INFO.NE.0 ) CALL ZLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF

         // ==== Clear out the trash, if necessary. ====

         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) CALL ZLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = DCMPLX( MAX( DBLE( MAX( 1, N ) ), DBLE( WORK( 1 ) ) ), RZERO )
      END IF

      // ==== End of ZHSEQR ====

      }
