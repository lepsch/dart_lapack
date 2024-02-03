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
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare ZLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0d0, 0.0d0 ), ONE = ( 1.0d0, 0.0d0 ) ;
      double             RZERO;
      const              RZERO = 0.0d0 ;
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
         INFO = -10
      } else if ( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO.NE.0 ) {

         // ==== Quick return in case of invalid argument. ====

         xerbla('ZHSEQR', -INFO );
         RETURN

      } else if ( N.EQ.0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         RETURN

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         zlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = DCMPLX( MAX( DBLE( WORK( 1 ) ), DBLE( MAX( 1, N ) ) ), RZERO )
         RETURN

      } else {

         // ==== copy eigenvalues isolated by ZGEBAL ====

         IF( ILO.GT.1 ) CALL ZCOPY( ILO-1, H, LDH+1, W, 1 )          IF( IHI.LT.N ) CALL ZCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 )

         // ==== Initialize Z, if requested ====

         IF( INITZ ) CALL ZLASET( 'A', N, N, ZERO, ONE, Z, LDZ )

         // ==== Quick return if possible ====

         if ( ILO.EQ.IHI ) {
            W( ILO ) = H( ILO, ILO )
            RETURN
         }

         // ==== ZLAHQR/ZLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'ZHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====

         if ( N.GT.NMIN ) {
            zlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         } else {

            // ==== Small matrix ====

            zlahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, INFO );

            if ( INFO.GT.0 ) {

               // ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
               // .    when ZLAHQR fails. ====

               KBOT = INFO

               if ( N.GE.NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call ZLAQR0 directly. ====

                  zlaqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from ZLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling ZLAQR0. ====

                  zlacpy('A', N, N, H, LDH, HL, NL );
                  HL( N+1, N ) = ZERO
                  zlaset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL )                   CALL ZLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT .OR. INFO.NE.0 ) CALL ZLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) CALL ZLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = DCMPLX( MAX( DBLE( MAX( 1, N ) ), DBLE( WORK( 1 ) ) ), RZERO )
      }

      // ==== End of ZHSEQR ====

      }
