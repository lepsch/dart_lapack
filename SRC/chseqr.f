      void chseqr(JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      // ..
      // .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

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
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RZERO;
      const              RZERO = 0.0 ;
      // ..
      // .. Local Arrays ..
      COMPLEX            HL( NL, NL ), WORKL( NL );
      // ..
      // .. Local Scalars ..
      int                KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- bool               LSAME;
      //- REAL               SROUNDUP_LWORK;
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

      WANTT = LSAME( JOB, 'S' );
      INITZ = LSAME( COMPZ, 'I' );
      WANTZ = INITZ || LSAME( COMPZ, 'V' );
      WORK( 1 ) = CMPLX( REAL( max( 1, N ) ), RZERO );
      LQUERY = LWORK == -1;

      INFO = 0;
      if ( !LSAME( JOB, 'E' ) && !WANTT ) {
         INFO = -1;
      } else if ( !LSAME( COMPZ, 'N' ) && !WANTZ ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -4;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -5;
      } else if ( LDH < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < max( 1, N ) ) ) {
         INFO = -10;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO != 0 ) {

         // ==== Quick return in case of invalid argument. ====

         xerbla('CHSEQR', -INFO );
         return;

      } else if ( N == 0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         return;

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         claqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = CMPLX( max( REAL( WORK( 1 ) ), REAL( max( 1, N ) ) ), RZERO );
         return;

      } else {

         // ==== copy eigenvalues isolated by CGEBAL ====

         if (ILO > 1) ccopy( ILO-1, H, LDH+1, W, 1 );
         IF( IHI < N ) ccopy( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 );

         // ==== Initialize Z, if requested ====

         if (INITZ) claset( 'A', N, N, ZERO, ONE, Z, LDZ );

         // ==== Quick return if possible ====

         if ( ILO == IHI ) {
            W( ILO ) = H( ILO, ILO );
            return;
         }

         // ==== CLAHQR/CLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'CHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK );
         NMIN = max( NTINY, NMIN );

         // ==== CLAQR0 for big matrices; CLAHQR for small ones ====

         if ( N > NMIN ) {
            claqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         } else {

            // ==== Small matrix ====

            clahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, INFO );

            if ( INFO > 0 ) {

               // ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds
               // .    when CLAHQR fails. ====

               KBOT = INFO;

               if ( N >= NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call CLAQR0 directly. ====

                  claqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from CLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling CLAQR0. ====

                  clacpy('A', N, N, H, LDH, HL, NL );
                  HL( N+1, N ) = ZERO;
                  claset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL );
                  claqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT || INFO != 0 ) CALL CLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         if( ( WANTT || INFO != 0 ) && N > 2 ) claset( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH );

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = CMPLX( max( REAL( max( 1, N ) ), REAL( WORK( 1 ) ) ), RZERO );
      }

      // ==== End of CHSEQR ====

      }
