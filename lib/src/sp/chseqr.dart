      void chseqr(final int JOB, final int COMPZ, final int N, final int ILO, final int IHI, final Matrix<double> H_, final int LDH, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final H = H_.dim();
  final Z = Z_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      Complex            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..


      // ==== Matrices of order NTINY or smaller must be processed by
      // .    CLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare CLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ilaenv(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      Complex            HL( NL, NL ), WORKL( NL );
      int                KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CLAHQR, CLAQR0, CLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL

      // ==== Decode and check the input parameters. ====

      WANTT = lsame( JOB, 'S' );
      INITZ = lsame( COMPZ, 'I' );
      WANTZ = INITZ || lsame( COMPZ, 'V' );
      WORK[1] = CMPLX( double( max( 1, N ) ), RZERO );
      LQUERY = LWORK == -1;

      INFO = 0;
      if ( !lsame( JOB, 'E' ) && !WANTT ) {
         INFO = -1;
      } else if ( !lsame( COMPZ, 'N' ) && !WANTZ ) {
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
         WORK[1] = CMPLX( max( REAL( WORK( 1 ) ), double( max( 1, N ) ) ), RZERO );
         return;

      } else {

         // ==== copy eigenvalues isolated by CGEBAL ====

         if (ILO > 1) ccopy( ILO-1, H, LDH+1, W, 1 );
         IF( IHI < N ) ccopy( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 );

         // ==== Initialize Z, if requested ====

         if (INITZ) claset( 'A', N, N, ZERO, ONE, Z, LDZ );

         // ==== Quick return if possible ====

         if ( ILO == IHI ) {
            W[ILO] = H( ILO, ILO );
            return;
         }

         // ==== CLAHQR/CLAQR0 crossover point ====

         NMIN = ilaenv( 12, 'CHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK );
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
                  HL[N+1][N] = ZERO;
                  claset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL );
                  claqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT || INFO != 0 ) CALL CLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         if( ( WANTT || INFO != 0 ) && N > 2 ) claset( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH );

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK[1] = CMPLX( max( REAL( max( 1, N ) ), double( WORK( 1 ) ) ), RZERO );
      }

      // ==== End of CHSEQR ====

      }
